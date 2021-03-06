//
// Created by Furukawa on 2016/05/16.
//
//=================================================================================//
//                                                                                 //
//  CT Reconstruct using Generic Algorithm                                         //
//                                                                                 //
//  Copyright (C) 2011-2016 Taichi-Furukawa                                        //
//                                                                                 //
//  This file is a portion of the ilab library. It is distributed under the MIT    //
//  License, available in the root of this distribution and at the following URL.  //
//  http://opensource.org/licenses/mit-license.php                                 //
//                                                                                 //
//=================================================================================//

#include "GACT.h"
#include "mtrand.h"
#include <iostream>

//class projection;

using namespace std;
using namespace ilab;


//=============================================================================
//class Individual
//=============================================================================
//constructor
Individual::Individual()
{

}
//initialize Gene randomly
void Individual::initialize(size_t x, size_t y) {
    random_device rd;
    mt19937_64 engine(rd());
    uniform_real_distribution<float> dist(0,random_maximum);
    gene.resize(x, y,1);
    for(auto&& g : gene.quantities()){
        g = dist(engine);
    }
    const auto intersect = [](float _x, float _y, float _r)
    {
        return _x * _x + _y * _y <= _r * _r;
    };
    const float center_x = static_cast<float>(gene.width()) / 2.0f;
    const float center_y = static_cast<float>(gene.height()) / 2.0f;
    for (size_t y = 0; y < gene.height(); y++) {
        for (size_t x = 0; x < gene.width(); x++){
            if (intersect(static_cast<float>(x) - center_x, static_cast<float>(y) - center_y, center_y)){
                gene.identity(x, y, 0) = ilab::blank_type::quantity;
            }else{
                gene.identity(x, y, 0) = ilab::blank_type::outside;
            }
        }
    }

    fitness = -FLT_MAX;
}
//既知情報としてなんらかのdistributionが与えられた際の初期個体生成
void Individual::initialize(size_t x, size_t y,distribution initial_dist) {
    gene.resize(x, y,1);
    for(auto&& g : gene.quantities()){
        g = 1.0;
    }
    const auto intersect = [](float _x, float _y, float _r)
    {
        return _x * _x + _y * _y <= _r * _r;
    };
    const float center_x = static_cast<float>(gene.width()) / 2.0f;
    const float center_y = static_cast<float>(gene.height()) / 2.0f;
    for (size_t j = 0; j < gene.height(); j++)
    {
        for (size_t i = 0; i < gene.width(); i++)
        {
            if (intersect(static_cast<float>(i) - center_x, static_cast<float>(j) - center_y, center_y))
            {
                gene.identity(i, j, 0) = ilab::blank_type::quantity;
                gene.quantity(i, j, 0) = abs(gaussian_noise(initial_dist.quantity(i, j, 0),0.1));
            }
            else
            {
                gene.identity(i, j, 0) = ilab::blank_type::outside;
            }
        }
    }

    fitness = -FLT_MAX;
}

bool Individual::operator<(const Individual& individual) const
{
    return fitness > individual.fitness;
}

float Individual::gaussian_noise(float mu, float sigma) {
    //対数正規分布だから名前変えろ
    random_device rd;
    mt19937_64 engine(rd());
    uniform_real_distribution<float> dist(0,1);
    float z = (float) (mu + sigma * sqrt(-2.0 * log(dist(engine))) * sin(2.0 * M_PI * dist(engine)));
    return exp(z);
}


//=============================================================================
//class GACT
//=============================================================================
GACT::GACT(ilab::projection& projections, string dist_data_path){

    p_data = projections;
    m_DimensionX = static_cast<int>(projections.height());
    m_DimensionY = static_cast<int>(projections.height());
    m_DimensionZ = static_cast<int>(projections.width());
    m_DimensionI = static_cast<int>(projections.counts());

    populationSize = 250;
    maxGeneration = 1500;
    crossover_pb = 1.0;
    mutation_pb = 0.5;
    tournamentSize = 4;
    bestFittness = 0;
    initial_image_path = dist_data_path;

    const auto weight_function = [](float _length)
    {
        //return 1.0 - _length;
        return 1.0f - 3.0f * _length * _length + 2.0f * _length * _length * _length;
    };
    projector.set_weight(weight_function);

}


distribution GACT::Evolution(){
    cout<<"Start evaluation"<<endl;
    distribution initial_dist(initial_image_path);
    if(initial_dist.empty()){
        cout<<"normal population method"<<endl;
        init_population();
    }else{
        cout<<"population method with Initial Distribution"<<endl;
        init_population_with_initial_dist(initial_dist);
    }
    projected_points =  projector.calculate_projected_points(population[0].gene,p_data.angles());

    random_maximum=3.0;

    cout<<"=====Evaluated "<< population.size() << " individuals====="<<endl;
    generation=1;
    ofstream logging;
    logging.open("evolution.csv",std::ios::out);
    logging<<"generation,fittness"<<endl;
    for(generation=1; generation <= maxGeneration; generation++){
        cout<<"===Generation:" << generation << "==="<<endl;
        //selection
        selection();
        //cout<<"end selection"<<endl;

        if(generation<1){
            int worth=0;
            double f = DBL_MAX;
            for(int i=0;i<population.size();i++){
                if (population[i].fitness<f){
                    f= population[i].fitness;
                    worth = i;
                }
            }
            population[worth] = elite;
        }
        crrossover();
        cout<<"end crossover"<<endl;
        mutate();
        cout<<"end mutation"<<endl;
        fittness();
        cout<<"end fittness"<<endl;
        best_individual();
        elite = bestIndividual;
        cout<<"Max : "<<fixed<<bestFittness<<endl<<endl;

        logging<<generation<<","<<bestFittness<<endl;
        if(generation%10 == 0) {
            save_individual(bestIndividual, generation);
        }
        cout<<endl;
    }
    division_each_angles(bestIndividual);
    return bestIndividual.gene;

}
//母集団から最も適応度の高い個体を探す
void GACT::best_individual(){
    double f = -DBL_MAX;
    for(auto&& ind:population){
        if (ind.fitness>f){
            f = ind.fitness;
            bestIndividual = ind;
            bestFittness = ind.fitness;
        }
    }
}
//個体の保存．命名のために世代数を渡す
void GACT::save_individual(Individual ind,int generation){
    for(int j = 0; j < ind.gene.width(); j++) {
        for (int k = 0; k < ind.gene.height(); k++){
            if(ind.gene.identity(static_cast<size_t>(j), static_cast<size_t>(k), 0)==ilab::blank_type::quantity) {
                ind.gene.quantity(static_cast<size_t>(j), static_cast<size_t>(k), 0) += 1;
            }else if(ind.gene.identity(static_cast<size_t>(j), static_cast<size_t>(k), 0)==ilab::blank_type::outside){
                ind.gene.quantity(static_cast<size_t>(j), static_cast<size_t>(k), 0) = 1;
            }

        }
    }
    distribution save_density;
    save_density.resize(m_DimensionX,m_DimensionY,m_DimensionZ);
    for (size_t y = 0; y < ind.gene.height(); y++) {
        for (size_t x = 0; x < ind.gene.width(); x++) {
            save_density.quantity(x,y,0) = ind.gene.quantity(x, y, 0);
        }
    }
    save_density.save("result/3d_density_gen" + to_string(generation));
}
//各角度における距離を算出して保存
void GACT::division_each_angles(Individual ind){
    double f=0.0;
    distribution b = ind.gene;
    ilab::projection reproject = projector.project(b,p_data.angles(),projected_points);

    vector<vector<double>> data(reproject.counts(),vector<double>(2));
    for(unsigned  int i=0;i<p_data.angles().size();i++) {
        data[i][0]=p_data.angle(i);
    }

    for(unsigned  int i=0;i<reproject.counts();i++){
        for(unsigned int x=0;x<reproject.height();x++) {
            if(reproject.identity(0,x,i)==ilab::blank_type::quantity){
                f+=abs(p_data.quantity(0, x, i) - reproject.quantity(0, x, i));
            }
        }
        data[i][1]=f;
        f=0.0;
    }
    FILE *fp = fopen( "div_eachangle.csv", "w" );
    if( fp != NULL ) {
        for(int i = 0; i < data.size(); i++ ) {
            fprintf( fp, "%f,%f\n", data[i][0], data[i][1]);
        }
        fclose( fp );
    }
}
//母集団を用意する
void GACT::init_population() {
    population.resize(populationSize);
    newGeneration.resize(populationSize);

    for(auto &&i : population){
        i = Individual();
        i.initialize(m_DimensionX,m_DimensionY);
        i.random_maximum=random_maximum;
    }
    save_individual(population[0],0);

}
//母集団を用意する．なんらかの既知情報が与えられてる時こちらを使う(overrideすればいいのに)
void GACT::init_population_with_initial_dist(distribution initial_dist) {
    population.resize(populationSize);
    newGeneration.resize(populationSize);

    if(initial_dist.width()!=m_DimensionX || initial_dist.height()!=m_DimensionY){
        cout<<"Size did not match"<<endl;
        init_population();
        return;
    }

    for(auto &&i : population){
        i = Individual();
        i.initialize(m_DimensionX,m_DimensionY,initial_dist);
        i.random_maximum=random_maximum;
    }
    save_individual(population[0],0);

}

//選択処理
void GACT::selection() {
    vector<Individual> tournament;
    tournament.resize(static_cast<size_t>(tournamentSize));

    for(unsigned i=0; i< newGeneration.size(); i++){
        for(unsigned j=0; j< tournament.size(); j++){
            tournament[j] = population[genrand_int(static_cast<int>(population.size()))];
        }
        sort(tournament.begin(),tournament.end());
        newGeneration[i] = tournament[0];

    }

    for(unsigned i = 0; i < population.size(); i++)
        population[i] = newGeneration[i];
}
//交叉処理
void GACT::crrossover() {

    int median;
    vector<Individual> child(2);
    median = static_cast<int>(population.size()) / 2;
    random_device rd; // obtain a random number from hardware
    mt19937 eng(rd()); // seed the generator
    uniform_int_distribution<> all_rand;

    for(int i = 0; i < median; i++) {
        if (genrand() < crossover_pb) {
            child[0] = population[i];
            child[1] = population[i + median];

            uniform_real_distribution<> rand_r(0,population[i].gene.width() / 2.0f);
            uniform_real_distribution<> rand_theta(0,2 * 3.141592f);
            double r = rand_r(eng);
            double theta = rand_theta(eng);
            int x1 = static_cast<int>(r*cos(theta)+r);
            int y1 = static_cast<int>(r*sin(theta)+r);

            r = rand_r(eng);
            theta = rand_theta(eng);
            int x2 = static_cast<int>(r*cos(theta)+r);
            int y2 = static_cast<int>(r*sin(theta)+r);

            int row_a,col_a,row_b,col_b;
            /* 正方形ないで一様ランダム
            uniform_int_distribution<> rand_row(0, static_cast<int>(child[0].gene.width())-1); // define the range 0_tablerows
            uniform_int_distribution<> rand_col(0, static_cast<int>(child[0].gene.height())-1); // define the range 0_tablerows
            row_a = rand_row(eng);
            col_a = rand_col(eng);
            row_b = rand_row(eng);
            col_b = rand_col(eng);
             */
            row_a = x1;
            col_a = y1;
            row_b = x2;
            col_b = y2;

            if (row_a > row_b){
                swap(row_a,row_b);
            }
            if (col_a > col_b){
                swap(col_a,col_b);
            }
            float v_0,v_1;
            for(int j=row_a;j<row_b;j++){
                for(int k=col_a;k<col_b;k++){
                    if(population[i].gene.identity(static_cast<size_t>(j), static_cast<size_t>(k), 0)==ilab::blank_type::quantity) {
                        v_0 = child[0].gene.quantity(static_cast<size_t>(j),static_cast<size_t>(k),0);
                        v_1 = child[1].gene.quantity(static_cast<size_t>(j),static_cast<size_t>(k),0);
                        //and swap
                        child[0].gene.quantity(static_cast<size_t>(j),static_cast<size_t>(k),0) = v_1;
                        child[1].gene.quantity(static_cast<size_t>(j),static_cast<size_t>(k),0) = v_0;
                    }
                }
            }
            //cout<<"swaped"<<endl;
            //cout<<child[0].gene<<endl<<" and "<<endl<<child[1].gene<<endl<<endl<<endl;

            population[i] = child[0];
            population[i + median] = child[1];

        }
    }
}
//突然変異処理
void GACT::mutate(){

    double probability;
    mt19937_64 engine(1);
    uniform_real_distribution<float> distribution(0,random_maximum);
    random_device rd; // obtain a random number from hardware
    mt19937 eng(rd()); // seed the generator

    //uniform_int_distribution<> rand_row(0, static_cast<int>(population[0].gene.width())-1); // define the range 0_tablerows
    //uniform_int_distribution<> rand_col(0, static_cast<int>(population[0].gene.height())-1); // define the range 0_tablerows
    uniform_real_distribution<> rand_r(0,population[0].gene.width() / 2.0f);
    uniform_real_distribution<> rand_theta(0,2 * 3.141592f);


    for(unsigned i = 0; i < population.size(); i++)
    {
        probability = genrand();
        if (probability < mutation_pb) {
            //int x = rand_row(eng);
            //int y = rand_col(eng);
            double r = rand_r(eng);
            double theta = rand_theta(eng);
            int x = static_cast<int>(r*cos(theta)+r);
            int y = static_cast<int>(r*sin(theta)+r);
            int mutation_range = 0;
            if(population[i].fitness>-800){
                mutation_range = 2;
            }else if(population[i].fitness>-1000){
                mutation_range = 6;
            }else{
                mutation_range = 10;
            }


            float rand_value = distribution(engine);

            int start_x = x-mutation_range;
            int start_y = y-mutation_range;
            if(start_x<0)
                start_x = 0;
            if(start_y<0)
                start_y = 0;

            int end_x = x+mutation_range;
            int end_y = y+mutation_range;
            if(end_x>population[i].gene.width()-1)
                end_x = static_cast<int>(population[i].gene.width())-1;
            if(end_y>population[i].gene.height()-1)
                end_y = static_cast<int>(population[i].gene.height())-1;

            for(int j = start_x; j < end_x; j++) {
                for (int k = start_y; k < end_y; k++){
                    if(population[i].gene.identity(static_cast<size_t>(j), static_cast<size_t>(k), 0)==ilab::blank_type::quantity) {
                        population[i].gene.quantity(static_cast<size_t>(j), static_cast<size_t>(k), 0) = rand_value;
                    }

                        //cout<<"mutate"<<endl;
                        //cout<<population[i].gene<<endl;
                }
            }
        }
    }


}
//母集団に含まれる個体に適応度を与える
void GACT::fittness(){
    for(auto&& ind:population){
        double f=0.0;
        distribution b = ind.gene;

        ilab::projection reproject = projector.project(b,p_data.angles(),projected_points);

        for(unsigned  int i=0;i<reproject.counts();i++){
            for(unsigned int x=0;x<reproject.height();x++) {
                if(reproject.identity(0,x,i)==ilab::blank_type::quantity){
                f += abs(p_data.quantity(0, x, i) - reproject.quantity(0, x, i));
                }
            }
        }

        ind.fitness = -1*(f/p_data.counts());
    }

}

