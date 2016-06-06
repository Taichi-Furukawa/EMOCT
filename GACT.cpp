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

class projection;

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
    uniform_real_distribution<float> dist(0,100);
    gene.resize(x, y,1);
    for(auto&& g : gene.quantities()){
        g = dist(engine);
    }

    fitness = -DBL_MAX;
}

bool Individual::operator<(const Individual& individual) const
{
    return fitness > individual.fitness;
}

//=============================================================================
//class GACT
//=============================================================================
GACT::GACT(ilab::projection& projections){

    p_data = projections;
    m_DimensionX = static_cast<int>(projections.height());
    m_DimensionY = static_cast<int>(projections.height());
    m_DimensionZ = static_cast<int>(projections.width());
    m_DimensionI = static_cast<int>(projections.counts());

    populationSize = 250;
    maxGeneration = 1000;
    crossover_pb = 1.0;
    mutation_pb = 0.5;
    tournamentSize = 4;
    bestFittness = 0;
}

distribution GACT::Evolution(){
    cout<<"Start evaluation"<<endl;
    init_population();
    cout<<"=====Evaluated "<< population.size() << " individuals====="<<endl;
    generation=1;

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
        cout<<"Max : "<<bestFittness<<endl<<endl;
        bestIndividual.gene.save("3d_density_gen" + generation);
        //cout<<"best individual"<<endl;
        //cout<<bestIndividual.gene<<endl;
        cout<<endl;
    }

    return bestIndividual.gene;

}

//generate first population

void GACT::best_individual(){
    double f = -DBL_MAX;
    for(auto&& ind:population){
        if (ind.fitness>f){
            bestIndividual = ind;
            bestFittness = ind.fitness;
        }
    }
}

void GACT::init_population() {
    population.resize(populationSize);
    newGeneration.resize(populationSize);

    for(auto &&i : population){
        i = Individual();
        i.initialize(m_DimensionX,m_DimensionY);
    }
}

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

void GACT::crrossover() {

    int median;
    vector<Individual> child(2);
    median = static_cast<int>(population.size()) / 2;
    random_device rd; // obtain a random number from hardware
    mt19937 eng(rd()); // seed the generator

    for(int i = 0; i < median; i++) {
        if (genrand() < crossover_pb) {
            child[0] = population[i];
            child[1] = population[i + median];
            //cout<<child[0].gene<<endl<<" and "<<endl<<child[1].gene<<endl;
            //Create a child Individuals
            //Create four random number in a range of Individual size
            uniform_int_distribution<> rand_row(0, static_cast<int>(child[0].gene.width())-1); // define the range 0_tablerows
            uniform_int_distribution<> rand_col(0, static_cast<int>(child[0].gene.height())-1); // define the range 0_tablerows
            int row,col,row_s,col_s;
            row = rand_row(eng);
            col = rand_col(eng);
            uniform_int_distribution<> rand_rowsize(1, static_cast<int>(child[0].gene.width())-row); // define the range 0_tablerows
            uniform_int_distribution<> rand_colsize(1, static_cast<int>(child[0].gene.height())-col); // define the range 0_tablerows
            row_s = rand_rowsize(eng);
            col_s = rand_colsize(eng);

            //Create sub-matrix following numbers
            //table<float> sub_mat_0(static_cast<size_t>(row_s),static_cast<size_t>(col_s));
            //table<float> sub_mat_1(static_cast<size_t>(row_s),static_cast<size_t>(col_s));

            float v_0,v_1;
            for(int j=row;j<row+row_s;j++){
                for(int k=col;k<col+col_s;k++){
                    v_0 = child[0].gene.quantity(static_cast<size_t>(j),static_cast<size_t>(k),0);
                    v_1 = child[1].gene.quantity(static_cast<size_t>(j),static_cast<size_t>(k),0);
                    //sub_mat_0.at(static_cast<size_t>(0+j),static_cast<size_t>(0+k)) = child[0].gene.at(static_cast<size_t>(j),static_cast<size_t>(k));
                    //sub_mat_1.at(static_cast<size_t>(0+j),static_cast<size_t>(0+k)) = child[1].gene.at(static_cast<size_t>(j),static_cast<size_t>(k));
                    //and swap
                    child[0].gene.quantity(static_cast<size_t>(j),static_cast<size_t>(k),0) = v_1;
                    child[1].gene.quantity(static_cast<size_t>(j),static_cast<size_t>(k),0) = v_0;
                }
            }
            //cout<<"swaped"<<endl;
            //cout<<child[0].gene<<endl<<" and "<<endl<<child[1].gene<<endl<<endl<<endl;

            population[i] = child[0];
            population[i + median] = child[1];

        }
    }
}

void GACT::mutate(){

    double probability;
    mt19937_64 engine(1);
    uniform_real_distribution<float> distribution(0,100);

    for(unsigned i = 0; i < population.size(); i++)
    {
        for(unsigned j = 0; j < population[i].gene.width(); j++) {
            for (unsigned k = 0; k < population[i].gene.height(); k++){
                probability = genrand();
                if (probability < mutation_pb) {
                    population[i].gene.quantity(j, k, 0) = distribution(engine);
                    //cout<<"mutate"<<endl;
                    //cout<<population[i].gene<<endl;
                }
            }
        }
    }

}

void GACT::fittness(){
    for(auto&& ind:population){
        double f=0.0;
        distribution b = ind.gene;

        const auto weight_function = [](float _length)
        {
            //return 1.0 - _length;
            return 1.0f - 3.0f * _length * _length + 2.0f * _length * _length * _length;
        };
        projector projector(weight_function);
        ilab::projection reproject = projector.project(b,p_data.angles());

        for(unsigned  int i=0;i<reproject.counts();i++){
            for(unsigned int x=0;x<reproject.height();x++) {

                f += abs(p_data.quantity(0, x, i) - reproject.quantity(0, x, i));
            }
        }

        ind.fitness = -1*(f/p_data.counts());
    }

}

