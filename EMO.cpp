//
// Created by Furukawa on 2017/05/19.
//

#include <iostream>
#include "EMO.h"

//=============================================================================
//class Individual
//=============================================================================
Individual::Individual(){
    gene = InverseDomain();
    fitness1 = FLT_MAX;
    fitness2 = FLT_MAX;
    rank = 0;
    local_distance = 0;

}

void Individual::initialize(size_t row, size_t colums){
    gene.resize(row,colums);
}

void Individual::initialize(distribution dist){
    gene = InverseDomain(dist);

}

void Individual::initialize(projection p_data){
    gene = InverseDomain(p_data);
    //std::mt19937 gen(1);
    //std::uniform_real_distribution<float> dist(-1000, 1000) ;

    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    //std::normal_distribution<float> dist(0.0, 10.0);結果をとったパラメータ
    std::normal_distribution<float> dist(0.0, 10.0);
    for(size_t i=0;i<gene.width();i++){
        for(size_t j=0;j<gene.height();j++){
            if(gene.infomation(i,j)==info_type::unknown){
                gene.quantity(i,j) = complex<float>(dist(engine),dist(engine));
            }
        }
    }
}
//=============================================================================
//class EMO
//=============================================================================
void make_p_data(projection &p ,size_t angle_count){
    distribution baseDist("experiment_data/no_object(196,196,1)-1.cfd");
    InverseDomain baseInv(baseDist);
    projection p_new(p.width(),baseDist.height(), angle_count);
    float dtheta = 180.0f / static_cast<float>(angle_count);
    vector<float> angles(angle_count);
    for (size_t i = 0; i < angles.size(); i++)
    {
        angles[i] = (float) ((dtheta * static_cast<float>(i)) * M_PI / 180.0);
    }
    p_new.angles() = angles;
    FFT<float> FFT;
    for (int i = 0; i < p_new.counts(); i++) {
        vector<complex<float>> freq(0);
        vector<float> time;
        for (int f=0;f<p_new.height();f++) {
            float r;
            float theta=p_new.angle(static_cast<size_t>(i));
            if (f>=p_new.height()/2) {
                r = p_new.height()/2-(f-p_new.height()/2);
                double x = r*sin(theta+M_PI)+((double)p_new.height()/2.0);
                double y = r*cos(theta+M_PI)+((double)p_new.height()/2.0);

                //cout<<"a "<<(int)x<<", "<<(int)y<<endl;
                complex<float> s = baseInv.quantity(static_cast<size_t>(x), static_cast<size_t>(y));
                freq.push_back(s);
            } else {
                r = f;
                double x = r*sin(theta)+((double)p_new.height()/2.0);
                double y = r*cos(theta)+((double)p_new.height()/2.0);
                //cout<<"b "<<static_cast<int>(x)<<", "<<static_cast<int>(y)<<endl;
                complex<float> s = baseInv.quantity(static_cast<size_t>(x), static_cast<size_t>(y));
                freq.push_back(s);
            }
        }

        FFT.inv(time,freq);
        for(int j=0;j<p_new.height();j++){
            p_new.quantity(0, static_cast<size_t>(j),static_cast<size_t>(i)) = time[j];
        }
    }
    p = p_new;
}

EMO::EMO(projection &projections) {
    p_data = projections;
    m_DimensionX = static_cast<int>(projections.height());
    m_DimensionY = static_cast<int>(projections.height());
    m_DimensionZ = static_cast<int>(projections.width());
    m_DimensionI = static_cast<int>(projections.counts());

    maxGeneration =100;
    generation = 0;
    searchSize = 100;
    archiveSize = 100;
    crossover_pb = 1.0;
    mutation_pb = 0.05;
    iterationOfGS = 25;

    bestFittness = 0;

    //InverseDomain visualize test-----------
    distribution dist("experiment_data/no_object(196,196,1)-1.cfd");
    InverseDomain invImg(dist);

    invImg.save_notshift("inverseimage.png");
    make_p_data(p_data,p_data.counts());
    cout<<""<<endl;
    InverseDomain inv_p_Img(p_data);
    inv_p_Img.save_notshift("inverse_p_data.png");

    FFT<float> FFT;
    vector<complex<float>> freq;
    for (int i = 0; i < p_data.counts(); i++) {
        vector<float> d(p_data.height());
        for (int j = 0; j < p_data.height(); ++j) {
            d[j] = p_data.quantity(0, static_cast<size_t >(j), static_cast<size_t>(i));
        }
        //dはpのi番目の投影データが格納された１次元配列
        FFT.fwd(freq, d);//fftする

        for (int f=0; f<freq.size();f++) {
            float r;
            float theta = p_data.angle(static_cast<size_t>(i));
            if (f>=freq.size()/2) {
                r = freq.size()/2-(f-freq.size()/2);
                double x = r*sin(theta+M_PI)+((double)p_data.height()/2.0);
                double y = r*cos(theta+M_PI)+((double)p_data.height()/2.0);
                invImg.quantity(static_cast<size_t>(x), static_cast<size_t>(y)) = freq[f];
                invImg.infomation(static_cast<size_t>(x), static_cast<size_t>(y)) = info_type::known;
            } else{
                r = f;
                double x = r*sin(theta)+((double)p_data.height()/2.0);
                double y = r*cos(theta)+((double)p_data.height()/2.0);
                invImg.quantity(static_cast<size_t>(x), static_cast<size_t>(y)) = freq[f];
                invImg.infomation(static_cast<size_t>(x), static_cast<size_t>(y)) = info_type::known;
            }
        }
    }
    invImg.save_notshift("inverseimage_with_p.png");
    if(invImg.TwoDimFFT().save("INVERSE")){
        cout<<"done"<<endl;
    }else{
        cout<<"faile"<<endl;
    }
    //InverseDomain visualize test-----------
    cout<<""<<endl;
}
void EMO::check_NaN(){
    for(int k=0;k<search_population.size();k++) {
        MatrixXcf child = search_population[k].gene.quantity2matrix();
        for (int x = 0; x < child.rows(); x++) {
            for (int y = 0; y < child.cols(); y++) {
                if (isnan(child(x, y).real()) || isnan(child(x, y).imag()) || isinf(child(x, y).real()) ||
                    isinf(child(x, y).imag())) {
                    cout << "base!! " << x << " " << y << " :" << child(x, y).real() << " "
                         << child(x, y).imag() << endl;
                    exit(0);
                }
            }
        }
    }
}

void EMO::evolution() {
    cout<<"Start evaluation"<<endl;
    init_population();
    vector<Individual> Rt;
    ofstream logging;
    logging.open("evolution.csv",std::ios::out);
    logging<<"generation,fittness"<<endl;
    for(generation=0;generation<maxGeneration;generation++){
        cout<<"===Generation:" << generation << "==="<<endl;

        fittness();
        for(auto ind:search_population){
            cout<<ind.fitness1+ind.fitness2<<endl;
        }
        cout<<"end fittness"<<endl;
        check_NaN();

        archive_population.clear();
        float distance_max = 0.0;
        vector<vector<Individual>> F;
        vector<Individual > P(search_population.size());
        P=search_population;
        for(int i=0;i<search_population.size();i++){
            for(int j=0;j<search_population.size();j++){
                if(i==j && distance_max<abs(search_population[i].fitness1-search_population[j].fitness1)+abs(search_population[i].fitness2-search_population[j].fitness2)){
                    distance_max = abs(search_population[i].fitness1-search_population[j].fitness1)+abs(search_population[i].fitness2-search_population[j].fitness2);
                }
            }
        }
        int r=0;
        vector<Individual> temp;
        while(P.size()!=0) {
            float best_distance = FLT_MAX;
            for (int i = 0; i < P.size(); i++) {
                if (best_distance > (P[i].fitness1 + P[i].fitness2)) {
                    best_distance = P[i].fitness1 + P[i].fitness2;
                }
            }
            vector<Individual>::iterator itr = P.begin();
            while (itr != P.end()) {
                if ((itr.base()->fitness1 + itr.base()->fitness2) == best_distance) {
                    itr.base()->rank = r;
                    temp.push_back(*itr.base());
                    itr = P.erase(itr);
                } else {
                    itr++;
                }
            }
            F.push_back(temp);
            temp.clear();
            r+=1;
        }

        archive_population.push_back(F[0][0]);
        F[0].erase(F[0].begin() + 0);
        vector<vector<float>> Fi;

        while(archive_population.size()<archiveSize/10) {
            for (int i = 0; i < F.size(); i++) {
                vector<float> I_num;
                for (int j = 0; j < F[i].size(); j++) {
                    int sum = 0;
                    for (int k = 0; k < archive_population.size(); k++) {
                        sum += pow(search_population.size() - archive_population.size(),
                                   2 * (distance_max - (abs(abs(archive_population[k].fitness1 - F[i][j].fitness1) -
                                                            abs(archive_population[k].fitness2 - F[i][j].fitness2))) /
                                                       distance_max));
                    }
                    I_num.push_back(F[i][j].rank + sum);
                }
                Fi.push_back(I_num);
            }

            float best_fi = FLT_MAX;
            int a = 0, b = 0;
            for (int i = 0; i < F.size(); i++) {
                for (int j = 0; j < F[i].size(); j++) {
                    if (Fi[i][j] < best_fi) {
                        best_fi = Fi[i][j];
                        a = i;
                        b = j;
                    }
                }
            }
            archive_population.push_back(F[a][b]);
            F[a].erase(F[a].begin() + b);
            Fi[a].erase(Fi[a].begin() + b);
        }

        Rt.insert(Rt.end(),search_population.begin(),search_population.end());
        for (int i = 0; i < F.size(); i++) {
            for (int j = 0; j < F[i].size(); j++) {
                Rt.push_back(F[i][j]);
            }
        }

        F.clear();
        Fi.clear();
        P.clear();
        r=0;
        while(Rt.size()!=0){
            float best_distance = FLT_MAX;
            for(int i=0;i<Rt.size();i++){
                if(best_distance>(Rt[i].fitness1+Rt[i].fitness2)){
                    best_distance = Rt[i].fitness1+Rt[i].fitness2;
                }
            }
            temp.clear();
            vector<Individual>::iterator itr = Rt.begin();
            int c = 0;
            while(itr != Rt.end()){
                if(c>Rt.size()-1){
                    break;
                }
                if((Rt[c].fitness1+Rt[c].fitness2) == best_distance){
                    Rt[c].rank = r;
                    temp.push_back(Rt[c]);
                    itr = Rt.erase(itr);
                }else{
                    itr++;
                }
                c+=1;
            }
            F.push_back(temp);
            r+=1;
        }

        for (int i = 0; i < F.size(); i++) {
            for (int j = 0; j < F[i].size(); j++) {

                archive_population.push_back(F[i][j]);
                if(archive_population.size()==archiveSize){
                    break;
                }
            }
        }
        selection();
        cout<<"end selection"<<endl;
        best_individual();
        logging<<generation<<","<<bestIndividual.fitness1+bestIndividual.fitness2<<endl;
        crrossover();
        cout<<"end crossover"<<endl;
        check_NaN();
        mutate();
        cout<<"end mutate"<<endl;

    }


}

void EMO::init_population() {
    search_population.resize(searchSize);
    for(auto &&i : search_population){
        i = Individual();
        i.initialize(p_data);
    }
}


void EMO::selection() {
    search_population.clear();

    vector<vector<Individual>> F;
    vector<Individual> Rt(archive_population.size());
    copy(archive_population.begin(),archive_population.end(),Rt.begin());
    int r=0;
    while(Rt.size()!=0){
        float best_distance = FLT_MAX;
        for(int i=0;i<Rt.size();i++){
            if(best_distance > (Rt[i].fitness1+Rt[i].fitness2)){
                best_distance = Rt[i].fitness1+Rt[i].fitness2;
            }
        }
        vector<Individual> temp;
        vector<Individual>::iterator itr = Rt.begin();
        while(itr != Rt.end()){
            if((itr.base()->fitness1+itr.base()->fitness2) == best_distance){
                itr.base()->rank = r;
                temp.push_back(*itr.base());
                itr = Rt.erase(itr);
            }else{
                itr++;
            }
        }
        F.push_back(temp);
        r+=1;
    }
    Rt.clear();
    for (int i = 0; i < F.size(); i++) {
        for (int j = 0; j < F[i].size(); j++) {
            if(j==0 || j==F[i].size()-1){
                F[i][j].local_distance = FLT_MAX;
            }else{
                F[i][j].local_distance = abs(F[i][j-1].fitness1-F[i][j+1].fitness1)+abs(F[i][j-1].fitness2-F[i][j+1].fitness2);
            }
            Rt.push_back(F[i][j]);
        }
    }

    vector<Individual> tournament;
    tournament.resize(static_cast<size_t>(2));
    for(unsigned i=0; i< searchSize; i++){
        for(unsigned j=0; j< tournament.size(); j++){
            tournament[j] = Rt[rand()%Rt.size()];
        }
        if(tournament[0].rank>tournament[1].rank){
            search_population.push_back(tournament[1]);
        }else if(tournament[0].rank<tournament[1].rank){
            search_population.push_back(tournament[0]);
        }else if(tournament[0].rank==tournament[1].rank){
            if(tournament[0].local_distance>tournament[1].local_distance){
                search_population.push_back(tournament[0]);
            }else{
                search_population.push_back(tournament[1]);
            }
        }
    }

}

void EMO::crrossover() {
    std::mt19937 engin;
    std::uniform_real_distribution<float> cross_dice(0.0f,1.0f) ;
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(int k=0;k<search_population.size();k++) {
        if (cross_dice(engin) < crossover_pb) {
            MatrixXcf child = search_population[k].gene.quantity2matrix();
            //cout<<"parent"<<child(child.rows()-1,child.cols()-1)<<endl;
            int cx = (int) child.cols() / 2;
            int cy = (int) child.rows() / 2;// 追記:rowとcol逆
            vector<vector<MatrixXcf>> parents_parts{{child.block(cx, cy, cx, cy)},
                                                    {child.block(0, 0, cx, cy)},
                                                    {child.block(0, cy, cx, cy)},
                                                    {child.block(cx, 0, cx, cy)}};

            for(int i=0;i<10;i++){
                std::uniform_int_distribution<int> parent_dice(0, static_cast<int>(search_population.size())-1);
                MatrixXcf parents = search_population[parent_dice(engin)].gene.quantity2matrix();
                parents_parts[0].push_back(parents.block(cx, cy, cx, cy));
                parents_parts[1].push_back(parents.block(0, 0, cx, cy));
                parents_parts[2].push_back(parents.block(0, cy, cx, cy));
                parents_parts[3].push_back(parents.block(cx, 0, cx, cy));
            }
            std::mt19937 gen;
            std::uniform_real_distribution<float> dist(0.5, 1.0);
            float weight=0;
            for(int i=0;i<parents_parts.size();i++){
                float wi = dist(gen);
                weight+=wi;
                MatrixXcf child_part = parents_parts[i][0];
                for(int x=0;x<child_part.rows();x++){
                    for(int y=0;y<child_part.cols();y++){
                        child_part(x,y) = complex<float>(child_part(x,y).real()*wi,child_part(x,y).imag()*wi);

                    }
                }
                float wj;
                MatrixXcf Gj(child_part.rows(), child_part.cols());
                MatrixXcf sum(child_part.rows(), child_part.cols());
                for (int x = 0; x < sum.rows(); x++) {
                    for (int J = 0; J < sum.cols(); J++) {
                        sum(x,J) = complex<float>(0,0);
                    }
                }
                for(int j=1;j<parents_parts[i].size();j++){
                    std::uniform_real_distribution<float> rand_wj(0.0f, 1.0f - weight);
                    Gj = parents_parts[i][j];
                    wj=rand_wj(gen);
                    weight+=wj;
                    for(int x=0;x<sum.rows();x++){
                        for(int y=0;y<sum.cols();y++){
                            float real;
                            float imag;
                            complex<float> s(sum(x,y).real(),sum(x,y).imag());
                            real = sum(x, y).real() + Gj(x, y).real() * wj;
                            imag = sum(x, y).imag() + Gj(x, y).imag() * wj;
                            sum(x, y) = std::complex<float>(real,imag);
                            if (isnan(sum(x, y).real()) || isnan(sum(x, y).imag()) || isinf(sum(x, y).real()) || isinf(sum(x, y).imag())) {
                                cout<<"NaN from sum!"<<x<<","<<y<<" "<<i<<","<<j<<":"<<s<<" + "<<Gj(x,y)<<" = "<<sum(x,y)<<endl;
                                exit(0);

                            }
                        }
                    }
                }
                for(int x=0;x<child_part.rows();x++){
                    for(int y=0;y<child_part.cols();y++){
                        //cout<<"base"<<I<<","<<J<<":"<<parents_parts[i][0](I,J)<<"then"<<child_part(I,J)<<" + "<<sum(I,J)<<" = "<<complex<float>(child_part(I,J).real()+sum(I,J).real(),child_part(I,J).imag()+sum(I,J).imag())<<endl;
                        child_part(x,y) = complex<float>(child_part(x,y).real()+sum(x,y).real(),child_part(x,y).imag()+sum(x,y).imag());
                    }
                }
                parents_parts[i][0] = child_part;
                weight=0;
            }

            child.block(cx, cy, cx, cy) = parents_parts[0][0];
            child.block(0, 0, cx, cy) = parents_parts[1][0];
            child.block(0, cy, cx, cy) = parents_parts[2][0];
            child.block(cx, 0, cx, cy) = parents_parts[3][0];
            //cout<<"child"<<child(child.rows()-1,child.cols()-1)<<endl;
            search_population[k].gene.matrix2quantity(child);
        }
    }
}

void EMO::mutate() {
    std::mt19937 engin;
    std::uniform_real_distribution<float> mutate_dice(0.0f,1.0f) ;
    for(auto&& ind:search_population){
        if(mutate_dice(engin)<mutation_pb) {
            MatrixXcf mat = ind.gene.quantity2matrix();
            int cx = (int) mat.cols() / 2;
            int cy = (int) mat.rows() / 2;// 追記:rowとcol逆
            vector<MatrixXcf> G;
            std::mt19937 gen;
            std::uniform_real_distribution<float> dist(0.5f, 1.0f - FLT_MIN);
            G.push_back(mat.block(cx, cy, cx, cy));
            G.push_back(mat.block(0, 0, cx, cy));
            G.push_back(mat.block(0, cy, cx, cy));
            G.push_back(mat.block(cx, 0, cx, cy));
            MatrixXcf Gj;
            float weight = 0;
            for (int i = 0; i < G.size(); i++) {
                float wi = dist(gen);
                weight+=wi;
                for(int x=0;x<G[i].rows();x++){
                    for(int y=0;y<G[i].cols();y++){
                        G[i](x,y) = complex<float>(G[i](x,y).real()*wi,G[i](x,y).imag()*wi);

                    }
                }
                float wj = 0.0f;
                MatrixXcf sum(G[i].rows(), G[i].cols());
                for (int x = 0; x < sum.rows(); x++) {
                    for (int y = 0; y < sum.cols(); y++) {
                        sum(x,y) = complex<float>(0,0);
                    }
                }
                for (int j = 0; j < G.size(); j++) {
                    if (j != i) {
                        std::uniform_real_distribution<float> rand_wj(0.0, 1.0f - weight);
                        Gj = G[j];
                        wj = rand_wj(gen);
                        weight+=wj;
                        for(int x=0;x<sum.rows();x++){
                            for(int y=0;y<sum.cols();y++){
                                float real;
                                float imag;
                                complex<float> s(sum(x,y).real(),sum(x,y).imag());
                                real = sum(x,y).real() + Gj(x,y).real() * wj;
                                imag = sum(x,y).imag() + Gj(x,y).imag() * wj;
                                sum(x,y) = std::complex<float>(real,imag);
                                if (isnan(sum(x,y).real()) || isnan(sum(x,y).imag()) || isinf(sum(x,y).real()) || isinf(sum(x,y).imag())) {
                                    cout<<"NaN from sum!"<<x<<","<<y<<" "<<i<<","<<j<<":"<<s<<" + "<<Gj(x,y)<<" = "<<sum(x,y)<<endl;
                                    exit(0);
                                }
                            }
                        }
                    }
                }
                for(int x=0;x<G[i].rows();x++){
                    for(int y=0;y<G[i].cols();y++){
                        G[i](x,y) = complex<float>(G[i](x,y).real()+sum(x,y).real(),G[i](x,y).imag()+sum(x,y).imag());
                    }
                }
                weight=0;
            }
            mat.block(cx, cy, cx, cy) = G[0];
            mat.block(0, 0, cx, cy) = G[1];
            mat.block(0, cy, cx, cy) = G[2];
            mat.block(cx, 0, cx, cy) = G[3];
            ind.gene.matrix2quantity(mat);
        }
    }
}

float e_out(InverseDomain g){
    const auto intersect = [](float _x, float _y, float _r)
    {
        return _x * _x + _y * _y <= _r * _r;
    };
    const float center_x = static_cast<float>(g.width()) / 2.0f;
    const float center_y = static_cast<float>(g.height()) / 2.0f;
    for (size_t y = 0; y < g.height(); y++) {
        for (size_t x = 0; x < g.width(); x++){
            if (intersect(static_cast<float>(x) - center_x, static_cast<float>(y) - center_y, center_y)){
                g.identity(x, y) = ilab::blank_type::quantity;
            }else{
                g.identity(x, y) = ilab::blank_type::outside;
            }
        }
    }
    float sum = 0.0f;
    for(size_t i=0;i<g.width();i++){
        for(size_t j=0;j<g.height();j++){
            if(g.identity(i,j)==ilab::blank_type::outside){
                sum+=abs(g.quantity(i,j));
            }
        }
    }

    return sum;
}

float e_image(InverseDomain g){
    float sum = 0.0f;
    for(size_t i=0;i<g.width();i++){
        for(size_t j=0;j<g.height();j++){
            if(g.identity(i,j)==ilab::blank_type::outside){
                sum+=abs(g.quantity(i,j).imag());
            }
        }
    }

    return sum;
}

float e_pos(InverseDomain g){
    float sum = 0.0f;
    for(size_t i=0;i<g.width();i++) {
        for (size_t j = 0; j < g.height(); j++) {
            if (g.identity(i, j) == ilab::blank_type::quantity && g.quantity(i, j).real() < 0) {
                sum += abs(g.quantity(i, j).real());
            }
        }
    }

    return sum;
}

float e_diff(InverseDomain Gafter,InverseDomain Gbefore){
    float sum = 0.0f;
    for(size_t i=0;i<Gafter.quantities().size();i++){
        sum+=abs(Gbefore.quantities()[i]-Gafter.quantities()[i]);
    }
    sum = sum/(Gbefore.width()*Gbefore.width());
    return sum;
}


void EMO::fittness() {
#ifdef _OPENMP
#pragma omp parallel for
#endif

    for(int k=0;k<search_population.size();k++){
        InverseDomain Gbefore = search_population[k].gene;
        //フーリエ逆変換してgを作る
        InverseDomain g = search_population[k].gene;
        MatrixXcf mat = g.quantity2matrix();
        g.shift(mat);
        g.matrix2quantity(mat);
        InverseDomain temp(g.width(),g.height());
        FFT<float> FFT;
        vector<complex<float>> time;
        for(int i=0;i<g.height();i++){
            vector<complex<float>> d(g.width());
            for(int j=0;j<g.width();j++){
                d[j]=g.quantity(static_cast<size_t >(j), static_cast<size_t>(i));
            }
            FFT.inv(time,d);
            for(int j=0;j<g.width();j++){
                temp.quantity(static_cast<size_t>(j), static_cast<size_t>(i))=time[j];
            }
        }
        for(int i=0;i<g.width();i++){
            vector<complex<float>> d(g.height());
            for(int j=0;j<g.height();j++){
                d[j]=temp.quantity(static_cast<size_t >(i), static_cast<size_t>(j));
            }
            FFT.inv(time,d);
            for(int j=0;j<g.height();j++){
                g.quantity(static_cast<size_t>(i), static_cast<size_t>(j))=time[j];
            }
        }


        search_population[k].fitness1 = e_out(g)+e_image(g)+e_pos(g);
        search_population[k]=gs_algorithm(search_population[k]);
        search_population[k].fitness2 = e_diff(search_population[k].gene,Gbefore);
    }
    /*
    distribution baseDist("experiment_data/no_object(196,196,1)_normalize.cfd");
    InverseDomain baseInv(baseDist);
    Individual baseInd;
    baseInd.gene = baseInv;

    InverseDomain Gbefore = baseInd.gene;

    //フーリエ逆変換してgを作る
    InverseDomain g = baseInd.gene;
    MatrixXcf mat = g.quantity2matrix();
    g.shift(mat);
    g.matrix2quantity(mat);
    InverseDomain temp(g.width(),g.height());
    FFT<float> FFT;
    vector<complex<float>> time;
    for(int i=0;i<g.height();i++){
        vector<complex<float>> d(g.width());
        for(int j=0;j<g.width();j++){
            d[j]=g.quantity(static_cast<size_t >(i), static_cast<size_t>(j));
        }
        FFT.inv(time,d);
        for(int j=0;j<g.width();j++){
            temp.quantity(static_cast<size_t>(i), static_cast<size_t>(j))=time[j];
        }
    }
    for(int i=0;i<g.width();i++){
        vector<complex<float>> d(g.height());
        for(int j=0;j<g.height();j++){
            d[j]=temp.quantity(static_cast<size_t >(j), static_cast<size_t>(i));
        }
        FFT.inv(time,d);
        for(int j=0;j<g.height();j++){
            g.quantity(static_cast<size_t>(j), static_cast<size_t>(i))=time[j];
        }
    }
    baseInd.fitness1 = e_out(g)+e_image(g)+e_pos(g);

    baseInd=gs_algorithm(baseInd);
    baseInd.fitness2 = e_diff(baseInd.gene,Gbefore);
    cout<<""<<endl;
     */
}

void EMO::best_individual() {
    vector<vector<Individual>> F;
    vector<Individual > P(search_population.size());
    P=search_population;
    int r=0;
    vector<Individual> temp;
    while(P.size()!=0) {
        float best_distance = FLT_MAX;
        for (int i = 0; i < P.size(); i++) {
            if (best_distance > (P[i].fitness1 + P[i].fitness2)) {
                best_distance = P[i].fitness1 + P[i].fitness2;
            }
        }
        vector<Individual>::iterator itr = P.begin();
        while (itr != P.end()) {
            if ((itr.base()->fitness1 + itr.base()->fitness2) == best_distance) {
                itr.base()->rank = r;
                temp.push_back(*itr.base());
                itr = P.erase(itr);
            } else {
                itr++;
            }
        }
        F.push_back(temp);
        temp.clear();
        r+=1;
    }
    for(int i=0;i<F[0].size();i++){
        save_individual(F[0][i], generation,i);
    }
    bestIndividual = F[0][0];
}

void EMO::save_individual(Individual ind, int gen,int label) {
    if(ind.gene.save_notshift("result/"+to_string(gen)+"gen_inverse"+to_string(label)+".png")){
        cout<<"save"<<endl;
    }else{
        cout<<"faile"<<endl;
    }
    if(ind.gene.TwoDimFFT().save("result/"+to_string(gen)+"gen_density"+to_string(label))){
        cout<<"save fittness="<<(ind.fitness1+ind.fitness2)<<endl;
    }else{
        cout<<"faile"<<endl;
    }

}

Individual EMO::gs_algorithm(Individual in) {
    Individual input = in;
    distribution test = input.gene.TwoDimFFT();
    InverseDomain tester(test);
    int k = 0;
    while (true){
        distribution dist = input.gene.TwoDimFFT();
        dist.identities() = input.gene.identities();

        const auto intersect = [](float _x, float _y, float _r)
        {
            return _x * _x + _y * _y <= _r * _r;
        };
        const float center_x = static_cast<float>(dist.width()) / 2.0f;
        const float center_y = static_cast<float>(dist.height()) / 2.0f;
        for (size_t y = 0; y < dist.height(); y++) {
            for (size_t x = 0; x < dist.width(); x++){
                if (intersect(static_cast<float>(x) - center_x, static_cast<float>(y) - center_y, center_y)){
                    dist.identity(x, y,0) = ilab::blank_type::quantity;
                }else{
                    dist.identity(x, y,0) = ilab::blank_type::outside;
                }
            }
        }

        for (int i = 0; i < dist.width(); i++) {
            for (int j = 0; j < dist.height(); j++) {
                if (dist.identity(static_cast<size_t>(i), static_cast<size_t>(j), 0) == ilab::blank_type::outside) {
                    dist.quantity(static_cast<size_t>(i), static_cast<size_t>(j), 0) = 0;
                }
            }
        }
        //dist.save("testDist"+to_string(k));
        InverseDomain inv(dist);
        //inv.save_notshift("testDist_inv"+to_string(k)+".png");
        if(k==iterationOfGS){
            input.gene = inv;
            return input;
        }

        FFT<float> FFT;
        vector<complex<float>> freq;
        for (int i = 0; i < p_data.counts(); i++) {
            vector<float> d(p_data.height());
            for (int j = 0; j < p_data.height(); ++j) {
                d[j] = p_data.quantity(0, static_cast<size_t >(j), static_cast<size_t>(i));
            }
            //dはpのi番目の投影データが格納された１次元配列
            FFT.fwd(freq, d);//fftする

            for (int f = 0; f < freq.size(); f++) {
                float r;
                float theta = p_data.angle(static_cast<size_t>(i));
                if (f>=freq.size()/2) {
                    r = freq.size()/2-(f-freq.size()/2);
                    double x = r*sin(theta+M_PI)+((double)freq.size()/2.0);
                    double y = r*cos(theta+M_PI)+((double)freq.size()/2.0);
                    inv.quantity(static_cast<size_t>(x), static_cast<size_t>(y)) = freq[f];
                    inv.infomation(static_cast<size_t>(x), static_cast<size_t>(y)) = info_type::known;
                    if (isnan(freq[f].real()) || isnan(freq[f].imag())) {
                        cout << "nan at gs " << freq[f] << endl;
                    }
                } else {
                    r = f;
                    double x = r*sin(theta)+((double)freq.size()/2.0);
                    double y = r*cos(theta)+((double)freq.size()/2.0);
                    inv.quantity(static_cast<size_t>(x), static_cast<size_t>(y)) = freq[f];
                    inv.infomation(static_cast<size_t>(x), static_cast<size_t>(y)) = info_type::known;
                    if (isnan(freq[f].real()) || isnan(freq[f].imag())) {
                        cout << "nan at gs " << freq[f] << endl;
                    }
                }
            }
        }
        k++;
        input.gene = inv;
    }
}