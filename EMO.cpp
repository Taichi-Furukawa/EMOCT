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

}

void Individual::initialize(size_t row, size_t colums){
    gene.resize(row,colums);
}

void Individual::initialize(distribution dist){
    gene = InverseDomain(dist);

}

void Individual::initialize(projection p_data){
    gene = InverseDomain(p_data);

    boost::random::random_device seed_gen;
    boost::random::mt19937 gen(seed_gen);
    boost::random::uniform_real_distribution<> dist;

    for(int i=0;i<gene.width();i++){
        for(int j=0;gene.height();j++){
            if(gene.infomation(static_cast<size_t>(i), static_cast<size_t>(j))==info_type::unknown){
                gene.quantity(static_cast<size_t>(i), static_cast<size_t>(j))=dist(gen);
            }
        }
    }


}
//=============================================================================
//class EMO
//=============================================================================
EMO::EMO(projection &projections) {
    p_data = projections;
    m_DimensionX = static_cast<int>(projections.height());
    m_DimensionY = static_cast<int>(projections.height());
    m_DimensionZ = static_cast<int>(projections.width());
    m_DimensionI = static_cast<int>(projections.counts());

    maxGeneration = 100;
    generation = 0;
    populationSize=100;
    archiveSize = 90;
    crossover_pb = 1.0;
    mutation_pb = 0.05;
    iterationOfGS = 25;

    bestFittness = 0;




    const auto weight_function = [](float _length)
    {
        return 1.0f - 3.0f * _length * _length + 2.0f * _length * _length * _length;
    };
    projector.set_weight(weight_function);

    distribution dist("art_dist/ART-DENSITY(197,197,70)_15.q");

    InverseDomain invImg(dist);
    InverseDomain inv_p_Img(p_data);
    invImg.save("inverseimage.png");
    inv_p_Img.save_notshift("inverse_p_data.png");
    if(inv_p_Img.TwoDimFFT().save("INVERSE")){
        cout<<"done"<<endl;
    }else{
        cout<<"faile"<<endl;
    }
}

void EMO::evolution() {

}

void EMO::init_population() {

}

void EMO::selection() {

}

void EMO::crrossover() {

}

void EMO::mutate() {

}
float e_diff(){

}

float e_out(distribution dist){
    float sum = 0.0;
    for(int i=0;i<dist.width();i++){
        for(int j=0;dist.height();j++){
            if(dist.identity(static_cast<size_t>(i), static_cast<size_t>(j),0)==ilab::blank_type::outside){
                sum += abs(dist.quantity(static_cast<size_t>(i), static_cast<size_t>(j),0));
            }
        }
    }
    return sum;
}

void EMO::fittness() {

}

void EMO::best_individual() {

}

void EMO::save_individual(Individual ind, int gen) {

}

InverseDomain EMO::gs_algorithm(Individual &in) {
    distribution dist = in.gene.TwoDimFFT();
    for(int i=0;i<dist.width();i++){
        for(int j=0;dist.height();j++){
            if(dist.identity(static_cast<size_t>(i), static_cast<size_t>(j),0)==ilab::blank_type::outside){
                dist.quantity(static_cast<size_t>(i), static_cast<size_t>(j),0)=0.0;
            }
        }
    }
    InverseDomain inv(dist);
    FFT<float> FFT;
    vector<complex<float>> freq;
    for(int i=0;i<p_data.counts();i++){
        vector<float> d(p_data.height());
        for (int j = 0; j < p_data.height(); ++j) {
            d[j]=p_data.quantity(0, static_cast<size_t >(j), static_cast<size_t>(i));
        }
        //dはpのi番目の投影データが格納された１次元配列
        FFT.fwd(freq,d);//fftする

        for(int f=0;f<freq.size();f++){
            float r;
            float theta = p_data.angle(static_cast<size_t>(i));
            if(f<freq.size()/2){

                r = freq.size()/2-f;
                int x = abs(static_cast<int>(r*cos(theta+M_PI)+freq.size()/2));
                int y = abs(static_cast<int>(r*sin(theta+M_PI)+freq.size()/2));
                inv.quantity(static_cast<size_t>(x), static_cast<size_t>(y))=freq[f];
                inv.infomation(static_cast<size_t>(x), static_cast<size_t>(y)) = info_type::known;
            }else{
                r = f-freq.size()/2;
                int x = static_cast<int>(r*cos(theta)+freq.size()/2);
                int y = static_cast<int>(r*sin(theta)+freq.size()/2);
                inv.quantity(static_cast<size_t>(x),static_cast<size_t>(y))=freq[f];
                inv.infomation(static_cast<size_t>(x),static_cast<size_t>(y)) = info_type::known;
            }
        }
    }
    return inv;
}
