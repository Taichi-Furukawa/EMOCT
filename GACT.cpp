//
// Created by Furukawa on 2016/05/16.
//

#include "GACT.h"
#include "mtrand.h"
#include <iostream>

using namespace std;
using namespace arch;


//=============================================================================
//class Individual
//=============================================================================
//constructor
Individual::Individual()
{

}
//initialize Gene randomly
void Individual::initialize(size_t row, size_t colums) {
    mt19937_64 engine(1);
    uniform_real_distribution<float> distribution(numeric_limits<float>::min(),numeric_limits<float>::max());
    gene.resize(row, colums);
    for(auto& g : gene){
        g = distribution(engine);
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
GACT::GACT(size_t popSize,int ngen,double cxpb, double mutpb){
    populationSize = popSize;
    maxGeneration = ngen;
    crossover_pb = cxpb;
    mutation_pb = mutpb;
    tournamentSize = 4;
    bestFittness = 0;
    
}

table<float> GACT::Evolution(){
    printf("Start evaluation");
    init_population();
    cout<<"=====Evaluated "<< population.size() << " individuals====="<<endl;

    for(generation; generation < maxGeneration; generation++){
        cout<<"===Generation:" << generation << "==="<<endl;
        //selection
        selection();
        crrossover();
        mutate();
        fittness();

    }

    return table<float>();
}

//generate first population
void GACT::init_population() {
    population.resize(populationSize);
    newGeneration.resize(populationSize);
    for(auto &i : population){
        i = Individual();
        i.initialize(100,100);
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
}

void GACT::crrossover() {

}

void GACT::mutate(){

}

void GACT::fittness(){

}

