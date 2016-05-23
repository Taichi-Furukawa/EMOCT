//
// Created by Furukawa on 2016/05/16.
//

#ifndef EMOCT_GACT_H
#define EMOCT_GACT_H
#include <array>
#include <cfloat>
#include <vector>
#include <random>
#include "libs/table.h"

using namespace std;
using namespace arch;


class Individual
{
public:
    //data
    table<float> gene;
    mutable double	fitness;

    //operator
    bool operator<(const Individual& individual) const;

    //method (constructor)
    Individual();
    void initialize(size_t row, size_t colums);


};


class GACT {
public:
    GACT(size_t popSize,int ngen,double cxpb, double mutpb);
    table<float> Evolution();
    void init_population();
    void selection();
    void crrossover();
    void mutate();
    void fittness();


private:
    int maxGeneration;
    int generation;
    size_t populationSize;
    int tournamentSize;
    double crossover_pb;
    double mutation_pb;

    double bestFittness;
    Individual bestIndividual;

    vector<Individual> population;
    vector<Individual> elites;
    vector<Individual> newGeneration;

};


#endif //EMOCT_GACT_H
