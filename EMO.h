//
// Created by Furukawa on 2017/05/19.
//

#ifndef EMOCT_EMO_H
#define EMOCT_EMO_H

#include <array>
#include <queue>
#include <cfloat>
#include <vector>
#include <random>
#include <boost/random.hpp>
#include <boost/random/random_device.hpp>
#include "CUDAMath.h"
#include "ilab.h"
#include "InverseDomain.h"


using namespace std;
using namespace ilab;


class Individual
{
public:
    //data
    InverseDomain gene;
    float fitness1;
    float fitness2;
    int rank;

    //property
    float random_maximum;

    //method (constructor)
    Individual();
    void initialize(size_t row, size_t colums);
    void initialize(distribution dist);
    void initialize(projection p_data);


};

class EMO {
public:
    EMO(projection &projections);
    void evolution();

private:
    void init_population();
    void selection();
    void crrossover();
    void mutate();
    void fittness();
    void best_individual();
    void save_individual(Individual ind,int gen);
    InverseDomain gs_algorithm(Individual &in);

    unsigned int		m_DimensionX;
    unsigned int		m_DimensionY;
    unsigned int		m_DimensionZ;
    unsigned int		m_DimensionI;
    projection p_data;
    ilab::projector projector;
    vector<distribution> projected_points;

    int maxGeneration;
    int generation;
    size_t populationSize;
    size_t archiveSize;
    double crossover_pb;
    double mutation_pb;
    double bestFittness;
    int iterationOfGS;
    Individual bestIndividual;
    vector<Individual> archive_population;
    vector<Individual> new_population;
    priority_queue<Individual> elites;
    Individual elite;

};


#endif //EMOCT_EMO_H
