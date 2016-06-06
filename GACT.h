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

#ifndef EMOCT_GACT_H
#define EMOCT_GACT_H
#include <array>
#include <queue>
#include <cfloat>
#include <vector>
#include <random>
#include "CUDAMath.h"
#include "ilab.h"

using namespace std;
using namespace ilab;


class Individual
{
public:
    //data
    distribution gene;
    mutable double	fitness;

    //operator
    bool operator<(const Individual& individual) const;

    //method (constructor)
    Individual();
    void initialize(size_t row, size_t colums);


};


class GACT {
public:
    GACT(ilab::projection &projections);
    distribution Evolution();

private:
    void init_population();
    void selection();
    void crrossover();
    void mutate();
    void fittness();
    void best_individual();

    unsigned int		m_DimensionX;
    unsigned int		m_DimensionY;
    unsigned int		m_DimensionZ;
    unsigned int		m_DimensionI;
    projection p_data;

    int maxGeneration;
    int generation;
    size_t populationSize;
    int tournamentSize;
    double crossover_pb;
    double mutation_pb;

    double bestFittness;
    Individual bestIndividual;

    vector<Individual> population;
    priority_queue<Individual> elites;
    Individual elite;
    vector<Individual> newGeneration;

};


#endif //EMOCT_GACT_H
