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
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/core/eigen.hpp>
#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>
#include <omp.h>


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
    float local_distance;
    int ID;

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
    void best_individual(vector<vector<Individual>> F);
    void save_individual(Individual ind,int gen,string label);
    Individual gs_algorithm(Individual in);
    void check_NaN();
    void save_pareto_set(vector<vector<Individual>> F);
    void combine(vector<Individual>,int gen);

    unsigned int		m_DimensionX;
    unsigned int		m_DimensionY;
    unsigned int		m_DimensionZ;
    unsigned int		m_DimensionI;
    projection p_data;

    int maxGeneration;
    int generation;
    size_t searchSize;
    size_t archiveSize;
    float crossover_pb;
    float mutation_pb;
    float bestFittness;
    int iterationOfGS;
    Individual bestIndividual;
    vector<Individual> archive_population;
    vector<Individual> search_population;
    priority_queue<Individual> elites;
    Individual elite;
};


#endif //EMOCT_EMO_H
