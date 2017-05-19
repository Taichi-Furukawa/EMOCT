//-----------------------------------------------------------------------------
// GA2K
//
// ISDL, Doshisha University Distributed Genetic Algorithm Software
// Version 1.4.2
//-----------------------------------------------------------------------------
// Parallel Genetic Algorithm Team of Intelligent Systems Design Laborato
// ry, Doshisha University.
//
// Copyright (C) 1999           Ikki Ohmukai,      All rights reserved.
// Copyright (C) 1999,2000,2001 Jun-ichi Yoshida,  All rights reserved.
// Copyright (C) 2001,2002,2003,2004,2005,2006  Tomoyuki Hiroyasu, All rights reserved.
// Copyright (C) 2001,2002,2003,2004,2005,2006  Mitsunori Miki,    All rights reserved.
// Copyright (C) 2002           Masaki Sano,       All rights reserved.
// Copyright (C) 2002,2003      Junji Sawada,      All rights reserved.
// Copyright (C) 2005,2006      Kengo Yoshii,      All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
// 3. The names of the authors should not be used to endorse or promote
//    products derived from this software without specific prior
//    written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE AUTHORS ``AS IS'' AND ANY EXPRESS OR
// IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT,
// INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// ga2k.cpp
//-----------------------------------------------------------------------------

#include <cstdio>
#include <cstdlib>
#include <cfloat>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <assert.h>
using namespace std;


#include "ga2k.h"
#include "randring.h"
#include "functions.h"
#include "interface.h"

///////////////////////////////////////////////////////////////////////////////
// global variable
///////////////////////////////////////////////////////////////////////////////

const char* function_name[] = {
	"Rastrigin",
	"Schwefel",
	"Griewank",
	"Rosenbrock",
	"Ridge",
	"Onemax",
};
const int function_num = sizeof(function_name)/sizeof(const char*);

const char* selection_name[] = {
	"Tournament",
	"Roulette"
};
const int selection_num = sizeof(selection_name)/sizeof(const char*);

ObjectiveProblem problem;



///////////////////////////////////////////////////////////////////////////////
// function
///////////////////////////////////////////////////////////////////////////////

//gray coding
static void gray_code(vector<char>& bitString, int startPoint, int bitLength, int num)
{
	int mask;
	bool f = false;
	
	for(mask = 1 << (sizeof(int) * 8 - 2); mask != 0; mask >>= 1)
	{
		if(f)
		{
			num ^= mask;
			if(num & mask) f = false;
		}
		else
			if(num & mask) f = true;
	}
	
	int i = startPoint + bitLength;
	while(startPoint < i){
		bitString[--i] = num & 1;
		num >>= 1;
	}
}


//gray decode
static int gray_decode(const vector<char>& bitString, int startPoint, int bitLength)
{
  int i, real, mask;
  //unsigned int mask;
	bool sw = false;

	for(i = real = 0; i < bitLength; i++)
	{
		if(bitString[startPoint + i])
			real |= 1 << (bitLength - i - 1);
	}
	
	for(mask = 1 << (bitLength - 1); mask != 0; mask >>= 1)
	{
	  //cerr << "mask : " << mask << endl;
	  if(mask==-1) exit(1);
	  if(sw)
		{
			if(real & mask) sw = false;
			real ^= mask;
		}
		else if(real & mask) sw = true;
	}
	return real;
}

//binary decode
static int binary_decode(const vector<char> bitString, int startPoint, int bitLength)
{
	int i, real;
	for(i = real = 0; i < bitLength; i++)
	{
		if(bitString[startPoint + i])
			real += 1 << (bitLength - i - 1);
	}

	return real;
}

//user definition problem
static double user_function(const vector<double>& real)
{
	vector<double> result;
	double fitness;

	result = problem.evaluate(real);
	fitness = problem.result_to_fitness(result);
	return fitness;
}


//=============================================================================
//class Gene
//=============================================================================

Gene::Gene()
{

}

Gene::~Gene()
{

}

void Gene::initialize(int geneLength)
{
	int i;

	bitString.resize(geneLength);

	for(i = 0; i < geneLength; i++)
	{
		bitString[i] = genrand_int(2);
	}
	fitness = -DBL_MAX;
}

void Gene::initialize(int geneLength, char value)
{
	int i;

	bitString.resize(geneLength);

	for(i = 0; i < geneLength; i++)
	{
		bitString[i] = value;
	}
	fitness = -DBL_MAX;
}

void Gene::initialize(int geneLength, const vector<char>& code)
{
	bitString.resize(geneLength);

	for(int i = 0; i < geneLength; i++)
	{
		bitString[i] = code[i];
	}
	fitness = -DBL_MAX;
}


void Gene::visualize()
{
	for(unsigned i = 0; i < bitString.size(); i++)
		cerr << static_cast<int>(bitString[i]);
	cerr << endl;
}

bool Gene::operator<(const Gene& gene) const
{
	return fitness > gene.fitness;
}

//=============================================================================
//class Island
//=============================================================================

Island::Island()
{

}

Island::~Island()
{

}

void Island::initialize(ifstream& fin)
{
	gene.resize(populationSize);
	migrant.resize(migrantSize);
	elite.resize(numElite);
	newGeneration.resize(populationSize);
	crossoverTable.resize(geneLength - 1);

	for(unsigned i = 0; i < gene.size(); i++)
		gene[i].initialize(geneLength);
	
	if(fin.is_open())
		setGeneFromFile(fin);
	for(unsigned i = 0; i < gene.size(); i++)
		fitness(gene[i]);

	for(unsigned i = 0; i < migrant.size(); i++)
		migrant[i].initialize(geneLength);

	for(unsigned i = 0; i < elite.size(); i++)
		elite[i].initialize(geneLength);

	sort(gene.begin(), gene.end());
	copy(gene.begin(), gene.begin() + elite.size(), elite.begin());
	
	for(unsigned i = 0; i < newGeneration.size(); i++)
		newGeneration[i].initialize(geneLength);

	for(int i = 0; i < geneLength - 1; i++)
		crossoverTable[i] = i + 1;
}

void Island::immigrate()
{
	if(migrantSize == 0)
		return;

	gene.resize(nativeSize);

	for(unsigned i = 0; i < migrant.size(); i++)
		gene.push_back(migrant[i]);

	//random_shuffle(gene.begin(), gene.end(), genrand_int);
	random_shuffle(gene.begin(), gene.end());
}

void Island::crossover()
{
	nPointCrossover();
}

void Island::nPointCrossover()
{
	int i, j, k;
	int median;
	vector<Gene> child(2);
	vector<int> crossoverPoint;
	median = gene.size() / 2;

	for(i = 0; i < median; i++)
	{
		if(genrand() < crossoverRate)
		{
		    random_shuffle(crossoverTable.begin(), crossoverTable.end(), genrand_int);
			
			crossoverPoint.resize(numCrossoverPoint);

			for(k = 0; k < numCrossoverPoint; k++)
				crossoverPoint[k] = crossoverTable[k];

			sort(crossoverPoint.begin(), crossoverPoint.end());

			child[0] = gene[i];
			child[1] = gene[i + median];

			for(unsigned k = 0; k < crossoverPoint.size(); k++)
			{
				for(j = crossoverPoint[k]; j < geneLength; j++)
				{
					if(k % 2 == 0)
					{
						child[0].bitString[j] = gene[i + median].bitString[j];
						child[1].bitString[j] = gene[i].bitString[j];
					}
					else
					{
						child[0].bitString[j] = gene[i].bitString[j];
						child[1].bitString[j] = gene[i + median].bitString[j];
					}
				}
			}
			gene[i] = child[0];
			gene[i + median] = child[1];
		}
	}
}


void Island::mutate()
{
	double probability;

	for(unsigned i = 0; i < gene.size(); i++)
	{
		for(unsigned j = 0; j < gene[i].bitString.size(); j++)
		{
			probability = genrand();
			if(probability < mutationRate)
				gene[i].bitString[j] ^= 1;
		}
	}
}

void Island::evaluate()
//���Τ�ɾ��
{
	for(unsigned i = 0; i < gene.size(); i++)
		fitness(gene[i]);
}


void Island::elitePreservation()
// ���꡼�Ȥ���¸
//��������Υ��꡼�Ȥȸ��ߤ�����������Τ򤢤碌�ơ����ξ�̤򥨥꡼�ȤȤ��롥
//�콸�ĥ������ϰ��Ū��popSize + numElite �ˤʤ롥
{
	vector<Gene> temp;
	temp = elite;
	
	sort(gene.begin(), gene.end());

	for(unsigned i = 0; i < elite.size(); i++)
		temp.push_back(gene[i]);

	sort(temp.begin(), temp.end());

	for(unsigned i = 0; i < elite.size(); i++){
		// ���꡼�ȸ��Τ򹹿�
		elite[i] = temp[i];
		// ���꡼���ᤷ
		gene.push_back(temp[i]);
	}

	bestFitness = elite[0].fitness;
}


void Island::setGeneFromFile(ifstream& fin)

{
	for(unsigned i = 0; i < gene.size(); ++i){
		vector<char>& bitString = gene[i].bitString;
		for(int j = 0; j < numVariable; ++j){
			double variable;
			fin >> variable;

			if(fin.fail())
				return;
			ᤷ
			if(variable < closedMin)
				variable = closedMin;
			else if(variable >= openMax)
				variable = openMax - delta;


			int num = static_cast<int>((variable - closedMin) / delta + 0.5);


			gray_code(bitString, j * variableLength, variableLength, num);
		}
	}
}


void Island::reproduce()
{
	if(selectionType == Tournament)
		tournamentSelection();
	else if(selectionType == Roulette)
		rouletteSelection();
	else{
		cerr << endl << "Invalid selection" << endl << "Exit" << endl;
		exit(1);
	}
	gene.resize(newGeneration.size());
	for(unsigned i = 0; i < gene.size(); i++)
		gene[i] = newGeneration[i];
}


//tournamentSize
void Island::tournamentSelection()
{
	vector<Gene> tournament;
	tournament.resize(tournamentSize);

	for(unsigned i=0; i< newGeneration.size(); i++){
		for(unsigned j=0; j< tournament.size(); j++){
			tournament[j] = gene[genrand_int(gene.size())];
		}
		sort(tournament.begin(),tournament.end());
		newGeneration[i] = tournament[0];
	}
}


void Island::rouletteSelection()
{


	double min_fit = DBL_MAX;
	for(unsigned i = 0; i < gene.size(); ++i)
		if(min_fit > gene[i].fitness)
			min_fit = gene[i].fitness;

	// make roulette
	roulette.resize(gene.size());
	roulette[0] = gene[0].fitness - min_fit;
	for(unsigned i = 1; i < gene.size(); ++i)
		roulette[i] = roulette[i-1] + gene[i].fitness - min_fit;
	double sum_of_fitness = roulette.back();
	assert(sum_of_fitness > 0);

	for(unsigned i = 0; i < newGeneration.size(); i++){
		// (0 <= r <= sum_of_fitness)
		double r = sum_of_fitness * genrand();

		vector<double>::iterator itr =
			std::lower_bound(roulette.begin(), roulette.end(), r);

		newGeneration[i] = gene[itr - roulette.begin()];
	}
}

void Island::fitness(const Gene& gene) const
{
	int		i;
	vector<double> real(numVariable);

	if(functionType == User){
		real = problem.decode(gene.bitString);
	}
	else{
		for(i = 0; i < numVariable; i++){
		  if(typeOfEncode == 0)
						real[i] = gray_decode(gene.bitString, i * variableLength, variableLength) * delta + closedMin;
		  
		
		  else
				real[i] = binary_decode(gene.bitString, i * variableLength, variableLength) * delta + closedMin;
		}
	}

	switch(functionType){
	case Rastrigin  : gene.fitness = rastrigin(real);        break;
	case Schwefel   : gene.fitness = schwefel(real);         break;
	case Griewank   : gene.fitness = griewank(real);         break;
	case Rosenbrock : gene.fitness = rosenbrock(real);       break;
	case Ridge      : gene.fitness = ridge(real);            break;
	case Onemax     : gene.fitness = onemax(gene.bitString); break;
	case User       : gene.fitness = user_function(real);    break;
	default :
		cerr << endl << "invalid function" << endl << "Exit" << endl;
		exit(1);
	}

	numEvaluate++;

// debug
//	copy(real.begin(), real.end(), ostream_iterator<double>(cerr, "\t"));
//	cerr << " \t" << gene.fitness << endl;
}

//=============================================================================
//class Globe
//=============================================================================

Globe::Globe()
{

}

Globe::~Globe()
{
	if(island_0.functionType == User)
		problem.finalize();
}

void Globe::initialize()
{
// Default settings

	fTotalPopSize = false;
	totalPopulationSize			= 400;

	// DGA parameters
	island_0.migrationRate		= 0.5;
	numIsland					= 40;
	migrationInterval			= 5;

	// parameters for chromosome
	island_0.geneLength			= 100;
	bestBitString.resize(100);
	island_0.numVariable		= 10;

	fInitialPopulation			= 0;
		// 0: random,  1: from file(design variables)

	// parameters for objective function
	island_0.functionType		= Rastrigin;
	island_0.openMax			= 5.12;
	island_0.closedMin			= -5.12;
	island_0.finValue			= -1E-10;

	// canonical GA parameters
	assert(totalPopulationSize % numIsland == 0);
	island_0.populationSize  	= totalPopulationSize / numIsland;// 10;
		// in each subpopulations
	island_0.numElite			= island_0.populationSize/2; //5;
	island_0.crossoverRate		= 1.0;
	island_0.mutationRate		= 1.0/island_0.geneLength; // 0.01;
	island_0.tournamentSize     = 4;
	island_0.selectionType		= Tournament; // 0: tournament,  1: roulette

	island_0.numCrossoverPoint	= 2;
	island_0.typeOfEncode		= 0;


	island_0.numEvaluate		= 0;
	maxGeneration		  		= 1000;


	maxTrial			= 1;
	sumMaxFitness		= 0;
	sumGeneration		= 0;
	numConvergent		= 0;
	sumTotalEvaluate	= 0;
	fStop				= true;
	fOutput				= false;

	reviseParameters();
}

void Globe::initialize(int argc, char *argv[])
{
	initialize();

	// command line
	for(int i = 1; i < argc; i++)
	{
		const char *p = argv[i] + 1;
		switch(*argv[i])
		{
		case 'b':
			island_0.geneLength = atoi(p);
			// change mutation rate
			island_0.mutationRate = 1.0/island_0.geneLength;
			break;

		case 'i':
			numIsland = atoi(p);
			break;

		case 'c':
			island_0.crossoverRate = atof(p);
			break;

		case 'f':
			island_0.functionType = static_cast<FunctionType>(atoi(p));
			if(atoi(p) == Rastrigin){
				island_0.closedMin = -5.12;
				island_0.openMax = 5.12;
			}
			else if(atoi(p) == Schwefel || atoi(p) == Griewank){
				island_0.closedMin = -512;
				island_0.openMax = 512;
			}
			else if(atoi(p) == Rosenbrock){
				island_0.closedMin = -2.048;
				island_0.openMax = 2.048;
			}
			else if(atoi(p) == Ridge){
				island_0.closedMin = -64;
				island_0.openMax = 64;
			}
			else if(atoi(p) == Onemax) //onemax
			{
				island_0.closedMin = 0.0;
				island_0.openMax = 0.0;
			}
			break;

		case 'I':
			migrationInterval = atoi(p);
			break;

		case 'm':
			island_0.mutationRate = atof(p);
			break;

		case 'M':
			island_0.migrationRate = atof(p);
			break;

		case 'o':
			fOutput = true;
			break;

		case 'p':
			island_0.populationSize = atoi(p);
//			island_0.numElite = island_0.populationSize / 2;
			break;

		case 'P':
			totalPopulationSize = atoi(p);
			fTotalPopSize = true;
			break;

		case 'z':
			maxTrial = atoi(p);
			break;

		case 'g':
			maxGeneration = atoi(p);
			break;

		case 'S':
			// 0: random,  1: from file(design variables)
			fInitialPopulation	= atoi(p);
			break;

		case 's':
			// 0: tournament,  1: roulette
			island_0.selectionType = static_cast<SelectionType>(atoi(p));
			break;

		case 'D':
			// number of variable
			island_0.numVariable = atoi(p);
			break;

		case 'e':
			// number of elites
			island_0.numElite = atoi(p);
			break;
		case 'F':
		    //set threshold
		    island_0.finValue = -1*atof(p);
		    break;
		}
	}

	reviseParameters();
	
	if(island_0.functionType == User){
		problem.initialize();
	}
}

void Globe::reviseParameters()
{
	assert(island_0.geneLength % island_0.numVariable == 0);
	assert(island_0.populationSize % 2 == 0);
	assert(totalPopulationSize % 2 == 0);
	assert(island_0.crossoverRate >= 0 && island_0.crossoverRate <= 1);
	assert(island_0.migrationRate >= 0 && island_0.migrationRate <= 1);
	assert(island_0.mutationRate >= 0 && island_0.mutationRate <= 1);
	assert(migrationInterval >= 1);
	assert(maxGeneration >= 1);
	assert(fInitialPopulation == 0 || fInitialPopulation == 1);
	assert(island_0.numVariable < 32);
	assert(island_0.numElite >= 1);
	
	island_0.variableLength = island_0.geneLength / island_0.numVariable;

	island_0.delta =
		(island_0.openMax - island_0.closedMin)
		/ ((double) (1U << island_0.variableLength));
}

void Globe::establish()
{
	if(fTotalPopSize)
	{
		island_0.populationSize = 10;
		island_0.numElite = island_0.populationSize / 2;
		numIsland = totalPopulationSize / island_0.populationSize;
	}
	island_0.migrantSize	= static_cast<int>(island_0.populationSize * island_0.migrationRate);
	island_0.nativeSize	= island_0.populationSize - island_0.migrantSize;

	populationSize	= island_0.populationSize;
	nativeSize		= island_0.nativeSize;
	migrantSize		= island_0.migrantSize;

	maxFitness		= -DBL_MAX;
	convGeneration	= maxGeneration;
	fConv			= false;
	totalEvaluate	= 0;

	island.resize(numIsland);

	ifstream fin;
	if(fInitialPopulation == 1){

		fin.open(INITIAL_POPULATION_FILE);
		if(!fin.is_open()){
			cerr << __FILE__ << ":(" << __LINE__ << ") can't open file "
				<< INITIAL_POPULATION_FILE << endl;
			exit(1);
		}
	}
	for(unsigned i = 0; i < island.size(); i++)
	{
		island[i] = island_0;
		island[i].initialize(fin);
	}

	fin.close();
}

void Globe::emigrate()
{
	if(migrantSize == 0 || island.size() < 2)
		return;

	RandomRing	destination(numIsland);


	destination.shuffle();

	for(unsigned i = 0; i < island.size(); i++)
	{

		random_shuffle( island[i].gene.begin(), island[i].gene.end());
		//random_shuffle( island[i].gene.begin(), island[i].gene.end(), genrand_int );
		

		for(int j = 0; j < migrantSize; j++)
			island[ destination[i] ].migrant[j] = island[i].gene[nativeSize + j];
	}
}


void Globe::execute()
{

    for(generation = 0; generation < maxGeneration; generation++)
	{

		if(generation % 10 == 0)
			cerr << "." << flush;
		
		// selection
		for(unsigned i = 0; i < island.size(); i++)
			island[i].reproduce();

		if(generation % migrationInterval == 0)
			emigrate();
	     
		for(unsigned i = 0; i < island.size(); i++)
		{
			if(generation % migrationInterval == 0 && numIsland > 1)
				island[i].immigrate();

			
			island[i].crossover();
			island[i].mutate();
			island[i].evaluate();
			island[i].elitePreservation();
			
			if(island[i].bestFitness > maxFitness){
			    maxFitness = island[i].bestFitness;
				bestBitString=island[i].elite[0].bitString;
			}
		}

		if(!fConv && maxFitness >= island_0.finValue){
		  	finalize(true);
		}

		sumFitness[generation] += maxFitness;
	}
	finalize(false);
	
}

void Globe::finalize(bool b)
{
	if(b)
	{
		maxFitness = 0.000000;
		sumMaxFitness += maxFitness;
		convGeneration = generation + 1;
		sumGeneration += convGeneration;
		numConvergent++;
		cerr << "*" << flush;
		fConv = true;
		
		if(fStop)
			for(; generation < maxGeneration; generation++)
				sumFitness[generation] += maxFitness;

		for(unsigned i = 0; i < island.size(); i++)
			totalEvaluate += island[i].numEvaluate;
		sumTotalEvaluate += totalEvaluate;
		
	}
	else if(!b && !fConv)
	{
		sumMaxFitness += maxFitness;
		for(unsigned i = 0; i < island.size(); i++)
			totalEvaluate += island[i].numEvaluate;
		sumTotalEvaluate += totalEvaluate;
	}
}


//�¹����˥ѥ�᡼����cerr�ǽ��Ϥ���
void Globe::showParameters() const
{
	cerr << "=========================================================" << endl;
	if(island_0.functionType != User)
		cerr << "Function       : " << function_name[island_0.functionType] << endl;
	else
		cerr << "Problem        : " << problem.problemName() << endl;
	cerr << "Chrom. Length  : " << island_0.geneLength << endl;
	if(island_0.functionType != User){
		cerr << "Num Variables  : " << island_0.numVariable
			<< " ( "  << island_0.variableLength << "bits/variable )" << endl;
		cerr << "Domain         : " << "[ " << island_0.closedMin << ", " << island_0.openMax << " ) (Resolution : " << island_0.delta  << ")" <<  endl;
	}
	else{
		cerr << "Num Variables  : " << problem.numVariable() << endl;
	}
	cerr << endl;

	cerr << "Population     : " << island_0.populationSize * numIsland << endl;
	cerr << "Num island     : " << numIsland << endl;
	cerr << "Pop/island     : " << island_0.populationSize << endl;
	cerr << "Elite/island   : " << island_0.numElite << endl;
	cerr << "init population: " << ((fInitialPopulation == 0) ? "random" : "from file") << endl;

	cerr << "Mig. Rate      : " << island_0.migrationRate << endl;
	cerr << "Mig. Interval  : " << migrationInterval << endl ;
	cerr << "Tournament size: " << island_0.tournamentSize << endl;
	cerr << "Selection      : "
		<< selection_name[island_0.selectionType] << endl;
		// 0: tournament,  1: roulette
	cerr << "Pc             : " << island_0.crossoverRate << endl;
	cerr << "Pm             : " << island_0.mutationRate << "  ( "
		<< island_0.mutationRate * island_0.geneLength << " / L )" << endl;
	cerr << endl;

	cerr << "Generations    : " << maxGeneration << endl;
	cerr << "Threshold      : " << island_0.finValue << endl;
	cerr << "Trials         : " << maxTrial << endl;
	cerr << "=========================================================" << endl;
	cerr << endl;
}


