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
// ga2k.h
//-----------------------------------------------------------------------------

#ifndef _GA2K_H_
#define _GA2K_H_


#include <vector>
#include <fstream>
#include "libs/table.h"
using namespace std;
using namespace arch;

#include "../mtrand.h"
#include "randring.h"

#define INITIAL_POPULATION_FILE "design_variable.ini"

///////////////////////////////////////////////////////////////////////////////
// global variable
///////////////////////////////////////////////////////////////////////////////

enum FunctionType {
	Rastrigin = 0,
	Schwefel = 1,
	Griewank = 2,
	Rosenbrock = 3,
	Ridge = 4,
	Onemax = 5,
	User = 99
};

enum SelectionType {
	Tournament = 0,
	Roulette = 1
};


///////////////////////////////////////////////////////////////////////////////
// class Gene
///////////////////////////////////////////////////////////////////////////////

class Gene
{
  public:
	//data
	table	bitString;
	mutable double	fitness;

	//operator
	bool operator<(const Gene& gene) const;

	//method (constructor)
	Gene();
	~Gene();
	void initialize(int geneLength);
	void initialize(int geneLength, char value);
	void initialize(int geneLength, const vector<char>& code);

	//method (operator)
	void visualize();
};


///////////////////////////////////////////////////////////////////////////////
// class Island
///////////////////////////////////////////////////////////////////////////////

class Island
{
  public:
	//data
	int		populationSize;

	double	crossoverRate;
	double	mutationRate;

	int		numElite;

	double	migrationRate;
	int		migrantSize;
	int		nativeSize;
	int 	tournamentSize;

	int		geneLength;
	int		variableLength;
	int		numVariable;
	double	openMax;
	double	closedMin;
	double	delta;
	double	finValue;

	FunctionType functionType;
	int		typeOfEncode;
	int		numCrossoverPoint;
	SelectionType selectionType;

	mutable int		numEvaluate;

	double	bestFitness;

	vector<double>	roulette;
	vector<int>		crossoverTable;

//gene
	vector<Gene>	gene;
	vector<Gene>	migrant;
	vector<Gene>	elite;
	vector<Gene>	newGeneration;

//method (constructor)
	Island();
	~Island();
	void initialize(std::ifstream& fin);

//method (operator)
	void immigrate();
	void crossover();
	void nPointCrossover();
	void mutate();
	void evaluate();
	void reproduce();
	void rouletteSelection();
	void tournamentSelection();
	void elitePreservation();
	void setGeneFromFile(std::ifstream& fin);


	void fitness(const Gene& gene) const;
};


///////////////////////////////////////////////////////////////////////////////
// class Globe
///////////////////////////////////////////////////////////////////////////////

class Globe
{
  public:
	//data
	vector<Island>	island;
	Island			island_0;
    vector<char> bestBitString;

	int					generation;
	int					maxGeneration;
	int					convGeneration;

	int					migrationInterval;
	int					numIsland;
	int					populationSize;
	int					totalPopulationSize;
	int					migrantSize;
	int					nativeSize;

	int 				fInitialPopulation;

	double				maxFitness;
	int					maxTrial;
	double				sumMaxFitness;
	int					sumGeneration;
	int					numConvergent;
	int					totalEvaluate;
	int					sumTotalEvaluate;

	vector<double>		sumFitness;

	Gene				elite;
	bool				fStop;
	bool				fConv;
	bool				fTotalPopSize;
	bool				fOutput;

	//method (constructor)
	Globe();
	virtual ~Globe();
	virtual void initialize();
	virtual void initialize(int argc, char *argv[]);
	virtual void reviseParameters();
	virtual void establish();

	//method (operator)
	virtual void emigrate();
	virtual void execute();
	virtual void showParameters() const;
	virtual void finalize(bool b);
};


#endif //_GA2K_H_
