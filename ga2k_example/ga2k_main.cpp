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
// ga2k_main.cpp
//-----------------------------------------------------------------------------

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
using namespace std;

#include "ga2k.h"


int main(int argc, char *argv[])
{
	int i, trial;
	Globe globe;

	if(argc < 2)
		globe.initialize();
	else
		globe.initialize(argc, argv);

	globe.sumFitness.resize(globe.maxGeneration);

	for(i = 0; i < globe.maxGeneration; i++)
			globe.sumFitness[i] = 0;


	for(trial = 0; trial < globe.maxTrial; trial++)
	{
	  
	  globe.establish();

		if(trial == 0)
			globe.showParameters();

		cerr << "trial " << trial + 1 << " |" << endl;
		globe.execute();

		cerr << endl << "fitness    : " << globe.maxFitness << endl;
		cerr << "gene       : ";
		for(i = 0; i < globe.island_0.geneLength; i++)
		  printf("%d",globe.bestBitString[i]);
		printf("\n");
		cerr << "generation : " << globe.convGeneration << endl;
		cerr << "evaluate   : " << globe.totalEvaluate << endl;
	}

	if(globe.fOutput)
		for(i = 0; i < globe.maxGeneration; i++)
			cout << globe.sumFitness[i] / globe.maxTrial << endl;
}
