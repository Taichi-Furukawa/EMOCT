//--------------------------------------------------------------------------
// mtrand.h
// GAPPA Parallel Distributed Environment Genetic Algorithm Library
// Mersenne Twister Random Number Generator
//--------------------------------------------------------------------------
// Copyright (C) 2000 Takashi Kawasaki, Tomoyuki Hiroyasu and Mitsunori Miki,
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright  notice
//    and this list of conditions.
// 2. Redistributions in binary form must reproduce the above  copyright notice
//    and this list of conditions in the documentation and/or other materials
//    provided with the distribution.
// 3. All advertising materials mentioning features or use of this software
//    must display the following acknowledgement:
//
//    This product includes software developed by Takashi Kawasaki,
//    Tomoyuki Hiroyasu and Mitsunori Miki except the functions of sgenrand()
//    and genrand().
//
// 4. The names of Takashi Kawasaki, Tomoyuki Hiroyasu and Mitsunori Miki
//    should not be used to endorse or promote products derived from this
//    software without specific prior written permission.
//--------------------------------------------------------------------------
// MT C++ Version Written by T.Kawasaki 10/04/1999
// C++ version enables to generate multiple random number sequence.
// First developed for the purpose of implementation of DEGA.
//--------------------------------------------------------------------------
// A C-program for MT19937: Real number version([0,1)-interval) (1998/4/6)
//   genrand() generates one pseudorandom real number (double)
// which is uniformly distributed on [0,1)-interval, for each
// call. sgenrand(seed) set initial values to the working area
// of 624 words. Before genrand(), sgenrand(seed) must be
// called once. (seed is any 32-bit integer except for 0).
// Integer generator is obtained by modifying two lines.
//   Coded by Takuji Nishimura, considering the suggestions by
// Topher Cooper and Marc Rieffel in July-Aug. 1997.
//--------------------------------------------------------------------------
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Library General Public
// License as published by the Free Software Foundation; either
// version 2 of the License, or (at your option) any later
// version.
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU Library General Public License for more details.
// You should have received a copy of the GNU Library General
// Public License along with this library; if not, write to the
// Free Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
// 02111-1307  USA
//--------------------------------------------------------------------------
// Copyright (C) 1997 Makoto Matsumoto and Takuji Nishimura.
// When you use this, send an email to: matumoto@math.keio.ac.jp
// with an appropriate reference to your work.
//--------------------------------------------------------------------------
// REFERENCE
// M. Matsumoto and T. Nishimura,
// "Mersenne Twister: A 623-Dimensionally Equidistributed Uniform
// Pseudo-Random Number Generator",
// ACM Transactions on Modeling and Computer Simulation,
// Vol. 8, No. 1, January 1998, pp 3--30.

#ifndef _MT_RAND_CPP_H_
#define _MT_RAND_CPP_H_

#ifdef __cplusplus

#define N_MTRAND_FACTOR 624
class MTRand
{
public:
	explicit MTRand(unsigned long seed);
	MTRand();
	void srand(unsigned long seed);
	void srandByTime();

	double genrand() const;
	int genrand_int(int n) const;

private:
	mutable unsigned long mt[N_MTRAND_FACTOR];
	mutable int mti;
};

#endif

#ifdef __cplusplus
extern "C" {
#endif
int seed_by_time();

// MT C Compatible Interface
double genrand();
int genrand_int(int n);
void sgenrand(unsigned long seed);
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
class ShuffledIntGenerator
{
private:
	int *ar;
	int sz;
	MTRand mt;

public:
	void srand(int n) {mt.srand(n);}
	void shuffle(int deepness = 1);
	void resize(int n);
	explicit ShuffledIntGenerator(int n) : ar(NULL) {resize(n);}
	ShuffledIntGenerator(int n, int seed) : ar(NULL) {srand(seed); resize(n);}
	~ShuffledIntGenerator() {delete[] ar;}
	inline const int *array() {return ar;}
	inline int operator[](int n) const {return ar[n];}
	inline int size() const {return sz;}
};
#endif

#endif // _MT_RAND_CPP_H_
