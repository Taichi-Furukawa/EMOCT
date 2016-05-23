//--------------------------------------------------------------------------
// mtrand.cpp
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
//  and this list of conditions.
// 2. Redistributions in binary form must reproduce the above  copyright notice
//  and this list of conditions in the documentation and/or other materials
//  provided with the distribution.
// 3. All advertising materials mentioning features or use of this software
//  must display the following acknowledgement:
//
//  This product includes software developed by Takashi Kawasaki,
//  Tomoyuki Hiroyasu and Mitsunori Miki except the functions of sgenrand()
//  and genrand().
//
// 4. The names of Takashi Kawasaki, Tomoyuki Hiroyasu and Mitsunori Miki
//  should not be used to endorse or promote products derived from this
//  software without specific prior written permission.
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

#if defined(_WIN32) && !defined(__CYGWIN__)
#define WIN32_PLATFORM 1
#include <windows.h>
#else
#define UNIX_PLATFORM 1
#include <sys/time.h>
#include <unistd.h>
#endif

#include "mtrand.h"

//   Period parameters
#define N N_MTRAND_FACTOR
#define M 397
#define MATRIX_A 0x9908b0df   /* constant vector a */
#define UPPER_MASK 0x80000000 /* most significant w-r bits */
#define LOWER_MASK 0x7fffffff /* least significant r bits */

// Tempering parameters
#define TEMPERING_MASK_B 0x9d2c5680
#define TEMPERING_MASK_C 0xefc60000
#define TEMPERING_SHIFT_U(y)  (y >> 11)
#define TEMPERING_SHIFT_S(y)  (y << 7)
#define TEMPERING_SHIFT_T(y)  (y << 15)
#define TEMPERING_SHIFT_L(y)  (y >> 18)

// This function makes number from time
int seed_by_time()
{
	static int co = 0;
	static int prev_value = 0;
	int value;

#if defined(WIN32_PLATFORM)
	LARGE_INTEGER tt;
	QueryPerformanceCounter(&tt);
	value = (int)tt.LowPart;

#else // UNIX
	struct timezone tz;
	struct timeval tv;
	gettimeofday(&tv, &tz);
	value = (int)tv.tv_usec;
#endif

	co++;
	if(value == prev_value) value += co;
	prev_value = value;
	return value;
}

/* initializing the array with a NONZERO seed */
void MTRand::srand(unsigned long seed)
{
	int i;
	for(i = 0; i < N; i++)
	{
		mt[i] = seed & 0xffff0000;
		seed = 69069 * seed + 1;
		mt[i] |= (seed & 0xffff0000) >> 16;
		seed = 69069 * seed + 1;
	}
	mti = N;
}

MTRand::MTRand(unsigned long seed)
{
	this->srand(seed);
}

void MTRand::srandByTime()
{
	this->srand((unsigned long)seed_by_time());
}

MTRand::MTRand()
{
	this->srand((unsigned long)seed_by_time());
}

double MTRand::genrand() const
{
	unsigned long y;
	const unsigned long mag01[2] = {0x0, MATRIX_A};
	/* mag01[x] = x * MATRIX_A  for x=0,1 */

	/* generate N words at one time */
	if (mti >= N)
	{
		int kk;

		for (kk = 0; kk < N - M; kk++) {
			y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
			mt[kk] = mt[kk + M] ^ (y >> 1) ^ mag01[y & 0x1];
		}
		for (;kk< N - 1; kk++) {
			y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
			mt[kk] = mt[kk + (M - N)] ^ (y >> 1) ^ mag01[y & 0x1];
		}
		y = (mt[N - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
		mt[N - 1] = mt[M - 1] ^ (y >> 1) ^ mag01[y & 0x1];

		mti = 0;
	}

	y = mt[mti++];
	y ^= TEMPERING_SHIFT_U(y);
	y ^= TEMPERING_SHIFT_S(y) & TEMPERING_MASK_B;
	y ^= TEMPERING_SHIFT_T(y) & TEMPERING_MASK_C;
	y ^= TEMPERING_SHIFT_L(y);

	return (double)y * 2.3283064365386963e-10; /* reals: [0,1)-interval */
	/* return y; */ /* for integer generation */
}

int MTRand::genrand_int(int size) const
{
	return (int)(genrand() * (double)size);
}

static MTRand mtrand_global;

void sgenrand(unsigned long seed)
{
	mtrand_global.srand(seed);
}

double genrand()
{
	return mtrand_global.genrand();
}

int genrand_int(int n)
{
	return mtrand_global.genrand_int(n);
}

void ShuffledIntGenerator::resize(int n)
{
	if(ar) delete[] ar;
	ar = new int[n];
	int i;
	for(i = 0; i < n; i++)
		ar[i] = i;

	sz = n;
	shuffle();
}

void ShuffledIntGenerator::shuffle(int deepness)
{
	int i;
	for(i = 0; i < deepness * sz; i++)
	{
		int& a = ar[mt.genrand_int(sz)];
		int& b = ar[mt.genrand_int(sz)];
		int tmp = a;
		a = b;
		b = tmp;
	}
}
