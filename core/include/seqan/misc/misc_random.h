// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================

#ifndef SEQAN_HEADER_MISC_RANDOM_H
#define SEQAN_HEADER_MISC_RANDOM_H

#include <cmath>
#include <cstdlib>
#include <ctime>

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
// Mersenne Twister Random Number Generator
// Implementation by Michael Brundage;
// with some modifications
// http://www.qbrundage.com/michaelb/pubs/essays/random_number_generation.html
//////////////////////////////////////////////////////////////////////////////


#define SEQAN_MERSENNE_MT_LEN			624
#define SEQAN_MERSENNE_MT_IA			397
#define SEQAN_MERSENNE_MT_IB			(SEQAN_MERSENNE_MT_LEN - SEQAN_MERSENNE_MT_IA)
#define SEQAN_MERSENNE_UPPER_MASK      0x80000000
#define SEQAN_MERSENNE_LOWER_MASK      0x7FFFFFFF
#define SEQAN_MERSENNE_MATRIX_A        0x9908B0DF
#define SEQAN_MERSENNE_TWIST(b,i,j)    ((b)[i] & SEQAN_MERSENNE_UPPER_MASK) | ((b)[j] & SEQAN_MERSENNE_LOWER_MASK)
#define SEQAN_MERSENNE_MAGIC(s)        (((s)&1)*SEQAN_MERSENNE_MATRIX_A)

//////////////////////////////////////////////////////////////////////////////
//forward declaration

inline unsigned mtRand(); 
inline void mtRandInit();

//////////////////////////////////////////////////////////////////////////////

template <typename T = void>
struct MersenneBuffer_
{
    static unsigned long buffer[SEQAN_MERSENNE_MT_LEN];
	static int index;
	static bool is_initialized;

};
//____________________________________________________________________________

template <typename T>
unsigned long MersenneBuffer_<T>::buffer[SEQAN_MERSENNE_MT_LEN];

//____________________________________________________________________________

template <typename T>
int MersenneBuffer_<T>::index = 0;

//____________________________________________________________________________

template <typename T>
bool MersenneBuffer_<T>::is_initialized = false;

//////////////////////////////////////////////////////////////////////////////
// NOTE: mtRandInit() must have been called at least once before mtRand() is used.
// Can also be called several times since it is protected against multiple initalizations. 

/**
.Function.mtRandInit:
..status:deprecated, use @Class.Rng@ and @Class.Pdf@ from $<seqan/random.h>$ instead
..summary:Initialize the Mersenne-Twister random number generator.
..cat:Random
..signature:mtRandInit()
..signature:mtRandInit(doSRand)
..param.doSRand:If $true$, the Mersenne-Twister is initialized with a random number from $std::rand$.
...type:nolink:bool
..remarks:@Function.mtRandInit@ must have been called at least once before @Function.mtRand@ is used.
..include:seqan/misc.h
*/

inline void 
mtRandInit(bool 
#ifndef SEQAN_NOSRAN
		   _doSRand
#endif
		   )
{
	// test whether mtRandInit was already initialized
	// return immediately if this is the case
	if (MersenneBuffer_<>::is_initialized) return;
	MersenneBuffer_<>::is_initialized = true;

#ifndef SEQAN_NOSRAN
	if (_doSRand)
		::std::srand((unsigned) ::std::time(0));
#endif

	int i;
	for (i = 0; i < SEQAN_MERSENNE_MT_LEN; i++)
		MersenneBuffer_<>::buffer[i] = ::std::rand();

	mtRand(); //pop the first number, since it is not as "random" as we like it
}

inline void 
mtRandInit()
{
	mtRandInit(true);
}

//////////////////////////////////////////////////////////////////////////////
// NOTE: mtRandInit() must be called once before mtRand() is used.


/**
.Function.mtRand:
..status:deprecated, use @Class.Rng@ and @Class.Pdf@ from $<seqan/random.h>$ instead
..summary:Return a Mersenne-Twister random number.
..cat:Random
..signature:mtRand()
..returns:A random number between 0 and $MaxValue<unsigned>::VALUE$.
...type:nolink:unsigned
..remarks:@Function.mtRandInit@ must have been called at least once before @Function.mtRand@ is used.
..see:Function.mtRandInit
..include:seqan/misc.h
*/

inline unsigned 
mtRand()
{
	unsigned long * b = MersenneBuffer_<>::buffer;
	int idx = MersenneBuffer_<>::index;
	unsigned long s;
	int i;
	
	if (idx == SEQAN_MERSENNE_MT_LEN*sizeof(unsigned long))
	{
		idx = 0;
		i = 0;
		for (; i < SEQAN_MERSENNE_MT_IB; i++) {
			s = SEQAN_MERSENNE_TWIST(b, i, i+1);
			b[i] = b[i + SEQAN_MERSENNE_MT_IA] ^ (s >> 1) ^ SEQAN_MERSENNE_MAGIC(s);
		}
		for (; i < SEQAN_MERSENNE_MT_LEN-1; i++) {
			s = SEQAN_MERSENNE_TWIST(b, i, i+1);
			b[i] = b[i - SEQAN_MERSENNE_MT_IB] ^ (s >> 1) ^ SEQAN_MERSENNE_MAGIC(s);
		}
        
		s = SEQAN_MERSENNE_TWIST(b, SEQAN_MERSENNE_MT_LEN-1, 0);
		b[SEQAN_MERSENNE_MT_LEN-1] = b[SEQAN_MERSENNE_MT_IA-1] ^ (s >> 1) ^ SEQAN_MERSENNE_MAGIC(s);
	}
	MersenneBuffer_<>::index = idx + sizeof(unsigned long);
	return *(unsigned long *)((unsigned char *)b + idx);
}

/**
.Function.mtRandDouble:
..cat:Random
..status:deprecated, use @Class.Rng@ and @Class.Pdf@ from $<seqan/random.h>$ instead
..summary:Return a random number between 0 and 1 using mtRand.
..see:Function.mtRand
..include:seqan/misc.h
*/
inline double
mtRandDouble()
{
    return static_cast<double>(mtRand()) / static_cast<double>(MaxValue<unsigned>::VALUE);
}

//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// Geometrical Distribution Heuristics
// from Numerical Recipes in C
//////////////////////////////////////////////////////////////////////////////


#define SEQAN_RG_IB1 1
#define SEQAN_RG_IB2 2
#define SEQAN_RG_IB5 16
#define SEQAN_RG_IB18 131072
#define SEQAN_RG_MASK ( SEQAN_RG_IB1 + SEQAN_RG_IB2 + SEQAN_RG_IB5 )

template <typename T>
inline T
geomRand()
{
	static unsigned long seed = ::std::rand();
	T value = 0;
	while ( true )
	{
		if( ( seed & SEQAN_RG_IB18 ) )
		{
			seed = ( ( seed ^ SEQAN_RG_MASK ) << 1 ) | SEQAN_RG_IB1;
			++value;
		}
		else 
		{
			seed <<= 1;
			break;
		}
	}
	return value;
}

//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// Normal Distribution Heuristics.
//
// Ported from Python.
//
// Kinderman and Monahan method. Reference: Kinderman,
// A.J. and Monahan, J.F., "Computer generation of random
// variables using the ratio of uniform deviates", ACM Trans
// Math Software, 3, (1977), pp257-260.
//////////////////////////////////////////////////////////////////////////////

static const double SEQAN_NV_MAGICCONST = 1.7155277699214135;  // == 4 * exp(-0.5)/sqrt(2.0)

inline double
normRand(double mean, double stddev)
{
    /* in Python:
        random = self.random
        while 1:
            u1 = random()
            u2 = 1.0 - random()
            z = NV_MAGICCONST*(u1-0.5)/u2
            zz = z*z/4.0
            if zz <= -_log(u2):
                break
        return mu + z*sigma
    */
    double z;
    while (true) {
        double u1 = mtRandDouble();
        double u2 = 1 - mtRandDouble();
        z = SEQAN_NV_MAGICCONST * (u1 - 0.5) / u2;
        double zz = z * z / 4.0;
        if (zz < -::std::log10(u2))
            break;
    }
    return mean + z * stddev;
}

//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// Lognormal Distribution Heuristics.

inline double
lognormRand(double mean, double stddev)
{
  return ::std::exp(normRand(mean, stddev));
}

//////////////////////////////////////////////////////////////////////////////

}  // namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
