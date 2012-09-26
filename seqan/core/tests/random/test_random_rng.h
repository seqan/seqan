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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Tests for the random number generation code in seqan/random.
// ==========================================================================

#ifndef TEST_RANDOM_TEST_RANDOM_RNG_H_
#define TEST_RANDOM_TEST_RANDOM_RNG_H_

// Test GetDefaultRang and defaultRng().
SEQAN_DEFINE_TEST(test_default_rng)
{
    using namespace seqan;

    // Test that calling the function works.
    typedef String<Dna> TTag;
    typedef typename GetDefaultRng<TTag>::Type TRng;
    TRng & rng = defaultRng(TTag());
    (void)rng;

    // Test that a reference is returned and the global state changes.
    typedef typename Value<TRng>::Type TValue;
    TValue x1 = pickRandomNumber(defaultRng(TTag()));
    TValue x2 = pickRandomNumber(defaultRng(TTag()));
    TValue x3 = pickRandomNumber(defaultRng(TTag()));
    SEQAN_ASSERT(x1 != x2 || x2 != x3);  // 3 times the same value is not probable!
}


// Construct MersenneTwister in all possible ways.
SEQAN_DEFINE_TEST(test_random_mt19937_constructors)
{
    using namespace seqan;

    {
        Rng<> mt;
    }
    {
        Rng<MersenneTwister> mt;
    }
    {
        Rng<MersenneTwister> mt(10);
    }
}

// Pick random numbers from the MT and make sure the same number is
// not returned twice in the first two picks.
SEQAN_DEFINE_TEST(test_random_mt19937_pick)
{
    using namespace seqan;

    Rng<MersenneTwister> mt(10);
    // Test pickRandomNumber().
    SEQAN_ASSERT_NEQ(pickRandomNumber(mt), pickRandomNumber(mt));
    // Test operator().
    SEQAN_ASSERT_NEQ(mt(), mt());
}

// Test metafunctiosn for mersenne twister RNG.
SEQAN_DEFINE_TEST(test_random_mt19937_metafunctions)
{
	using namespace seqan;
	
	typedef Value<Rng<MersenneTwister> >::Type TValue;
	
	TValue m = MinValue<Rng<MersenneTwister> >::VALUE;
	SEQAN_ASSERT_EQ(MinValue<TValue>::VALUE, m);
	TValue M = MaxValue<Rng<MersenneTwister> >::VALUE;
	SEQAN_ASSERT_EQ(MaxValue<TValue>::VALUE, M);
}

// Construct RngFunctor specialization in all possible ways.
SEQAN_DEFINE_TEST(test_random_rng_functor_constructors)
{
    using namespace seqan;
    
	typedef Rng<MersenneTwister> TMersenneTwister;
	typedef Pdf<Uniform<int> > TUniformPdf;
    
    {
        TMersenneTwister mt;
        TUniformPdf uniformPdf(0, 10);
        
        Rng<RngFunctor<TMersenneTwister, TUniformPdf> > rng(mt, uniformPdf);
    }
}

// Pick random number from RngFunctor and make sure it is the same as when
// directly using a MT and a Pdf.
SEQAN_DEFINE_TEST(test_random_rng_functor_pick)
{
    using namespace seqan;
    
    const int SEED = 10;

    typedef Rng<MersenneTwister> TMersenneTwister;
	typedef Pdf<Uniform<int> > TUniformPdf;
    typedef Rng<RngFunctor<TMersenneTwister, TUniformPdf> > TRngFunctor;

    // Compute by using the raw MT and uniform Pdf.
    String<int> rawInts;
    {
	    TMersenneTwister mt(SEED);
    	TUniformPdf uniform(10, 100);
        
        for (int i = 0; i < 100; ++i)
            appendValue(rawInts, pickRandomNumber(mt, uniform));
    }
    
    // Use RngFunctor with pickRandomNumber() and check equality.
    {
	    TMersenneTwister mt(SEED);
    	TUniformPdf uniform(10, 100);
	    TRngFunctor rngFunctor(mt, uniform);
        
        for (int i = 0; i < 100; ++i) {
            SEQAN_ASSERT_EQ_MSG(unsigned(rawInts[i]), pickRandomNumber(rngFunctor), "i = %d", i);
        }
    }
    
    // Use RngFunctor with operator() and check equality.
    {
	    TMersenneTwister mt(SEED);
    	TUniformPdf uniform(10, 100);
	    TRngFunctor rngFunctor(mt, uniform);
        
        for (int i = 0; i < 100; ++i) {
            SEQAN_ASSERT_EQ_MSG(rawInts[i], rngFunctor(), "i = %d", i);
        }
    }
}

#endif  // TEST_RANDOM_TEST_RANDOM_RNG_H_
