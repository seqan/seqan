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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================

#ifndef TESTS_MODIFIER_TEST_MODIFIER_FUNCTORS_H_
#define TESTS_MODIFIER_TEST_MODIFIER_FUNCTORS_H_

#include <seqan/basic.h>
#include <seqan/modifier.h>

using namespace seqan;


SEQAN_DEFINE_TEST(test_modifier_functors_functor_upcase) {
    FunctorUpcase<char> func;

    SEQAN_ASSERT_EQ('A', func('a'));
    SEQAN_ASSERT_EQ('A', func('A'));
    SEQAN_ASSERT_EQ('!', func('!'));
}


SEQAN_DEFINE_TEST(test_modifier_functors_functor_lowcase) {
    FunctorLowcase<char> func;

    SEQAN_ASSERT_EQ('a', func('a'));
    SEQAN_ASSERT_EQ('a', func('A'));
    SEQAN_ASSERT_EQ('!', func('!'));
}


SEQAN_DEFINE_TEST(test_modifier_functors_dna_complement) {
    {
        FunctorComplement<Dna> func;

        SEQAN_ASSERT_EQ(Dna('C'), func(Dna('G')));
        SEQAN_ASSERT_EQ(Dna('G'), func(Dna('C')));
        SEQAN_ASSERT_EQ(Dna('A'), func(Dna('T')));
        SEQAN_ASSERT_EQ(Dna('T'), func(Dna('A')));
    }
    {
        FunctorComplement<Dna5> func;

        SEQAN_ASSERT_EQ(Dna5('C'), func(Dna5('G')));
        SEQAN_ASSERT_EQ(Dna5('G'), func(Dna5('C')));
        SEQAN_ASSERT_EQ(Dna5('A'), func(Dna5('T')));
        SEQAN_ASSERT_EQ(Dna5('T'), func(Dna5('A')));
        SEQAN_ASSERT_EQ(Dna5('N'), func(Dna5('N')));
    }
}

#endif  // TESTS_MODIFIER_TEST_MODIFIER_FUNCTORS_H_
