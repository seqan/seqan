// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
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
// Tests for the control structure part of the metaprogramming library.
// ==========================================================================

#ifndef SEQAN_TESTS_BASIC_TEST_BASIC_METAPROGRAMMING_CONTROL_H_
#define SEQAN_TESTS_BASIC_TEST_BASIC_METAPROGRAMMING_CONTROL_H_

SEQAN_DEFINE_TEST(test_basic_metaprogramming_control_if)
{
    using namespace seqan;

    bool b = IfC<true, False, True>::Type::VALUE;
    SEQAN_ASSERT_EQ(b, false);
    b = IfC<false, False, True>::Type::VALUE;
    SEQAN_ASSERT_EQ(b, true);
}

// Support code for test_basic_metaprogramming_control_switch.

struct Nothing_ {};

int switchTest(Nothing_ const &) { return -1; }
int switchTest(seqan::False const &) { return 0; }
int switchTest(seqan::True const &) { return 1; }
int switchTest(seqan::NilCase const &) { return 2; }

template <int X>
struct SwitchTest
{
    typedef typename seqan::Switch<
        X,
        seqan::Case<-1, Nothing_,
        seqan::Case<0, seqan::False,
        seqan::Case<1, seqan::True
        > > > >::Type Type;
};

SEQAN_DEFINE_TEST(test_basic_metaprogramming_control_switch)
{
    using namespace seqan;

    typedef typename SwitchTest<-1>::Type T1;
    SEQAN_ASSERT_EQ(switchTest(T1()), -1);

    typedef typename SwitchTest<0>::Type T2;
    SEQAN_ASSERT_EQ(switchTest(T2()), 0);

    typedef typename SwitchTest<1>::Type T3;
    SEQAN_ASSERT_EQ(switchTest(T3()), 1);

    typedef typename SwitchTest<2>::Type T4;
    SEQAN_ASSERT_EQ(switchTest(T4()), 2);
}

// Support code for test_basic_metaprogramming_control_loop.

struct LoopTestWorker_
{
    template <typename TArg>
    static inline void body(TArg & arg, int I)
    {
        arg[I - 1] += I;
    }
};

SEQAN_DEFINE_TEST(test_basic_metaprogramming_control_loop)
{
    using namespace seqan;

    int i[5] = { 1, 2, 3, 4, 5 };
    Loop<LoopTestWorker_, 5>::run(i);

    SEQAN_ASSERT_EQ(i[0], 2);
    SEQAN_ASSERT_EQ(i[1], 4);
    SEQAN_ASSERT_EQ(i[2], 6);
    SEQAN_ASSERT_EQ(i[3], 8);
    SEQAN_ASSERT_EQ(i[4], 10);
}

SEQAN_DEFINE_TEST(test_basic_metaprogramming_control_loop_reverse)
{
    using namespace seqan;

    int i[5] = { 1, 2, 3, 4, 5 };
    LoopReverse<LoopTestWorker_, 5>::run(i);

    SEQAN_ASSERT_EQ(i[0], 2);
    SEQAN_ASSERT_EQ(i[1], 4);
    SEQAN_ASSERT_EQ(i[2], 6);
    SEQAN_ASSERT_EQ(i[3], 8);
    SEQAN_ASSERT_EQ(i[4], 10);
}

#endif  // #ifndef SEQAN_TESTS_BASIC_TEST_BASIC_METAPROGRAMMING_CONTROL_H_
