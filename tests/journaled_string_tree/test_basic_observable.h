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
// Author: Rene Rahn <rene.rahn@fu-berlin.de>
// ==========================================================================
// Tests for the observable pattern.
// ==========================================================================

#ifndef TESTS_JOURNALED_STRING_TREE_TEST_BASIC_OBSERVABLE_H_
#define TESTS_JOURNALED_STRING_TREE_TEST_BASIC_OBSERVABLE_H_

#include <seqan/basic.h>
#include <seqan/journaled_string_tree.h>

using namespace seqan;

namespace test_basic_observer

{
struct Foo_
{

    unsigned i = 0;

    Foo_() : i(1)
    {}

    Foo_(unsigned num) : i(num)
    {}
};

struct Bar_
{

    unsigned i = 0;

    Bar_() : i(2)
    {}

    Bar_(unsigned num) : i(num)
    {}
};

inline void notify(Foo_ & f, Nothing const & /*tag*/)
{
    f.i +=2;
}

inline void notify(Bar_ & b, Nothing const & /*tag*/)
{
    b.i += 3;
}

}

SEQAN_DEFINE_TEST(test_basic_observable_length)
{
    {  // Length 2.
        ObserverList<test_basic_observer::Bar_, test_basic_observer::Foo_> subject;
        SEQAN_ASSERT(length(subject) == 2u);
    }

    {  // Length 3.
        ObserverList<test_basic_observer::Bar_, test_basic_observer::Bar_, test_basic_observer::Foo_> subject;
        SEQAN_ASSERT(length(subject) == 3);
    }

    {  // Length 0.
        ObserverList<> subject;
        SEQAN_ASSERT(length(subject) == 0u);
    }

}

SEQAN_DEFINE_TEST(test_basic_observable_set_observer)
{

    ObserverList<test_basic_observer::Bar_, test_basic_observer::Foo_, test_basic_observer::Bar_> subject;
    test_basic_observer::Bar_ bar(4);
    setObserver<0>(subject, bar);
    setObserver<1>(subject, test_basic_observer::Foo_(6));
    setObserver<2>(subject, test_basic_observer::Bar_(8));

    SEQAN_ASSERT(std::get<0>(subject).i == 4);
    SEQAN_ASSERT(std::get<1>(subject).i == 6);
    SEQAN_ASSERT(std::get<2>(subject).i == 8);
}

SEQAN_DEFINE_TEST(test_basic_observable_notify)
{
    ObserverList<test_basic_observer::Bar_, test_basic_observer::Foo_, test_basic_observer::Bar_> subject;
    test_basic_observer::Bar_ bar(4);
    setObserver<0>(subject, bar);
    setObserver<1>(subject, test_basic_observer::Foo_(6));
    setObserver<2>(subject, test_basic_observer::Bar_(8));

    notify(subject, Nothing());

    SEQAN_ASSERT(std::get<0>(subject).i == 7);
    SEQAN_ASSERT(std::get<1>(subject).i == 8);
    SEQAN_ASSERT(std::get<2>(subject).i == 11);
}

SEQAN_DEFINE_TEST(test_basic_observable_make_observer_list)
{
    auto list = makeObserverList(test_basic_observer::Bar_(4), test_basic_observer::Foo_(6), test_basic_observer::Bar_(5));

    SEQAN_ASSERT(std::get<0>(list).i == 4);
    SEQAN_ASSERT(std::get<1>(list).i == 6);
    SEQAN_ASSERT(std::get<2>(list).i == 5);
}

#endif // TESTS_JOURNALED_STRING_TREE_TEST_BASIC_OBSERVABLE_H_
