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
// Author: Andreas Gogol-Doering <andreas.doering@mdc-berlin.de>
// ==========================================================================

#ifndef TESTS_BASIC_TEST_BASIC_PROXY_H_
#define TESTS_BASIC_TEST_BASIC_PROXY_H_

#include <sstream>

// TODO(holtgrew): Use the helper struct from test construct/destruct.

SEQAN_DEFINE_TEST(test_basic_proxy_iterator)
{
    using namespace seqan;

    int i1[] = {10, 20, 30};
    int * pi1 = i1;
    Proxy<IteratorProxy<int *> > px(pi1);
    SEQAN_ASSERT_EQ(px, 10);

//assign
    px = 11;
    SEQAN_ASSERT_EQ(i1[0], 11);

    int i2 = 12;
    px = i2;
    SEQAN_ASSERT_EQ(i1[0], 12);

    int * pi2 = i1 + 1;
    Proxy<IteratorProxy<int *> > px2(pi2);
    px = px2;
    SEQAN_ASSERT_EQ(i1[0], 20);
    SEQAN_ASSERT_EQ(px, 20);

//copy ctor
    Proxy<IteratorProxy<int *> > px3(px2);
    SEQAN_ASSERT_EQ(px3, 20);


//assign
    char s1[100] = "";
    char * it1 = s1;
    Proxy<IteratorProxy<char *> > px4(it1);

    assign(px4, 'X');
    SEQAN_ASSERT_EQ(px4, 'X');

    char c1 = 'a';
    assign(px4, c1);
    SEQAN_ASSERT_EQ(px4, 'a');

   SEQAN_ASSERT_FAIL("Move me to other tests in this header!");
}

SEQAN_DEFINE_TEST(test_basic_proxy_iterator_constructors)
{
    using namespace seqan;

    typedef Proxy<IteratorProxy<char *> > TProxy;

    // Construct from iterator.
    char data[10] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    char * it = &data[3];
    TProxy proxy(it);
    int x = proxy;
    SEQAN_ASSERT_EQ(3, x);

    // Copy constructor;
    TProxy proxy2(proxy);
    int y = proxy2;
    SEQAN_ASSERT_EQ(3, y);
}

SEQAN_DEFINE_TEST(test_basic_proxy_iterator_assign)
{
    using namespace seqan;

    typedef Proxy<IteratorProxy<int *> > TProxy;

    int data[3] = { 0, 1, 2 };
    int * it = &data[1];

    TProxy proxy(it);
    assign(proxy, 4);

    SEQAN_ASSERT_EQ(data[0], 0);
    SEQAN_ASSERT_EQ(data[1], 4);
    SEQAN_ASSERT_EQ(data[2], 2);
}

// TODO(holtgrew): Test for set() after writing it.
// TODO(holtgrew): Test for move() after writing it.

SEQAN_DEFINE_TEST(test_basic_proxy_iterator_getValue)
{
    using namespace seqan;

    typedef Proxy<IteratorProxy<int *> > TProxy;

    int data[3] = { 0, 1, 2 };
    int * it = &data[1];

    TProxy proxy(it);

    SEQAN_ASSERT_EQ(getValue(proxy), 1);
}

SEQAN_DEFINE_TEST(test_basic_proxy_iterator_comparators)
{
    using namespace seqan;

    typedef Proxy<IteratorProxy<int *> > TProxy;

    int data[4] = { 0, 1, 2, 1 };
    int * it = &data[1];
    int * it2 = &data[2];
    int * it3 = &data[3];

    TProxy proxy(it);
    TProxy proxy2(it2);
    TProxy proxy3(it3);

    SEQAN_ASSERT(proxy == 1);
    SEQAN_ASSERT(1 == proxy);
    SEQAN_ASSERT(proxy != 2);
    SEQAN_ASSERT(2 != proxy);
    SEQAN_ASSERT(proxy <= 1);
    SEQAN_ASSERT(1 <= proxy);
    SEQAN_ASSERT(proxy >= 1);
    SEQAN_ASSERT(2 >= proxy);
    SEQAN_ASSERT(proxy < 2);
    SEQAN_ASSERT(0 < proxy);
    SEQAN_ASSERT(proxy > 0);
    SEQAN_ASSERT(2 > proxy);

    SEQAN_ASSERT(proxy == proxy3);
    SEQAN_ASSERT(proxy != proxy2);
    SEQAN_ASSERT(proxy < proxy2);
    SEQAN_ASSERT(proxy <= proxy2);
    SEQAN_ASSERT(proxy2 > proxy);
    SEQAN_ASSERT(proxy2 >= proxy);
}

SEQAN_DEFINE_TEST(test_basic_proxy_iterator_stream_read)
{
    using namespace seqan;

    std::stringstream ss;
    ss << 33;

    typedef Proxy<IteratorProxy<int *> > TProxy;

    int x = 0;
    TProxy proxy(&x);

    ss >> proxy;

    int y = proxy;
    SEQAN_ASSERT_EQ(y, 33);
}

SEQAN_DEFINE_TEST(test_basic_proxy_iterator_stream_write)
{
    using namespace seqan;

    std::stringstream ss;
    int x = 33;

    typedef Proxy<IteratorProxy<int *> > TProxy;
    TProxy proxy(&x);

    ss << proxy;

    SEQAN_ASSERT(ss.str() == "33");
}

#endif  // #ifndef TESTS_BASIC_TEST_BASIC_PROXY_H_
