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

#ifndef TESTS_BASIC_TEST_BASIC_ITERATOR_ZIP_H_
#define TESTS_BASIC_TEST_BASIC_ITERATOR_ZIP_H_

#include <tuple>
#include <vector>
#include <string>

#include <seqan/basic.h>
#include <seqan/sequence.h>

// --------------------------------------------------------------------------
// Tests for Zip Iterator
// --------------------------------------------------------------------------

SEQAN_DEFINE_TEST(test_basic_iterator_zip_metafunctions)
{
    using namespace seqan;

    // Iterator.
    {
        using TFirst = typename Iterator<std::vector<int> >::Type;
        using TSecond = typename Iterator<std::string>::Type;
        using TIterator = Iter<std::tuple<TFirst, TSecond>, ZipIterator>;

        bool b = IsSameType<typename Difference<TIterator>::Type, typename Difference<TFirst>::Type>::VALUE;
        SEQAN_ASSERT(b);
        b = IsSameType<typename Position<TIterator>::Type, typename Position<TFirst>::Type>::VALUE;
        SEQAN_ASSERT(b);
        b = IsSameType<typename Size<TIterator>::Type, typename Size<TFirst>::Type>::VALUE;
        SEQAN_ASSERT(b);

        b = IsSameType<typename Value<TIterator>::Type, std::tuple<int, char> >::VALUE;
        SEQAN_ASSERT(b);
        b = IsSameType<typename GetValue<TIterator>::Type, std::tuple<const int&, const char&> >::VALUE;
        SEQAN_ASSERT(b);
        b = IsSameType<typename Reference<TIterator>::Type, std::tuple<int&, char&> >::VALUE;
        SEQAN_ASSERT(b);

        b = IsSameType<typename Container<TIterator>::Type, std::tuple<std::vector<int>, std::string> >::VALUE;
        SEQAN_ASSERT(b);
    }
    // Const-Iterators.
    {
        using TFirst = typename Iterator<std::vector<int> >::Type;
        using TSecond = typename Iterator<std::string>::Type;
        using TIterator = Iter<std::tuple<TFirst const, TSecond const>, ZipIterator>;

        bool b = IsSameType<typename Difference<TIterator>::Type, typename Difference<TFirst const>::Type>::VALUE;
        SEQAN_ASSERT(b);
        b = IsSameType<typename Position<TIterator>::Type, typename Position<TFirst const>::Type>::VALUE;
        SEQAN_ASSERT(b);
        b = IsSameType<typename Size<TIterator>::Type, typename Size<TFirst const>::Type>::VALUE;
        SEQAN_ASSERT(b);

        // TODO(holtgrew): This is inconsistent
        b = IsSameType<typename Value<TIterator>::Type, std::tuple<int, char> >::VALUE;
        SEQAN_ASSERT(b);
        b = IsSameType<typename GetValue<TIterator>::Type, std::tuple<int const &, char const &> >::VALUE;
        SEQAN_ASSERT(b);
        b = IsSameType<typename Reference<TIterator>::Type, std::tuple<int &, char &> >::VALUE;
        SEQAN_ASSERT(b);

        b = IsSameType<typename Container<TIterator>::Type, std::tuple<std::vector<int>, std::string> >::VALUE;
        SEQAN_ASSERT(b);
    }

    // Non-const-Iterators with const containers.
    {
        using TFirst = typename Iterator<std::vector<int> const>::Type;
        using TSecond = typename Iterator<std::string const>::Type;
        using TIterator = Iter<std::tuple<TFirst, TSecond>, ZipIterator>;

        bool b = IsSameType<typename Difference<TIterator>::Type, typename Difference<TFirst>::Type>::VALUE;
        SEQAN_ASSERT(b);
        b = IsSameType<typename Position<TIterator>::Type, typename Position<TFirst>::Type>::VALUE;
        SEQAN_ASSERT(b);
        b = IsSameType<typename Size<TIterator>::Type, typename Size<TFirst>::Type>::VALUE;
        SEQAN_ASSERT(b);

        // TODO(holtgrew): This is inconsistent
        b = IsSameType<typename Value<TIterator>::Type, std::tuple<int, char> >::VALUE;
        SEQAN_ASSERT(b);
        b = IsSameType<typename GetValue<TIterator>::Type, std::tuple<int const &, char const &> >::VALUE;
        SEQAN_ASSERT(b);
        b = IsSameType<typename Reference<TIterator>::Type, std::tuple<int const &, char const &> >::VALUE;
        SEQAN_ASSERT(b);

        b = IsSameType<typename Container<TIterator>::Type, std::tuple<std::vector<int> const, std::string const> >::VALUE;
        SEQAN_ASSERT(b);
    }
}

SEQAN_DEFINE_TEST(test_basic_iterator_zip_constructors)
{
    using namespace seqan;

    using TFirst = typename Iterator<std::vector<int>, Standard>::Type;
    using TSecond = typename Iterator<std::string, Standard>::Type;

    using TIterator = Iter<std::tuple<TFirst, TSecond>, ZipIterator>;

    std::vector<int> vec{1, 2, 3};
    std::string str("abc");

    // Construct from list of iterators.
    TIterator it(begin(vec, Standard()), begin(str, Standard()));

    SEQAN_ASSERT(std::get<0>(it.dataIter) == begin(vec, Standard()));
    SEQAN_ASSERT(std::get<1>(it.dataIter) == begin(str, Standard()));
}

SEQAN_DEFINE_TEST(test_basic_iterator_zip_make_zip_iterator)
{
    using namespace seqan;

    using TFirst = typename Iterator<std::vector<int>, Standard>::Type;
    using TSecond = typename Iterator<std::string, Standard>::Type;

    using TIterator = Iter<std::tuple<TFirst, TSecond>, ZipIterator>;

    std::vector<int> vec{1, 2, 3};
    std::string str("abc");

    TIterator it = makeZipIterator(begin(vec, Standard()), begin(str, Standard()));

    SEQAN_ASSERT_EQ(*std::get<0>(it.dataIter), 1);
    SEQAN_ASSERT_EQ(*std::get<1>(it.dataIter), 'a');
}

SEQAN_DEFINE_TEST(test_basic_iterator_zip_transport)
{
    using namespace seqan;

    using TFirst = typename Iterator<std::vector<int>, Standard>::Type;
    using TSecond = typename Iterator<std::string, Standard>::Type;

    using TIterator = Iter<std::tuple<TFirst, TSecond>, ZipIterator>;

    std::vector<int> vec{1, 2, 3};
    std::string str("abc");

    TIterator it(begin(vec, Standard()), begin(str, Standard()));

    // assign()
    {
        TIterator it2;
        assign(it2, it);
        SEQAN_ASSERT(it == it2);
    }
    // set()
    {
        TIterator it2;
        seqan::set(it2, it);
        SEQAN_ASSERT(it == it2);
    }
    // move()
    {
        TIterator it2;
        move(it2, it);
        SEQAN_ASSERT(it == it2);
    }
}

SEQAN_DEFINE_TEST(test_basic_iterator_zip_transport_value)
{
    using namespace seqan;

    // assignValue()
    {
        std::vector<CDStruct> vec;
        vec.resize(3);
        resetCDStructStatics();
        std::string str("abc");

        auto it = makeZipIterator(begin(vec, Standard()), begin(str, Standard()));
        assignValue(it, std::forward_as_tuple(vec[1], str[1]));

        SEQAN_ASSERT_EQ(std::get<0>(it.dataIter)->copiedFrom, -1);
        SEQAN_ASSERT_EQ(std::get<0>(it.dataIter)->movedFrom, -1);
        SEQAN_ASSERT_EQ(std::get<0>(it.dataIter)->setFrom, -1);
        SEQAN_ASSERT_EQ(std::get<0>(it.dataIter)->assignedFrom, vec[1].id);

        SEQAN_ASSERT_EQ(CDStruct::lastOther, &vec[1]);
        SEQAN_ASSERT_EQ(CDStruct::defaultConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::copyConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::moveConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::moves, 0);
        SEQAN_ASSERT_EQ(CDStruct::destructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::assignments, 1);
    }
    // moveValue()
    {
        std::vector<CDStruct> vec;
        vec.resize(3);
        resetCDStructStatics();
        std::string str("abc");

        auto it = makeZipIterator(begin(vec, Standard()), begin(str, Standard()));
        moveValue(it, std::forward_as_tuple(vec[1], str[1]));

        SEQAN_ASSERT_EQ(std::get<0>(it.dataIter)->copiedFrom, -1);
        SEQAN_ASSERT_EQ(std::get<0>(it.dataIter)->movedFrom, vec[1].id);
        SEQAN_ASSERT_EQ(std::get<0>(it.dataIter)->setFrom, -1);
        SEQAN_ASSERT_EQ(std::get<0>(it.dataIter)->assignedFrom, -1);

        SEQAN_ASSERT_EQ(CDStruct::lastOther, &vec[1]);
        SEQAN_ASSERT_EQ(CDStruct::defaultConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::copyConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::moveConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::moves, 1);
        SEQAN_ASSERT_EQ(CDStruct::destructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::assignments, 0);
        SEQAN_ASSERT_EQ(CDStruct::sets, 0);
    }
    // There is no setValue() with iterator semantics.
}

SEQAN_DEFINE_TEST(test_basic_iterator_zip_movement)
{
    using namespace seqan;

    std::vector<int> vec{0,1,2,3,4,5,6,7,8,9};
    std::string str = "abcdefghij";

    // goNext/operator++
    {
        auto it = makeZipIterator(begin(vec, Standard()) + 4, begin(str, Standard()) + 4);
        goNext(it);
        SEQAN_ASSERT_EQ(std::get<0>(*it), vec[5]);
        SEQAN_ASSERT_EQ(std::get<1>(*it), str[5]);

        it = makeZipIterator(begin(vec, Standard()) + 4, begin(str, Standard()) + 4);
        it++;
        SEQAN_ASSERT_EQ(std::get<0>(*it), vec[5]);
        SEQAN_ASSERT_EQ(std::get<1>(*it), str[5]);

        it = makeZipIterator(begin(vec, Standard()) + 4, begin(str, Standard()) + 4);
        ++it;
        SEQAN_ASSERT_EQ(std::get<0>(*it), vec[5]);
        SEQAN_ASSERT_EQ(std::get<1>(*it), str[5]);
    }
    // goPrevious/operator--
    {
        auto it = makeZipIterator(begin(vec, Standard()) + 4, begin(str, Standard()) + 4);
        goPrevious(it);
        SEQAN_ASSERT_EQ(std::get<0>(*it), vec[3]);
        SEQAN_ASSERT_EQ(std::get<1>(*it), str[3]);

        it = makeZipIterator(begin(vec, Standard()) + 4, begin(str, Standard()) + 4);
        it--;
        SEQAN_ASSERT_EQ(std::get<0>(*it), vec[3]);
        SEQAN_ASSERT_EQ(std::get<1>(*it), str[3]);

        it = makeZipIterator(begin(vec, Standard()) + 4, begin(str, Standard()) + 4);
        --it;
        SEQAN_ASSERT_EQ(std::get<0>(*it), vec[3]);
        SEQAN_ASSERT_EQ(std::get<1>(*it), str[3]);
    }
    // goFurther/operator+=/operator-=
    {
        auto it = makeZipIterator(begin(vec, Standard()) + 4, begin(str, Standard()) + 4);
        goFurther(it, 2);
        SEQAN_ASSERT_EQ(std::get<0>(*it), vec[6]);
        SEQAN_ASSERT_EQ(std::get<1>(*it), str[6]);

        it = makeZipIterator(begin(vec, Standard()) + 4, begin(str, Standard()) + 4);
        goFurther(it, -2);
        SEQAN_ASSERT_EQ(std::get<0>(*it), vec[2]);
        SEQAN_ASSERT_EQ(std::get<1>(*it), str[2]);
    }
}

SEQAN_DEFINE_TEST(test_basic_iterator_zip_arithmetics)
{
    using namespace seqan;

    std::vector<int> vec{0,1,2,3,4,5,6,7,8,9};
    std::string str = "abcdefghij";

    auto it = makeZipIterator(begin(vec, Standard()) + 4, begin(str, Standard()) + 4);
    auto it2 = makeZipIterator(begin(vec, Standard()) + 6, begin(str, Standard()) + 6);
    SEQAN_ASSERT_EQ(&std::get<0>(*it2), &vec[6]);
    SEQAN_ASSERT_EQ(&std::get<1>(*it2), &str[6]);
    
    it2 = it - 2;
    SEQAN_ASSERT_EQ(&std::get<0>(*it2), &vec[2]);
    SEQAN_ASSERT_EQ(&std::get<1>(*it2), &str[2]);
    
    it2 = it + 2;
    SEQAN_ASSERT_EQ(it2 - it, 2);
}

#endif  // #ifndef TESTS_BASIC_TEST_BASIC_ITERATOR_ZIP_H_
