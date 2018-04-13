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
// Tests for the STL adaptions of seqan.
// ==========================================================================

// TODO(holtgrew): Split into a file for STL strings and one for STL lists?  Writing templatized tests do not make too much sense, I guess, becauses lists and strings are so dissimilar.

#ifndef TEST_SEQUENCE_TEST_SEQUENCE_STD_ADAPTIONS_H_
#define TEST_SEQUENCE_TEST_SEQUENCE_STD_ADAPTIONS_H_


// Tests the return types and existence of the metafunctions for STL vectors.
SEQAN_DEFINE_TEST(test_sequence_adaptions_metafunctions_std_vector)
{
    using namespace seqan;

    typedef int TElement;
    typedef std::vector<TElement> TVector;
    typedef TVector const TConstVector;

    // Test IsContiguous<>::VALUE
    {
        bool b = IsContiguous<TVector>::VALUE;
        SEQAN_ASSERT(b);
        b = IsContiguous<TConstVector>::VALUE;
        SEQAN_ASSERT(b);
    }
    // Test Value<>::VALUE
    {
        typedef Value<TVector>::Type TValue;
        bool b = IsSameType<TValue, TElement>::VALUE;
        SEQAN_ASSERT(b);
        typedef Value<TConstVector>::Type TConstValue;
        b = IsSameType<TConstValue, TElement>::VALUE;
        SEQAN_ASSERT(b);
    }
    // Test GetValue<>::VALUE
    {
        typedef GetValue<TVector>::Type TGetValue;
        bool b = IsSameType<TGetValue, TElement const &>::VALUE;
        SEQAN_ASSERT(b);
        typedef GetValue<TConstVector>::Type TConstGetValue;
        b = IsSameType<TConstGetValue, TElement const &>::VALUE;
        SEQAN_ASSERT(b);
    }
    // Test GetReference<>::VALUE
    {
        typedef Reference<TVector>::Type TReference;
        bool b = IsSameType<TReference, TElement &>::VALUE;
        SEQAN_ASSERT(b);
        typedef Reference<TConstVector>::Type TConstReference;
        b = IsSameType<TConstReference, TElement const &>::VALUE;
        SEQAN_ASSERT(b);
    }
    // Test Iterator<, Rooted>::VALUE
    {
        typedef Iterator<TVector, Rooted>::Type TIterator;
        typedef Iter<TVector, AdaptorIterator<Iter<TVector, StdIteratorAdaptor> > > TExpected;
        bool b = IsSameType<TIterator, TExpected>::VALUE;
        SEQAN_ASSERT(b);
        typedef Iterator<TConstVector, Rooted>::Type TConstIterator;
        typedef Iter<TConstVector, AdaptorIterator<Iter<TConstVector, StdIteratorAdaptor> > > TExpectedConst;
        b = IsSameType<TConstIterator, TExpectedConst>::VALUE;
        SEQAN_ASSERT(b);
    }
    // Test Iterator<, Standard>::VALUE
    {
        typedef Iterator<TVector, Standard>::Type TIterator;
        typedef Iter<TVector, StdIteratorAdaptor> TExpected;
        bool b = IsSameType<TIterator, TExpected>::VALUE;
        SEQAN_ASSERT(b);
        typedef Iterator<TConstVector, Standard>::Type TConstIterator;
        typedef Iter<TConstVector, StdIteratorAdaptor> TExpectedConst;
        b = IsSameType<TConstIterator, TExpectedConst>::VALUE;
        SEQAN_ASSERT(b);
    }
    // Test Position<>::VALUE
    {
        typedef Position<TVector>::Type TPosition;
        bool b = IsSameType<TPosition, TVector::size_type>::VALUE;
        SEQAN_ASSERT(b);
        typedef Position<TConstVector>::Type TConstPosition;
        b = IsSameType<TConstPosition, TVector::size_type>::VALUE;
        SEQAN_ASSERT(b);
    }
    // Test Size<>::VALUE
    {
        typedef Size<TVector>::Type TPosition;
        bool b = IsSameType<TPosition, TVector::size_type>::VALUE;
        SEQAN_ASSERT(b);
        typedef Size<TConstVector>::Type TConstPosition;
        b = IsSameType<TConstPosition, TVector::size_type>::VALUE;
        SEQAN_ASSERT(b);
    }
}


// Test iterators for STL strings.
SEQAN_DEFINE_TEST(test_sequence_adaptions_iterators_std_vector)
{
    using namespace seqan;

    // Test const iterator.
    {
        std::vector<int> const vec(2,100);
        //typedef Iterator<std::vector<int> const>::Type TIterator;

        std::vector<int> vecCopy;
        vecCopy.resize(2);
        copy(vec.begin(), vec.begin()+2, vecCopy.begin());

        SEQAN_ASSERT_EQ(vecCopy[0],vec[0]);
        SEQAN_ASSERT_EQ(vecCopy[1],vec[1]);

    }

    // Test non-const iterator.
    {
        std::vector<int> vec(2,100);
        //typedef Iterator<std::vector<int> >::Type TIterator;

        std::vector<int> vecCopy;
        vecCopy.resize(2);
        copy(vec.begin(), vec.begin()+2, vecCopy.begin());

        SEQAN_ASSERT_EQ(vecCopy[0],vec[0]);
        SEQAN_ASSERT_EQ(vecCopy[1],vec[1]);
    }
}


// Tests for the basic sequence functions for STL strings,
// e.g. value(), front(), back().
SEQAN_DEFINE_TEST(test_sequence_adaptions_sequence_interface_std_vector)
{
    using namespace seqan;

    std::vector<int> vec(2,100);

    // value(str, i), getValue(str, i)
    SEQAN_ASSERT_EQ(value(vec, 0), 100);
    SEQAN_ASSERT_EQ(value(vec, 1), 100);
     // front(), back()
    SEQAN_ASSERT_EQ(front(vec),100);
    SEQAN_ASSERT_EQ(back(vec),100);

    // length()
    SEQAN_ASSERT_EQ(length(vec), 2u);

    // TODO(holtgrew): Anything else missing? Probably...
}


// Tests for the memory allocation and reservation related functions
// for STL strings.
SEQAN_DEFINE_TEST(test_sequence_adaptions_sequence_memory_std_vector)
{
    using namespace seqan;

    // Test resize function -- resize down.
    {
        std::vector<int> vec(8,100);
        resize(vec, 5);

        SEQAN_ASSERT_EQ(vec.size(),5u);
    }

    // Test resize function -- resize up.
    {
         std::vector<int>  vec(2,100);
        resize(vec, 6);

          SEQAN_ASSERT_EQ(vec.size(),6u);
    }

    // Tests reserve function.
    {
        std::vector<int> vec;
        reserve(vec, 10);
        SEQAN_ASSERT_GEQ(capacity(vec), 10u);
    }
    {
        std::vector<int> vec;
        reserve(vec, 10, Generous());
        SEQAN_ASSERT_GEQ(capacity(vec), 10u);
    }
    // We cannot guarantee that the behaviour is supported by the STL
    // implementation with the Exact() tag.
    {
        std::vector<int> vec;
        reserve(vec, 10, Exact());
        SEQAN_ASSERT_EQ(capacity(vec), 10u);
    }
    // test first replace
    {
        std::vector<int> vec_target(5,100);
        std::vector<int> vec_source(3,10);
        std::vector<int> vec_source2(2,20);
        typename Position< std::vector<int> >::Type pos_begin = 3;
        typename Position< std::vector<int> >::Type pos_end = 4;

        // replace with insertion
        replace(vec_target,pos_begin,pos_end,vec_source);

        SEQAN_ASSERT_EQ(vec_target[2],100);
        SEQAN_ASSERT_EQ(vec_target[3],10);
        SEQAN_ASSERT_EQ(vec_target[5],10);
        SEQAN_ASSERT_EQ(vec_target[6],100);

        // replace with shrinking
        pos_begin = 2;
        pos_end = 5;
        replace(vec_target,pos_begin,pos_end,vec_source2);
        SEQAN_ASSERT_EQ(vec_target[1],100);
        SEQAN_ASSERT_EQ(vec_target[2],20);
        SEQAN_ASSERT_EQ(vec_target[3],20);
        SEQAN_ASSERT_EQ(vec_target[4],10);


    }
    // test replace with limits
    {
        std::vector<int> vec_target(8,100);
        std::vector<int> vec_source(6,10);
        typename Position< std::vector<int> >::Type pos_begin = 4;
        typename Position< std::vector<int> >::Type pos_end = 8;


        // replace with insertion
        typename Size< std::vector<int> >::Type limit = 5;

        replace(vec_target,pos_begin,pos_end,vec_source,limit);

        SEQAN_ASSERT_EQ(vec_target[3],100);
        SEQAN_ASSERT_EQ(vec_target[4],10);
        SEQAN_ASSERT_EQ(length(vec_target),9u);
    }
}


// **************************************************************************
// Tests the return types and existence of the metafunctions for STL strings.
SEQAN_DEFINE_TEST(test_sequence_adaptions_metafunctions_std_string)
{
    using namespace seqan;

    typedef int TElement;
    typedef std::basic_string<TElement> TString;
    typedef TString const TConstString;

    // Test IsContiguous<>::VALUE
    {
        bool b = IsContiguous<TString>::VALUE;
        SEQAN_ASSERT(b);
        b = IsContiguous<TConstString>::VALUE;
        SEQAN_ASSERT(b);
    }
    // Test Value<>::VALUE
    {
        typedef Value<TString>::Type TValue;
        bool b = IsSameType<TValue, TElement>::VALUE;
        SEQAN_ASSERT(b);
        typedef Value<TConstString>::Type TConstValue;
        b = IsSameType<TConstValue, TElement>::VALUE;
        SEQAN_ASSERT(b);
    }
    // Test GetValue<>::VALUE
    {
        typedef GetValue<TString>::Type TGetValue;
        bool b = IsSameType<TGetValue, TElement const &>::VALUE;
        SEQAN_ASSERT(b);
        typedef GetValue<TConstString>::Type TConstGetValue;
        b = IsSameType<TConstGetValue, TElement const &>::VALUE;
        SEQAN_ASSERT(b);
    }
    // Test Reference<>::VALUE
    {
        typedef Reference<TString>::Type TReference;
        bool b = IsSameType<TReference, TElement &>::VALUE;
        SEQAN_ASSERT(b);
        typedef Reference<TConstString>::Type TConstReference;
        b = IsSameType<TConstReference, TElement const &>::VALUE;
        SEQAN_ASSERT(b);
    }
    // Test Iterator<, Rooted>::VALUE
    {
        typedef Iterator<TString, Rooted>::Type TIterator;
        typedef Iter<TString, AdaptorIterator<Iterator<TString, Standard>::Type> > TExpected;
        bool b = IsSameType<TIterator, TExpected>::VALUE;
        SEQAN_ASSERT(b);
        typedef Iterator<TConstString, Rooted>::Type TConstIterator;
        typedef Iter<TConstString, AdaptorIterator<Iterator<TConstString, Standard>::Type> > TExpectedConst;
        b = IsSameType<TConstIterator, TExpectedConst>::VALUE;
        SEQAN_ASSERT(b);
    }
    // Test Iterator<, Standard>::VALUE
    {
        typedef Iterator<TString, Standard>::Type TIterator;
        typedef Iter<TString, StdIteratorAdaptor> TExpected;
        bool b = IsSameType<TIterator, TExpected>::VALUE;
        SEQAN_ASSERT(b);
        typedef Iterator<TConstString, Standard>::Type TConstIterator;
        typedef Iter<TString const, StdIteratorAdaptor> TExpectedConst;
        b = IsSameType<TConstIterator, TExpectedConst>::VALUE;
        SEQAN_ASSERT(b);
    }
    // Test Position<>::VALUE
    {
        typedef Position<TString>::Type TPosition;
        bool b = IsSameType<TPosition, TString::size_type>::VALUE;
        SEQAN_ASSERT(b);
        typedef Position<TConstString>::Type TConstPosition;
        b = IsSameType<TConstPosition, TString::size_type>::VALUE;
        SEQAN_ASSERT(b);
    }
    // Test Size<>::VALUE
    {
        typedef Size<TString>::Type TPosition;
        bool b = IsSameType<TPosition, TString::size_type>::VALUE;
        SEQAN_ASSERT(b);
        typedef Size<TConstString>::Type TConstPosition;
        b = IsSameType<TConstPosition, TString::size_type>::VALUE;
        SEQAN_ASSERT(b);
    }
}


// Test iterators for STL strings.
SEQAN_DEFINE_TEST(test_sequence_adaptions_iterators_std_string)
{
    using namespace seqan;

    // Test const iterator.
    {
        std::string const str = "Unimportant contents.";
        typedef Iterator<std::string const>::Type TIterator;

        std::string strCopy;
        for (TIterator it = begin(str, Standard()); it != end(str, Standard()); ++it)
            appendValue(strCopy, value(it));

        SEQAN_ASSERT_EQ(str, strCopy);
    }

    // Test non-const iterator.
    {
        std::string str = "Unimportant contents.";
        typedef Iterator<std::string>::Type TIterator;

        std::string strCopy;
        for (TIterator it = begin(str, Standard()); it != end(str, Standard()); ++it)
            appendValue(strCopy, value(it));

        SEQAN_ASSERT_EQ(str, strCopy);
    }
}


// Tests for the basic sequence functions for STL strings,
// e.g. value(), front(), back().
SEQAN_DEFINE_TEST(test_sequence_adaptions_sequence_interface_std_string)
{
    using namespace seqan;

    std::string str = "Hello World!";

    // value(str, i), getValue(str, i)
    SEQAN_ASSERT_EQ(value(str, 0), 'H');
    SEQAN_ASSERT_EQ(value(str, 4), 'o');
    SEQAN_ASSERT_EQ(getValue(str, 0), 'H');
    SEQAN_ASSERT_EQ(getValue(str, 4), 'o');

    // front(), back()
    SEQAN_ASSERT_EQ(front(str), 'H');
    SEQAN_ASSERT_EQ(back(str), '!');

    // length()
    SEQAN_ASSERT_EQ(length(str), 12u);

    // TODO(holtgrew): Anything else missing? Probably...
}


// Tests for the memory allocation and reservation related functions
// for STL strings.
SEQAN_DEFINE_TEST(test_sequence_adaptions_sequence_memory_std_string)
{
    using namespace seqan;

    // Test resize function -- resize down.
    {
        std::string str = "Hello world!";
        resize(str, 5);
        SEQAN_ASSERT_EQ(str, "Hello");
    }

    // Test resize function -- resize up.
    {
        std::string str = "12345";
        resize(str, 6);
        // The following gives an assertion in positional setValue() if not resized properly.
        str[5] = '6';
        SEQAN_ASSERT_EQ(str, "123456");
    }

    // Tests reserve function.
    {
        std::string str;
        reserve(str, 10);
        SEQAN_ASSERT_GEQ(capacity(str), 10u);
    }
    {
        std::string str;
        reserve(str, 10, Generous());
        SEQAN_ASSERT_GEQ(capacity(str), 10u);
    }
    // We cannot guarantee that the behaviour is supported by the STL
    // implementation with the Exact() tag.
//    {
//        std::string str;
//        reserve(str, 10, Exact());
//        SEQAN_ASSERT_EQ(capacity(str), 10u);
//    }
}


// *******************************************************************
// Tests the return types and existence of the metafunctions for STL lists.
SEQAN_DEFINE_TEST(test_sequence_adaptions_metafunctions_std_list)
{
    using namespace seqan;

    typedef int TElement;
    typedef std::list<TElement> TList;
    typedef TList const TConstList;

    // Test IsContiguous<>::VALUE
    {
        bool b = IsContiguous<TList>::VALUE;
        SEQAN_ASSERT_NOT(b);
        b = IsContiguous<TConstList>::VALUE;
        SEQAN_ASSERT_NOT(b);
    }
    // Test Value<>::VALUE
    {
        typedef Value<TList>::Type TValue;
        bool b = IsSameType<TValue, TElement>::VALUE;
        SEQAN_ASSERT(b);
        typedef Value<TConstList>::Type TConstValue;
        b = IsSameType<TConstValue, TElement>::VALUE;
        SEQAN_ASSERT(b);
    }
    // Test GetValue<>::VALUE
    {
        typedef GetValue<TList>::Type TGetValue;
        bool b = IsSameType<TGetValue, TElement const &>::VALUE;
        SEQAN_ASSERT(b);
        typedef GetValue<TConstList>::Type TConstGetValue;
        b = IsSameType<TConstGetValue, TElement const &>::VALUE;
        SEQAN_ASSERT(b);
    }
    // Test Reference<>::VALUE
    {
        typedef Reference<TList>::Type TReference;
        bool b = IsSameType<TReference, TElement &>::VALUE;
        SEQAN_ASSERT(b);
        typedef Reference<TConstList>::Type TConstReference;
        b = IsSameType<TConstReference, TElement const &>::VALUE;
        SEQAN_ASSERT(b);
    }
    // Test Iterator<, Rooted>::VALUE
    {
        typedef Iterator<TList, Rooted>::Type TIterator;
        typedef Iter<TList, AdaptorIterator<Iter<TList, StdIteratorAdaptor> > > TExpected;
        bool b = IsSameType<TIterator, TExpected>::VALUE;
        SEQAN_ASSERT(b);
        typedef Iterator<TConstList, Rooted>::Type TConstIterator;
        typedef Iter<TConstList, AdaptorIterator<Iter<TConstList, StdIteratorAdaptor> > > TExpectedConst;
        b = IsSameType<TConstIterator, TExpectedConst>::VALUE;
        SEQAN_ASSERT(b);
    }
    // Test Iterator<, Standard>::VALUE
    {
        typedef Iterator<TList, Standard>::Type TIterator;
        typedef Iter<TList, StdIteratorAdaptor> TExpected;
        bool b = IsSameType<TIterator, TExpected>::VALUE;
        SEQAN_ASSERT(b);
        typedef Iterator<TConstList, Standard>::Type TConstIterator;
        typedef Iter<TConstList, StdIteratorAdaptor> TExpectedConst;
        b = IsSameType<TConstIterator, TExpectedConst>::VALUE;
        SEQAN_ASSERT(b);
    }
    // Test Position<>::VALUE
    {
        typedef Position<TList>::Type TPosition;
        bool b = IsSameType<TPosition, TList::size_type>::VALUE;
        SEQAN_ASSERT(b);
        typedef Position<TConstList>::Type TConstPosition;
        b = IsSameType<TConstPosition, TList::size_type>::VALUE;
        SEQAN_ASSERT(b);
    }
    // Test Size<>::VALUE
    {
        typedef Size<TList>::Type TPosition;
        bool b = IsSameType<TPosition, TList::size_type>::VALUE;
        SEQAN_ASSERT(b);
        typedef Size<TConstList>::Type TConstPosition;
        b = IsSameType<TConstPosition, TList::size_type>::VALUE;
        SEQAN_ASSERT(b);
    }
}


// Test iterators for STL lists.
SEQAN_DEFINE_TEST(test_sequence_adaptions_iterators_std_list)
{
    using namespace seqan;

    typedef int TElement;

    // Test Standard, non-const iterators.
    {
        typedef std::list<TElement> TList;
        typedef Iterator<TList, Standard>::Type TIterator;

        TList list;
        appendValue(list, 1);
        appendValue(list, 2);
        appendValue(list, 3);

        // The command sequence in the following is a bit arbitrary
        // but should robustly test that the iterators work correctly.
        TIterator it = seqan::begin(list);
        SEQAN_ASSERT_EQ(1, *it);
        ++it;
        SEQAN_ASSERT_EQ(2, *it);
        --it;
        SEQAN_ASSERT_EQ(1, *it);
        TIterator itEnd = seqan::end(list);
        SEQAN_ASSERT_NOT(it == itEnd);
        ++it;
        ++it;
        ++it;
        SEQAN_ASSERT(it == itEnd);

        // The following does not apply to const iterators.
        it = seqan::begin(list);
        *seqan::begin(list) = 4;
        SEQAN_ASSERT_EQ(4, *it);
    }
    // Test Standard, const iterators.
    {
        typedef std::list<TElement> TList;
        typedef Iterator<TList const, Standard>::Type TIterator;

        TList mutableList;
        appendValue(mutableList, 1);
        appendValue(mutableList, 2);
        appendValue(mutableList, 3);

        TList const & list = mutableList;

        // The command sequence in the following is a bit arbitrary
        // but should robustly test that the iterators work correctly.
        TIterator it = seqan::begin(list);
        SEQAN_ASSERT_EQ(1, *it);
        ++it;
        SEQAN_ASSERT_EQ(2, *it);
        --it;
        SEQAN_ASSERT_EQ(1, *it);
        TIterator itEnd = seqan::end(list);
        SEQAN_ASSERT_NOT(it == itEnd);
        ++it;
        ++it;
        ++it;
        SEQAN_ASSERT(it == itEnd);
    }
    // Test Rooted, non-const iterators.
    {
        typedef std::list<TElement> TList;
        typedef Iterator<TList, Rooted>::Type TIterator;

        TList list;
        appendValue(list, 1);
        appendValue(list, 2);
        appendValue(list, 3);

        // The command sequence in the following is a bit arbitrary
        // but should robustly test that the iterators work correctly.
        TIterator it = seqan::begin(list);
        SEQAN_ASSERT_EQ(1, *it);
        ++it;
        SEQAN_ASSERT_EQ(2, *it);
        --it;
        SEQAN_ASSERT_EQ(1, *it);
        TIterator itEnd = seqan::end(list);
        SEQAN_ASSERT_NOT(it == itEnd);
        ++it;
        ++it;
        ++it;
        SEQAN_ASSERT(it == itEnd);

        // The following does not apply to const iterators.
        it = seqan::begin(list);
        *seqan::begin(list) = 4;
        SEQAN_ASSERT_EQ(4, *it);
    }
    // Test Rooted, const iterators.
    {
        typedef std::list<TElement> TList;
        typedef Iterator<TList const, Rooted>::Type TIterator;

        TList mutableList;
        appendValue(mutableList, 1);
        appendValue(mutableList, 2);
        appendValue(mutableList, 3);

        TList const & list = mutableList;

        // The command sequence in the following is a bit arbitrary
        // but should robustly test that the iterators work correctly.
        TIterator it = seqan::begin(list);
        SEQAN_ASSERT_EQ(1, *it);
        ++it;
        SEQAN_ASSERT_EQ(2, *it);
        --it;
        SEQAN_ASSERT_EQ(1, *it);
        TIterator itEnd = seqan::end(list);
        SEQAN_ASSERT_NOT(it == itEnd);
        ++it;
        ++it;
        ++it;
        SEQAN_ASSERT(it == itEnd);
    }
}


// Test the basic sequence interface implemented for STL list, e.g. front() and back().
SEQAN_DEFINE_TEST(test_sequence_adaptions_sequence_interface_std_list)
{
    using namespace seqan;

    typedef int TElement;

    // Test with non-const container.
    {
        typedef std::list<TElement> TList;
        typedef Iterator<TList>::Type TIterator;

        // Prepare list...
        TList list;
        appendValue(list, 1);
        appendValue(list, 2);
        appendValue(list, 3);

        // Test reading front and back.
        SEQAN_ASSERT_EQ(1, front(list));
        SEQAN_ASSERT_EQ(3, back(list));

        // Test assigning to front and back.
        front(list) = -1;
        back(list) = -3;

        TIterator it = seqan::begin(list);
        SEQAN_ASSERT_EQ(-1, *it);
        ++it;
        SEQAN_ASSERT_EQ(2, *it);
        ++it;
        SEQAN_ASSERT_EQ(-3, *it);

        // Test appending and prepending values.
        prependValue(list, 42);
        appendValue(list, 43);
        SEQAN_ASSERT_EQ(42, front(list));
        SEQAN_ASSERT_EQ(43, back(list));

        // Test length().
        SEQAN_ASSERT_EQ(5u, length(list));

        // Test clear().
        {
            TList listCopy(list);
            SEQAN_ASSERT_EQ(5u, length(listCopy));
            clear(listCopy);
            SEQAN_ASSERT_EQ(0u, length(listCopy));
        }
    }
    // Test with const container.
    {
        typedef std::list<TElement> TList;
        //typedef Iterator<TList>::Type TIterator;

        // Prepare list...
        TList mutableList;
        appendValue(mutableList, 1);
        appendValue(mutableList, 2);
        appendValue(mutableList, 3);

        TList const & list = mutableList;

        // Test reading front and back.
        SEQAN_ASSERT_EQ(1, front(list));
        SEQAN_ASSERT_EQ(3, back(list));

        // Test length().
        SEQAN_ASSERT_EQ(3u, length(list));

        // Test reserve().
        reserve(mutableList, 4);
        SEQAN_ASSERT_EQ(3u, length(list));
        SEQAN_ASSERT_EQ(3u, length(mutableList));
        SEQAN_ASSERT_EQ(3u, capacity(list));
        SEQAN_ASSERT_EQ(3u, capacity(mutableList));
    }
}

// Tests the return types and existence of the metafunctions for STL arrays.
SEQAN_DEFINE_TEST(test_sequence_adaptions_metafunctions_std_array)
{
    using namespace seqan;

    typedef int TElement;
    typedef std::array<TElement, 1> TArray;
    typedef TArray const TConstArray;

    // Test IsContiguous<>::VALUE
    {
        bool b = IsContiguous<TArray>::VALUE;
        SEQAN_ASSERT(b);
        b = IsContiguous<TConstArray>::VALUE;
        SEQAN_ASSERT(b);
    }
    // Test Value<>::VALUE
    {
        typedef Value<TArray>::Type TValue;
        bool b = IsSameType<TValue, TElement>::VALUE;
        SEQAN_ASSERT(b);
        typedef Value<TConstArray>::Type TConstValue;
        b = IsSameType<TConstValue, TElement>::VALUE;
        SEQAN_ASSERT(b);
    }
    // Test GetValue<>::VALUE
    {
        typedef GetValue<TArray>::Type TGetValue;
        bool b = IsSameType<TGetValue, TElement const &>::VALUE;
        SEQAN_ASSERT(b);
        typedef GetValue<TConstArray>::Type TConstGetValue;
        b = IsSameType<TConstGetValue, TElement const &>::VALUE;
        SEQAN_ASSERT(b);
    }
    // Test GetReference<>::VALUE
    {
        typedef Reference<TArray>::Type TReference;
        bool b = IsSameType<TReference, TElement &>::VALUE;
        SEQAN_ASSERT(b);
        typedef Reference<TConstArray>::Type TConstReference;
        b = IsSameType<TConstReference, TElement const &>::VALUE;
        SEQAN_ASSERT(b);
    }
    // Test Iterator<, Rooted>::VALUE
    {
        typedef Iterator<TArray, Rooted>::Type TIterator;
        typedef Iter<TArray, AdaptorIterator<Iter<TArray, StdIteratorAdaptor> > > TExpected;
        bool b = IsSameType<TIterator, TExpected>::VALUE;
        SEQAN_ASSERT(b);
        typedef Iterator<TConstArray, Rooted>::Type TConstIterator;
        typedef Iter<TConstArray, AdaptorIterator<Iter<TConstArray, StdIteratorAdaptor> > > TExpectedConst;
        b = IsSameType<TConstIterator, TExpectedConst>::VALUE;
        SEQAN_ASSERT(b);
    }
    // Test Iterator<, Standard>::VALUE
    {
        typedef Iterator<TArray, Standard>::Type TIterator;
        typedef Iter<TArray, StdIteratorAdaptor> TExpected;
        bool b = IsSameType<TIterator, TExpected>::VALUE;
        SEQAN_ASSERT(b);
        typedef Iterator<TConstArray, Standard>::Type TConstIterator;
        typedef Iter<TConstArray, StdIteratorAdaptor> TExpectedConst;
        b = IsSameType<TConstIterator, TExpectedConst>::VALUE;
        SEQAN_ASSERT(b);
    }
    // Test Position<>::VALUE
    {
        typedef Position<TArray>::Type TPosition;
        bool b = IsSameType<TPosition, TArray::size_type>::VALUE;
        SEQAN_ASSERT(b);
        typedef Position<TConstArray>::Type TConstPosition;
        b = IsSameType<TConstPosition, TArray::size_type>::VALUE;
        SEQAN_ASSERT(b);
    }
    // Test Size<>::VALUE
    {
        typedef Size<TArray>::Type TPosition;
        bool b = IsSameType<TPosition, TArray::size_type>::VALUE;
        SEQAN_ASSERT(b);
        typedef Size<TConstArray>::Type TConstPosition;
        b = IsSameType<TConstPosition, TArray::size_type>::VALUE;
        SEQAN_ASSERT(b);
    }
}


// Test iterators for STL strings.
SEQAN_DEFINE_TEST(test_sequence_adaptions_iterators_std_array)
{
    using namespace seqan;

    // Test const iterator.
    {
        std::array<int, 3> const vec = { {100, 101, 102} };
        //typedef Iterator<std::array<int> const>::Type TIterator;

        std::array<int, 2>  vecCopy;
        std::copy(vec.begin(), vec.begin()+2, vecCopy.begin());

        SEQAN_ASSERT_EQ(vecCopy[0],vec[0]);
        SEQAN_ASSERT_EQ(vecCopy[1],vec[1]);

    }

    // Test non-const iterator.
    {
        std::array<int, 3> vec = { {100, 101, 102} };
        //typedef Iterator<std::array<int> >::Type TIterator;

        std::array<int, 2> vecCopy;
        std::copy(vec.begin(), vec.begin()+2, vecCopy.begin());

        SEQAN_ASSERT_EQ(vecCopy[0],vec[0]);
        SEQAN_ASSERT_EQ(vecCopy[1],vec[1]);
    }
}


// Tests for the basic sequence functions for STL strings,
// e.g. value(), front(), back().
SEQAN_DEFINE_TEST(test_sequence_adaptions_sequence_interface_std_array)
{
    using namespace seqan;

    std::array<int, 2> vec = { {100, 101} };

    // value(str, i), getValue(str, i)
    SEQAN_ASSERT_EQ(value(vec, 0), 100);
    SEQAN_ASSERT_EQ(value(vec, 1), 101);
    // front(), back()
    SEQAN_ASSERT_EQ(front(vec),100);
    SEQAN_ASSERT_EQ(back(vec),101);

    // length()
    SEQAN_ASSERT_EQ(length(vec), 2u);

    // TODO(holtgrew): Anything else missing? Probably...
}

#endif  // TEST_SEQUENCE_TEST_SEQUENCE_STD_ADAPTIONS_H_
