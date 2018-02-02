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
// Tests for iterators on the journaled string.
// ==========================================================================

#ifndef TEST_SEQUENCE_JOURNALED_TEST_SEQUENCE_JOURNALED_ITERATOR_H_
#define TEST_SEQUENCE_JOURNALED_TEST_SEQUENCE_JOURNALED_ITERATOR_H_

#include <cstdlib>
#include <sstream>
#include <string>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/sequence_journaled.h>

using namespace seqan;

// Test difference operator: operator-()
template <typename THostSpec, typename TStringJournalSpec, typename TBuffSpec>
void testJournaledStringIteratorDifference(THostSpec const &, TStringJournalSpec const &, TBuffSpec const &)
{
    typedef String<char, Journaled<THostSpec, TStringJournalSpec, TBuffSpec> > TJournaledString;

    CharString charStr = "test";

    // Case 1: Difference between begin and end.
    {
        TJournaledString journaledString;
        setHost(journaledString, charStr);

        typedef typename Iterator<TJournaledString, Standard>::Type TIterator;
        TIterator itBegin = begin(journaledString, Standard());
        TIterator itEnd = end(journaledString, Standard());
        SEQAN_ASSERT_EQ(4, itEnd - itBegin);
        SEQAN_ASSERT_EQ(-4, itBegin - itEnd);
    }

    // Case 2: Difference between non-begin and end.
    {
        TJournaledString journaledString;
        setHost(journaledString, charStr);

        typedef typename Iterator<TJournaledString, Standard>::Type TIterator;
        TIterator it = begin(journaledString, Standard());
        ++it; ++it;
        TIterator itEnd = end(journaledString, Standard());
        SEQAN_ASSERT_EQ(2, itEnd - it);
        SEQAN_ASSERT_EQ(-2, it - itEnd);
    }

    // Case 3: Slightly more complicated, between begin and end.
    {
        TJournaledString journaledString;
        setHost(journaledString, charStr);
        insert(journaledString, 2, CharString("XX"));

        typedef typename Iterator<TJournaledString, Standard>::Type TIterator;
        TIterator itBegin = begin(journaledString, Standard());
        TIterator itEnd = end(journaledString, Standard());
        SEQAN_ASSERT_EQ(6, itEnd - itBegin);
        SEQAN_ASSERT_EQ(-6, itBegin - itEnd);
    }

    // Case 4: Slightly more complicated, between non-begin and end.
    {
        TJournaledString journaledString;
        setHost(journaledString, charStr);
        insert(journaledString, 2, CharString("XX"));

        typedef typename Iterator<TJournaledString, Standard>::Type TIterator;
        TIterator it = begin(journaledString, Standard());
        ++it; ++it;
        TIterator itEnd = end(journaledString, Standard());
        SEQAN_ASSERT_EQ(4, itEnd - it);
        SEQAN_ASSERT_EQ(-4, it - itEnd);
    }
}

// Test sum operator: operator+()
template <typename THostSpec, typename TStringJournalSpec, typename TBuffSpec>
void testJournaledStringIteratorSum(THostSpec const &, TStringJournalSpec const &, TBuffSpec const &)
{
    typedef String<char, Journaled<THostSpec, TStringJournalSpec, TBuffSpec> > TJournaledString;
    typedef typename Iterator<TJournaledString, Standard>::Type TIterator;

    CharString charStr = "test";
    TJournaledString journaledString;

    setHost(journaledString, charStr);

    // Case 1: Go from beginning to end.
    {
        TIterator it = begin(journaledString, Standard());
        it += 4;
        TIterator itEnd = end(journaledString, Standard());
        SEQAN_ASSERT(it == itEnd);
    }

    // Case 2: Go in between beginning and end.
    {
        TIterator it = begin(journaledString, Standard());
        ++it; ++it;
        TIterator it2 = begin(journaledString, Standard());
        it2 += 2;
        SEQAN_ASSERT(it == it2);
    }

    // Case 3: Go over end of journal string.
    {
        TIterator it = begin(journaledString, Standard());
        it += 5;
        TIterator itEnd = end(journaledString, Standard());
        SEQAN_ASSERT(it == itEnd);
        it += 10;
        SEQAN_ASSERT(it == itEnd);
    }
}

// Test relation operators on iterators.
template <typename THostSpec, typename TStringJournalSpec, typename TBuffSpec>
void testJournaledStringIteratorRelations(THostSpec const &, TStringJournalSpec const &, TBuffSpec const &)
{
    typedef String<char, Journaled<THostSpec, TStringJournalSpec, TBuffSpec> > TJournaledString;
    typedef typename Iterator<TJournaledString, Standard>::Type TIterator;

    CharString charStr = "test";
    TJournaledString journaledString;
    setHost(journaledString, charStr);

    TIterator it1 = begin(journaledString, Standard());
    TIterator it2 = begin(journaledString, Standard());
    ++it2; ++it2;
    TIterator it3 = end(journaledString, Standard());

    SEQAN_ASSERT(it1 < it2);
    SEQAN_ASSERT(it1 <= it2);
    SEQAN_ASSERT(it2 < it3);
    SEQAN_ASSERT(it2 <= it3);

    SEQAN_ASSERT(it1 <= it1);
    SEQAN_ASSERT(it2 <= it2);
    SEQAN_ASSERT(it3 <= it3);

    SEQAN_ASSERT(it2 > it1);
    SEQAN_ASSERT(it2 >= it1);
    SEQAN_ASSERT(it3 > it2);
    SEQAN_ASSERT(it3 >= it2);

    SEQAN_ASSERT(it1 >= it1);
    SEQAN_ASSERT(it2 >= it2);
    SEQAN_ASSERT(it3 >= it3);
}

template <typename THostSpec, typename TStringJournalSpec, typename TBuffSpec>
void testJournaledStringIteratorDecrement(THostSpec const &, TStringJournalSpec const &, TBuffSpec const &)
{
    typedef String<char, Journaled<THostSpec, TStringJournalSpec, TBuffSpec> > TJournalString;
    typedef typename Iterator<TJournalString, Standard>::Type TIterator;
    typedef typename Host<TJournalString>::Type THost;

    {
        TJournalString journal("aaaa");

        insert(journal,2,"xx"); // aaxxaa
        erase(journal,4,5); // aaxxa
        insert(journal,0,"y"); //yaaxxa

        SEQAN_ASSERT_EQ(journal, "yaaxxa");

        {
            TIterator it = end(journal);
            --it;
            SEQAN_ASSERT_EQ(value(it), 'a');
            --it;
            SEQAN_ASSERT_EQ(value(it), 'x');
            --it;
            SEQAN_ASSERT_EQ(value(it), 'x');
            --it;
            SEQAN_ASSERT_EQ(value(it), 'a');
            --it;
            SEQAN_ASSERT_EQ(value(it), 'a');
            --it;
            SEQAN_ASSERT_EQ(value(it), 'y');
        }

        {
            TIterator it = end(journal);
            it--;
            SEQAN_ASSERT_EQ(value(it), 'a');
            it--;
            SEQAN_ASSERT_EQ(value(it), 'x');
            it--;
            SEQAN_ASSERT_EQ(value(it), 'x');
            it--;
            SEQAN_ASSERT_EQ(value(it), 'a');
            it--;
            SEQAN_ASSERT_EQ(value(it), 'a');
            it--;
            SEQAN_ASSERT_EQ(value(it), 'y');
        }

        {
            TIterator it = end(journal)-2;

            SEQAN_ASSERT_EQ(value(it), 'x');
            it -= 2;
            SEQAN_ASSERT_EQ(value(it), 'a');
            it -= 1;
            SEQAN_ASSERT_EQ(value(it), 'a');
            it -= 5;
            SEQAN_ASSERT_EQ(value(it), 'y');
        }
    }

    {
        THost src = "aaaaaaa";
        TJournalString journal(src);

        SEQAN_ASSERT_EQ(host(journal), src);

        insert(journal,2,"xyxyx"); // aaxyxyxaaaaa
        erase(journal,8,10); // aaxyxyxaaa

        {
            TIterator it = end(journal);
            --it;
            SEQAN_ASSERT_EQ(value(it), 'a');
            --it;
            SEQAN_ASSERT_EQ(value(it), 'a');
            --it;
            SEQAN_ASSERT_EQ(value(it), 'a');
            --it;
            SEQAN_ASSERT_EQ(value(it), 'x');
            --it;
            SEQAN_ASSERT_EQ(value(it), 'y');
            --it;
            SEQAN_ASSERT_EQ(value(it), 'x');
            --it;
            SEQAN_ASSERT_EQ(value(it), 'y');
            --it;
            SEQAN_ASSERT_EQ(value(it), 'x');
            --it;
            SEQAN_ASSERT_EQ(value(it), 'a');
            --it;
            SEQAN_ASSERT_EQ(value(it), 'a');
        }

        {
            TIterator it = end(journal);
            it--;
            SEQAN_ASSERT_EQ(value(it), 'a');
            it--;
            SEQAN_ASSERT_EQ(value(it), 'a');
            it--;
            SEQAN_ASSERT_EQ(value(it), 'a');
            it--;
            SEQAN_ASSERT_EQ(value(it), 'x');
            it--;
            SEQAN_ASSERT_EQ(value(it), 'y');
            it--;
            SEQAN_ASSERT_EQ(value(it), 'x');
            it--;
            SEQAN_ASSERT_EQ(value(it), 'y');
            it--;
            SEQAN_ASSERT_EQ(value(it), 'x');
            it--;
            SEQAN_ASSERT_EQ(value(it), 'a');
            it--;
            SEQAN_ASSERT_EQ(value(it), 'a');
        }

        {
            TIterator it = end(journal)-3;

            SEQAN_ASSERT_EQ(value(it), 'a');
            it -= 1;
            SEQAN_ASSERT_EQ(value(it), 'x');
            it -= 3;
            SEQAN_ASSERT_EQ(value(it), 'y');
            it -= 2;
            SEQAN_ASSERT_EQ(value(it), 'a');
            it -= 20;
            SEQAN_ASSERT_EQ(value(it), 'a');

            it = end(journal) - 20;
            SEQAN_ASSERT(it == TIterator(begin(journal)));
        }
    }
}

template <typename THostSpec, typename TStringJournalSpec, typename TBuffSpec>
void testJournaledStringIteratorSetPosition(THostSpec const &, TStringJournalSpec const &, TBuffSpec const &)
{
    typedef String<char, Journaled<THostSpec, TStringJournalSpec, TBuffSpec> > TJournalString;
    typedef typename Iterator<TJournalString, Standard>::Type TIterator;
    typedef typename Host<TJournalString>::Type THost;

    {
        THost hostSeq = "aacac";
        TJournalString journal(hostSeq);

        insert(journal,2,"bb"); // aabbcac
        erase(journal,5,6); // aaxxcc
        insert(journal,0,"d"); //daabbcc

        TIterator journalIt = begin(journal);
        SEQAN_ASSERT_EQ(value(journalIt),'d');
        setPosition(journalIt, 3);
        SEQAN_ASSERT_EQ(value(journalIt),'b');
        ++journalIt;
        SEQAN_ASSERT_EQ(value(journalIt),'b');
        ++journalIt;
        SEQAN_ASSERT_EQ(value(journalIt),'c');
        ++journalIt;
        SEQAN_ASSERT_EQ(value(journalIt),'c');
    }

    {
        THost hostSeq = "aacac";
        TJournalString journal(hostSeq);

        insert(journal,2,"bb"); // aabbcac
        erase(journal,5,6); // aaxxcc
        insert(journal,0,"d"); //daabbcc

        TIterator journalIt = begin(journal);
        SEQAN_ASSERT_EQ(value(journalIt),'d');
        setPosition(journalIt, 5);
        SEQAN_ASSERT_EQ(value(journalIt),'c');
        --journalIt;
        SEQAN_ASSERT_EQ(value(journalIt),'b');
        --journalIt;
        SEQAN_ASSERT_EQ(value(journalIt),'b');
        --journalIt;
        SEQAN_ASSERT_EQ(value(journalIt),'a');
    }
}

template <typename THostSpec, typename TStringJournalSpec, typename TBuffSpec>
void testJournaledStringIteratorPosition(THostSpec const &, TStringJournalSpec const &, TBuffSpec const &)
{
    typedef String<char, Journaled<THostSpec, TStringJournalSpec, TBuffSpec> > TJournalString;
    typedef typename Iterator<TJournalString, Standard>::Type TIterator;
    typedef typename Host<TJournalString>::Type THost;
    typedef typename Position<TIterator>::Type TPos;

    {
        THost hostSeq = "aacac";
        TJournalString journal(hostSeq);

        insert(journal,2,"bb"); // aabbcac
        erase(journal,5,6); // aaxxcc
        insert(journal,0,"d"); //daabbcc

        TIterator journalIt = begin(journal);
        SEQAN_ASSERT_EQ(position(journalIt), static_cast<TPos>(0));
        setPosition(journalIt, 3);
        SEQAN_ASSERT_EQ(position(journalIt), static_cast<TPos>(3));
        ++journalIt;
        SEQAN_ASSERT_EQ(position(journalIt), static_cast<TPos>(4));
        journalIt += 2;
        SEQAN_ASSERT_EQ(position(journalIt), static_cast<TPos>(6));
        journalIt -= 5;
        SEQAN_ASSERT_EQ(position(journalIt), static_cast<TPos>(1));
    }
}

template <typename THostSpec, typename TStringJournalSpec, typename TBuffSpec>
void testJournaledStringIteratorAtBegin(THostSpec const &, TStringJournalSpec const &, TBuffSpec const &)
{
    typedef String<char, Journaled<THostSpec, TStringJournalSpec, TBuffSpec> > TJournalString;
    typedef typename Iterator<TJournalString, Rooted>::Type TIterator;
    typedef typename Host<TJournalString>::Type THost;

    {
        THost hostSeq = "aacac";
        TJournalString journal(hostSeq);

        insert(journal,2,"bb"); // aabbcac
        erase(journal,5,6); // aaxxcc
        insert(journal,0,"d"); //daabbcc

        TIterator journalIt = begin(journal, Rooted());
        SEQAN_ASSERT(atBegin(journalIt));
        ++journalIt;
        SEQAN_ASSERT_NOT(atBegin(journalIt));
    }
}

template <typename THostSpec, typename TStringJournalSpec, typename TBuffSpec>
void testJournaledStringIteratorAtEnd(THostSpec const &, TStringJournalSpec const &, TBuffSpec const &)
{
    typedef String<char, Journaled<THostSpec, TStringJournalSpec, TBuffSpec> > TJournalString;
    typedef typename Iterator<TJournalString, Rooted>::Type TIterator;
    typedef typename Host<TJournalString>::Type THost;

    {
        THost hostSeq = "aacac";
        TJournalString journal(hostSeq);

        insert(journal,2,"bb"); // aabbcac
        erase(journal,5,6); // aaxxcc
        insert(journal,0,"d"); //daabbcc

        TIterator journalIt = end(journal, Rooted());
        SEQAN_ASSERT(atEnd(journalIt));
        --journalIt;
        SEQAN_ASSERT_NOT(atEnd(journalIt));
    }
}

template <typename THostSpec, typename TStringJournalSpec, typename TBuffSpec>
void testJournaledStringIteratorGoBegin(THostSpec const &, TStringJournalSpec const &, TBuffSpec const &)
{
    typedef String<char, Journaled<THostSpec, TStringJournalSpec, TBuffSpec> > TJournalString;
    typedef typename Iterator<TJournalString, Rooted>::Type TIterator;
    typedef typename Host<TJournalString>::Type THost;

    {
        THost hostSeq = "aacac";
        TJournalString journal(hostSeq);

        insert(journal,2,"bb"); // aabbcac
        erase(journal,5,6); // aaxxcc
        insert(journal,0,"d"); //daabbcc

        TIterator journalIt = end(journal, Rooted());
        SEQAN_ASSERT_NOT(atBegin(journalIt));
        goBegin(journalIt);
        SEQAN_ASSERT(atBegin(journalIt));
    }
}

template <typename THostSpec, typename TStringJournalSpec, typename TBuffSpec>
void testJournaledStringIteratorGoEnd(THostSpec const &, TStringJournalSpec const &, TBuffSpec const &)
{
    typedef String<char, Journaled<THostSpec, TStringJournalSpec, TBuffSpec> > TJournalString;
    typedef typename Iterator<TJournalString, Rooted>::Type TIterator;
    typedef typename Host<TJournalString>::Type THost;

    {
        THost hostSeq = "aacac";
        TJournalString journal(hostSeq);

        insert(journal,2,"bb"); // aabbcac
        erase(journal,5,6); // aaxxcc
        insert(journal,0,"d"); //daabbcc

        TIterator journalIt = begin(journal, Rooted());
        SEQAN_ASSERT_NOT(atEnd(journalIt));
        goEnd(journalIt);
        SEQAN_ASSERT(atEnd(journalIt));
    }
}

template <typename THostSpec, typename TStringJournalSpec, typename TBuffSpec>
void testJournaledStringIteratorContainer(THostSpec const &, TStringJournalSpec const &, TBuffSpec const &)
{
    typedef String<char, Journaled<THostSpec, TStringJournalSpec, TBuffSpec> > TJournalString;
    typedef typename Iterator<TJournalString, Standard>::Type TIterator;
    typedef typename Host<TJournalString>::Type THost;

    {
        THost hostSeq = "aacac";
        TJournalString journal(hostSeq);

        insert(journal,2,"bb"); // aabbcac
        erase(journal,5,6); // aaxxcc
        insert(journal,0,"d"); //daabbcc

        TIterator journalIt = begin(journal, Standard());
        SEQAN_ASSERT_EQ(container(journalIt), journal);
    }
}

// Tag: SortedArray()

SEQAN_DEFINE_TEST(test_sequence_journaled_sorted_array_iterator_difference)
{
    testJournaledStringIteratorDifference(Alloc<void>(), SortedArray(), Alloc<void>());
    testJournaledStringIteratorDifference(Alloc<Nothing>(), SortedArray(), Alloc<void>());
}

SEQAN_DEFINE_TEST(test_sequence_journaled_sorted_array_iterator_sum)
{
    testJournaledStringIteratorSum(Alloc<void>(), SortedArray(), Alloc<void>());
    testJournaledStringIteratorSum(Alloc<Nothing>(), SortedArray(), Alloc<void>());
}

SEQAN_DEFINE_TEST(test_sequence_journaled_sorted_array_iterator_relations)
{
    testJournaledStringIteratorRelations(Alloc<void>(), SortedArray(), Alloc<void>());
    testJournaledStringIteratorRelations(Alloc<Nothing>(), SortedArray(), Alloc<void>());
}

SEQAN_DEFINE_TEST(test_sequence_journaled_sorted_array_iterator_decrement)
{
    testJournaledStringIteratorDecrement(Alloc<void>(), SortedArray(), Alloc<void>());
    testJournaledStringIteratorDecrement(Alloc<Nothing>(), SortedArray(), Alloc<void>());
}

SEQAN_DEFINE_TEST(test_sequence_journaled_sorted_array_iterator_set_position)
{
    testJournaledStringIteratorSetPosition(Alloc<void>(), SortedArray(), Alloc<void>());
    testJournaledStringIteratorSetPosition(Alloc<Nothing>(), SortedArray(), Alloc<void>());
}

SEQAN_DEFINE_TEST(test_sequence_journaled_sorted_array_iterator_position)
{
    testJournaledStringIteratorPosition(Alloc<void>(), SortedArray(), Alloc<void>());
    testJournaledStringIteratorPosition(Alloc<Nothing>(), SortedArray(), Alloc<void>());
}

SEQAN_DEFINE_TEST(test_sequence_journaled_sorted_array_iterator_rooted_at_begin)
{
    testJournaledStringIteratorAtBegin(Alloc<void>(), SortedArray(), Alloc<void>());
    testJournaledStringIteratorAtBegin(Alloc<Nothing>(), SortedArray(), Alloc<void>());
}

SEQAN_DEFINE_TEST(test_sequence_journaled_sorted_array_iterator_rooted_at_end)
{
    testJournaledStringIteratorAtEnd(Alloc<void>(), SortedArray(), Alloc<void>());
    testJournaledStringIteratorAtEnd(Alloc<Nothing>(), SortedArray(), Alloc<void>());
}

SEQAN_DEFINE_TEST(test_sequence_journaled_sorted_array_iterator_rooted_go_begin)
{
    testJournaledStringIteratorGoBegin(Alloc<void>(), SortedArray(), Alloc<void>());
    testJournaledStringIteratorGoBegin(Alloc<Nothing>(), SortedArray(), Alloc<void>());
}

SEQAN_DEFINE_TEST(test_sequence_journaled_sorted_array_iterator_rooted_go_end)
{
    testJournaledStringIteratorGoEnd(Alloc<void>(), SortedArray(), Alloc<void>());
    testJournaledStringIteratorGoEnd(Alloc<Nothing>(), SortedArray(), Alloc<void>());
}

SEQAN_DEFINE_TEST(test_sequence_journaled_sorted_array_iterator_rooted_container)
{
    testJournaledStringIteratorContainer(Alloc<void>(), SortedArray(), Alloc<void>());
    testJournaledStringIteratorContainer(Alloc<Nothing>(), SortedArray(), Alloc<void>());
}

#endif  // TEST_SEQUENCE_JOURNALED_TEST_SEQUENCE_JOURNALED_ITERATOR_H_
