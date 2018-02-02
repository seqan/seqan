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
// Tests for the sequence_journaled module.  Note that we only test journaled
// strings here.  In a perfect world, we would have tests for each atomic
// part of the module but for the moment, this has to suffice.  Instead of
// testing each of the many corner cases for the Journal Entries data
// structures, we rely on "fuzzying", i.e. executing a random set of
// operations and hope that all corner cases occur there.
// ==========================================================================

#ifndef TEST_SEQUENCE_JOURNALED_TEST_SEQUENCE_JOURNALED_H_
#define TEST_SEQUENCE_JOURNALED_TEST_SEQUENCE_JOURNALED_H_

#include <cstdlib>
#include <sstream>
#include <string>
#include <random>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/sequence_journaled.h>

using namespace seqan;

// Test assign(), operator=()
template <typename TStringJournalSpec>
void testJournaledStringAssign(TStringJournalSpec const &)
{
    typedef String<char, Journaled<Alloc<void>, TStringJournalSpec> > TJournaledString;

    CharString charStr = "test";
    TJournaledString journaledString(charStr);

    CharString charStr2 = "not a test!";

    // Test assignment operator with other string type.
    {
        journaledString = charStr2;

        std::stringstream tmp1;
        tmp1 << journaledString;
        SEQAN_ASSERT_EQ("not a test!", tmp1.str());

        std::stringstream tmp2;
        tmp2 << host(journaledString);
        SEQAN_ASSERT_EQ("test", tmp2.str());
    }

    // Test assign() with other string.
    {
        assign(journaledString, charStr2);

        std::stringstream tmp1;
        tmp1 << journaledString;
        SEQAN_ASSERT_EQ("not a test!", tmp1.str());

        std::stringstream tmp2;
        tmp2 << host(journaledString);
        SEQAN_ASSERT_EQ("test", tmp2.str());
    }

    TJournaledString journaledString2(charStr2);

    // Test assignment operator with same journaled string type.
    {
        // We need to copy over all information from the source to the target.
        journaledString = journaledString2;

        std::stringstream tmp1;
        tmp1 << journaledString;
        SEQAN_ASSERT_EQ("not a test!", tmp1.str());

        std::stringstream tmp2;
        tmp2 << host(journaledString);
        SEQAN_ASSERT_EQ("not a test!", tmp2.str());
    }

    // Test assign() with same journaled string type.
    {
        assign(journaledString, journaledString2);

        std::stringstream tmp1;
        tmp1 << journaledString;
        SEQAN_ASSERT_EQ("not a test!", tmp1.str());

        std::stringstream tmp2;
        tmp2 << host(journaledString);
        SEQAN_ASSERT_EQ("not a test!", tmp2.str());
    }
}

// Test set()
template <typename TStringJournalSpec>
void testJournaledStringSet(TStringJournalSpec const &)
{
    typedef String<char, Journaled<Alloc<void>, TStringJournalSpec> > TJournaledString;

    CharString charStr = "test";
    TJournaledString journaledString(charStr);

    CharString charStr2 = "not a test!";
    TJournaledString journaledString2(charStr2);

    // Test set() with other string.
    {
        set(journaledString, charStr2);

        std::stringstream tmp1;
        tmp1 << journaledString;
        SEQAN_ASSERT_EQ("not a test!", tmp1.str());

        std::stringstream tmp2;
        tmp2 << host(journaledString);
        SEQAN_ASSERT_EQ("test", tmp2.str());
    }

    // Test set() with same journaled string type.
    {
        set(journaledString, journaledString2);

        std::stringstream tmp1;
        tmp1 << journaledString;
        SEQAN_ASSERT_EQ("not a test!", tmp1.str());

        std::stringstream tmp2;
        tmp2 << host(journaledString);
        SEQAN_ASSERT_EQ("not a test!", tmp2.str());
    }
}

template <typename TJournalSpec>
void testJournaledStringReset(TJournalSpec const & /*spec*/)
{
    CharString hostStr = "test";

    String<char, Journaled<Alloc<>, TJournalSpec> > journaledString(hostStr);

    insert(journaledString, 4, "XXL");
    erase(journaledString, 1, 2);

    std::stringstream test;
    test << journaledString;

    SEQAN_ASSERT_EQ(test.str(), "tstXXL");
    reset(journaledString);

    test.str("");
    test.clear();
    test << journaledString;

    SEQAN_ASSERT_EQ(test.str(), "");
    SEQAN_ASSERT_EQ(empty(journaledString._holder), true);
    SEQAN_ASSERT_EQ(empty(journaledString._insertionBuffer), true);
    SEQAN_ASSERT_EQ(length(journaledString), 0u);
}


// Test setHost(), host().
template <typename TStringJournalSpec>
void testJournaledStringHost(TStringJournalSpec const &)
{
    CharString charStr = "test";
    String<char, Journaled<Alloc<void>, TStringJournalSpec> > journaledString;

    setHost(journaledString, charStr);
    SEQAN_ASSERT_EQ(&charStr, &host(journaledString));
    SEQAN_ASSERT_EQ(charStr, host(journaledString));
}


// Test clear().
template <typename TStringJournalSpec>
void testJournaledStringClear(TStringJournalSpec const &)
{
    typedef String<char, Journaled<Alloc<void>, TStringJournalSpec> > TJournaledString;
    typedef typename Size<TJournaledString>::Type TSize;

    CharString charStr = "test";
    TJournaledString journaledString(charStr);

    insert(journaledString, 4, "!");
    {
        std::stringstream tmp;
        tmp << journaledString;
        SEQAN_ASSERT_EQ("test!", tmp.str());
    }

    clear(journaledString);
    {
        std::stringstream tmp;
        tmp << journaledString;
        SEQAN_ASSERT_EQ("test", tmp.str());
        SEQAN_ASSERT_EQ(static_cast<TSize>(4), length(journaledString));
    }
}

// Test empty().
template <typename TStringJournalSpec>
void testJournaledStringEmpty(TStringJournalSpec const &)
{
    typedef String<char, Journaled<Alloc<void>, TStringJournalSpec> > TJournaledString;

    CharString charStr = "test";

    {  // Empty string.
        TJournaledString journaledString;
        SEQAN_ASSERT_EQ(empty(journaledString), true);
    }

    {  // Host set.
        TJournaledString journaledString;
        setHost(journaledString, charStr);
        SEQAN_ASSERT_EQ(empty(journaledString), false);
    }

    {  // Host is set but empty;
        TJournaledString journaledString;
        create(journaledString._holder);
        SEQAN_ASSERT_EQ(empty(journaledString), true);
    }

    {  // Host is set and empty but insertion buffer is filled;
        TJournaledString journaledString;
        CharString otherString = "Other String";
        assign(journaledString, otherString);
        SEQAN_ASSERT_EQ(empty(journaledString), false);
    }
}


// Test erase() with position only
template <typename TStringJournalSpec>
void testJournaledStringErasePosition(TStringJournalSpec const &)
{
    CharString charStr = "test";
    String<char, Journaled<Alloc<void>, TStringJournalSpec> > journaledString(charStr);

    erase(journaledString, 1);
    std::stringstream tmp;
    tmp << journaledString;
    SEQAN_ASSERT_EQ("tst", tmp.str());
}


// Test erase() with begin/end parameters.
template <typename TStringJournalSpec>
void testJournaledStringEraseBeginEnd(TStringJournalSpec const &)
{
    CharString charStr = "test";
    String<char, Journaled<Alloc<void>, TStringJournalSpec> > journaledString(charStr);

    erase(journaledString, 1, 3);
    std::stringstream tmp;
    tmp << journaledString;
    SEQAN_ASSERT_EQ("tt", tmp.str());
}


// Test insert().
template <typename TStringJournalSpec>
void testJournaledStringInsert(TStringJournalSpec const &)
{
    CharString charStr = "test";
    String<char, Journaled<Alloc<void>, TStringJournalSpec> > journaledString(charStr);

    insert(journaledString, 1, "!!");
    std::stringstream tmp;
    tmp << journaledString;
    SEQAN_ASSERT_EQ("t!!est", tmp.str());
}


// Test insertValue().
template <typename TStringJournalSpec>
void testJournaledStringInsertValue(TStringJournalSpec const &)
{
    CharString charStr = "test";
    String<char, Journaled<Alloc<void>, TStringJournalSpec> > journaledString(charStr);

    insert(journaledString, 2, "X");
    std::stringstream tmp;
    tmp << journaledString;
    SEQAN_ASSERT_EQ("teXst", tmp.str());
}


// Test assignValue().
template <typename TStringJournalSpec>
void testJournaledStringAssignValue(TStringJournalSpec const &)
{
    CharString charStr = "test";
    String<char, Journaled<Alloc<void>, TStringJournalSpec> > journaledString(charStr);

    assignValue(journaledString, 2, 'X');
    std::stringstream tmp;
    tmp << journaledString;
    SEQAN_ASSERT_EQ("teXt", tmp.str());
}


// Test operator[]().
template <typename TStringJournalSpec>
void testJournaledStringSubscriptOperator(TStringJournalSpec const &)
{
    CharString charStr = "test";
    String<char, Journaled<Alloc<void>, TStringJournalSpec> > journaledString(charStr);

    SEQAN_ASSERT_EQ(journaledString[0], 't');
    SEQAN_ASSERT_EQ(journaledString[3], 't');

    // static_cast<Nothing>(journaledString[2]);
    journaledString[2] = 'X';

    SEQAN_ASSERT_EQ(journaledString[0], 't');
    SEQAN_ASSERT_EQ(journaledString[2], 'X');
    SEQAN_ASSERT_EQ(journaledString[3], 't');

    std::stringstream tmp;
    tmp << journaledString;
    SEQAN_ASSERT_EQ("teXt", tmp.str());
}


// Test assignInfix().
template <typename TStringJournalSpec>
void testJournaledStringAssignInfix(TStringJournalSpec const &)
{
    CharString charStr = "test";
    String<char, Journaled<Alloc<void>, TStringJournalSpec> > journaledString(charStr);

    assignInfix(journaledString, 2, 4, "quick brown fox");
    std::stringstream tmp;
    tmp << journaledString;
    SEQAN_ASSERT_EQ("tequick brown fox", tmp.str());
}


// Test length().
template <typename TStringJournalSpec>
void testJournaledStringLength(TStringJournalSpec const &)
{
    CharString charStr = "test";
    String<char, Journaled<Alloc<void>, TStringJournalSpec> > journaledString(charStr);

    assignInfix(journaledString, 2, 3, "XX");
    SEQAN_ASSERT_EQ(length("teXXt"), length(journaledString));
}

// Test flatten().
template <typename TStringJournalSpec>
void testJournalStringFlatten(TStringJournalSpec const &)
{
    typedef String<char, Journaled<Alloc<void>, TStringJournalSpec> > TJournaledString;

    {  // Test insertion at beginning.
        CharString charStr = "test";
        TJournaledString journaledString(charStr);
        insert(journaledString, 0, "XX");
        SEQAN_ASSERT_EQ(host(journaledString), "test");
        SEQAN_ASSERT_EQ(journaledString, "XXtest");
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).virtualPosition, 0u);
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).physicalPosition, 0u);
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).length, 2u);
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).segmentSource, SOURCE_PATCH);

        flatten(journaledString);

        SEQAN_ASSERT_EQ(host(journaledString), "XXtest");
        SEQAN_ASSERT_EQ(journaledString, "XXtest");
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).virtualPosition, 0u);
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).physicalPosition, 0u);
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).length, length("XXtest"));
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).segmentSource, SOURCE_ORIGINAL);
    }

    {  // Test insertion in middle part.
        CharString charStr = "test";
        TJournaledString journaledString(charStr);
        insert(journaledString, 2, "XX");
        SEQAN_ASSERT_EQ(host(journaledString), "test");
        SEQAN_ASSERT_EQ(journaledString, "teXXst");
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).virtualPosition, 0u);
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).physicalPosition, 0u);
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).length, 2u);
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).segmentSource, SOURCE_ORIGINAL);

        flatten(journaledString);

        SEQAN_ASSERT_EQ(host(journaledString), "teXXst");
        SEQAN_ASSERT_EQ(journaledString, "teXXst");
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).virtualPosition, 0u);
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).physicalPosition, 0u);
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).length, length("teXXst"));
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).segmentSource, SOURCE_ORIGINAL);
    }

    {  // Test insertion at end part.
        CharString charStr = "test";
        TJournaledString journaledString(charStr);
        insert(journaledString, 4, "XX");
        SEQAN_ASSERT_EQ(host(journaledString), "test");
        SEQAN_ASSERT_EQ(journaledString, "testXX");
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).virtualPosition, 0u);
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).physicalPosition, 0u);
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).length, 4u);
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).segmentSource, SOURCE_ORIGINAL);

        flatten(journaledString);

        SEQAN_ASSERT_EQ(host(journaledString), "testXX");
        SEQAN_ASSERT_EQ(journaledString, "testXX");
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).virtualPosition, 0u);
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).physicalPosition, 0u);
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).length, length("testXX"));
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).segmentSource, SOURCE_ORIGINAL);
    }

    {  // Test deletion at beginning.
        CharString charStr = "test";
        TJournaledString journaledString(charStr);
        erase(journaledString, 0, 2);
        SEQAN_ASSERT_EQ(host(journaledString), "test");
        SEQAN_ASSERT_EQ(journaledString, "st");
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).virtualPosition, 0u);
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).physicalPosition, 2u);
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).length, 2u);
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).segmentSource, SOURCE_ORIGINAL);

        flatten(journaledString);

        SEQAN_ASSERT_EQ(host(journaledString), "st");
        SEQAN_ASSERT_EQ(journaledString, "st");
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).virtualPosition, 0u);
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).physicalPosition, 0u);
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).length, length("st"));
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).segmentSource, SOURCE_ORIGINAL);
    }

    {  // Test deletion in middle part.
        CharString charStr = "test";
        TJournaledString journaledString(charStr);
        erase(journaledString, 1, 3);
        SEQAN_ASSERT_EQ(host(journaledString), "test");
        SEQAN_ASSERT_EQ(journaledString, "tt");
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).virtualPosition, 0u);
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).physicalPosition, 0u);
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).length, 1u);
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).segmentSource, SOURCE_ORIGINAL);

        flatten(journaledString);

        SEQAN_ASSERT_EQ(host(journaledString), "tt");
        SEQAN_ASSERT_EQ(journaledString, "tt");
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).virtualPosition, 0u);
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).physicalPosition, 0u);
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).length, length("tt"));
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).segmentSource, SOURCE_ORIGINAL);
    }

    {  // Test deletion at end part.
        CharString charStr = "test";
        TJournaledString journaledString(charStr);
        erase(journaledString, 2, 4);
        SEQAN_ASSERT_EQ(host(journaledString), "test");
        SEQAN_ASSERT_EQ(journaledString, "te");
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).virtualPosition, 0u);
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).physicalPosition, 0u);
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).length, 2u);
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).segmentSource, SOURCE_ORIGINAL);

        flatten(journaledString);

        SEQAN_ASSERT_EQ(host(journaledString), "te");
        SEQAN_ASSERT_EQ(journaledString, "te");
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).virtualPosition, 0u);
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).physicalPosition, 0u);
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).length, length("te"));
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).segmentSource, SOURCE_ORIGINAL);
    }

    {  // Test insertion/deletion at beginning.
        CharString charStr = "test";
        TJournaledString journaledString(charStr);
        insert(journaledString, 0, "XX");
        erase(journaledString, 2, 4);
        SEQAN_ASSERT_EQ(host(journaledString), "test");
        SEQAN_ASSERT_EQ(journaledString, "XXst");
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).virtualPosition, 0u);
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).physicalPosition, 0u);
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).length, 2u);
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).segmentSource, SOURCE_PATCH);

        flatten(journaledString);

        SEQAN_ASSERT_EQ(host(journaledString), "XXst");
        SEQAN_ASSERT_EQ(journaledString, "XXst");
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).virtualPosition, 0u);
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).physicalPosition, 0u);
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).length, length("XXst"));
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).segmentSource, SOURCE_ORIGINAL);
    }

    {  // Test insertion/deletion in middle part.
        CharString charStr = "test";
        TJournaledString journaledString(charStr);
        erase(journaledString, 1, 3);
        insert(journaledString, 1, "XX");
        SEQAN_ASSERT_EQ(host(journaledString), "test");
        SEQAN_ASSERT_EQ(journaledString, "tXXt");
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).virtualPosition, 0u);
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).physicalPosition, 0u);
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).length, 1u);
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).segmentSource, SOURCE_ORIGINAL);

        flatten(journaledString);

        SEQAN_ASSERT_EQ(host(journaledString), "tXXt");
        SEQAN_ASSERT_EQ(journaledString, "tXXt");
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).virtualPosition, 0u);
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).physicalPosition, 0u);
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).length, length("tXXt"));
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).segmentSource, SOURCE_ORIGINAL);
    }

    {  // Test insertion/deletion at end part.
        CharString charStr = "test";
        TJournaledString journaledString(charStr);
        insert(journaledString, 4, "XX");
        erase(journaledString, 2, 4);
        SEQAN_ASSERT_EQ(host(journaledString), "test");
        SEQAN_ASSERT_EQ(journaledString, "teXX");
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).virtualPosition, 0u);
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).physicalPosition, 0u);
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).length, 2u);
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).segmentSource, SOURCE_ORIGINAL);

        flatten(journaledString);

        SEQAN_ASSERT_EQ(host(journaledString), "teXX");
        SEQAN_ASSERT_EQ(journaledString, "teXX");
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).virtualPosition, 0u);
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).physicalPosition, 0u);
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).length, length("teXX"));
        SEQAN_ASSERT_EQ(value(begin(journaledString._journalEntries)).segmentSource, SOURCE_ORIGINAL);
    }
}


// Test conversion of virtual to host position.
template <typename TStringJournalSpec>
void testJournaledStringVirtualToHostPosition()
{
    {  // Test right-based projection.
        CharString testSrc = "leftright";
        //  0123   45678
        //--left---right
        //zzl--tyyyri--txx
        //012  345678  901
        String<char, Journaled<Alloc<> > > jString(testSrc);

        insert(jString, 9, "xx");
        erase(jString, 6, 8);
        insert(jString, 4, "yyy");
        erase(jString, 1, 3);
        insert(jString, 0, "zz");

        SEQAN_ASSERT_EQ(virtualToHostPosition(jString, 0u), 0u);
        SEQAN_ASSERT_EQ(virtualToHostPosition(jString, 1u), 0u);
        SEQAN_ASSERT_EQ(virtualToHostPosition(jString, 2u), 0u);
        SEQAN_ASSERT_EQ(virtualToHostPosition(jString, 3u), 3u);
        SEQAN_ASSERT_EQ(virtualToHostPosition(jString, 4u), 4u);
        SEQAN_ASSERT_EQ(virtualToHostPosition(jString, 5u), 4u);
        SEQAN_ASSERT_EQ(virtualToHostPosition(jString, 6u), 4u);
        SEQAN_ASSERT_EQ(virtualToHostPosition(jString, 7u), 4u);
        SEQAN_ASSERT_EQ(virtualToHostPosition(jString, 8u), 5u);
        SEQAN_ASSERT_EQ(virtualToHostPosition(jString, 9u), 8u);
        SEQAN_ASSERT_EQ(virtualToHostPosition(jString, 10u), 9u);
        SEQAN_ASSERT_EQ(virtualToHostPosition(jString, 11u), 9u);
        SEQAN_ASSERT_EQ(virtualToHostPosition(jString, 20u), 9u);
    }
    {
        CharString charStr = "ABCDE";
        String<char, Journaled<Alloc<void>, TStringJournalSpec> > journaledString(charStr);
        insert(journaledString, 1, "XX");
        SEQAN_ASSERT_EQ(0u, virtualToHostPosition(journaledString, 0u));
        SEQAN_ASSERT_EQ(1u, virtualToHostPosition(journaledString, 1u));
        SEQAN_ASSERT_EQ(1u, virtualToHostPosition(journaledString, 2u));
        SEQAN_ASSERT_EQ(1u, virtualToHostPosition(journaledString, 3u));
        SEQAN_ASSERT_EQ(2u, virtualToHostPosition(journaledString, 4u));
        SEQAN_ASSERT_EQ(3u, virtualToHostPosition(journaledString, 5u));
        SEQAN_ASSERT_EQ(4u, virtualToHostPosition(journaledString, 6u));
        SEQAN_ASSERT_EQ(5u, virtualToHostPosition(journaledString, 7u));
    }
}

template <typename TJournaledStringSpec>
void testJournaledStringHostToVirtualPosition()
{
    /*
     * write concrete test scenario
     * Given a JournalString with insertions and deletions
     * if the source position of ...
     */
    { // test 1: insertion
                            //01234567
        CharString hostStr = "ACGTACGT";
        String<char, Journaled<Alloc<>, TJournaledStringSpec> > journaledStr(hostStr);

        insert(journaledStr, 2, "TTT");
                        //ACTTTGTACGT
        // journaled str: 01234567890
        // host str:      01---234567
        //              //AC   GTACGT

        SEQAN_ASSERT_EQ(0u, hostToVirtualPosition(journaledStr, 0));
        SEQAN_ASSERT_EQ(1u, hostToVirtualPosition(journaledStr, 1));
        SEQAN_ASSERT_EQ(5u, hostToVirtualPosition(journaledStr, 2));
        SEQAN_ASSERT_EQ(6u, hostToVirtualPosition(journaledStr, 3));
        SEQAN_ASSERT_EQ(7u, hostToVirtualPosition(journaledStr, 4));
        SEQAN_ASSERT_EQ(8u, hostToVirtualPosition(journaledStr, 5));
        SEQAN_ASSERT_EQ(9u, hostToVirtualPosition(journaledStr, 6));
        SEQAN_ASSERT_EQ(10u, hostToVirtualPosition(journaledStr, 7));
    }

    { // test 2: deletion
                            //01234567
        CharString hostStr = "ACGTACGT";
        String<char, Journaled<Alloc<>, TJournaledStringSpec> > journaledStr(hostStr);

        erase(journaledStr, 2, 5);
                        //AC---CGT
        // journaled str: 01222234
        // host str:      01234567
        //              //ACGTACGT

        SEQAN_ASSERT_EQ(0u, hostToVirtualPosition(journaledStr, 0));
        SEQAN_ASSERT_EQ(1u, hostToVirtualPosition(journaledStr, 1));
        SEQAN_ASSERT_EQ(2u, hostToVirtualPosition(journaledStr, 2));
        SEQAN_ASSERT_EQ(2u, hostToVirtualPosition(journaledStr, 3));
        SEQAN_ASSERT_EQ(2u, hostToVirtualPosition(journaledStr, 4));
        SEQAN_ASSERT_EQ(2u, hostToVirtualPosition(journaledStr, 5));
        SEQAN_ASSERT_EQ(3u, hostToVirtualPosition(journaledStr, 6));
        SEQAN_ASSERT_EQ(4u, hostToVirtualPosition(journaledStr, 7));
    }

    { // test 3: mixed
                            //01234567
        CharString hostStr = "ACGTACGT";
        String<char, Journaled<Alloc<>, TJournaledStringSpec> > journaledStr(hostStr);

        erase(journaledStr, 2, 5);
        insert(journaledStr, 2, "TTT");
        insert(journaledStr, 5, "AAA");
        insert(journaledStr, 11, "CC");

                        //ACTTT---AAACGTCC
        // journaled str: 0123455556789012
        // host str:      01---234---567--
        //              //AC---GTA---CGT--

        SEQAN_ASSERT_EQ(0u, hostToVirtualPosition(journaledStr, 0));
        SEQAN_ASSERT_EQ(1u, hostToVirtualPosition(journaledStr, 1));
        SEQAN_ASSERT_EQ(8u, hostToVirtualPosition(journaledStr, 2));
        SEQAN_ASSERT_EQ(8u, hostToVirtualPosition(journaledStr, 3));
        SEQAN_ASSERT_EQ(8u, hostToVirtualPosition(journaledStr, 4));
        SEQAN_ASSERT_EQ(8u, hostToVirtualPosition(journaledStr, 5));
        SEQAN_ASSERT_EQ(9u, hostToVirtualPosition(journaledStr, 6));
        SEQAN_ASSERT_EQ(10u, hostToVirtualPosition(journaledStr, 7));
        SEQAN_ASSERT_EQ(13u, hostToVirtualPosition(journaledStr, 8));
        SEQAN_ASSERT_EQ(13u, hostToVirtualPosition(journaledStr, 20));
    }

    { // test 4: empty journal string
                            //01234567
        CharString hostStr = "ACGTACGT";
        String<char, Journaled<Alloc<>, TJournaledStringSpec> > journaledStr(hostStr);

                        //ACGTACGT
        // journaled str: 01234567
        // host str:      01234567
        //              //ACGTACGT

        SEQAN_ASSERT_EQ(0u, hostToVirtualPosition(journaledStr, 0));
        SEQAN_ASSERT_EQ(1u, hostToVirtualPosition(journaledStr, 1));
        SEQAN_ASSERT_EQ(2u, hostToVirtualPosition(journaledStr, 2));
        SEQAN_ASSERT_EQ(3u, hostToVirtualPosition(journaledStr, 3));
        SEQAN_ASSERT_EQ(4u, hostToVirtualPosition(journaledStr, 4));
        SEQAN_ASSERT_EQ(5u, hostToVirtualPosition(journaledStr, 5));
        SEQAN_ASSERT_EQ(6u, hostToVirtualPosition(journaledStr, 6));
        SEQAN_ASSERT_EQ(7u, hostToVirtualPosition(journaledStr, 7));
    }
}


template <typename TStringJournalSpec>
void testJournaledStringCopyConstructor(TStringJournalSpec const &)
{
    CharString charStr = "test";
    typedef String<char, Journaled<Alloc<void>, TStringJournalSpec> > TJournaledString;
    TJournaledString * journaledStringPtr = new TJournaledString(charStr);
    TJournaledString journaledString(*journaledStringPtr);
    delete journaledStringPtr;

    insert(journaledString, 1, "XX");

    std::stringstream ss;
    std::stringstream hss;

    ss << journaledString;
    hss << host(journaledString);

    SEQAN_ASSERT_EQ("tXXest", ss.str());
    SEQAN_ASSERT_EQ("test", hss.str());
}

template <typename TStringJournalSpec>
void testJournaledStringBeginEndIteratorStandard(TStringJournalSpec const &)
{
    CharString charStr = "test";
    typedef String<char, Journaled<Alloc<void>, TStringJournalSpec> > TJournaledString;
    TJournaledString journaledString(charStr);
    insert(journaledString, 2, "X");

    typedef typename Iterator<TJournaledString, Standard>::Type TIterator;

    {  // Test pre-increment iteration.
        CharString buffer;
        for (TIterator it = begin(journaledString, Standard()), itend = end(journaledString, Standard()); it != itend; ++it)
            appendValue(buffer, *it);
        SEQAN_ASSERT_EQ("teXst", buffer);
    }
    {  // Test post-increment iteration.
        CharString buffer;
        for (TIterator it = begin(journaledString, Standard()), itend = end(journaledString, Standard()); it != itend; it++)
            appendValue(buffer, *it);
        SEQAN_ASSERT_EQ("teXst", buffer);
    }
}

template <typename TStringJournalSpec>
void testJournaledStringBeginEndIteratorRooted(TStringJournalSpec const &)
{
    CharString charStr = "test";
    typedef String<char, Journaled<Alloc<void>, TStringJournalSpec> > TJournaledString;
    TJournaledString journaledString(charStr);
    insert(journaledString, 2, "X");

    typedef typename Iterator<TJournaledString, Rooted>::Type TIterator;

    {  // Test pre-increment iteration.
        CharString buffer;
        for (TIterator it = begin(journaledString, Rooted()), itend = end(journaledString, Rooted()); it != itend; ++it)
            appendValue(buffer, *it);
        SEQAN_ASSERT_EQ("teXst", buffer);
    }
    {  // Test post-increment iteration.
        CharString buffer;
        for (TIterator it = begin(journaledString, Rooted()), itend = end(journaledString, Rooted()); it != itend; it++)
            appendValue(buffer, *it);
        SEQAN_ASSERT_EQ("teXst", buffer);
    }
}


template <typename TStringJournalSpec>
void testJournaledStringBeginEndConstIterator(TStringJournalSpec const &)
{
    CharString charStr = "test";
    String<char, Journaled<Alloc<void>, TStringJournalSpec> > journaledString(charStr);
    insert(journaledString, 2, "X");

    typedef String<char, Journaled<Alloc<void>, TStringJournalSpec> > TJournaledString;

    // Test with non-const iterator.
    {
        typedef typename Iterator<TJournaledString, Standard>::Type TIterator;
        CharString buffer;
        for (TIterator it = begin(journaledString, Standard()), itend = end(journaledString, Standard()); it != itend; ++it)
            appendValue(buffer, *it);
        SEQAN_ASSERT_EQ("teXst", buffer);
    }
    // Test with const iterator.
    {
        typedef typename Iterator<TJournaledString const, Standard>::Type TIterator;
        String<char, Journaled<Alloc<void>, TStringJournalSpec> > const & constSJ = journaledString;
        CharString buffer;
        for (TIterator it = begin(constSJ, Standard()), itend = end(constSJ, Standard()); it != itend; ++it)
            appendValue(buffer, *it);
        SEQAN_ASSERT_EQ("teXst", buffer);
    }
    // Test comparison of const and non-const iterators.
//     {
//         typedef typename Iterator<String<char, Journaled<Alloc<void>, TStringJournalSpec>, Standard>::Type TNonConstIterator;
//         typedef typename Iterator<String<char, Journaled<Alloc<void>, TStringJournalSpec> const, Standard>::Type TIterator;
//         String<char, Journaled<Alloc<void>, TStringJournalSpec> const & constSJ = journaledString;

//         SEQAN_ASSERT(begin(journaledString, Standard()) == begin(constSJ, Standard()));
//     }
}


template <typename TStringJournalSpec>
void testJournaledStringSubscriptOperatorRandomized(TStringJournalSpec const &)
{
    using namespace seqan;

    const unsigned ASSIGN_COUNT = 2;
    const unsigned LENGTH = 1000;
    const unsigned SEED = 42;
    std::mt19937 rng(SEED);

#define RAND_CHAR() (static_cast<char>(std::uniform_int_distribution<int>('A', 'Z')(rng)))

    // Build random reference and host string.
    String<char> string;
    reserve(string, LENGTH);
    String<char> host;
    reserve(host, LENGTH);
    for (unsigned i = 0; i < LENGTH; ++i) {
        char c = RAND_CHAR();
        appendValue(string, c);
        appendValue(host, c);
    }

    // Create journal string over host and randomly assign infixes in
    // both to fill the tree.
    String<char, Journaled<Alloc<void>, TStringJournalSpec> > journaledString(host);
    for (unsigned i = 0; i < ASSIGN_COUNT; ++i) {
        unsigned begin = 0;
        unsigned end = 0;
        while (begin == end) {
            begin = std::uniform_int_distribution<int>(0, length(string) - 1)(rng);
            end = std::uniform_int_distribution<int>(0, length(string))(rng);
        }
        if (begin > end)
            std::swap(begin, end);
        unsigned len = end - begin;
        String<char> buffer;
        reserve(buffer, len);
        for (unsigned i = 0; i < len; ++i)
            appendValue(buffer, RAND_CHAR());
        replace(string, begin, end, buffer);
        assignInfix(journaledString, begin, end, buffer);
    }

#undef RAND_CHAR

    SEQAN_ASSERT_EQ(length(journaledString), length(string));

    std::stringstream tmp;
    tmp << journaledString;
    CharString str2(tmp.str());
    SEQAN_ASSERT_EQ(str2, string);

    // Now, test the subscript operator.
    for (unsigned i = 0; i < length(journaledString); ++i) {
        SEQAN_ASSERT_EQ_MSG(journaledString[i], string[i], "i = %d", i);
    }
}


// Perform a random insert/edit/delete of the sequence journal.
template <typename TStringJournalSpec>
void testJournaledStringFuzzying(TStringJournalSpec const &)
{
    using namespace seqan;

    const unsigned INITIAL_LENGTH = 100;
    const unsigned NUM_CHANGES = 100;
    const unsigned MAX_INSERT = 100;

    const unsigned SEED = 42;
    std::mt19937 rng(SEED);

#define RAND_CHAR() (static_cast<char>(std::uniform_int_distribution<int>('A', 'Z')(rng)))

    // Build random reference and host string.
    String<char> string;
    reserve(string, INITIAL_LENGTH);
    String<char> host;
    reserve(host, INITIAL_LENGTH);
    for (unsigned i = 0; i < INITIAL_LENGTH; ++i) {
        char c = RAND_CHAR();
        appendValue(string, c);
        appendValue(host, c);
    }

    // Output of initial sequences.
//     std::cout << "reference = " << string << std::endl;
//     std::cout << "host = " << host << std::endl;

    typedef String<char, Journaled<Alloc<void>, TStringJournalSpec> > TJournaledString;

    // Construct sequence journal on host.
    TJournaledString journaledString(host);

//     unsigned nextId = 0;
//     std::cerr << "digraph {" << std::endl;

    // We will use a string stream to test the string result of tmp.
    {
        std::stringstream tmp;
        tmp << journaledString;
//         SEQAN_ASSERT_EQ(string, tmp.str());
//         std::cout << "string = " << string << std::endl;
//         std::cout << "jrnld  = " << tmp.str() << std::endl;
//         std::cout << "  tree = " << journaledString._journalEntries << std::endl;
//         std::cout << "  orig = " << value(journaledString._host) << std::endl;
//         std::cout << "  buff = " << journaledString._insertionBuffer << std::endl;
//         journalTreeToDot(std::cerr, nextId, journaledString._journalEntries);
    }

    size_t expectedLength = length(string);

    for (unsigned i = 0; i < NUM_CHANGES; ++i) {
//         std::cout << "i == " << i << std::endl;
        unsigned changeType = std::uniform_int_distribution<int>(0, 2)(rng);
        if (changeType == 0) {  // edit
            if (length(string) == 0)
                continue;
            unsigned begin = 0;
            unsigned end = 0;
            while (begin == end) {
                begin = std::uniform_int_distribution<int>(0, length(string) - 1)(rng);
                end = std::uniform_int_distribution<int>(0, length(string))(rng);
            }
            if (begin > end)
                std::swap(begin, end);
            unsigned len = end - begin;
            String<char> buffer;
            reserve(buffer, len);
            for (unsigned i = 0; i < len; ++i)
                appendValue(buffer, RAND_CHAR());
            // Perform insert.
//             std::cout << "assignInfix(journaledString, " << begin << ", " << end << ", \"" << buffer << "\")" << std::endl;
//             std::cout << "assignInfix(journaledString, " << begin << ", " << end << ", buffer, len(buffer) == " << length(buffer) << ")" << std::endl;
            replace(string, begin, end, buffer);
//             std::cout << "pre assign infix " << length(journaledString) << std::endl;
            assignInfix(journaledString, begin, end, buffer);
//             std::cout << "post assign infix " << length(journaledString) << std::endl;
        } else if (changeType == 1) {  // insert
            unsigned begin = 0;
            unsigned len = 0;
            while (len == 0) {
                if (length(string) == 0)
                    begin = 0;
                else
                    begin = std::uniform_int_distribution<int>(0, length(string) - 1)(rng);
                len = std::uniform_int_distribution<int>(0, MAX_INSERT)(rng);
            }
            expectedLength += len;
            String<char> buffer;
            reserve(buffer, len);
            for (unsigned i = 0; i < len; ++i)
                appendValue(buffer, RAND_CHAR());
            // Perform insert.
//             std::cout << "insert(journaledString, " << begin << ", \"" << buffer << "\")" << std::endl;
//             std::cout << "insert(journaledString, " << begin << ", buffer)" << std::endl;
            replace(string, begin, begin, buffer);
            insert(journaledString, begin, buffer);
        } else if (changeType == 2) {  // delete
            if (length(string) == 0)
                continue;
            //std::cerr << "length(string) == " << length(string) << std::endl;
            //std::cerr << "string == " << string << std::endl;
            //std::cerr << "journal string== " << journaledString << std::endl;
            //std::cerr << "journal string== " << journaledString._journalEntries << std::endl;
            unsigned begin = 0;
            unsigned end = 0;
            while (begin == end) {
                begin = std::uniform_int_distribution<int>(0, length(string) - 1)(rng);
                end = std::uniform_int_distribution<int>(0, length(string))(rng);
            }
            if (begin > end)
                std::swap(begin, end);
            expectedLength -= (end - begin);
            // Perform erase.
//             std::stringstream tmp;
//             tmp << journaledString;
//             std::cout << ",---" << std::endl;
//             std::cout << "| string = " << string << std::endl;
//             std::cout << "| jrnld  = " << tmp.str() << std::endl;
//             std::cout << "|   tree = " << journaledString._journalEntries << std::endl;
//             std::cout << "|   orig = " << value(journaledString._host) << std::endl;
//             std::cout << "|   buff = " << journaledString._insertionBuffer << std::endl;
//             std::cout << "erase(journaledString, " << begin << ", " << end << ")" << std::endl;
//             std::cout << "`---" << std::endl;
//             std::cout << journaledString._journalEntries << std::endl;
            erase(string, begin, end);
            erase(journaledString, begin, end);
//             std::cout << journaledString._journalEntries << std::endl;
        } else {
            SEQAN_ASSERT_FAIL("Invalid change type.");
        }

        {
            // Check via stream operator<< into stringstream.
            std::stringstream tmp;
            tmp << journaledString;
            SEQAN_ASSERT_EQ(expectedLength, length(tmp.str()));
            SEQAN_ASSERT_EQ(expectedLength, length(string));
            SEQAN_ASSERT_EQ(string, tmp.str());
//             std::cout << "string = " << string << std::endl << "tmp.str() = " << tmp.str() << std::endl;
            // Check via iterator on the journal string.
            std::string buffer;
            reserve(buffer, length(journaledString));
            typedef typename Iterator<TJournaledString, Standard>::Type TIterator;
//             std::cout << journaledString._journalEntries << std::endl;
            for (TIterator it = begin(journaledString), itend = end(journaledString, Standard()); it != itend; ++it) {
                appendValue(buffer, *it);
            }
//             std::cout << "buffer = " << buffer << std::endl;
            SEQAN_ASSERT_EQ(expectedLength, length(buffer));
            SEQAN_ASSERT_EQ(buffer, tmp.str());
        }

        {
            // Check operator+ and operator+= on sequence journal iterators;
            typedef typename Iterator<TJournaledString, Standard>::Type TJournaledStringIterator;
            typedef typename Iterator<String<char>, Standard>::Type TCharStringIterator;

            TJournaledStringIterator sjIt = begin(journaledString);
            TCharStringIterator csIt = begin(string);
            size_t remaining = length(journaledString);

            while (remaining > 1) {
                SEQAN_ASSERT(csIt != end(string) - 1);
                size_t len = std::uniform_int_distribution<int>(0, length(remaining))(rng);
                remaining -= len;
                if (remaining == 0)
                    break;
                SEQAN_ASSERT_EQ(*(sjIt + len), *(csIt + len));
                sjIt += len;
                csIt += len;
                SEQAN_ASSERT_EQ(*sjIt, *csIt);
            }
        }

        {
            // Check operator- on sequence journal iterators;
            typedef typename Iterator<TJournaledString, Standard>::Type TJournaledStringIterator;
            typedef typename Iterator<String<char>, Standard>::Type TCharStringIterator;

            TJournaledStringIterator sjIt = begin(journaledString);
            TJournaledStringIterator sjIt2 = sjIt;
            TCharStringIterator csIt = begin(string);
            TCharStringIterator csIt2 = csIt;
            size_t remaining = length(journaledString);

            while (remaining > 1) {
                SEQAN_ASSERT(csIt != end(string) - 1);
                size_t len = std::uniform_int_distribution<int>(0, length(remaining))(rng);
                remaining -= len;
                if (remaining == 0)
                    break;
                SEQAN_ASSERT_EQ(*(sjIt + len), *(csIt + len));
                sjIt += len;
                csIt += len;
                SEQAN_ASSERT_EQ(*sjIt, *csIt);

                SEQAN_ASSERT_EQ(csIt2 - csIt, sjIt2 - sjIt);
                SEQAN_ASSERT_EQ(csIt - csIt2, sjIt - sjIt2);
                csIt2 = csIt;
                sjIt2 = sjIt;
            }
        }
    }

#undef RAND_CHAR
}


template <typename TStringJournalSpec>
void testJournaledStringSegmentsReadOnly(TStringJournalSpec const &)
{
    typedef String<char, Journaled<Alloc<void>, TStringJournalSpec> > TJournaledString;
    typedef typename Prefix<TJournaledString>::Type TPrefix;
    typedef typename Suffix<TJournaledString>::Type TSuffix;
    typedef typename Infix<TJournaledString>::Type TInfix;

    CharString charStr = "test";
    TJournaledString journaledString;
    setHost(journaledString, charStr);
    insert(journaledString, 2, "XX");

    // Prefixes.
    {
        TPrefix prefix1 = prefix(journaledString, 3);
        SEQAN_ASSERT(prefix1 == CharString("teX"));
        TPrefix prefix2(journaledString, 3);
        SEQAN_ASSERT(prefix2 == CharString("teX"));
    }
    // Suffixes.
    {
        TSuffix suffix1 = suffix(journaledString, 3);
        SEQAN_ASSERT(suffix1 == CharString("Xst"));
        TSuffix suffix2(journaledString, 3);
        SEQAN_ASSERT(suffix2 == CharString("Xst"));

    }
    // Infixes.
    {
        TInfix infix1 = infix(journaledString, 1, 5);
        SEQAN_ASSERT(infix1 == CharString("eXXs"));
        TInfix infix2(journaledString, 1, 5);
        SEQAN_ASSERT(infix2 == CharString("eXXs"));
    }
}


template <typename TStringJournalSpec>
void testJournaledStringReplace(TStringJournalSpec const &)
{
    typedef String<char, Journaled<Alloc<void>, TStringJournalSpec> > TJournaledString;

    CharString charStr = "test";

    // Replace prefix.
    {
        TJournaledString journaledString(charStr);
        insert(journaledString, 2, "XX");
        replace(journaledString, 0, 3, "ABCD");
        SEQAN_ASSERT_EQ(journaledString, "ABCDXst");
        SEQAN_ASSERT_EQ(charStr, "test");
    }

    // Replace suffix.
    {
        TJournaledString journaledString(charStr);
        insert(journaledString, 2, "XX");
        replace(journaledString, 3, 6, "ABCD");
        SEQAN_ASSERT_EQ(journaledString, "teXABCD");
        SEQAN_ASSERT_EQ(charStr, "test");
    }

    // Infixes.
    {
        TJournaledString journaledString(charStr);
        insert(journaledString, 2, "XX");
        replace(journaledString, 1, 5, "ABCD");
        SEQAN_ASSERT_EQ(journaledString, "tABCDt");
        SEQAN_ASSERT_EQ(charStr, "test");
    }
}

// Tag: SortedArray()

SEQAN_DEFINE_TEST(test_sequence_journaled_sorted_array_assign)
{
    testJournaledStringAssign(SortedArray());
}

SEQAN_DEFINE_TEST(test_sequence_journaled_sorted_array_set)
{
    testJournaledStringSet(SortedArray());
}

SEQAN_DEFINE_TEST(test_sequence_journaled_sorted_array_host)
{
    testJournaledStringHost(SortedArray());
}

SEQAN_DEFINE_TEST(test_sequence_journaled_sorted_array_clear)
{
    testJournaledStringClear(SortedArray());
}

SEQAN_DEFINE_TEST(test_sequence_journaled_sorted_array_empty)
{
    testJournaledStringEmpty(SortedArray());
}

SEQAN_DEFINE_TEST(test_sequence_journaled_sorted_array_erase_position)
{
    testJournaledStringErasePosition(SortedArray());
}

SEQAN_DEFINE_TEST(test_sequence_journaled_sorted_array_erase_begin_end)
{
    testJournaledStringEraseBeginEnd(SortedArray());
}

SEQAN_DEFINE_TEST(test_sequence_journaled_sorted_array_insert)
{
    testJournaledStringInsert(SortedArray());
}

SEQAN_DEFINE_TEST(test_sequence_journaled_sorted_array_insert_value)
{
    testJournaledStringInsertValue(SortedArray());
}

SEQAN_DEFINE_TEST(test_sequence_journaled_sorted_array_assign_value)
{
    testJournaledStringAssignValue(SortedArray());
}

SEQAN_DEFINE_TEST(test_sequence_journaled_sorted_array_subscript_operator)
{
    testJournaledStringSubscriptOperator(SortedArray());
}

SEQAN_DEFINE_TEST(test_sequence_journaled_sorted_array_assign_infix)
{
    testJournaledStringAssignInfix(SortedArray());
}

SEQAN_DEFINE_TEST(test_sequence_journaled_sorted_array_length)
{
    testJournaledStringLength(SortedArray());
}

SEQAN_DEFINE_TEST(test_sequence_journaled_sorted_array_virtual_to_host_position)
{
    testJournaledStringVirtualToHostPosition<SortedArray>();
}

SEQAN_DEFINE_TEST(test_sequence_journaled_sorted_array_host_to_virtual_position)
{
    testJournaledStringHostToVirtualPosition<SortedArray>();
}

SEQAN_DEFINE_TEST(test_sequence_journaled_sorted_array_reset)
{
    testJournaledStringReset(seqan::SortedArray());
}

SEQAN_DEFINE_TEST(test_sequence_journaled_sorted_array_copy_constructor)
{
    testJournaledStringCopyConstructor(SortedArray());
}

SEQAN_DEFINE_TEST(test_sequence_journaled_sorted_array_begin_end_iterator)
{
    testJournaledStringBeginEndIteratorStandard(SortedArray());
    testJournaledStringBeginEndIteratorRooted(SortedArray());
}

SEQAN_DEFINE_TEST(test_sequence_journaled_sorted_array_begin_end_const_iterator)
{
    testJournaledStringBeginEndConstIterator(SortedArray());
}

SEQAN_DEFINE_TEST(test_sequence_journaled_sorted_array_subscript_operator_randomized)
{
    testJournaledStringSubscriptOperatorRandomized(SortedArray());
}

SEQAN_DEFINE_TEST(test_sequence_journaled_sorted_array_fuzzying)
{
    testJournaledStringFuzzying(SortedArray());
}

SEQAN_DEFINE_TEST(test_sequence_journaled_sorted_array_segments_read_only)
{
    testJournaledStringSegmentsReadOnly(SortedArray());
}

SEQAN_DEFINE_TEST(test_sequence_journaled_sorted_array_segments_read_write)
{
    testJournaledStringReplace(SortedArray());
}

SEQAN_DEFINE_TEST(test_sequence_journaled_sorted_array_flatten)
{
    testJournalStringFlatten(SortedArray());
}

#endif  // TEST_SEQUENCE_JOURNALED_TEST_SEQUENCE_JOURNALED_H_
