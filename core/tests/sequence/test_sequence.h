// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2012, Knut Reinert, FU Berlin
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
// Author: Jochen Singer <jochen.singer@fu-berlin.de>
// ==========================================================================
// This file contains functions to test the functionality of the sequence
// module.
// ==========================================================================

#ifndef CORE_TESTS_SEQUENCE_TEST_SEQUENCE_H_
#define CORE_TESTS_SEQUENCE_TEST_SEQUENCE_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

// --------------------------------------------------------------------------
// CountingChar is used to test sequences of non simple data types.
// --------------------------------------------------------------------------

struct CountingChar
{
    char value;                     // value of the object
    static unsigned numConstruct;   // number of constructor calls
    static unsigned numDeconstruct; // number of destructor calls

    CountingChar()
    {
        numConstruct += 1;
    }

    CountingChar(char const & value) : value(value)
    {
        numConstruct += 1;
    }

    CountingChar(CountingChar const & other) : value(other.value)
    {
        numConstruct += 1;
    }

    ~CountingChar()
    {
        numDeconstruct += 1;
    }

    static void clear()
    {
        numConstruct = 0;
        numDeconstruct = 0;
    }

bool operator==(CountingChar const & other) const
    {
        return value == other.value;
    }

    bool operator>(CountingChar const & other) const
    {
        return value > other.value;
    }

    bool operator<(CountingChar const & other) const
    {
        return value < other.value;
    }
};

template <typename TStream>
inline TStream & operator<<(TStream & stream, CountingChar const & countingChar)
{
    stream << countingChar.value;

    return stream;
}

template <typename TStream, typename TSpec>
inline TStream & operator<<(TStream & stream, seqan::String<CountingChar, TSpec> const & string)
{
    for (unsigned i = 0; i < length(string); ++i)
        stream << string[i];

    return stream;
}

unsigned CountingChar::numConstruct = 0;
unsigned CountingChar::numDeconstruct = 0;

// --------------------------------------------------------------------------
// Generic String Tests
// --------------------------------------------------------------------------

// Test whether sequences are copy constructible.
template <typename TString>
void testSequenceCopyConstructible(TString & /*Tag*/)
{
    using namespace seqan;

    TString string1 = "ACGCTAGCAT";
    TString string2(string1);
    SEQAN_ASSERT_EQ(string1, string2);
}

// Test whether sequences are default constructible.
template <typename TString>
void testSequenceDefaultConstructible(TString & /*Tag*/)
{
    using namespace seqan;

    TString string;
    SEQAN_ASSERT_EQ(begin(string), end(string));
}

// Test operator<().
template <typename TString>
void testSequenceLess(TString & /*Tag*/)
{
    using namespace seqan;

    // Nothing is smaller than something.
    {
        TString string1;
        TString string2 = "A";
        SEQAN_ASSERT_EQ(string1 < string2, true);
    }

    // Equal is not smaller
    {
        TString string1 = "A";
        TString string2 = "A";
        SEQAN_ASSERT_EQ(string1 < string2, false);
    }

    // Sequences with lex. smaller characters and equal length are smaller.
    {
        TString string1 = "A";
        TString string2 = "C";
        SEQAN_ASSERT_EQ(string1 < string2, true);
    }

    // Sequences with lex. smaller characters but larger length are smaller.
    {
        TString string1 = "AA";
        TString string2 = "C";
        SEQAN_ASSERT_EQ(string1 < string2, true);
    }

    // Sequences with equal characters but smaller length are smaller.
    {
        TString string1 = "AA";
        TString string2 = "AAA";
        SEQAN_ASSERT_EQ(string1 < string2, true);
    }
}

// Test operator<=().
template <typename TString>
void testSequenceLessEqual(TString & /*Tag*/)
{
    using namespace seqan;

    // Nothing is equal to nothing.
    {
        TString string1;
        TString string2;
        SEQAN_ASSERT_EQ(string1 <= string2, true);
    }

    // Nothing is smaller than something.
    {
        TString string1;
        TString string2 = "A";
        SEQAN_ASSERT_EQ(string1 <= string2, true);
    }

    // Sequences of equal characters and length are equal.
    {
        TString string1 = "A";
        TString string2 = "A";
        SEQAN_ASSERT_EQ(string1 <= string2, true);
    }

    // Sequences with lex. smaller characters and equal length are smaller.
    {
        TString string1 = "A";
        TString string2 = "C";
        SEQAN_ASSERT_EQ(string1 <= string2, true);
    }

    // Sequences with lex. smaller characters but larger length are smaller.
    {
        TString string1 = "AA";
        TString string2 = "C";
        SEQAN_ASSERT_EQ(string1 <= string2, true);
    }

    // Sequences with equal characters but smaller length are smaller.
    {
        TString string1 = "AA";
        TString string2 = "AAA";
        SEQAN_ASSERT_EQ(string1 <= string2, true);
    }
}

// Test operator>().
template <typename TString>
void testSequenceGreater(TString & /*Tag*/)
{
    using namespace seqan;

    // Something is greater than nothing.
    {
           TString string1 = "A";
        TString string2;
        SEQAN_ASSERT_EQ(string1 > string2, true);
    }

    // Equal is not greater
    {
        TString string1 = "A";
        TString string2 = "A";
        SEQAN_ASSERT_EQ(string1 > string2, false);
    }

    // Sequences with lex. greater characters and equal length are greater.
    {
        TString string1 = "C";
        TString string2 = "A";
        SEQAN_ASSERT_EQ(string1 > string2, true);
    }

    // Sequences with equal characters but larger length are larger.
    {
        TString string1 = "AA";
        TString string2 = "A";
        SEQAN_ASSERT_EQ(string1 > string2, true);
    }

    // Sequences with lex. greater characters are greater.
    {
        TString string1 = "C";
        TString string2 = "AA";
        SEQAN_ASSERT_EQ(string1 > string2, true);
    }
}

// Test operator>=().
template <typename TString>
void testSequenceGreaterEqual(TString & /*Tag*/)
{
    using namespace seqan;

    // Nothing is equal to nothing.
    {
        TString string1;
        TString string2;
        SEQAN_ASSERT_EQ(string1 >= string2, true);
    }

    // Something is greater than nothing.
    {
        TString string1 = "A";
        TString string2;
        SEQAN_ASSERT_EQ(string1 >= string2, true);
    }

    // Sequences of equal characters and length are equal.
    {
        TString string1 = "A";
        TString string2 = "A";
        SEQAN_ASSERT_EQ(string1 >= string2, true);
    }

    // Sequences with lex. greater characters and equal length are greater.
    {
        TString string1 = "C";
        TString string2 = "A";
        SEQAN_ASSERT_EQ(string1 >= string2, true);
    }

    // Sequences with equal characters but larger length are greater.
    {
        TString string1 = "AA";
        TString string2 = "A";
        SEQAN_ASSERT_EQ(string1 >= string2, true);
    }

    // Sequences with lex. greater characters are greater.
    {
        TString string1 = "C";
        TString string2 = "AA";
        SEQAN_ASSERT_EQ(string1 > string2, true);
    }
}

// Test operator==().
template <typename TString>
void testSequenceEqual(TString & /*Tag*/)
{
    using namespace seqan;

    // Nothing is equal to nothing.
    {
        TString string1;
        TString string2;
        SEQAN_ASSERT_EQ(string1 == string2, true);
    }

    // Something is greater than nothing.
    {
        TString string1 = "A";
        TString string2;
        SEQAN_ASSERT_EQ(string1 == string2, false);
    }
    {
        TString string1;
        TString string2 = "A";
        SEQAN_ASSERT_EQ(string1 == string2, false);
    }

    // Sequences of equal characters and length are equal.
    {
        TString string1 = "A";
        TString string2 = "A";
        SEQAN_ASSERT_EQ(string1 == string2, true);
    }

    // Sequences of equal characters but different length are not equal.
    {
        TString string1 = "A";
        TString string2 = "AA";
        SEQAN_ASSERT_EQ(string1 == string2, false);
    }
    {
        TString string1 = "AA";
        TString string2 = "A";
        SEQAN_ASSERT_EQ(string1 == string2, false);
    }

    // Sequences of different characters are not equal.
    {
        TString string1 = "A";
        TString string2 = "C";
        SEQAN_ASSERT_EQ(string1 == string2, false);
    }
    {
        TString string1 = "C";
        TString string2 = "A";
        SEQAN_ASSERT_EQ(string1 == string2, false);
    }
}

// Test operator!=().
template <typename TString>
void testSequenceUnequal(TString & /*Tag*/)
{
    using namespace seqan;

    // Nothing is equal to nothing.
    {
        TString string1;
            TString string2;
            SEQAN_ASSERT_EQ(string1 != string2, false);
    }

    // Something is greater than nothing.
    {
        TString string1 = "A";
        TString string2;
        SEQAN_ASSERT_EQ(string1 != string2, true);
    }
    {
        TString string1;
        TString string2 = "A";
        SEQAN_ASSERT_EQ(string1 != string2, true);
    }

    // Sequences of equal characters and length are equal.
    {
        TString string1 = "A";
        TString string2 = "A";
        SEQAN_ASSERT_EQ(string1 != string2, false);
    }

    // Sequences of equal characters but different length are not equal.
    {
        TString string1 = "A";
        TString string2 = "AA";
        SEQAN_ASSERT_EQ(string1 != string2, true);
    }
    {
        TString string1 = "AA";
        TString string2 = "A";
        SEQAN_ASSERT_EQ(string1 != string2, true);
    }

    // Sequences of different characters are not equal.
    {
        TString string1 = "A";
        TString string2 = "C";
        SEQAN_ASSERT_EQ(string1 != string2, true);
    }
    {
        TString string1 = "C";
        TString string2 = "A";
        SEQAN_ASSERT_EQ(string1 != string2, true);
    }
}

// Test of append().
template <typename TString>
void testSequenceAppend(TString & /*Tag*/)
{
    using namespace seqan;

    // Test the append function on empty strings
    {
        TString string1 = "";
        TString string2 = "";
        append(string1, string2);
        SEQAN_ASSERT_EQ(string1, "");
    }

    // Test the append function on one empty and one non empty string
    {
        TString string1 = "ACGTACGTACGT";
        TString string2 = "";
        append(string1, string2);
        SEQAN_ASSERT_EQ(string1, "ACGTACGTACGT");
    }

    // Test the append function on one empty and one non empty string
    {
        TString string1 = "";
        TString string2 = "TTGGATTAACC";
        append(string1, string2);
        SEQAN_ASSERT_EQ(string1, "TTGGATTAACC");
    }

    // Test the append function on two non empty strings.
    {
        TString string1 = "ACGTACGTACGT";
        TString string2 = "TTGGATTAACCC";
        append(string1, string2);
        SEQAN_ASSERT_EQ(string1, "ACGTACGTACGTTTGGATTAACCC");
    }
}

// Test of appendValue().
template <typename TString>
void testSequenceAppendValue(TString & /*Tag*/)
{
    using namespace seqan;

    typedef typename Value<TString>::Type TValue;

    TString string = "";

    // Test the appendValue function
    TValue value = 'A';
    appendValue(string, value);
    SEQAN_ASSERT_EQ(string, "A");

    // Test the appendValue function
    appendValue(string, 'A');
    SEQAN_ASSERT_EQ(string, "AA");
}

// Test whether sequences are assignable.
template <typename TString>
void testSequenceAssignable(TString & /*Tag*/)
{
    using namespace seqan;

    // Test on an empty string.
    TString string1 = "";

    // Test assign() function.
    {
        TString string2;
        assign(string2, string1);
        SEQAN_ASSERT_EQ(string1, string2);
    }

    // Test operator=().
    {
        TString string2;  // Separate definition and assignment on purpose.
        string2 = string1;
        SEQAN_ASSERT_EQ(string1, string2);
    }

    // Test the basic concept on a non empty string.
    string1 = "ACGTACGTACGT";

    // Test the assign() function.
    {
        TString string2;
        assign(string2, string1);
        SEQAN_ASSERT_EQ(string1, string2);
    }

    // Test operator=().
    {
        TString string2;  // Separate definition and assignment on purpose.
        string2 = string1;
        SEQAN_ASSERT_EQ(string1, string2);
    }
}

// Test of assignValue().
template <typename TString>
void testSequenceAssignValue(TString & /*Tag*/)
{
    using namespace seqan;

    typedef typename Value<TString>::Type TValue;

    TString string = "";

    resize(string, 2);
    assignValue(string, 1, TValue('G'));
    SEQAN_ASSERT_EQ(string[1], TValue('G'));

}

// We need two back() tests, since back() returns a reference or a copy
// (depending on TString). We check whether we can modify the reference.
// Test of back() for non const strings.
template <typename TString>
void testSequenceBack(TString & /*Tag*/)
{
    using namespace seqan;

    typedef typename Value<TString>::Type TValue;
    TString string = "ACGT";

    // val is a reference in contrast to the const version of back().
    TValue & val = back(string);
    val = 'A';
    SEQAN_ASSERT_EQ(val, string[length(string) - 1]);
}

// Test of back() for const strings.
template <typename TString>
void testSequenceBack(TString const & /*Tag*/)
{
    using namespace seqan;

    typedef typename Value<TString>::Type TValue;
    TString const string = "ACGT";

    // val is not a reference in contrast to the non const version of back().
    TValue val = back(string);
    SEQAN_ASSERT_EQ(val, string[length(string) - 1]);
}

// Test of begin().
template <typename TString>
void testSequenceBegin(TString & /*Tag*/)
{
    using namespace seqan;

    typedef typename Value<TString>::Type TValue;

    TString string = "ACGT";
    SEQAN_ASSERT_EQ(*begin(string), TValue('A'));
    SEQAN_ASSERT_EQ(*begin(string, Standard()), TValue('A'));
    SEQAN_ASSERT_EQ(*begin(string, Rooted()), TValue('A'));
}

// Test of beginPosition().
template <typename TString>
void testSequenceBeginPosition(TString & /*Tag*/)
{
    using namespace seqan;

    // Test on an empty string.
    TString string1 = "";
    SEQAN_ASSERT_EQ(beginPosition(string1), 0u);

    // Test on a non empty string.
    TString string2 = "ACGT";
    SEQAN_ASSERT_EQ(beginPosition(string2), 0u);
}

// Test of capacity().
template <typename TString>
void testSequenceCapacity(TString & /*Tag*/)
{
    using namespace seqan;

    // Test on an empty string.
    TString string1 = "";
    SEQAN_ASSERT_GEQ(capacity(string1), length(string1));

    // Test on a non empty string.
    TString string2 = "ACGTACGTACGT";
    SEQAN_ASSERT_GEQ(capacity(string2), length(string2));
}

// Test of clear().
template <typename TString>
void testSequenceClear(TString & /*Tag*/)
{
    using namespace seqan;

    // Test on an empty string.
    TString string = "";
    clear(string);
    SEQAN_ASSERT_EQ(begin(string), end(string));
    SEQAN_ASSERT_EQ(capacity(string), 0u);

    // Test on a non empty string.
    string = "ACGTACGTACGT";
    clear(string);
    SEQAN_ASSERT_EQ(string, TString());
}

// Test of end().
template <typename TString>
void testSequenceEnd(TString & /*Tag*/)
{
    using namespace seqan;

    typedef typename Value<TString>::Type TValue;

    TString string = "ACGT";
    typename Iterator<TString>::Type iter = end(string);
    --iter;
    SEQAN_ASSERT_EQ(*iter, TValue('T'));

    typename Iterator<TString, Standard>::Type standardIter = end(string, Standard());
    --standardIter;
    SEQAN_ASSERT_EQ(*standardIter, TValue('T'));

    typename Iterator<TString, Rooted>::Type rootedIter = end(string, Rooted());
    --rootedIter;
    SEQAN_ASSERT_EQ(*rootedIter, TValue('T'));
}

// Test of endPosition().
template <typename TString>
void testSequenceEndPosition(TString & /*Tag*/)
{
    using namespace seqan;

    // Test on an empty string.
    TString string1 = "";
    SEQAN_ASSERT_EQ(endPosition(string1), length(string1));

    // Test on a non empty string.
    TString string2 = "ACGT";
    SEQAN_ASSERT_EQ(endPosition(string2), length(string2));
}

// Test of erase().
template <typename TString>
void testSequenceErase(TString & /*Tag*/)
{
    using namespace seqan;

    // Test on an empty string.
    TString string = "";

    // TODO (singer): erasing an empty string is not possible.
    // Error message: std::bad_alloc.
    //erase(string, 0);
    //SEQAN_ASSERT_EQ(string, "");

    //erase(string1, 0, 0);
    //SEQAN_ASSERT_EQ(string, "");

    // Test on a non empty string.
    string = "ACGTACGTACGT";
    erase(string, 1);
    SEQAN_ASSERT_EQ(string, "AGTACGTACGT");

    erase(string, 2, 5);
    SEQAN_ASSERT_EQ(string, "AGGTACGT");
}

// Test of eraseBack().
template <typename TString>
void testSequenceEraseBack(TString & /*Tag*/)
{
    using namespace seqan;

    // Test on an empty string.
    TString string = "";

    // TODO (singer): eraseBack() on an empty string is not possible.
    // Error message: "String must have more than 0 characters in eraseBack()!".
    // If erase() should work on empty strings than eraseBack() as well ???
    // eraseBack(string);
    // SEQAN_ASSERT_EQ(string, "");

    // Test on a non empty string.
    string = "ACGTACGTACGT";
    eraseBack(string);
    SEQAN_ASSERT_EQ(string, "ACGTACGTACG");

    string = "A";
    eraseBack(string);
    SEQAN_ASSERT_EQ(string, "");
}

// Test of front() for non const strings.
template <typename TString>
void testSequenceFront(TString & /*Tag*/)
{
    using namespace seqan;

    typedef typename Value<TString>::Type TValue;
    TString string = "ACGT";

    // val is a reference in contrast to the const version of front()
    TValue & val = front(string);
    val = 'A';
    SEQAN_ASSERT_EQ(val, string[0]);
}

// Test of front() for const strings.
template <typename TString>
void testSequenceFront(TString const & /*Tag*/)
{
    using namespace seqan;

    typedef typename Value<TString>::Type TValue;
    TString string = "ACGT";

    // val is not a reference in contrast to the non const version of front()
    TValue val = front(string);
    SEQAN_ASSERT_EQ(val, string[0]);
}

// Test of getValue().
template <typename TString>
void testSequenceGetValue(TString & /*Tag*/)
{
    using namespace seqan;

    typedef typename Value<TString>::Type TValue;

     // In contrast to value(), getValue() does not return a reference but a copy.
     // We test this using the variable value_.
    TString string = "ACGT";
    TValue dummy_ = 'T';
    TValue & value_ = dummy_;
    SEQAN_ASSERT_EQ(value_, TValue('T'));

    value_ = getValue(string, 0);
    SEQAN_ASSERT_EQ(value_, TValue('A'));

    value_ = 'T';
    SEQAN_ASSERT_EQ(getValue(string, 0), TValue('A'));
}

// Test of insert().
template <typename TString>
void testSequenceInsert(TString & /*Tag*/)
{
    using namespace seqan;

    //typedef typename Value<TString>::Type TValue;

    // Test of inserting an empty string.
    TString string1 = "";
    TString string2 = "";
    insert(string1, 0, string2);
    SEQAN_ASSERT_EQ(string1, "");

    // Test of inserting an string.
    string1 = "A";
    string2 = "ACGT";
    insert(string1, 0, string2);
    SEQAN_ASSERT_EQ(string1, "ACGTA");
}

// Test of insertValue().
template <typename TString>
void testSequenceInsertValue(TString & /*Tag*/)
{
    using namespace seqan;

    //typedef typename Value<TString>::Type TValue;

    // Test of inserting into an empty string.
    TString string = "";
    insert(string, 0, 'A');
    SEQAN_ASSERT_EQ(string, "A");

    // Test of inserting into a non empty string.
    insert(string, 0, 'C');
    SEQAN_ASSERT_EQ(string, "CA");
}

// Test of iter().
template <typename TString>
void testSequenceIter(TString & /*Tag*/)
{
    using namespace seqan;

    typedef typename Iterator<TString>::Type TIterator;
    typedef typename Iterator<TString, Standard>::Type TStandardIterator;
    typedef typename Iterator<TString, Rooted>::Type TRootedIterator;

    // Test on an empty string.
    {
        TString string;
        TIterator iterator = iter(string, 0);
        TStandardIterator standardIterator = iter(string, 0);
        TRootedIterator rootedIterator = iter(string, 0);
        SEQAN_ASSERT_EQ(iterator, begin(string));
        SEQAN_ASSERT_EQ(standardIterator, begin(string));
        SEQAN_ASSERT_EQ(rootedIterator, begin(string));
    }

    // Test on a non empty string.
    {
        TString string = "A";
        TIterator iterator = iter(string, 0);
        TStandardIterator standardIterator = iter(string, 0);
        TRootedIterator rootedIterator = iter(string, 0);
        SEQAN_ASSERT_EQ(getValue(iterator), getValue(string, 0));
        SEQAN_ASSERT_EQ(getValue(standardIterator), getValue(string, 0));
        SEQAN_ASSERT_EQ(getValue(rootedIterator), getValue(string, 0));
    }

    // Test on a non empty string.
    {
        TString string = "ACGT";
        TIterator iterator = iter(string, 3);
        TStandardIterator standardIterator = iter(string, 3);
        TRootedIterator rootedIterator = iter(string, 3);
        SEQAN_ASSERT_EQ(getValue(iterator), getValue(string, 3));
        SEQAN_ASSERT_EQ(getValue(standardIterator), getValue(string, 3));
        SEQAN_ASSERT_EQ(getValue(rootedIterator), getValue(string, 3));
    }
}

// Test of length().
template <typename TString>
void testSequenceLength(TString & /*Tag*/)
{
    using namespace seqan;

    // Test on an empty string.
    TString string1;
    SEQAN_ASSERT_EQ(length(string1), 0u);

    // Test on a non empty string.
    TString string2 = "CGTACGTATC";
    SEQAN_ASSERT_EQ(length(string2), 10u);
}

// Test of value().
template <typename TString>
void testSequenceMoveValue(TString & /*Tag*/)
{
    using namespace seqan;

    TString string = "";

    resize(string, 2);
    moveValue(string, 1, 'G');
    SEQAN_ASSERT_EQ(string[1], 'G');
}

// Test of replace().
template <typename TString>
void testSequenceReplace(TString & /*Tag*/)
{
    using namespace seqan;

    TString string1 = "";
    TString string2 = "";

    // This is problematic according to the documentation.
    // 0 can be a position or an iterator causing compiler errors.
    replace(string1, 0, 0, string2);
    SEQAN_ASSERT_EQ(string1, "");

    string1 = "ACACACAC";
    string2 = "";

    replace(string1, 4, 4, string2);
    SEQAN_ASSERT_EQ(string1, "ACACACAC");

    string1 = "";
    string2 = "GTGTGTGT";

    replace(string1, 0, 0, string2);
    SEQAN_ASSERT_EQ(string1, "GTGTGTGT");

    string1 = "ACACACAC";
    string2 = "GTGTGTGT";

    replace(string1, 4, 4, string2);
    SEQAN_ASSERT_EQ(string1, "ACACGTGTGTGTACAC");
}

// Test of reserve().
template <typename TString>
void testSequenceReserve(TString & /*Tag*/)
{
    using namespace seqan;

    TString string = "";

    reserve(string, 0);
    SEQAN_ASSERT_EQ(capacity(string), 0u);

    reserve(string, 1000);
    SEQAN_ASSERT_GEQ(capacity(string), 10u);

    // If the the new capacity is smaller than the current one
    // the new capacity must be larger or equal to the current length.
    reserve(string, 1);
    SEQAN_ASSERT_GEQ(capacity(string), length(string));
}

// Test of resize().
template <typename TString>
void testSequenceResize(TString & /*Tag*/)
{
    using namespace seqan;

    typedef typename Value<TString>::Type TValue;

    TString string = "";

    resize(string, 0);
    SEQAN_ASSERT_EQ(length(string), 0u);

    resize(string, 10);
    SEQAN_ASSERT_EQ(length(string), 10u);
    // TODO (singer): resize should initialize newly allocated memory,
    // which it does not at the moment!
    //SEQAN_ASSERT_EQ(string[0], TValue());
    //SEQAN_ASSERT_EQ(string[0], TValue());

    resize(string, 0);
    SEQAN_ASSERT_EQ(length(string), 0u);

    resize(string, 10, TValue('C'));
    SEQAN_ASSERT_EQ(string[0], TValue('C'));
    SEQAN_ASSERT_EQ(string[0], TValue('C'));
}

// Test of swap().
template <typename TString>
void testSequenceSwap(TString & /*Tag*/)
{
    using namespace seqan;

    TString string1 = "";
    TString string2 = "";
    TString string3 = string1;
    TString string4 = string2;

    swap(string1, string2);
    SEQAN_ASSERT_EQ(string1, string4);
    SEQAN_ASSERT_EQ(string2, string3);

    string1 = "ACGT";
    string2 = "";
    string3 = string1;
    string4 = string2;

    swap(string1, string2);
    SEQAN_ASSERT_EQ(string1, string4);
    SEQAN_ASSERT_EQ(string2, string3);

    string1 = "";
    string2 = "ACGT";
    string3 = string1;
    string4 = string2;

    swap(string1, string2);
    SEQAN_ASSERT_EQ(string1, string4);
    SEQAN_ASSERT_EQ(string2, string3);

    string1 = "ACAC";
    string2 = "GTGT";
    string3 = string1;
    string4 = string2;

    swap(string1, string2);
    SEQAN_ASSERT_EQ(string1, string4);
    SEQAN_ASSERT_EQ(string2, string3);
}

// Test of value().
template <typename TString>
void testSequenceValue(TString & /*Tag*/)
{
    using namespace seqan;

    typedef typename Value<TString>::Type TValue;

    // In contrast to getValue(), value() does not return a copy but a reference.
    // We test this using the variable value_.
    TString string = "ACAC";
    TValue & value_ = value(string, 0);
    SEQAN_ASSERT_EQ(value_, 'A');

    value_ = 'G';
    SEQAN_ASSERT_EQ(value_, 'G');
    SEQAN_ASSERT_EQ(string, "GCAC");
}

// --------------------------------------------------------------------------
// Testing Alloc Strings With Simple Types
// --------------------------------------------------------------------------

// Test whether sequences are copy constructible.
SEQAN_DEFINE_TEST(test_sequence_alloc_string_dna_copy_constructible)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag;
    testSequenceCopyConstructible(tag);

    String<Dna, Alloc<> > const constTag;
    testSequenceCopyConstructible(constTag);
}

// Test whether sequences are default constructible.
SEQAN_DEFINE_TEST(test_sequence_alloc_string_dna_default_constructible)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag;
    testSequenceDefaultConstructible(tag);

    String<Dna, Alloc<> > const constTag;
    testSequenceDefaultConstructible(constTag);
}

// Test operator<()
SEQAN_DEFINE_TEST(test_sequence_alloc_string_dna_less)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag;
    testSequenceLess(tag);

    String<Dna, Alloc<> > const constTag;
    testSequenceLess(constTag);
}

// Test operator<=()
SEQAN_DEFINE_TEST(test_sequence_alloc_string_dna_less_equal)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag;
    testSequenceLessEqual(tag);

    String<Dna, Alloc<> > const constTag;
    testSequenceLessEqual(constTag);
}

// Test operator>()
SEQAN_DEFINE_TEST(test_sequence_alloc_string_dna_greater)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag;
    testSequenceGreater(tag);

    String<Dna, Alloc<> > const constTag;
    testSequenceGreater(constTag);
}

// Test operator>=()
SEQAN_DEFINE_TEST(test_sequence_alloc_string_dna_greater_equal)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag;
    testSequenceGreaterEqual(tag);

    String<Dna, Alloc<> > const constTag;
    testSequenceGreaterEqual(constTag);
}

// Test operator==()
SEQAN_DEFINE_TEST(test_sequence_alloc_string_dna_equal)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag;
    testSequenceEqual(tag);

    String<Dna, Alloc<> > const constTag;
    testSequenceEqual(constTag);
}

// Test operator!=()
SEQAN_DEFINE_TEST(test_sequence_alloc_string_dna_unequal)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag;
    testSequenceUnequal(tag);

    String<Dna, Alloc<> > const constTag;
    testSequenceUnequal(constTag);
}

// Test of append().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_dna_append)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag;
    testSequenceAppend(tag);
}

// Test of appendValue().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_dna_append_value)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag;
    testSequenceAppendValue(tag);
}

// Test whether sequences are assignable.
SEQAN_DEFINE_TEST(test_sequence_alloc_string_dna_assignable)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag;
    testSequenceAssignable(tag);
}

// Test of assignValue().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_dna_assign_value)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag;
    testSequenceAssignValue(tag);
}

// Test of back().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_dna_back)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag;
    testSequenceBack(tag);

    String<Dna, Alloc<> > const constTag;
    testSequenceBack(constTag);
}

// Test of begin().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_dna_begin)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag;
    testSequenceBegin(tag);

    String<Dna, Alloc<> > const constTag;
    testSequenceBegin(constTag);
}

// Test of beginPosition().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_dna_begin_position)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag;
    testSequenceBeginPosition(tag);

    String<Dna, Alloc<> > const constTag;
    testSequenceBeginPosition(constTag);
}

// Test of capacity().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_dna_capacity)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag;
    testSequenceCapacity(tag);

    String<Dna, Alloc<> > const constTag;
    testSequenceCapacity(constTag);
}

// Test of clear().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_dna_clear)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag;
    testSequenceClear(tag);
}

// Test of end().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_dna_end)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag;
    testSequenceEnd(tag);

    String<Dna, Alloc<> > const constTag;
    testSequenceEnd(constTag);
}

// Test of endPosition().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_dna_end_position)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag;
    testSequenceEndPosition(tag);

    String<Dna, Alloc<> > const constTag;
    testSequenceEndPosition(constTag);
}

// Test of erase().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_dna_erase)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag;
    testSequenceErase(tag);
}

// Test of eraseBack().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_dna_erase_back)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag;
    testSequenceEraseBack(tag);
}

// Test of front().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_dna_front)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag;
    testSequenceFront(tag);
}

// Test of getValue().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_dna_get_value)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag;
    testSequenceGetValue(tag);

    String<Dna, Alloc<> > const constTag;
    testSequenceGetValue(constTag);
}

// Test of insert().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_dna_insert)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag;
    testSequenceInsert(tag);
}

// Test of insertValue().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_dna_insert_value)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag;
    testSequenceInsertValue(tag);
}

// Test of iter().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_dna_iter)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag;
    testSequenceIter(tag);

    String<Dna, Alloc<> > const constTag;
    testSequenceIter(constTag);
}

// Test of length().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_dna_length)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag;
    testSequenceLength(tag);

    String<Dna, Alloc<> > const constTag;
    testSequenceLength(constTag);
}

// Test of moveValue().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_dna_move_value)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag;
    testSequenceAssignValue(tag);
}

// Test of replace().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_dna_replace)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag;
    testSequenceReplace(tag);
}

// Test of reserve().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_dna_reserve)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag;
    testSequenceReserve(tag);
}

// Test of resize().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_dna_resize)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag;
    testSequenceResize(tag);
}

// Test of swap().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_dna_swap)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag;
    testSequenceSwap(tag);
}

// Test of value().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_dna_value)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag;
    testSequenceValue(tag);
}

// --------------------------------------------------------------------------
// Testing Alloc Strings With Non Simple Types
// --------------------------------------------------------------------------

// Test whether sequences are copy constructible.
SEQAN_DEFINE_TEST(test_sequence_alloc_string_counting_char_copy_constructible)
{
    using namespace seqan;

    CountingChar::clear();

    String<CountingChar, Alloc<> > tag;
    testSequenceCopyConstructible(tag);

    String<CountingChar, Alloc<> > const constTag;
    testSequenceCopyConstructible(constTag);

    SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
    SEQAN_ASSERT_GT(CountingChar::numConstruct, 0u);
}

// Test whether sequences are default constructible.
SEQAN_DEFINE_TEST(test_sequence_alloc_string_counting_char_default_constructible)
{
    using namespace seqan;

    CountingChar::clear();

    String<Dna, Alloc<> > tag;
    testSequenceDefaultConstructible(tag);

    String<CountingChar, Alloc<> > const constTag;
    testSequenceDefaultConstructible(constTag);

    SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
}

// Test operator<()
SEQAN_DEFINE_TEST(test_sequence_alloc_string_counting_char_less)
{
    using namespace seqan;

    CountingChar::clear();

    String<Dna, Alloc<> > tag;
    testSequenceLess(tag);

    String<CountingChar, Alloc<> > const constTag;
    testSequenceLess(constTag);

    SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
}

// Test operator<=()
SEQAN_DEFINE_TEST(test_sequence_alloc_string_counting_char_less_equal)
{
    using namespace seqan;

    CountingChar::clear();

    String<Dna, Alloc<> > tag;
    testSequenceLessEqual(tag);

    String<CountingChar, Alloc<> > const constTag;
    testSequenceLessEqual(constTag);

    SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
}

// Test operator>()
SEQAN_DEFINE_TEST(test_sequence_alloc_string_counting_char_greater)
{
    using namespace seqan;

    CountingChar::clear();

    String<Dna, Alloc<> > tag;
    testSequenceGreater(tag);

    String<CountingChar, Alloc<> > const constTag;
    testSequenceGreater(constTag);

    SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
}

// Test operator>=()
SEQAN_DEFINE_TEST(test_sequence_alloc_string_counting_char_greater_equal)
{
    using namespace seqan;

    CountingChar::clear();

    String<Dna, Alloc<> > tag;
    testSequenceGreaterEqual(tag);

    String<CountingChar, Alloc<> > const constTag;
    testSequenceGreaterEqual(constTag);

    SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
}

// Test operator==()
SEQAN_DEFINE_TEST(test_sequence_alloc_string_counting_char_equal)
{
    using namespace seqan;

    CountingChar::clear();

    String<Dna, Alloc<> > tag;
    testSequenceEqual(tag);

    String<CountingChar, Alloc<> > const constTag;
    testSequenceEqual(constTag);

    SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
}

// Test operator!=()
SEQAN_DEFINE_TEST(test_sequence_alloc_string_counting_char_unequal)
{
    using namespace seqan;

    CountingChar::clear();

    String<Dna, Alloc<> > tag;
    testSequenceUnequal(tag);

    String<CountingChar, Alloc<> > const constTag;
    testSequenceUnequal(constTag);

    SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
}


// Test of append().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_counting_char_append)
{
    using namespace seqan;

    CountingChar::clear();

    String<CountingChar, Alloc<> > tag;
    testSequenceAppend(tag);

    SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
    SEQAN_ASSERT_GT(CountingChar::numConstruct, 0u);
}

// Test of appendValue().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_counting_char_append_value)
{
    using namespace seqan;

    CountingChar::clear();

    String<CountingChar, Alloc<> > tag;
    testSequenceAppendValue(tag);

    SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
    SEQAN_ASSERT_GT(CountingChar::numConstruct, 0u);
}

// Test whether sequences are assignable.
SEQAN_DEFINE_TEST(test_sequence_alloc_string_counting_char_assignable)
{
    using namespace seqan;

    CountingChar::clear();

    String<CountingChar, Alloc<> > tag;
    testSequenceAssignable(tag);

    SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
    SEQAN_ASSERT_GT(CountingChar::numConstruct, 0u);
}

// Test of assignValue().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_counting_char_assign_value)
{
    using namespace seqan;

    CountingChar::clear();

    String<CountingChar, Alloc<> > tag;
    testSequenceAssignValue(tag);

    SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
    SEQAN_ASSERT_GT(CountingChar::numConstruct, 0u);
}

// Test of back().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_counting_char_back)
{
    using namespace seqan;

    CountingChar::clear();

    String<CountingChar, Alloc<> > tag;
    testSequenceBack(tag);

    String<CountingChar, Alloc<> > const constTag;
    testSequenceBack(constTag);

    SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
    SEQAN_ASSERT_GT(CountingChar::numConstruct, 0u);
}

// Test of begin().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_counting_char_begin)
{
    using namespace seqan;

    CountingChar::clear();

    String<CountingChar, Alloc<> > tag;
    testSequenceBegin(tag);

    String<CountingChar, Alloc<> > const constTag;
    testSequenceBegin(constTag);

    SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
    SEQAN_ASSERT_GT(CountingChar::numConstruct, 0u);
}

// Test of beginPosition().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_counting_char_begin_position)
{
    using namespace seqan;

    CountingChar::clear();

    String<CountingChar, Alloc<> > tag;
    testSequenceBeginPosition(tag);

    String<CountingChar, Alloc<> > const constTag;
    testSequenceBeginPosition(constTag);

    SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
    SEQAN_ASSERT_GT(CountingChar::numConstruct, 0u);
}

// Test of capacity().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_counting_char_capacity)
{
    using namespace seqan;

    CountingChar::clear();

    String<CountingChar, Alloc<> > tag;
    testSequenceCapacity(tag);

    String<CountingChar, Alloc<> > const constTag;
    testSequenceCapacity(constTag);

    SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
    SEQAN_ASSERT_GT(CountingChar::numConstruct, 0u);
}

// Test of clear().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_counting_char_clear)
{
    using namespace seqan;

    CountingChar::clear();

    String<CountingChar, Alloc<> > tag;
    testSequenceClear(tag);

    SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
    SEQAN_ASSERT_GT(CountingChar::numConstruct, 0u);
}

// Test of end().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_counting_char_end)
{
    using namespace seqan;

    CountingChar::clear();

    String<CountingChar, Alloc<> > tag;
    testSequenceEnd(tag);

    String<CountingChar, Alloc<> > const constTag;
    testSequenceEnd(constTag);

    SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
    SEQAN_ASSERT_GT(CountingChar::numConstruct, 0u);
}

// Test of front().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_counting_char_front)
{
    using namespace seqan;

    CountingChar::clear();

    String<CountingChar, Alloc<> > tag;
    testSequenceFront(tag);

    String<CountingChar, Alloc<> > const constTag;
    testSequenceFront(constTag);

    SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
    SEQAN_ASSERT_GT(CountingChar::numConstruct, 0u);
}

// Test of endPosition().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_counting_char_end_position)
{
    using namespace seqan;

    CountingChar::clear();

    String<CountingChar, Alloc<> > tag;
    testSequenceEndPosition(tag);

    String<CountingChar, Alloc<> > const constTag;
    testSequenceEndPosition(constTag);

    SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
    SEQAN_ASSERT_GT(CountingChar::numConstruct, 0u);
}

// Test of eraseBack().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_counting_char_erase_back)
{
    using namespace seqan;

    CountingChar::clear();

    String<CountingChar, Alloc<> > tag;
    testSequenceEraseBack(tag);

    SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
    SEQAN_ASSERT_GT(CountingChar::numConstruct, 0u);
}

// Test of erase().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_counting_char_erase)
{
    using namespace seqan;

    CountingChar::clear();

    String<CountingChar, Alloc<> > tag;
    testSequenceErase(tag);

    SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
    SEQAN_ASSERT_GT(CountingChar::numConstruct, 0u);
}

// Test of getValue().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_counting_char_get_value)
{
    using namespace seqan;

    CountingChar::clear();

    String<CountingChar, Alloc<> > tag;
    testSequenceGetValue(tag);

    String<CountingChar, Alloc<> > const constTag;
    testSequenceGetValue(constTag);

    SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
    SEQAN_ASSERT_GT(CountingChar::numConstruct, 0u);
}

// Test of insert().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_counting_char_insert)
{
    using namespace seqan;

    CountingChar::clear();

    String<CountingChar, Alloc<> > tag;
    testSequenceInsert(tag);

    SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
    SEQAN_ASSERT_GT(CountingChar::numConstruct, 0u);
}

// Test of insertValue().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_counting_char_insert_value)
{
    using namespace seqan;

    CountingChar::clear();

    String<CountingChar, Alloc<> > tag;
    testSequenceInsertValue(tag);

    SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
    SEQAN_ASSERT_GT(CountingChar::numConstruct, 0u);
}

// Test of iter().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_counting_char_iter)
{
    using namespace seqan;

    CountingChar::clear();

    String<CountingChar, Alloc<> > tag;
    testSequenceIter(tag);

    SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
    SEQAN_ASSERT_GT(CountingChar::numConstruct, 0u);
}

// Test length().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_counting_char_length)
{
    using namespace seqan;

    CountingChar::clear();

    String<CountingChar, Alloc<> > tag;
    testSequenceLength(tag);

    String<CountingChar, Alloc<> > const constTag;
    testSequenceLength(constTag);

    SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
    SEQAN_ASSERT_GT(CountingChar::numConstruct, 0u);
}

// Test of moveValue().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_counting_char_move_value)
{
    using namespace seqan;

    CountingChar::clear();

    String<CountingChar, Alloc<> > tag;
    testSequenceMoveValue(tag);

    SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
    SEQAN_ASSERT_GT(CountingChar::numConstruct, 0u);
}

// Test of replace().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_counting_char_replace)
{
    using namespace seqan;

    CountingChar::clear();

    String<CountingChar, Alloc<> > tag;
    testSequenceReplace(tag);

    SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
    SEQAN_ASSERT_GT(CountingChar::numConstruct, 0u);
}

// Test of reserve().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_counting_char_reserve)
{
    using namespace seqan;

    CountingChar::clear();

    String<CountingChar, Alloc<> > tag;
    testSequenceReserve(tag);

    SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
}

// Test of resize().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_counting_char_resize)
{
    using namespace seqan;

    CountingChar::clear();

    String<CountingChar, Alloc<> > tag;
    testSequenceResize(tag);

    SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
    SEQAN_ASSERT_GT(CountingChar::numConstruct, 0u);
}

// Test of swap().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_counting_char_swap)
{
    using namespace seqan;

    CountingChar::clear();

    String<CountingChar, Alloc<> > tag;
    testSequenceSwap(tag);

    SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
    SEQAN_ASSERT_GT(CountingChar::numConstruct, 0u);
}

// Test of value().
SEQAN_DEFINE_TEST(test_sequence_alloc_string_counting_char_value)
{
    using namespace seqan;

    CountingChar::clear();

    String<CountingChar, Alloc<> > tag;
    testSequenceValue(tag);

    SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
    SEQAN_ASSERT_GT(CountingChar::numConstruct, 0u);
}
#endif  // CORE_TESTS_SEQUENCE_TEST_SEQUENCE_H_
