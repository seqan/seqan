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

#ifndef CORE_TESTS_SEQUENCE_TEST_INFIX_H_
#define CORE_TESTS_SEQUENCE_TEST_INFIX_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include "test_sequence.h"

// --------------------------------------------------------------------------
// Generic String Tests
// --------------------------------------------------------------------------

// Test whether PrefixSegments are copy constructible.
template <typename TString>
void testPrefixConstructible(TString & /*Tag*/)
{
    using namespace seqan;

    TString string = "ACGTACGTACGT";

    // Test prefix() without start and end
    {
        Segment<TString, PrefixSegment> prefixSeg(string);
        SEQAN_ASSERT_EQ(prefixSeg, string);
    }

    // Test prefix() with start and end
    {
        Segment<TString, PrefixSegment> prefixSeg(string, 6);
        SEQAN_ASSERT_EQ(prefixSeg, "ACGTAC");
    }

    // Test prefix() with an infix as host with start and end
    {
        Segment<TString, PrefixSegment> prefixSeg1(string, 6);
        Segment<Segment<TString, PrefixSegment>, PrefixSegment > prefixSeg2(prefixSeg1, 2);
        SEQAN_ASSERT_EQ(prefixSeg2, "AC");
    }
}

// Test whether InfixSegments are copy constructible.
template <typename TString>
void testInfixConstructible(TString & /*Tag*/)
{
    using namespace seqan;

    TString string = "ACGTACGTACGT";

    // Test infix() without start and end
    {
        Segment<TString, InfixSegment> infixSeg(string);
        SEQAN_ASSERT_EQ(infixSeg, string);
    }

    // Test infix() with start and end
    {
        Segment<TString, InfixSegment> infixSeg(string, 4, 10);
        SEQAN_ASSERT_EQ(infixSeg, "ACGTAC");
    }

    // Test infix() with an infix as host with start and end
    {
        Segment<TString, InfixSegment> infixSeg(string, 4, 10);
        Segment<Segment<TString, InfixSegment> > infixSeg2(infixSeg, 1, 3);
        SEQAN_ASSERT_EQ(infixSeg2, "CG");
    }
}

// Test whether SuffixSegments are copy constructible.
template <typename TString>
void testSuffixConstructible(TString & /*Tag*/)
{
    using namespace seqan;

    TString string = "ACGTACGTACGT";

    // Test suffix() without start and end
    {
        Segment<TString, SuffixSegment> suffixSeg(string);
        SEQAN_ASSERT_EQ(suffixSeg, string);
    }

    // Test suffix() with start and end
    {
        Segment<TString, SuffixSegment> suffixSeg(string, 6);
        SEQAN_ASSERT_EQ(suffixSeg, "GTACGT");
    }

    // Test suffix() with an infix as host with start and end
    {
        Segment<TString, SuffixSegment> suffixSeg1(string, 6);
        Segment<Segment<TString, SuffixSegment>, SuffixSegment > suffixSeg2(suffixSeg1, 4);
        SEQAN_ASSERT_EQ(suffixSeg2, "GT");
    }
}

// TODO (singer): the constructor of segments has problems if 0 is one of the
// positions because 0 could be an iterator or a position.
// I would propose that segments are positional only.
// Test whether PrefixSegments are copy constructible.
template <typename TString>
void testPrefixCopyConstructible(TString & /*Tag*/)
{
    using namespace seqan;

    TString string = "ACGTACGTACGT";

    // Test copy constructible on empty template
    {
        Segment<TString, PrefixSegment> prefixSeg(string);
        Segment<TString, PrefixSegment> prefixSegCopy(prefixSeg);
        SEQAN_ASSERT_EQ(prefixSeg, prefixSegCopy);
    }

    // Test copy constructible on non empty template
    {
        Segment<TString, PrefixSegment> prefixSeg(string, 6);
        Segment<TString, PrefixSegment> prefixSegCopy(prefixSeg);
        SEQAN_ASSERT_EQ(prefixSeg, "ACGTAC");
        SEQAN_ASSERT_EQ(prefixSeg, prefixSegCopy);
    }
}

// TODO (singer): the constructor of segments has problems if 0 is one of the
// positions because 0 could be an iterator or a position.
// I would propose that segments are positional only.
// Test whether InfixSegments are copy constructible.
template <typename TString>
void testInfixCopyConstructible(TString & /*Tag*/)
{
    using namespace seqan;

    TString string = "ACGTACGTACGT";

    // Test copy constructible on empty template
    {
        Segment<TString, InfixSegment> infixSeg(string);
        Segment<TString, InfixSegment> infixSegCopy(infixSeg);
        SEQAN_ASSERT_EQ(infixSeg, infixSegCopy);
    }

   // Test copy constructible on non empty template
    {
        Segment<TString, InfixSegment> infixSeg(string, 4, 10);
        Segment<TString, InfixSegment> infixSegCopy(infixSeg);
        SEQAN_ASSERT_EQ(infixSeg, "ACGTAC");
        SEQAN_ASSERT_EQ(infixSeg, infixSegCopy);
    }
}

// TODO (singer): the constructor of segments has problems if 0 is one of the
// positions because 0 could be an iterator or a position.
// I would propose that segments are positional only.
// Test whether SuffixSegments are copy constructible.
template <typename TString>
void testSuffixCopyConstructible(TString & /*Tag*/)
{
    using namespace seqan;

    TString string = "ACGTACGTACGT";

    // Test copy constructible on empty template
    {
        Segment<TString, SuffixSegment> suffixSeg(string);
        Segment<TString, SuffixSegment> suffixSegCopy(suffixSeg);
        SEQAN_ASSERT_EQ(suffixSeg, suffixSegCopy);
    }

    // Test copy constructible on non empty template
    {
        Segment<TString, SuffixSegment> suffixSeg(string, 6);
        Segment<TString, SuffixSegment> suffixSegCopy(suffixSeg);
        SEQAN_ASSERT_EQ(suffixSeg, "GTACGT");
        SEQAN_ASSERT_EQ(suffixSeg, suffixSegCopy);
    }
}

// Test whether PrefixSegments are default constructible.
template <typename TString>
void testPrefixDefaultConstructible(TString & /*Tag*/)
{
    // TODO (singer): There is no description what happens
    // if we apply setHost(). Through experimentation it became
    // obvious that the start and end of the segment are set to the
    // begin of the host. Why???
    // A workaround is to use setEnd().
    // setHost() should be symmetric to construction from string.
    using namespace seqan;

    Segment<TString, PrefixSegment> prefixSeg;
    TString string = "ACGTACGTACGT";

    setHost(prefixSeg, string);
    setEnd(prefixSeg, end(string));

    SEQAN_ASSERT_EQ(begin(prefixSeg), begin(string));
    SEQAN_ASSERT_EQ(end(prefixSeg), end(string));
    SEQAN_ASSERT_EQ(prefixSeg, string);
}

// Test whether InfixSegments are default constructible.
template <typename TString>
void testInfixDefaultConstructible(TString & /*Tag*/)
{
    // TODO (singer): There is no description what happens
    // if we apply setHost(). Through experimentation it became
    // obvious that the start and end of the segment are set to the
    // begin of the host. Why???
    // A workaround is to use setEnd().
    // setHost() should be symmetric to construction from string.
    using namespace seqan;

    Segment<TString, InfixSegment> infixSeg;
    TString string = "ACGTACGTACGT";

    setHost(infixSeg, string);
    setEnd(infixSeg, end(string));

    SEQAN_ASSERT_EQ(begin(infixSeg), begin(string));
    SEQAN_ASSERT_EQ(end(infixSeg), end(string));
    SEQAN_ASSERT_EQ(infixSeg, string);
}

// Test whether SuffixSegments are default constructible.
template <typename TString>
void testSuffixDefaultConstructible(TString & /*Tag*/)
{
    // TODO (singer): There is no description what happens
    // if we apply setHost(). Through experimentation it became
    // obvious that the start and end of the segment are set to the
    // begin of the host. Why???
    // A workaround is to use setEnd().
    // setHost() should be symmetric to construction from string.
    using namespace seqan;

    Segment<TString, SuffixSegment> suffixSeg;
    TString string = "ACGTACGTACGT";

    setHost(suffixSeg, string);
    setEnd(suffixSeg, end(string));

    SEQAN_ASSERT_EQ(begin(suffixSeg), begin(string));
    SEQAN_ASSERT_EQ(end(suffixSeg), end(string));
    SEQAN_ASSERT_EQ(suffixSeg, string);
}

// Test operator<().
template <typename TString>
void testSegmentLess(TString & /*Tag*/)
{
    using namespace seqan;

    typedef Segment<TString, PrefixSegment> TPrefix;
    typedef Segment<TString, InfixSegment> TInfix;
    typedef Segment<TString, SuffixSegment> TSuffix;

    TString string = "ACGTACGTAAAA";

    {
        TPrefix prefixSeg1(string, 4);
        TPrefix prefixSeg2(string, 8);
        TInfix infixSeg(string, 5, 9);
        TSuffix suffixSeg(string, 3);

        SEQAN_ASSERT_EQ(prefixSeg1 < prefixSeg2, true);
        SEQAN_ASSERT_EQ(prefixSeg1 < infixSeg, true);
        SEQAN_ASSERT_EQ(prefixSeg1 < suffixSeg, true);
    }

    {
        TInfix infixSeg1(string, 4, 8);
        TInfix infixSeg2(string, 5, 9);
        TPrefix prefixSeg(string, 5);
        TSuffix suffixSeg(string, 3);

        SEQAN_ASSERT_EQ(infixSeg1 < infixSeg2, true);
        SEQAN_ASSERT_EQ(infixSeg1 < prefixSeg, true);
        SEQAN_ASSERT_EQ(infixSeg1 < suffixSeg, true);
    }

    {
        TSuffix suffixSeg1(string, 8);
        TSuffix suffixSeg2(string, 5);
        TPrefix prefixSeg(string, 5);
        TInfix infixSeg(string, 4, 8);

        SEQAN_ASSERT_EQ(suffixSeg1 < suffixSeg2, true);
        SEQAN_ASSERT_EQ(suffixSeg1 < prefixSeg, true);
        SEQAN_ASSERT_EQ(suffixSeg1 < infixSeg, true);
    }
}

// Test operator<=().
template <typename TString>
void testSegmentLessEqual(TString & /*Tag*/)
{
    using namespace seqan;

    typedef Segment<TString, PrefixSegment> TPrefix;
    typedef Segment<TString, InfixSegment> TInfix;
    typedef Segment<TString, SuffixSegment> TSuffix;

    TString string = "ACGTACGTAAAA";

    {
        TPrefix prefixSeg1(string, 4);
        TPrefix prefixSeg2(string, 8);
        TInfix infixSeg(string, 4, 8);
        TSuffix suffixSeg(string, 4);

        SEQAN_ASSERT_EQ(prefixSeg1 <= prefixSeg2, true);
        SEQAN_ASSERT_EQ(prefixSeg1 <= infixSeg, true);
        SEQAN_ASSERT_EQ(prefixSeg1 <= suffixSeg, true);
    }

    {
        TInfix infixSeg1(string, 4, 8);
        TInfix infixSeg2(string, 0, 4);
        TPrefix prefixSeg(string, 5);
        TSuffix suffixSeg(string, 3);

        SEQAN_ASSERT_EQ(infixSeg1 <= infixSeg2, true);
        SEQAN_ASSERT_EQ(infixSeg1 <= prefixSeg, true);
        SEQAN_ASSERT_EQ(infixSeg1 <= suffixSeg, true);
    }

    {
        TSuffix suffixSeg1(string, 8);
        TSuffix suffixSeg2(string, 5);
        TPrefix prefixSeg(string, 4);
        TInfix infixSeg(string, 4, 8);

        SEQAN_ASSERT_EQ(suffixSeg1 <= suffixSeg2, true);
        SEQAN_ASSERT_EQ(suffixSeg1 <= prefixSeg, true);
        SEQAN_ASSERT_EQ(suffixSeg1 <= infixSeg, true);
    }
}

// Test operator>().
template <typename TString>
void testSegmentGreater(TString & /*Tag*/)
{
    using namespace seqan;

    typedef Segment<TString, PrefixSegment> TPrefix;
    typedef Segment<TString, InfixSegment> TInfix;
    typedef Segment<TString, SuffixSegment> TSuffix;

    TString string = "ACGTACGTAAAA";

    {
        TPrefix prefixSeg1(string, 5);
        TPrefix prefixSeg2(string, 4);
        TInfix infixSeg(string, 4, 8);
        TSuffix suffixSeg(string, 8);

        SEQAN_ASSERT_EQ(prefixSeg1 > prefixSeg2, true);
        SEQAN_ASSERT_EQ(prefixSeg1 > infixSeg, true);
        SEQAN_ASSERT_EQ(prefixSeg1 > suffixSeg, true);
    }

    {
        TInfix infixSeg1(string, 0, 12);
        TInfix infixSeg2(string, 0, 4);
        TPrefix prefixSeg(string, 5);
        TSuffix suffixSeg(string, 8);

        SEQAN_ASSERT_EQ(infixSeg1 > infixSeg2, true);
        SEQAN_ASSERT_EQ(infixSeg1 > prefixSeg, true);
        SEQAN_ASSERT_EQ(infixSeg1 > suffixSeg, true);
    }

    {
        TSuffix suffixSeg1(string, 4);
        TSuffix suffixSeg2(string, 9);
        TPrefix prefixSeg(string, 4);
        TInfix infixSeg(string, 4, 8);

        SEQAN_ASSERT_EQ(suffixSeg1 > suffixSeg2, true);
        SEQAN_ASSERT_EQ(suffixSeg1 > prefixSeg, true);
        SEQAN_ASSERT_EQ(suffixSeg1 > infixSeg, true);
    }
}

// Test operator>=().
template <typename TString>
void testSegmentGreaterEqual(TString & /*Tag*/)
{
    using namespace seqan;

    typedef Segment<TString, PrefixSegment> TPrefix;
    typedef Segment<TString, InfixSegment> TInfix;
    typedef Segment<TString, SuffixSegment> TSuffix;

    TString string = "ACGTACGTAAAA";

    {
        TPrefix prefixSeg1(string, 4);
        TPrefix prefixSeg2(string, 4);
        TInfix infixSeg(string, 4, 8);
        TSuffix suffixSeg(string, 8);

        SEQAN_ASSERT_EQ(prefixSeg1 >= prefixSeg2, true);
        SEQAN_ASSERT_EQ(prefixSeg1 >= infixSeg, true);
        SEQAN_ASSERT_EQ(prefixSeg1 >= suffixSeg, true);
    }

    {
        TInfix infixSeg1(string, 0, 12);
        TInfix infixSeg2(string, 0, 4);
        TPrefix prefixSeg(string, 5);
        TSuffix suffixSeg(string, 9);

        SEQAN_ASSERT_EQ(infixSeg1 >= infixSeg2, true);
        SEQAN_ASSERT_EQ(infixSeg1 >= prefixSeg, true);
        SEQAN_ASSERT_EQ(infixSeg1 >= suffixSeg, true);
    }

    {
        TSuffix suffixSeg1(string, 4);
        TSuffix suffixSeg2(string, 9);
        TPrefix prefixSeg(string, 4);
        TInfix infixSeg(string, 4, 8);

        SEQAN_ASSERT_EQ(suffixSeg1 >= suffixSeg2, true);
        SEQAN_ASSERT_EQ(suffixSeg1 >= prefixSeg, true);
        SEQAN_ASSERT_EQ(suffixSeg1 >= infixSeg, true);
    }
}

// Test operator==().
template <typename TString>
void testSegmentEqual(TString & /*Tag*/)
{
    using namespace seqan;

    typedef Segment<TString, PrefixSegment> TPrefix;
    typedef Segment<TString, InfixSegment> TInfix;
    typedef Segment<TString, SuffixSegment> TSuffix;

    TString string = "ACGTACGT";

    {
        TPrefix prefixSeg1(string, 4);
        TPrefix prefixSeg2(string, 4);
        TInfix infixSeg(string, 4, 8);
        TSuffix suffixSeg(string, 4);

        SEQAN_ASSERT_EQ(prefixSeg1 == prefixSeg2, true);
        SEQAN_ASSERT_EQ(prefixSeg1 == infixSeg, true);
        SEQAN_ASSERT_EQ(prefixSeg1 == suffixSeg, true);
    }

    {
        TInfix infixSeg1(string, 0, 4);
        TInfix infixSeg2(string, 4, 8);
        TPrefix prefixSeg(string, 4);
        TSuffix suffixSeg(string, 4);

        SEQAN_ASSERT_EQ(infixSeg1 == infixSeg2, true);
        SEQAN_ASSERT_EQ(infixSeg1 == prefixSeg, true);
        SEQAN_ASSERT_EQ(infixSeg1 == suffixSeg, true);
    }

    {
        TSuffix suffixSeg1(string, 4);
        TSuffix suffixSeg2(string, 4);
        TPrefix prefixSeg(string, 4);
        TInfix infixSeg(string, 4, 8);

        SEQAN_ASSERT_EQ(suffixSeg1 == suffixSeg2, true);
        SEQAN_ASSERT_EQ(suffixSeg1 == prefixSeg, true);
        SEQAN_ASSERT_EQ(suffixSeg1 == infixSeg, true);
    }
}

// Test operator!=().
template <typename TString>
void testSegmentUnequal(TString & /*Tag*/)
{
    using namespace seqan;

    typedef Segment<TString, PrefixSegment> TPrefix;
    typedef Segment<TString, InfixSegment> TInfix;
    typedef Segment<TString, SuffixSegment> TSuffix;

    TString string = "ACGTACGT";

    {
        TPrefix prefixSeg1(string, 4);
        TPrefix prefixSeg2(string, 5);
        TInfix infixSeg(string, 4, 7);
        TSuffix suffixSeg(string, 3);

        SEQAN_ASSERT_EQ(prefixSeg1 != prefixSeg2, true);
        SEQAN_ASSERT_EQ(prefixSeg1 != infixSeg, true);
        SEQAN_ASSERT_EQ(prefixSeg1 != suffixSeg, true);
    }

    {
        TInfix infixSeg1(string, 1, 4);
        TInfix infixSeg2(string, 4, 8);
        TPrefix prefixSeg(string, 4);
        TSuffix suffixSeg(string, 4);

        SEQAN_ASSERT_EQ(infixSeg1 != infixSeg2, true);
        SEQAN_ASSERT_EQ(infixSeg1 != prefixSeg, true);
        SEQAN_ASSERT_EQ(infixSeg1 != suffixSeg, true);
    }

    {
        TSuffix suffixSeg1(string, 3);
        TSuffix suffixSeg2(string, 4);
        TPrefix prefixSeg(string, 4);
        TInfix infixSeg(string, 4, 8);

        SEQAN_ASSERT_EQ(suffixSeg1 != suffixSeg2, true);
        SEQAN_ASSERT_EQ(suffixSeg1 != prefixSeg, true);
        SEQAN_ASSERT_EQ(suffixSeg1 != infixSeg, true);
    }
}

// Test whether sequences are assignable.
template <typename TString>
void testSegmentAssignable(TString & /*Tag*/)
{
    using namespace seqan;

    typedef Segment<TString, PrefixSegment> TPrefix;
    typedef Segment<TString, InfixSegment> TInfix;
    typedef Segment<TString, SuffixSegment> TSuffix;

    TString string = "ACGTACGT";

    // Test assign() function.
    {
        // Test with default segment (the segment has no host).
        {
            TPrefix prefixSeg1;
            TPrefix prefixSeg2(string, 5);
            assign(prefixSeg1, prefixSeg2);
            SEQAN_ASSERT_EQ(prefixSeg1, prefixSeg2);
        }
        {
            TInfix infixSeg1;
            TInfix infixSeg2(string, 0, 5);
            assign(infixSeg1, infixSeg2);
            SEQAN_ASSERT_EQ(infixSeg1, infixSeg2);

            TInfix infixSeg3;
            TPrefix prefixSeg(string, 5);
            assign(infixSeg3, prefixSeg);
            SEQAN_ASSERT_EQ(infixSeg3, prefixSeg);

            TInfix infixSeg4;
            TSuffix suffixSeg(string, 3);
            assign(infixSeg4, suffixSeg);
            SEQAN_ASSERT_EQ(infixSeg4, suffixSeg);
        }
        {
            TSuffix suffixSeg1;
            TSuffix suffixSeg2(string, 5);
            assign(suffixSeg1, suffixSeg2);
            SEQAN_ASSERT_EQ(suffixSeg1, suffixSeg2);
        }
    }
    {
        // Test with initialized segment (the segment has a host).
        // In contrast to above the content of the host is changed.
        // TODO (singer): the described behaviour is not desired.
        {
            TPrefix prefixSeg1(string, 4);
            TPrefix prefixSeg2(string, 4);
            assign(prefixSeg1, prefixSeg2);
            SEQAN_ASSERT_EQ(prefixSeg1, prefixSeg2);
        }
        {
            TInfix infixSeg1(string, 0, 4);
            TInfix infixSeg2(string, 4, 8);
            assign(infixSeg1, infixSeg2);
            SEQAN_ASSERT_EQ(infixSeg1, infixSeg2);
        }
        {
            TSuffix suffixSeg1(string, 4);
            TSuffix suffixSeg2(string, 4);
            assign(suffixSeg1, suffixSeg2);
            SEQAN_ASSERT_EQ(suffixSeg1, suffixSeg2);
        }
    }

    // Test operator=()
    {
        // Test with default segment (the segment has no host).
        {
            TPrefix prefixSeg1;
            TPrefix prefixSeg2(string, 5);
            prefixSeg1 = prefixSeg2;
            SEQAN_ASSERT_EQ(prefixSeg1, prefixSeg2);
        }
        {
            TInfix infixSeg1;
            TInfix infixSeg2(string, 0, 5);
            infixSeg1 = infixSeg2;
            SEQAN_ASSERT_EQ(infixSeg1, infixSeg2);

            TInfix infixSeg3;
            TPrefix prefixSeg(string, 5);
            infixSeg3 = prefixSeg;
            SEQAN_ASSERT_EQ(infixSeg3, prefixSeg);

            TInfix infixSeg4;
            TSuffix suffixSeg(string, 3);
            infixSeg4 = suffixSeg;
            SEQAN_ASSERT_EQ(infixSeg4, suffixSeg);
        }
        {
            TSuffix suffixSeg1;
            TSuffix suffixSeg2(string, 5);
            suffixSeg1 = suffixSeg2;
            SEQAN_ASSERT_EQ(suffixSeg1, suffixSeg2);
        }
    }
    {
        // Test with initialized segment (the segment has a host).
        // In contrast to above the content of the host is changed.
        // TODO (singer): the described behaviour is not desired.
        {
            TPrefix prefixSeg1(string, 4);
            TPrefix prefixSeg2(string, 4);
            prefixSeg1 = prefixSeg2;
            SEQAN_ASSERT_EQ(prefixSeg1, prefixSeg2);
        }
        {
            TInfix infixSeg1(string, 0, 4);
            TInfix infixSeg2(string, 4, 8);
            infixSeg1 = infixSeg2;
            SEQAN_ASSERT_EQ(infixSeg1, infixSeg2);
        }
        {
            TSuffix suffixSeg1(string, 4);
            TSuffix suffixSeg2(string, 4);
            suffixSeg1 = suffixSeg2;
            SEQAN_ASSERT_EQ(suffixSeg1, suffixSeg2);
        }
    }
}

// Test of assignValue().
template <typename TString>
void testSegmentAssignValue(TString & /*Tag*/)
{
    using namespace seqan;

    typedef typename Value<TString>::Type TValue;
    typedef Segment<TString, PrefixSegment> TPrefix;
    typedef Segment<TString, InfixSegment> TInfix;
    typedef Segment<TString, SuffixSegment> TSuffix;

    TString string = "ACGTACGTACGT";

    TPrefix prefixSeg(string, 4);
    TInfix infixSeg(string, 4, 8);
    TSuffix suffixSeg(string, 8);

    // TODO (singer): assignValue for segments not in docu.
    assignValue(prefixSeg, 0, TValue('G'));
    assignValue(infixSeg, 0, TValue('G'));
    assignValue(suffixSeg, 0, TValue('G'));
    SEQAN_ASSERT_EQ(string[0], TValue('G'));
    SEQAN_ASSERT_EQ(string[4], TValue('G'));
    SEQAN_ASSERT_EQ(string[8], TValue('G'));
}

// We need two back() tests, since back() returns a reference or a copy.
// We check whether we can modify the reference.
// Test of back() for non const strings.
template <typename TString>
void testSegmentBack(TString & /*Tag*/)
{
    using namespace seqan;

    typedef typename Value<TString>::Type TValue;
    typedef Segment<TString, PrefixSegment> TPrefix;
    typedef Segment<TString, InfixSegment> TInfix;
    typedef Segment<TString, SuffixSegment> TSuffix;

    TString string = "AAAAAAAAAA";

    TPrefix prefixSeg(string, 3);
    TInfix infixSeg(string, 3, 8);
    TSuffix suffixSeg(string, 8);

    // val is a reference in contrast to the const version of back().
    TValue & preVal = back(prefixSeg);
    preVal = 'C';
    TValue & infVal = back(infixSeg);
    infVal = 'C';
    TValue & suffVal = back(suffixSeg);
    suffVal = 'C';
    SEQAN_ASSERT_EQ(string, "AACAAAACAC");
}

// Test of back() for const strings.
template <typename TString>
void testSegmentBack(TString const & /*Tag*/)
{
    using namespace seqan;

    typedef typename Value<TString>::Type TValue;
    typedef Segment<TString const, PrefixSegment> TPrefix;
    typedef Segment<TString const, InfixSegment> TInfix;
    typedef Segment<TString const, SuffixSegment> TSuffix;

    TString string = "AACAAAAGAT";

    TPrefix prefixSeg(string, 3);
    TInfix infixSeg(string, 3, 8);
    TSuffix suffixSeg(string, 8);

    // val is a reference in contrast to the const version of back().
    TValue preVal = back(prefixSeg);
    SEQAN_ASSERT_EQ(preVal, 'C');
    TValue infVal = back(infixSeg);
    SEQAN_ASSERT_EQ(infVal, 'G');
    TValue suffVal = back(suffixSeg);
    SEQAN_ASSERT_EQ(suffVal, 'T');
}

// Test of begin().
template <typename TString>
void testSegmentBegin(TString & /*Tag*/)
{
    using namespace seqan;

    typedef typename Value<TString>::Type TValue;
    typedef Segment<TString, PrefixSegment> TPrefix;
    typedef Segment<TString, InfixSegment> TInfix;
    typedef Segment<TString, SuffixSegment> TSuffix;

    TString string = "CAAGAAAATA";
    TPrefix prefixSeg(string, 3);
    TInfix infixSeg(string, 3, 8);
    TSuffix suffixSeg(string, 8);

    SEQAN_ASSERT_EQ(*begin(prefixSeg), TValue('C'));
    SEQAN_ASSERT_EQ(*begin(prefixSeg, Standard()), TValue('C'));
    SEQAN_ASSERT_EQ(*begin(prefixSeg, Rooted()), TValue('C'));

    SEQAN_ASSERT_EQ(*begin(infixSeg), TValue('G'));
    SEQAN_ASSERT_EQ(*begin(infixSeg, Standard()), TValue('G'));
    SEQAN_ASSERT_EQ(*begin(infixSeg, Rooted()), TValue('G'));

    SEQAN_ASSERT_EQ(*begin(suffixSeg), TValue('T'));
    SEQAN_ASSERT_EQ(*begin(suffixSeg, Standard()), TValue('T'));
    SEQAN_ASSERT_EQ(*begin(suffixSeg, Rooted()), TValue('T'));
}

// Test of beginPosition().
template <typename TString>
void testSegmentBeginPosition(TString & /*Tag*/)
{
    using namespace seqan;
    typedef Segment<TString, PrefixSegment> TPrefix;
    typedef Segment<TString, InfixSegment> TInfix;
    typedef Segment<TString, SuffixSegment> TSuffix;

    TString string = "CAAGAAAATA";
    TPrefix prefixSeg(string, 3);
    TInfix infixSeg(string, 3, 8);
    TSuffix suffixSeg(string, 8);

    SEQAN_ASSERT_EQ(beginPosition(prefixSeg), 0u);
    SEQAN_ASSERT_EQ(beginPosition(infixSeg), 3u);
    SEQAN_ASSERT_EQ(beginPosition(suffixSeg), 8u);
}

// Test of clear().
// template <typename TString>
// void testSegmentClear(TString & /*Tag*/)
// {
//     using namespace seqan;
//     typedef Segment<TString, PrefixSegment> TPrefix;
//     typedef Segment<TString, InfixSegment> TInfix;
//     typedef Segment<TString, SuffixSegment> TSuffix;
// 
//     // Since clear does not have to free anything it is only
//     // necessary that the code compiles.
//     // TODO (singer): at the moment clear erases parts of the host.
//     TString string = "CAAGAAAATA";
//     TPrefix prefixSeg(string, 3);
//     clear(prefixSeg);
//     TInfix infixSeg(string, 3, 7);
//     clear(infixSeg);
//     TSuffix suffixSeg(string, 1);
//     clear(suffixSeg);
// }

// Test of end().
template <typename TString>
void testSegmentEnd(TString & /*Tag*/)
{
    using namespace seqan;

    typedef typename Value<TString>::Type TValue;
    typedef Segment<TString, PrefixSegment> TPrefix;
    typedef Segment<TString, InfixSegment> TInfix;
    typedef Segment<TString, SuffixSegment> TSuffix;

    TString string = "AACAAAAGAT";
    TPrefix prefixSeg(string, 3);
    TInfix infixSeg(string, 3, 8);
    TSuffix suffixSeg(string, 8);

    typename Iterator<TPrefix, Standard>::Type preIter = end(prefixSeg);
    typename Iterator<TPrefix, Standard>::Type preStandardIter = end(prefixSeg, Standard());
    typename Iterator<TPrefix, Rooted>::Type preRootedIter = end(prefixSeg, Rooted());
    --preIter;
    --preStandardIter;
    --preRootedIter;
    SEQAN_ASSERT_EQ(*preIter, TValue('C'));
    SEQAN_ASSERT_EQ(*preStandardIter, TValue('C'));
    SEQAN_ASSERT_EQ(*preRootedIter, TValue('C'));

    typename Iterator<TInfix, Standard>::Type infIter = end(infixSeg);
    typename Iterator<TInfix, Standard>::Type infStandardIter = end(infixSeg, Standard());
    typename Iterator<TInfix, Rooted>::Type infRootedIter = end(infixSeg, Rooted());
    --infIter;
    --infStandardIter;
    --infRootedIter;
    SEQAN_ASSERT_EQ(*infIter, TValue('G'));
    SEQAN_ASSERT_EQ(*infStandardIter, TValue('G'));
    SEQAN_ASSERT_EQ(*infRootedIter, TValue('G'));

    typename Iterator<TSuffix, Standard>::Type sufixIter = end(suffixSeg);
    typename Iterator<TSuffix, Standard>::Type suffStandardIter = end(suffixSeg, Standard());
    typename Iterator<TSuffix, Rooted>::Type suffRootedIter = end(suffixSeg, Rooted());
    --sufixIter;
    --suffStandardIter;
    --suffRootedIter;
    SEQAN_ASSERT_EQ(*sufixIter, TValue('T'));
    SEQAN_ASSERT_EQ(*suffStandardIter, TValue('T'));
    SEQAN_ASSERT_EQ(*suffRootedIter, TValue('T'));
}

// Test of endPosition().
template <typename TString>
void testSegmentEndPosition(TString & /*Tag*/)
{
    using namespace seqan;
    typedef Segment<TString, PrefixSegment> TPrefix;
    typedef Segment<TString, InfixSegment> TInfix;
    typedef Segment<TString, SuffixSegment> TSuffix;

    TString string = "CAAGAAAATA";
    TPrefix prefixSeg(string, 3);
    TInfix infixSeg(string, 3, 8);
    TSuffix suffixSeg(string, 8);

    SEQAN_ASSERT_EQ(endPosition(prefixSeg), 3u);
    SEQAN_ASSERT_EQ(endPosition(infixSeg), 8u);
    SEQAN_ASSERT_EQ(endPosition(suffixSeg), 10u);
}

// We need two front() tests, since back() returns a reference or a copy.
// We check whether we can modify the reference.
// Test of front() for non const strings.
template <typename TString>
void testSegmentFront(TString & /*Tag*/)
{
    using namespace seqan;

    typedef typename Value<TString>::Type TValue;
    typedef Segment<TString, PrefixSegment> TPrefix;
    typedef Segment<TString, InfixSegment> TInfix;
    typedef Segment<TString, SuffixSegment> TSuffix;

    TString string = "AAAAAAAAAA";

    TPrefix prefixSeg(string, 3);
    TInfix infixSeg(string, 3, 8);
    TSuffix suffixSeg(string, 8);

    // val is a reference in contrast to the const version of back().
    TValue & preVal = front(prefixSeg);
    preVal = 'C';
    TValue & infVal = front(infixSeg);
    infVal = 'C';
    TValue & suffVal = front(suffixSeg);
    suffVal = 'C';
    SEQAN_ASSERT_EQ(string, "CAACAAAACA");
}

// Test of back() for const strings.
template <typename TString>
void testSegmentFront(TString const & /*Tag*/)
{
    using namespace seqan;

    typedef typename Value<TString>::Type TValue;
    typedef Segment<TString const, PrefixSegment> TPrefix;
    typedef Segment<TString const, InfixSegment> TInfix;
    typedef Segment<TString const, SuffixSegment> TSuffix;

    TString string = "CAAGAAAATA";

    TPrefix prefixSeg(string, 3);
    TInfix infixSeg(string, 3, 8);
    TSuffix suffixSeg(string, 8);

    // val is a reference in contrast to the const version of back().
    TValue preVal = front(prefixSeg);
    SEQAN_ASSERT_EQ(preVal, 'C');
    TValue infVal = front(infixSeg);
    SEQAN_ASSERT_EQ(infVal, 'G');
    TValue suffVal = front(suffixSeg);
    SEQAN_ASSERT_EQ(suffVal, 'T');
}

// Test of getValue().
template <typename TString>
void testSegmentGetValue(TString & /*Tag*/)
{
    using namespace seqan;

    typedef typename Value<TString>::Type TValue;
    typedef Segment<TString const, PrefixSegment> TPrefix;
    typedef Segment<TString const, InfixSegment> TInfix;
    typedef Segment<TString const, SuffixSegment> TSuffix;

    TString string = "CAAGAAAATA";

    TPrefix prefixSeg(string, 3);
    TInfix infixSeg(string, 3, 8);
    TSuffix suffixSeg(string, 8);

     // In contrast to value(), getValue() does not return a reference but a copy.
    TValue preValue = getValue(prefixSeg, 0);
    TValue infValue = getValue(infixSeg, 0);
    TValue suffValue = getValue(suffixSeg, 0);

    SEQAN_ASSERT_EQ(preValue, TValue('C'));
    SEQAN_ASSERT_EQ(infValue, TValue('G'));
    SEQAN_ASSERT_EQ(suffValue, TValue('T'));
}

// Test of iter().
template <typename TString>
void testSegmentIter(TString & /*Tag*/)
{
    using namespace seqan;

    //typedef typename Value<TString>::Type TValue;

    TString string = "ACGTACGTACGT";

    {
        typedef Segment<TString, PrefixSegment> TPrefix;
        typedef typename Iterator<TPrefix>::Type TIterator;
        typedef typename Iterator<TPrefix, Standard>::Type TStandardIterator;
        typedef typename Iterator<TPrefix, Rooted>::Type TRootedIterator;

        TPrefix prefixSeg(string, 4);

        TIterator iterator = iter(prefixSeg, 2);
        TStandardIterator standardIterator = iter(prefixSeg, 2);
        TRootedIterator rootedIterator = iter(prefixSeg, 2);
        SEQAN_ASSERT_EQ(getValue(iterator), getValue(string, 2));
        SEQAN_ASSERT_EQ(getValue(standardIterator), getValue(string, 2));
        SEQAN_ASSERT_EQ(getValue(rootedIterator), getValue(string, 2));
    }
    {
        typedef Segment<TString, InfixSegment> TInfix;
        typedef typename Iterator<TInfix>::Type TIterator;
        typedef typename Iterator<TInfix, Standard>::Type TStandardIterator;
        typedef typename Iterator<TInfix, Rooted>::Type TRootedIterator;

        TInfix infixSeg(string, 4, 8);

        TIterator iterator = iter(infixSeg, 2);
        TStandardIterator standardIterator = iter(infixSeg, 2);
        TRootedIterator rootedIterator = iter(infixSeg, 2);
        SEQAN_ASSERT_EQ(getValue(iterator), getValue(string, 6));
        SEQAN_ASSERT_EQ(getValue(standardIterator), getValue(string, 6));
        SEQAN_ASSERT_EQ(getValue(rootedIterator), getValue(string, 6));
    }
    {
        typedef Segment<TString, SuffixSegment> TSuffix;
        typedef typename Iterator<TSuffix>::Type TIterator;
        typedef typename Iterator<TSuffix, Standard>::Type TStandardIterator;
        typedef typename Iterator<TSuffix, Rooted>::Type TRootedIterator;

        TSuffix suffixSeg(string, 4);

        TIterator iterator = iter(suffixSeg, 2);
        TStandardIterator standardIterator = iter(suffixSeg, 2);
        TRootedIterator rootedIterator = iter(suffixSeg, 2);
        SEQAN_ASSERT_EQ(getValue(iterator), getValue(string, 6));
        SEQAN_ASSERT_EQ(getValue(standardIterator), getValue(string, 6));
        SEQAN_ASSERT_EQ(getValue(rootedIterator), getValue(string, 6));
    }
}

// Test of length().
template <typename TString>
void testSegmentLength(TString & /*Tag*/)
{
    using namespace seqan;

    //typedef typename Value<TString>::Type TValue;
    typedef Segment<TString, PrefixSegment> TPrefix;
    typedef Segment<TString, InfixSegment> TInfix;
    typedef Segment<TString, SuffixSegment> TSuffix;

    TString string = "ACGTACGTACGT";
    TPrefix prefixSeg(string, 4);
    TInfix infixSeg(string, 4, 8);
    TSuffix suffixSeg(string, 8);
    SEQAN_ASSERT_EQ(length(infixSeg), 4u);
    SEQAN_ASSERT_EQ(length(infixSeg), 4u);
    SEQAN_ASSERT_EQ(length(infixSeg), 4u);
}

// Test of value().
template <typename TString>
void testSegmentMoveValue(TString & /*Tag*/)
{
    using namespace seqan;

    //typedef typename Value<TString>::Type TValue;
    typedef Segment<TString, PrefixSegment> TPrefix;
    typedef Segment<TString, InfixSegment> TInfix;
    typedef Segment<TString, SuffixSegment> TSuffix;

    TString string = "AAAAAAAAAAAA";
    TPrefix prefixSeg(string, 4);
    TInfix infixSeg(string, 4, 8);
    TSuffix suffixSeg(string, 8);

    moveValue(prefixSeg, 3, 'C');
    moveValue(infixSeg, 1, 'C');
    moveValue(suffixSeg, 3, 'C');
    SEQAN_ASSERT_EQ(string, "AAACACAAAAAC");
}

// Test of value().
template <typename TString>
void testSegmentValue(TString & /*Tag*/)
{
    using namespace seqan;

    typedef typename Value<TString>::Type TValue;
    typedef Segment<TString, PrefixSegment> TPrefix;
    typedef Segment<TString, InfixSegment> TInfix;
    typedef Segment<TString, SuffixSegment> TSuffix;

    TString string = "CAAGAAAATA";
    TPrefix prefixSeg(string, 3);
    TInfix infixSeg(string, 3, 8);
    TSuffix suffixSeg(string, 8);

    // In contrast to getValue(), value() does not return a copy but a reference.
    // We test this using the variable value_.
    TValue & preValue = value(prefixSeg, 0);
    SEQAN_ASSERT_EQ(preValue, 'C');
    TValue & infValue = value(infixSeg, 0);
    SEQAN_ASSERT_EQ(infValue, 'G');
    TValue & suffValue = value(suffixSeg, 0);
    SEQAN_ASSERT_EQ(suffValue, 'T');

    preValue = 'A';
    infValue = 'A';
    suffValue = 'A';
    SEQAN_ASSERT_EQ(string, "AAAAAAAAAA");
}

// --------------------------------------------------------------------------
// Testing Alloc Strings With Simple Types
// --------------------------------------------------------------------------

// Test whether sequences are copy constructible.
SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_constructible)
{
    using namespace seqan;

    {
        String<Dna, Alloc<> > tag;
        testInfixConstructible(tag);
        testPrefixConstructible(tag);
        testSuffixConstructible(tag);
    }
//     {
//         String<Dna, Block<1024> > tag;
//         testInfixConstructible(tag);
//         testPrefixConstructible(tag);
//         testSuffixConstructible(tag);
//     }
//     {
//         String<Dna, CStyle> tag;
//         testInfixConstructible(tag);
//         testPrefixConstructible(tag);
//         testSuffixConstructible(tag);
//     } 
//     {
//         String<Dna, External<> > tag;
//         testInfixConstructible(tag);
//         testPrefixConstructible(tag);
//         testSuffixConstructible(tag);
//     } 
//     {
// //         String<Dna, MMap<> > tag;
// //         testInfixConstructible(tag);
// //         testPrefixConstructible(tag);
// //         testSuffixConstructible(tag);
//     } 
//     {
// //        String<Dna, Packed<> > tag;
// //        testInfixConstructible(tag);
// //        testPrefixConstructible(tag);
// //        testSuffixConstructible(tag);
//     } 

    String<Dna, Alloc<> > const constTag;
    testInfixConstructible(constTag);
    testPrefixConstructible(constTag);
    testSuffixConstructible(constTag);
}

// Test whether sequences are copy constructible.
SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_copy_constructible)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag;
    testInfixCopyConstructible(tag);
    testPrefixCopyConstructible(tag);
    testSuffixCopyConstructible(tag);

    String<Dna, Alloc<> > const constTag;
    testInfixCopyConstructible(constTag);
    testPrefixCopyConstructible(constTag);
    testSuffixCopyConstructible(constTag);
}

// Test whether sequences are default constructible.
SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_default_constructible)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag;
    testPrefixDefaultConstructible(tag);
    testInfixDefaultConstructible(tag);
    testSuffixDefaultConstructible(tag);

    String<Dna, Alloc<> > const constTag;
    testPrefixDefaultConstructible(constTag);
    testInfixDefaultConstructible(constTag);
    testSuffixDefaultConstructible(constTag);
}

// Test operator<()
SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_less)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag;
    testSegmentLess(tag);

    String<Dna, Alloc<> > const constTag;
    testSegmentLess(constTag);
}

// Test operator<=()
SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_less_equal)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag;
    testSegmentLessEqual(tag);

    String<Dna, Alloc<> > const constTag;
    testSegmentLessEqual(constTag);
}

// Test operator>()
SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_greater)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag;
    testSegmentGreater(tag);

    String<Dna, Alloc<> > const constTag;
    testSegmentGreater(constTag);
}

// Test operator>=()
SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_greater_equal)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag;
    testSegmentGreaterEqual(tag);

    String<Dna, Alloc<> > const constTag;
    testSegmentGreaterEqual(constTag);
}

// Test operator==()
SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_equal)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag;
    testSegmentEqual(tag);

    String<Dna, Alloc<> > const constTag;
    testSegmentEqual(constTag);
}

// Test operator!=()
SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_unequal)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag;
    testSegmentUnequal(tag);

    String<Dna, Alloc<> > const constTag;
    testSegmentUnequal(constTag);
}

// Test whether sequences are assignable.
SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_assignable)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag;
    testSegmentAssignable(tag);
}

// Test of assignValue().
SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_assign_value)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag;
    testSegmentAssignValue(tag);
}

// Test of back().
SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_back)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag;
    testSegmentBack(tag);

    String<Dna, Alloc<> > const constTag;
    testSegmentBack(constTag);
}

// Test of begin().
SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_begin)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag;
    testSegmentBegin(tag);

    String<Dna, Alloc<> > const constTag;
    testSegmentBegin(constTag);
}

// Test of beginPosition().
SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_begin_position)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag;
    testSegmentBeginPosition(tag);

    String<Dna, Alloc<> > const constTag;
    //testSegmentBeginPosition(constTag);
}

// // Test of clear().
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_clear)
// {
//     using namespace seqan;
// 
//     String<Dna, Alloc<> > tag;
//     testSegmentClear(tag);
// }

// Test of end().
SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_end)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag;
    testSegmentEnd(tag);

    String<Dna, Alloc<> > const constTag;
    testSegmentEnd(constTag);
}

// Test of endPosition().
SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_end_position)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag;
    testSegmentEndPosition(tag);

    String<Dna, Alloc<> > const constTag;
    testSegmentEndPosition(constTag);
}

// Test of front().
SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_front)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag;
    testSegmentFront(tag);

    String<Dna, Alloc<> > const consTag;
    testSegmentFront(consTag);
}

// Test of getValue().
SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_get_value)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag;
    testSegmentGetValue(tag);

    String<Dna, Alloc<> > const constTag;
    testSegmentGetValue(constTag);
}

// Test of iter().
SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_iter)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag;
    testSegmentIter(tag);

    String<Dna, Alloc<> > const constTag;
    testSegmentIter(constTag);
}

// Test of length().
SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_length)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag;
    testSegmentLength(tag);

    String<Dna, Alloc<> > const constTag;
    testSegmentLength(constTag);
}

// Test of moveValue().
SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_move_value)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag;
    testSegmentMoveValue(tag);
}

// Test of value().
SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_value)
{
    using namespace seqan;

    String<Dna, Alloc<> > tag;
    testSegmentValue(tag);
}

// TODO (singer): add tests to ensure that the underlying string of a segment is not modified. 

#endif // CORE_TESTS_SEQUENCE_TEST_INFIX_H_

