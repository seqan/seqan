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
// Unit tests for the journaled string set.
// ==========================================================================

#ifndef TESTS_JOURNALED_SET_TEST_JOURNALED_SET_H_
#define TESTS_JOURNALED_SET_TEST_JOURNALED_SET_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/journaled_set.h>

using namespace seqan;

SEQAN_DEFINE_TEST(test_journaled_set_constructor)
{
    typedef String<char, Journaled<Alloc<>, SortedArray > > TJString;

// ----------------------------------------------------------------------------
// Default Constructor.
// ----------------------------------------------------------------------------

    {  // Test default constructor with some string.
        StringSet<CharString, Owner<JournaledSet> > jSet;
        SEQAN_ASSERT_EQ(empty(jSet.strings), true);
        SEQAN_ASSERT_EQ(length(jSet.limits), 1u);
        SEQAN_ASSERT_EQ(jSet.limitsValid, true);
        SEQAN_ASSERT_EQ(empty(jSet._globalRefHolder), true);
    }

    {  // Test default constructor with some string.
        StringSet<TJString, Owner<JournaledSet> > jSet;
        SEQAN_ASSERT_EQ(empty(jSet.strings), true);
        SEQAN_ASSERT_EQ(length(jSet.limits), 1u);
        SEQAN_ASSERT_EQ(jSet.limitsValid, true);
        SEQAN_ASSERT_EQ(empty(jSet._globalRefHolder), true);
    }

// ----------------------------------------------------------------------------
// Copy Constructor.
// ----------------------------------------------------------------------------

    CharString ref = "This is the reference sequence";
    {  // Test copy constructor with non const.
        StringSet<TJString, Owner<JournaledSet> > jSet;
        appendValue(jSet.strings, TJString());
        appendValue(jSet.strings, TJString());
        setValue(jSet._globalRefHolder, ref);
        StringSet<TJString, Owner<JournaledSet> > jSet2(jSet);

        SEQAN_ASSERT_EQ(length(jSet2.strings), 2u);
        SEQAN_ASSERT_EQ(length(jSet2.limits), 1u);
        SEQAN_ASSERT_EQ(jSet2.limitsValid, true);
        SEQAN_ASSERT_EQ(&value(jSet2._globalRefHolder), &ref);
        SEQAN_ASSERT_EQ(&value(jSet2._globalRefHolder), &value(jSet._globalRefHolder));
    }

    {  // Test copy constructor with non const.
        StringSet<TJString, Owner<JournaledSet> > jSet;
        appendValue(jSet.strings, TJString());
        appendValue(jSet.strings, TJString());
        setValue(jSet._globalRefHolder, ref);
        const StringSet<TJString, Owner<JournaledSet> > jSet2(jSet);

        SEQAN_ASSERT_EQ(length(jSet2.strings), 2u);
        SEQAN_ASSERT_EQ(length(jSet2.limits), 1u);
        SEQAN_ASSERT_EQ(jSet2.limitsValid, true);
        SEQAN_ASSERT_EQ(&value(jSet2._globalRefHolder), &ref);
        SEQAN_ASSERT_EQ(&value(jSet2._globalRefHolder), &value(jSet._globalRefHolder));
    }
}

SEQAN_DEFINE_TEST(test_journaled_set_assign)
{
    typedef String<char, Journaled<Alloc<>, SortedArray > > TJString;

    CharString ref = "This is the reference sequence";

    // Test assign with non const lvalue.
    StringSet<TJString, Owner<JournaledSet> > jSet;
    appendValue(jSet.strings, TJString());
    appendValue(jSet.strings, TJString());
    setValue(jSet._globalRefHolder, ref);
    StringSet<TJString, Owner<JournaledSet> > jSet2;

    jSet2 = jSet;

    SEQAN_ASSERT_EQ(length(jSet2.strings), 2u);
    SEQAN_ASSERT_EQ(length(jSet2.limits), 1u);
    SEQAN_ASSERT_EQ(jSet2.limitsValid, true);
    SEQAN_ASSERT_EQ(&value(jSet2._globalRefHolder), &ref);
    SEQAN_ASSERT_EQ(&value(jSet2._globalRefHolder), &value(jSet._globalRefHolder));

    // Test assign with const lvalue.
    const StringSet<TJString, Owner<JournaledSet> > jSetConst(jSet);
    StringSet<TJString, Owner<JournaledSet> > jSet3;
    jSet3 = jSetConst;

    SEQAN_ASSERT_EQ(length(jSet3.strings), 2u);
    SEQAN_ASSERT_EQ(length(jSet3.limits), 1u);
    SEQAN_ASSERT_EQ(jSet3.limitsValid, true);
    SEQAN_ASSERT_EQ(&value(jSet3._globalRefHolder), &ref);
    SEQAN_ASSERT_EQ(&value(jSet3._globalRefHolder), &value(jSet._globalRefHolder));
}

SEQAN_DEFINE_TEST(test_journaled_set_set_host)
{
    typedef String<char, Journaled<Alloc<>, SortedArray > > TJString;

    CharString ref = "This is the reference sequence";

    StringSet<TJString, Owner<JournaledSet> > jSet;
    setHost(jSet, ref);

    SEQAN_ASSERT_EQ(&value(jSet._globalRefHolder), &ref);
    SEQAN_ASSERT_EQ(&value(jSet._globalRefHolder), &value(jSet._globalRefHolder));
}

SEQAN_DEFINE_TEST(test_journaled_set_host)
{
    typedef String<char, Journaled<Alloc<>, SortedArray > > TJString;

    CharString ref = "This is the reference sequence";

    StringSet<TJString, Owner<JournaledSet> > jSet;
    setHost(jSet, ref);

    SEQAN_ASSERT_EQ(&host(jSet), &ref);
    SEQAN_ASSERT_EQ(host(jSet), "This is the reference sequence");

    const StringSet<TJString, Owner<JournaledSet> > jSetConst(jSet);
    SEQAN_ASSERT_EQ(&host(jSetConst), &ref);
    SEQAN_ASSERT_EQ(host(jSetConst), "This is the reference sequence");
}

SEQAN_DEFINE_TEST(test_journaled_set_resize)
{
    typedef String<char, Journaled<Alloc<>, SortedArray > > TJString;

    CharString ref = "This is the reference sequence";

    {  // Resize without fill value.
        StringSet<TJString, Owner<JournaledSet> > jSet;
        setHost(jSet, ref);

        resize(jSet, 10, Exact());

        SEQAN_ASSERT_EQ(length(jSet), 10u);

        resize(jSet, 128, Generous());

        SEQAN_ASSERT_EQ(length(jSet), 128u);
    }

    {  // Resize with fill value.
        TJString fillValue;
        StringSet<TJString, Owner<JournaledSet> > jSet;
        setHost(jSet, ref);
        setHost(fillValue, host(jSet));
        insert(fillValue, 7, "xxx");
        erase(fillValue, 15, 19);

        resize(jSet, 10, fillValue, Exact());
        insert(fillValue, 0, "xxx");

        SEQAN_ASSERT_EQ(length(jSet), 10u);
        for (unsigned i = 0; i < 10u; ++i)
        {
            SEQAN_ASSERT_EQ(&host(jSet[i]), &ref);
            SEQAN_ASSERT_EQ(jSet[i], "This isxxx the rence sequence");
        }

        resize(jSet, 128, fillValue, Generous());
        SEQAN_ASSERT_EQ(length(jSet), 128u);
        for (unsigned i = 0; i < 10u; ++i)
        {
            SEQAN_ASSERT_EQ(&host(jSet[i]), &ref);
            SEQAN_ASSERT_EQ(jSet[i], "This isxxx the rence sequence");
        }
        for (unsigned i = 10; i < 128u; ++i)
        {
            SEQAN_ASSERT_EQ(&host(jSet[i]), &ref);
            SEQAN_ASSERT_EQ(jSet[i], "xxxThis isxxx the rence sequence");
        }
    }
}

SEQAN_DEFINE_TEST(test_journaled_set_append_value)
{
    typedef String<char, Journaled<Alloc<>, SortedArray > > TJString;

    CharString ref = "This is the reference sequence";

    StringSet<TJString, Owner<JournaledSet> > jSet;
    setHost(jSet, ref);

    // Assign empty string.
    TJString tmpEmpty;
    appendValue(jSet, tmpEmpty);

    // Assign non-empty string.
    TJString tmp(ref);
    TJString tmp2(ref);
    insert(tmp2, 7, "xxx");
    erase(tmp2, 15, 19);

    appendValue(jSet, tmp);
    appendValue(jSet, tmp2);

    // Assign different string;
    CharString  tmpChar = "Different string type";
    appendValue(jSet, tmpChar);

    SEQAN_ASSERT_EQ(length(jSet.strings), 4u);
    SEQAN_ASSERT_EQ(&value(jSet._globalRefHolder), &ref);
    SEQAN_ASSERT_EQ(&value(jSet._globalRefHolder), &value(jSet._globalRefHolder));
    SEQAN_ASSERT_EQ(empty(jSet.strings[0]), true);
    SEQAN_ASSERT_EQ(jSet.strings[1], ref);
    SEQAN_ASSERT_EQ(&host(jSet.strings[1]), &ref);
    SEQAN_ASSERT_EQ(jSet.strings[2], "This isxxx the rence sequence");
    SEQAN_ASSERT_EQ(&host(jSet.strings[2]), &ref);
    SEQAN_ASSERT_EQ(jSet.strings[3], "Different string type");
    SEQAN_ASSERT_NEQ(&host(jSet.strings[3]), &ref);
}

SEQAN_DEFINE_TEST(test_journaled_set_assign_value)
{
    typedef String<char, Journaled<Alloc<>, SortedArray > > TJString;

    CharString ref = "This is the reference sequence";

    StringSet<TJString, Owner<JournaledSet> > jSet;
    setHost(jSet, ref);
    resize(jSet, 4, Exact());

    // Assign empty string.
    TJString tmpEmpty;
    assignValue(jSet, 0, tmpEmpty);

    // Assign non-empty string.
    TJString tmp(ref);
    TJString tmp2(ref);
    insert(tmp2, 7, "xxx");
    erase(tmp2, 15, 19);

    assignValue(jSet, 1, tmp);
    assignValue(jSet, 2, tmp2);

    // Assign different string;
    CharString  tmpChar = "Different string type";
    assignValue(jSet, 3, tmpChar);

    SEQAN_ASSERT_EQ(length(jSet.strings), 4u);
    SEQAN_ASSERT_EQ(&value(jSet._globalRefHolder), &ref);
    SEQAN_ASSERT_EQ(&value(jSet._globalRefHolder), &value(jSet._globalRefHolder));
    SEQAN_ASSERT_EQ(empty(jSet.strings[0]), true);
    SEQAN_ASSERT_EQ(jSet.strings[1], ref);
    SEQAN_ASSERT_EQ(&host(jSet.strings[1]), &ref);
    SEQAN_ASSERT_EQ(jSet.strings[2], "This isxxx the rence sequence");
    SEQAN_ASSERT_EQ(&host(jSet.strings[2]), &ref);
    SEQAN_ASSERT_EQ(jSet.strings[3], "Different string type");
    SEQAN_ASSERT_NEQ(&host(jSet.strings[3]), &ref);
}

SEQAN_DEFINE_TEST(test_journaled_set_value)
{
    typedef String<char, Journaled<Alloc<>, SortedArray > > TJString;

    CharString ref = "This is the reference sequence";

    StringSet<TJString, Owner<JournaledSet> > jSet;
    setHost(jSet, ref);
    resize(jSet, 4, Exact());

    // Assign empty string.
    TJString tmpEmpty;
    assignValue(jSet, 0, tmpEmpty);

    // Assign non-empty string.
    TJString tmp(ref);
    TJString tmp2(ref);
    insert(tmp2, 7, "xxx");
    erase(tmp2, 15, 19);

    assignValue(jSet, 1, tmp);
    assignValue(jSet, 2, tmp2);

    // Assign different string;
    CharString  tmpChar = "Different string type";
    assignValue(jSet, 3, tmpChar);

    SEQAN_ASSERT_EQ(length(jSet.strings), 4u);
    SEQAN_ASSERT_EQ(&value(jSet._globalRefHolder), &ref);
    SEQAN_ASSERT_EQ(&value(jSet._globalRefHolder), &value(jSet._globalRefHolder));
    SEQAN_ASSERT_EQ(empty(value(jSet, 0)), true);
    SEQAN_ASSERT_EQ(value(jSet, 1), ref);
    SEQAN_ASSERT_EQ(&host(jSet.strings[1]), &ref);
    SEQAN_ASSERT_EQ(value(jSet, 2), "This isxxx the rence sequence");
    SEQAN_ASSERT_EQ(&host(value(jSet, 2)), &ref);
    SEQAN_ASSERT_EQ(value(jSet, 3), "Different string type");
    SEQAN_ASSERT_NEQ(&host(value(jSet, 3)), &ref);
}

SEQAN_DEFINE_TEST(test_journaled_set_empty)
{
    typedef String<char, Journaled<Alloc<>, SortedArray > > TJString;

    CharString ref = "This is the reference sequence";

    StringSet<TJString, Owner<JournaledSet> > jSet;
    setHost(jSet, ref);

    SEQAN_ASSERT_EQ(empty(jSet), true);
    resize(jSet, 4, Exact());
    SEQAN_ASSERT_EQ(empty(jSet), false);
}

SEQAN_DEFINE_TEST(test_journaled_set_clear)
{
    typedef String<char, Journaled<Alloc<>, SortedArray > > TJString;

    CharString ref = "This is the reference sequence";

    StringSet<TJString, Owner<JournaledSet> > jSet;
    setHost(jSet, ref);

    SEQAN_ASSERT_EQ(empty(jSet), true);
    resize(jSet, 4, Exact());
    SEQAN_ASSERT_EQ(empty(jSet), false);
    clear(jSet);
    SEQAN_ASSERT_EQ(empty(jSet), true);
}

SEQAN_DEFINE_TEST(test_journaled_set_iterator)
{
    typedef String<char, Journaled<Alloc<>, SortedArray > > TJString;
    typedef Iterator<StringSet<TJString, Owner<JournaledSet> > >::Type TIter;

    CharString ref = "This is the reference sequence";

    StringSet<TJString, Owner<JournaledSet> > jSet;
    setHost(jSet, ref);
    resize(jSet, 4, Exact());

    // Assign empty string.
    TJString tmpEmpty;
    assignValue(jSet, 0, tmpEmpty);

    // Assign non-empty string.
    TJString tmp(ref);
    TJString tmp2(ref);
    insert(tmp2, 7, "xxx");
    erase(tmp2, 15, 19);

    assignValue(jSet, 1, tmp);
    assignValue(jSet, 2, tmp2);

    // Assign different string;
    CharString  tmpChar = "Different string type";
    assignValue(jSet, 3, tmpChar);

    TIter it = begin(jSet, Standard());
    TIter itEnd = end(jSet, Standard());

    SEQAN_ASSERT_EQ(empty(*it), true);
    ++it;
    SEQAN_ASSERT_EQ(*it, ref);
    SEQAN_ASSERT_EQ(&host(*it), &ref);
    it += 2;
    SEQAN_ASSERT_EQ(*it, "Different string type");
    SEQAN_ASSERT_NEQ(&host(*it), &ref);
    SEQAN_ASSERT((it + 1) == itEnd);

    it--;
    SEQAN_ASSERT_EQ(*it, "This isxxx the rence sequence");
    SEQAN_ASSERT_EQ(&host(*it), &ref);
}

#endif  // TESTS_JOURNALED_SET_TEST_JOURNALED_SET_H_
