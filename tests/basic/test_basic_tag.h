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
// Tests for the SeqAn tag header.  Most of the tests here just check that
// the give tags, structs etc. are there by instantiating them.
// ==========================================================================

#ifndef TEST_BASIC_TEST_BASIC_TAG_H_
#define TEST_BASIC_TEST_BASIC_TAG_H_

#include <seqan/basic.h>

// Helper type for Tag<> definition below.
struct Foo_;

// Check the Tag<T> struct template is there.
SEQAN_DEFINE_TEST(test_basic_tag_tag_struct)
{
    typedef Tag<Foo_> Foo;
}

// Check that basic tags Default and Nothing are there.
SEQAN_DEFINE_TEST(test_basic_tag_tag_basic_tags)
{
    Nothing nothing;
    (void)nothing;
    Default default_;
    (void)default_;
}

// Check that the Move tag is there.
// TODO(holtgrew): Does it actually belong here?
SEQAN_DEFINE_TEST(test_basic_tag_move)
{
    Move move;
    (void)move;
}

// Misc tags, that probably do not belong here but closer to where they are used.
SEQAN_DEFINE_TEST(test_basic_tag_misc_tags1)
{
    MinimalCtor tag1;
    (void)tag1;
    NonMinimalCtor tag2;
    (void)tag2;
}

// Misc tags, that probably do not belong here but closer to where they are used.
SEQAN_DEFINE_TEST(test_basic_tag_misc_tags2)
{
    GoEnd tag3;
    (void)tag3;
}

// Test for the TagList and TagSelector constructs and LENGTH on such lists.
SEQAN_DEFINE_TEST(test_basic_tag_tag_list_selector)
{
    typedef TagList<Nothing, TagList<Default> > TList;
    unsigned len = LENGTH<TList>::VALUE;
    SEQAN_ASSERT_EQ(len, 2u);
    bool b = IsSameType<typename TagSelector<TList>::Base, TagSelector<TagList<Default> > >::Type::VALUE;
    SEQAN_ASSERT(b);
}

// The DotDrawing tag does probably not belong here.
SEQAN_DEFINE_TEST(test_basic_tag_misc_tags3)
{
    DotDrawing tag;
    (void)tag;
}

// Various alignment-related tags do probably not belong here.
SEQAN_DEFINE_TEST(test_basic_tag_misc_tags4)
{
    HammingDistance tag1;
    (void)tag1;
    LevenshteinDistance tag2;
    (void)tag2;
    EditDistance tag3;
    (void)tag3;
    NeedlemanWunsch tag4;
    (void)tag4;
    BandedNeedlemanWunsch tag5;
    (void)tag5;
    Gotoh tag6;
    (void)tag6;
    BandedGotoh tag7;
    (void)tag7;
    MyersBitVector tag8;
    (void)tag8;
    MyersHirschberg tag9;
    (void)tag9;
    Hirschberg tag10;
    (void)tag10;
    Lcs tag11;
    (void)tag11;
    SmithWaterman tag12;
    (void)tag12;
    BandedSmithWaterman tag13;
    (void)tag13;
    SmithWatermanClump tag14;
    (void)tag14;
    WatermanEggert tag15;
    (void)tag15;
    BandedSmithWatermanClump tag16;
    (void)tag16;
    BandedWatermanEggert tag17;
    (void)tag17;
    Blat tag18;
    (void)tag18;
}

#endif  // TEST_BASIC_TEST_BASIC_TAG_H_
