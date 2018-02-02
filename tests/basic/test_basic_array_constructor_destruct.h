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
// Tests for fundamental tags.
// ==========================================================================

#ifndef SEQAN_TESTS_BASIC_TEST_BASIC_FUNDAMENTAL_TAGS_H_
#define SEQAN_TESTS_BASIC_TEST_BASIC_FUNDAMENTAL_TAGS_H_

// Test for the Tag<T> template.
struct SomeTag_;
typedef seqan::Tag<SomeTag_> SomeTag;

SEQAN_DEFINE_TEST(test_basic_fundamental_tags_tag)
{
    SomeTag instance;
    (void)instance;
}

// Test that predefined tags are there.
SEQAN_DEFINE_TEST(test_basic_fundamental_tags_tags)
{
    using namespace seqan;

#define SEQAN_TAG_INSTANCE_TEST(x) x instance ## x; (void) instance ## x

    SEQAN_TAG_INSTANCE_TEST(Default);
    SEQAN_TAG_INSTANCE_TEST(Nothing);
    SEQAN_TAG_INSTANCE_TEST(Move);
    SEQAN_TAG_INSTANCE_TEST(MinimalCtor);
    // The following should not go into basic at all.
    SEQAN_TAG_INSTANCE_TEST(DotDrawing);
    SEQAN_TAG_INSTANCE_TEST(HammingDistance);
    SEQAN_TAG_INSTANCE_TEST(LevenshteinDistance);
    SEQAN_TAG_INSTANCE_TEST(EditDistance);
    SEQAN_TAG_INSTANCE_TEST(NeedlemanWunsch);
    SEQAN_TAG_INSTANCE_TEST(BandedNeedlemanWunsch);
    SEQAN_TAG_INSTANCE_TEST(Gotoh);
    SEQAN_TAG_INSTANCE_TEST(BandedGotoh);
    SEQAN_TAG_INSTANCE_TEST(MyersBitVector);
    SEQAN_TAG_INSTANCE_TEST(MyersHirschberg);
    SEQAN_TAG_INSTANCE_TEST(Hirschberg);
    SEQAN_TAG_INSTANCE_TEST(Lcs);
    SEQAN_TAG_INSTANCE_TEST(SmithWaterman);
    SEQAN_TAG_INSTANCE_TEST(BandedSmithWaterman);
    SEQAN_TAG_INSTANCE_TEST(SmithWatermanClump);
    SEQAN_TAG_INSTANCE_TEST(WatermanEggert);
    SEQAN_TAG_INSTANCE_TEST(BandedSmithWatermanClump);
    SEQAN_TAG_INSTANCE_TEST(BandedWatermanEggert);
    SEQAN_TAG_INSTANCE_TEST(Blat);

#undef SEQAN_TAG_INSTANCE_TEST
}

SEQAN_DEFINE_TEST(test_basic_fundamental_tags_tag_list)
{
    using namespace seqan;

    typedef TagList<void, TagList<void, TagList<void> > > TTagList;

    SEQAN_ASSERT(+(SameType_<typename TTagList::Type, void>::VALUE));
}

SEQAN_DEFINE_TEST(test_basic_fundamental_tags_tag_selector)
{
    using namespace seqan;

}

SEQAN_DEFINE_TEST(test_basic_fundamental_tags_length_tag_list)
{
    using namespace seqan;

    typedef TagList<>                                           TList0;
    typedef TagList<void, void>                                 TList1;
    typedef TagList<void, TagList<void, void> >                 TList2;
    typedef TagList<void, TagList<void, TagList<void, void> > > TList3;

    // Using unary-plus trick to create rvalue.
    SEQAN_ASSERT_EQ(+LENGTH<void>::VALUE, 0);
    SEQAN_ASSERT_EQ(+LENGTH<TList0>::VALUE, 1);
    SEQAN_ASSERT_EQ(+LENGTH<TList1>::VALUE, 1);
    SEQAN_ASSERT_EQ(+LENGTH<TList2>::VALUE, 2);
    SEQAN_ASSERT_EQ(+LENGTH<TList3>::VALUE, 3);
}

#endif  // #ifndef SEQAN_TESTS_BASIC_TEST_BASIC_FUNDAMENTAL_TAGS_H_
