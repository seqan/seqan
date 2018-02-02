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

#ifndef SEQAN_TESTS_ALIGN_TEST_ALIGN_FRAGMENT_H_
#define SEQAN_TESTS_ALIGN_TEST_ALIGN_FRAGMENT_H_

#include <seqan/align.h>

SEQAN_DEFINE_TEST(test_align_fragment)
{
    using namespace seqan;

    typedef StringSet<CharString, Dependent<> > CharStringSet;
    typedef    Id<CharStringSet>::Type TId;
    typedef    Size<CharStringSet>::Type TSize;

    CharStringSet str;
    CharString str0("annual");
    assignValueById(str, str0);
    CharString str1("anneal");
    assignValueById(str, str1);

    // Fragment: SeqId1, Begin1, SeqId2, Begin2, Length of Fragment
    Fragment<> f(0,4,1,4,2);
    SEQAN_ASSERT_EQ(f.seqId1, 0u);
    SEQAN_ASSERT_EQ(f.begin1, 4u);
    SEQAN_ASSERT_EQ(f.seqId2, 1u);
    SEQAN_ASSERT_EQ(f.begin2, 4u);
    SEQAN_ASSERT_EQ(f.len, 2u);
    SEQAN_ASSERT_EQ(fragmentBegin(f, 0), 4u);
    SEQAN_ASSERT_EQ(fragmentBegin(f, 1), 4u);
    SEQAN_ASSERT_EQ(fragmentLength(f, 0), 2u);
    SEQAN_ASSERT_EQ(fragmentLength(f, 1), 2u);
    SEQAN_ASSERT_EQ(sequenceId(f, 0), 0u);
    SEQAN_ASSERT_EQ(sequenceId(f, 1), 1u);
    SEQAN_ASSERT_EQ(label(f, str, 0), "al");
    SEQAN_ASSERT_EQ(label(f, str, 1), "al");

    TId id2;
    TSize pos2;
    getProjectedPosition(f, 0, 5, id2, pos2);
    SEQAN_ASSERT_EQ(pos2, 5u);
    SEQAN_ASSERT_EQ(id2, 1u);
    getProjectedPosition(f, 1, 5, id2, pos2);
    SEQAN_ASSERT_EQ(pos2, 5u);
    SEQAN_ASSERT_EQ(id2, 0u);

    // Fragment: SeqId1, Begin1, SeqId2, Begin2, Length of Fragment
    Fragment<> f2(0, 0, 1, 4, 1);
    SEQAN_ASSERT_EQ(f2.seqId1, 0u);
    SEQAN_ASSERT_EQ(f2.begin1, 0u);
    SEQAN_ASSERT_EQ(f2.seqId2, 1u);
    SEQAN_ASSERT_EQ(f2.begin2, 4u);
    SEQAN_ASSERT_EQ(f2.len, 1u);
    SEQAN_ASSERT_EQ(fragmentBegin(f2, 0), 0u);
    SEQAN_ASSERT_EQ(fragmentBegin(f2, 1), 4u);
    SEQAN_ASSERT_EQ(fragmentLength(f2, 0), 1u);
    SEQAN_ASSERT_EQ(fragmentLength(f2, 1), 1u);
    SEQAN_ASSERT_EQ(fragmentLength(f2), 1u);
    SEQAN_ASSERT_EQ(label(f2, str, 0), "a");
    SEQAN_ASSERT_EQ(label(f2, str, 1), "a");
    getProjectedPosition(f2, 0, 0, id2, pos2);
    SEQAN_ASSERT_EQ(pos2, 4u);
    SEQAN_ASSERT_EQ(id2, 1u);
    getProjectedPosition(f2, 1, 4, id2, pos2);
    SEQAN_ASSERT_EQ(pos2, 0u);
    SEQAN_ASSERT_EQ(id2, 0u);
}

SEQAN_DEFINE_TEST(test_align_reversable_fragment)
{
    using namespace seqan;

    typedef StringSet<CharString, Dependent<> > CharStringSet;
    typedef    Id<CharStringSet>::Type TId;
    typedef    Size<CharStringSet>::Type TSize;

    CharStringSet str;
    CharString str0("annual");
    assignValueById(str, str0);
    CharString str1("anneal");
    assignValueById(str, str1);

    // Reversable Fragment
    typedef Fragment<Size<Fragment<> >::Type, ExactReversableFragment<> > TRevFrag;
    TRevFrag fRev(0, 4, 1, 4, 2);

    SEQAN_ASSERT_EQ(fRev.seqId1, 0u);
    SEQAN_ASSERT_EQ(fRev.begin1, 4u);
    SEQAN_ASSERT_EQ(fRev.seqId2, 1u);
    SEQAN_ASSERT_EQ(fRev.begin2, 4u);
    SEQAN_ASSERT_EQ(fRev.len, 2u);
    SEQAN_ASSERT_EQ(fragmentBegin(fRev, 0), 4u);
    SEQAN_ASSERT_EQ(fragmentBegin(fRev, 1), 4u);
    SEQAN_ASSERT_EQ(fragmentLength(fRev, 0), 2u);
    SEQAN_ASSERT_EQ(fragmentLength(fRev, 1), 2u);
    SEQAN_ASSERT_EQ(sequenceId(fRev, 0), 0u);
    SEQAN_ASSERT_EQ(sequenceId(fRev, 1), 1u);
    SEQAN_ASSERT_EQ(label(fRev, str, 0), "al");
    SEQAN_ASSERT_EQ(label(fRev, str, 1), "al");
    SEQAN_ASSERT(!isReversed(fRev));

    TId id2;
    TSize pos2;
    getProjectedPosition(fRev, 0, 5, id2, pos2);
    SEQAN_ASSERT_EQ(pos2, 5);
    SEQAN_ASSERT_EQ(id2, 1u);
    getProjectedPosition(fRev, 1, 5, id2, pos2);
    SEQAN_ASSERT_EQ(pos2, 5);
    SEQAN_ASSERT_EQ(id2, 0);

    TRevFrag fRev2(0, 4, 1, 4, 2, true);

    SEQAN_ASSERT_EQ(fRev2.seqId1, 0u);
    SEQAN_ASSERT_EQ(fRev2.begin1, 4u);
    SEQAN_ASSERT_EQ(fRev2.seqId2, 1u);
    SEQAN_ASSERT_EQ(fRev2.begin2, 4u);
    SEQAN_ASSERT_EQ(fRev2.len, 2u);
    SEQAN_ASSERT_EQ(fragmentBegin(fRev2, 0), 4u);
    SEQAN_ASSERT_EQ(fragmentBegin(fRev2, 1), 4u);
    SEQAN_ASSERT_EQ(fragmentLength(fRev2, 0), 2u);
    SEQAN_ASSERT_EQ(fragmentLength(fRev2, 1), 2u);
    SEQAN_ASSERT_EQ(sequenceId(fRev2, 0), 0u);
    SEQAN_ASSERT_EQ(sequenceId(fRev2, 1), 1u);
    SEQAN_ASSERT_EQ(label(fRev2, str, 0), "al");
    SEQAN_ASSERT_EQ(label(fRev2, str, 1), "al");
    SEQAN_ASSERT(isReversed(fRev2));
    getProjectedPosition(fRev2, 0, 5, id2, pos2);
    SEQAN_ASSERT_EQ(pos2, 4u);
    SEQAN_ASSERT_EQ(id2, 1u);
    getProjectedPosition(fRev2, 1, 5, id2, pos2);
    SEQAN_ASSERT_EQ(pos2, 4u);
    SEQAN_ASSERT_EQ(id2, 0u);
}

#endif  // #ifndef SEQAN_TESTS_ALIGN_TEST_ALIGN_FRAGMENT_H_
