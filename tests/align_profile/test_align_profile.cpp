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

#include <seqan/basic.h>
#include <seqan/file.h>

#include <seqan/align_profile.h>

// A simple test for the ProfileChar type.
SEQAN_DEFINE_TEST(test_align_profile_profile_test)
{
    typedef seqan::ProfileChar<seqan::Dna, int> TDnaProfile;
    seqan::String<TDnaProfile> profile = "CGAT";

    SEQAN_ASSERT_EQ(length(profile), 4u);

    SEQAN_ASSERT_EQ(profile[0].count[0], 0);
    SEQAN_ASSERT_EQ(profile[0].count[1], 1);
    SEQAN_ASSERT_EQ(profile[0].count[2], 0);
    SEQAN_ASSERT_EQ(profile[0].count[3], 0);
    SEQAN_ASSERT_EQ(profile[0].count[4], 0);

    SEQAN_ASSERT_EQ(profile[1].count[0], 0);
    SEQAN_ASSERT_EQ(profile[1].count[1], 0);
    SEQAN_ASSERT_EQ(profile[1].count[2], 1);
    SEQAN_ASSERT_EQ(profile[1].count[3], 0);
    SEQAN_ASSERT_EQ(profile[1].count[4], 0);

    SEQAN_ASSERT_EQ(profile[2].count[0], 1);
    SEQAN_ASSERT_EQ(profile[2].count[1], 0);
    SEQAN_ASSERT_EQ(profile[2].count[2], 0);
    SEQAN_ASSERT_EQ(profile[2].count[3], 0);
    SEQAN_ASSERT_EQ(profile[2].count[4], 0);

    SEQAN_ASSERT_EQ(profile[3].count[0], 0);
    SEQAN_ASSERT_EQ(profile[3].count[1], 0);
    SEQAN_ASSERT_EQ(profile[3].count[2], 0);
    SEQAN_ASSERT_EQ(profile[3].count[3], 1);
    SEQAN_ASSERT_EQ(profile[3].count[4], 0);
}

// Align a profile with a sequence.
SEQAN_DEFINE_TEST(test_align_profile_align_profile_sequence)
{
    typedef seqan::ProfileChar<seqan::Dna, int> TDnaProfile;
    typedef seqan::String<TDnaProfile> TProfileString;

    TProfileString profile = "CGAT";
    seqan::DnaString seq = "CGGAAT";

    seqan::Gaps<TProfileString> gapsH(profile);
    seqan::Gaps<seqan::DnaString> gapsV(seq);

    seqan::Score<int, seqan::ProfileSeqScore> sScheme(profile);

    int val = globalAlignment(gapsH, gapsV, sScheme, seqan::NeedlemanWunsch());
    SEQAN_ASSERT_EQ(val, -2097152);

    SEQAN_ASSERT_EQ(length(gapsV), 6u);
    for (unsigned i = 0; i < 6u; ++i)
        SEQAN_ASSERT_NOT_MSG(isGap(gapsV, i), "i == %u", i);

    SEQAN_ASSERT_EQ(length(gapsH), 6u);
    SEQAN_ASSERT_NOT(isGap(gapsH, 0));
    SEQAN_ASSERT(isGap(gapsH, 1));
    SEQAN_ASSERT_NOT(isGap(gapsH, 2));
    SEQAN_ASSERT(isGap(gapsH, 3));
    SEQAN_ASSERT_NOT(isGap(gapsH, 4));
    SEQAN_ASSERT_NOT(isGap(gapsH, 5));
}


// Align a profile with a sequence.
SEQAN_DEFINE_TEST(test_align_profile_align_profile_sequence_frac)
{
    typedef seqan::ProfileChar<seqan::Dna, int> TDnaProfile;
    typedef seqan::String<TDnaProfile> TProfileString;

    TProfileString profile = "CGAT";
    seqan::DnaString seq = "CGGAAT";

    seqan::Gaps<TProfileString> gapsH(profile);
    seqan::Gaps<seqan::DnaString> gapsV(seq);

    seqan::Score<int, seqan::ProfileSeqFracScore> sScheme(profile);

    int val = globalAlignment(gapsH, gapsV, sScheme, seqan::NeedlemanWunsch());
    SEQAN_ASSERT_EQ(val, -2097152);

    SEQAN_ASSERT_EQ(length(gapsV), 6u);
    for (unsigned i = 0; i < 6u; ++i)
        SEQAN_ASSERT_NOT_MSG(isGap(gapsV, i), "i == %u", i);

    SEQAN_ASSERT_EQ(length(gapsH), 6u);
    SEQAN_ASSERT_NOT(isGap(gapsH, 0));
    SEQAN_ASSERT(isGap(gapsH, 1));
    SEQAN_ASSERT_NOT(isGap(gapsH, 2));
    SEQAN_ASSERT(isGap(gapsH, 3));
    SEQAN_ASSERT_NOT(isGap(gapsH, 4));
    SEQAN_ASSERT_NOT(isGap(gapsH, 5));
}

// Call addToProfile() on a profile and a sequence.
//
// This will first align the sequence to the profile and then update the profile with the alignment information.
SEQAN_DEFINE_TEST(test_align_profile_add_to_profile)
{
    typedef seqan::ProfileChar<seqan::Dna, int> TDnaProfile;
    typedef seqan::String<TDnaProfile> TProfileString;

    TProfileString profile = "CGAT";
    seqan::DnaString seq = "CGGAAT";

    addToProfile(profile, seq);

    SEQAN_ASSERT_EQ(length(profile), 6u);

    SEQAN_ASSERT_EQ(profile[0].count[0], 0);
    SEQAN_ASSERT_EQ(profile[0].count[1], 2);
    SEQAN_ASSERT_EQ(profile[0].count[2], 0);
    SEQAN_ASSERT_EQ(profile[0].count[3], 0);
    SEQAN_ASSERT_EQ(profile[0].count[4], 0);

    SEQAN_ASSERT_EQ(profile[1].count[0], 0);
    SEQAN_ASSERT_EQ(profile[1].count[1], 0);
    SEQAN_ASSERT_EQ(profile[1].count[2], 2);
    SEQAN_ASSERT_EQ(profile[1].count[3], 0);
    SEQAN_ASSERT_EQ(profile[1].count[4], 0);

    SEQAN_ASSERT_EQ(profile[2].count[0], 0);
    SEQAN_ASSERT_EQ(profile[2].count[1], 0);
    SEQAN_ASSERT_EQ(profile[2].count[2], 1);
    SEQAN_ASSERT_EQ(profile[2].count[3], 0);
    SEQAN_ASSERT_EQ(profile[2].count[4], 1);

    SEQAN_ASSERT_EQ(profile[3].count[0], 1);
    SEQAN_ASSERT_EQ(profile[3].count[1], 0);
    SEQAN_ASSERT_EQ(profile[3].count[2], 0);
    SEQAN_ASSERT_EQ(profile[3].count[3], 0);
    SEQAN_ASSERT_EQ(profile[3].count[4], 1);

    SEQAN_ASSERT_EQ(profile[4].count[0], 2);
    SEQAN_ASSERT_EQ(profile[4].count[1], 0);
    SEQAN_ASSERT_EQ(profile[4].count[2], 0);
    SEQAN_ASSERT_EQ(profile[4].count[3], 0);
    SEQAN_ASSERT_EQ(profile[4].count[4], 0);

    SEQAN_ASSERT_EQ(profile[5].count[0], 0);
    SEQAN_ASSERT_EQ(profile[5].count[1], 0);
    SEQAN_ASSERT_EQ(profile[5].count[2], 0);
    SEQAN_ASSERT_EQ(profile[5].count[3], 2);
    SEQAN_ASSERT_EQ(profile[5].count[4], 0);
}

// Call addToProfile() three times.
SEQAN_DEFINE_TEST(test_align_profile_add_to_profile_multiple)
{
    typedef seqan::ProfileChar<seqan::Dna, int> TDnaProfile;
    typedef seqan::String<TDnaProfile> TProfileString;

    TProfileString profile = "CGAT";
    seqan::DnaString seq1 = "CGGAAT";
    seqan::DnaString seq2 = "CGGAT";
    seqan::DnaString seq3 = "AGAAT";

    addToProfile(profile, seq1);
    addToProfile(profile, seq2);
    addToProfile(profile, seq3);

    SEQAN_ASSERT_EQ(length(profile), 6u);

    SEQAN_ASSERT_EQ(profile[0].count[0], 1);
    SEQAN_ASSERT_EQ(profile[0].count[1], 3);
    SEQAN_ASSERT_EQ(profile[0].count[2], 0);
    SEQAN_ASSERT_EQ(profile[0].count[3], 0);
    SEQAN_ASSERT_EQ(profile[0].count[4], 0);

    SEQAN_ASSERT_EQ(profile[1].count[0], 0);
    SEQAN_ASSERT_EQ(profile[1].count[1], 0);
    SEQAN_ASSERT_EQ(profile[1].count[2], 4);
    SEQAN_ASSERT_EQ(profile[1].count[3], 0);
    SEQAN_ASSERT_EQ(profile[1].count[4], 0);

    SEQAN_ASSERT_EQ(profile[2].count[0], 1);
    SEQAN_ASSERT_EQ(profile[2].count[1], 0);
    SEQAN_ASSERT_EQ(profile[2].count[2], 2);
    SEQAN_ASSERT_EQ(profile[2].count[3], 0);
    SEQAN_ASSERT_EQ(profile[2].count[4], 1);

    SEQAN_ASSERT_EQ(profile[3].count[0], 1);
    SEQAN_ASSERT_EQ(profile[3].count[1], 0);
    SEQAN_ASSERT_EQ(profile[3].count[2], 0);
    SEQAN_ASSERT_EQ(profile[3].count[3], 0);
    SEQAN_ASSERT_EQ(profile[3].count[4], 3);

    SEQAN_ASSERT_EQ(profile[4].count[0], 4);
    SEQAN_ASSERT_EQ(profile[4].count[1], 0);
    SEQAN_ASSERT_EQ(profile[4].count[2], 0);
    SEQAN_ASSERT_EQ(profile[4].count[3], 0);
    SEQAN_ASSERT_EQ(profile[4].count[4], 0);

    SEQAN_ASSERT_EQ(profile[5].count[0], 0);
    SEQAN_ASSERT_EQ(profile[5].count[1], 0);
    SEQAN_ASSERT_EQ(profile[5].count[2], 0);
    SEQAN_ASSERT_EQ(profile[5].count[3], 4);
    SEQAN_ASSERT_EQ(profile[5].count[4], 0);
}

// Call addToProfile() on a profile and a sequence.
//
// This will first align the sequence to the profile and then update the profile with the alignment information.
SEQAN_DEFINE_TEST(test_align_profile_add_to_profile_banded)
{
    typedef seqan::ProfileChar<seqan::Dna, int> TDnaProfile;
    typedef seqan::String<TDnaProfile> TProfileString;

    TProfileString profile = "CGAT";
    seqan::DnaString seq = "CGGAAT";

    addToProfile(profile, seq, -3, 3);

    SEQAN_ASSERT_EQ(length(profile), 6u);

    SEQAN_ASSERT_EQ(profile[0].count[0], 0);
    SEQAN_ASSERT_EQ(profile[0].count[1], 2);
    SEQAN_ASSERT_EQ(profile[0].count[2], 0);
    SEQAN_ASSERT_EQ(profile[0].count[3], 0);
    SEQAN_ASSERT_EQ(profile[0].count[4], 0);

    SEQAN_ASSERT_EQ(profile[1].count[0], 0);
    SEQAN_ASSERT_EQ(profile[1].count[1], 0);
    SEQAN_ASSERT_EQ(profile[1].count[2], 2);
    SEQAN_ASSERT_EQ(profile[1].count[3], 0);
    SEQAN_ASSERT_EQ(profile[1].count[4], 0);

    SEQAN_ASSERT_EQ(profile[2].count[0], 0);
    SEQAN_ASSERT_EQ(profile[2].count[1], 0);
    SEQAN_ASSERT_EQ(profile[2].count[2], 1);
    SEQAN_ASSERT_EQ(profile[2].count[3], 0);
    SEQAN_ASSERT_EQ(profile[2].count[4], 1);

    SEQAN_ASSERT_EQ(profile[3].count[0], 1);
    SEQAN_ASSERT_EQ(profile[3].count[1], 0);
    SEQAN_ASSERT_EQ(profile[3].count[2], 0);
    SEQAN_ASSERT_EQ(profile[3].count[3], 0);
    SEQAN_ASSERT_EQ(profile[3].count[4], 1);

    SEQAN_ASSERT_EQ(profile[4].count[0], 2);
    SEQAN_ASSERT_EQ(profile[4].count[1], 0);
    SEQAN_ASSERT_EQ(profile[4].count[2], 0);
    SEQAN_ASSERT_EQ(profile[4].count[3], 0);
    SEQAN_ASSERT_EQ(profile[4].count[4], 0);

    SEQAN_ASSERT_EQ(profile[5].count[0], 0);
    SEQAN_ASSERT_EQ(profile[5].count[1], 0);
    SEQAN_ASSERT_EQ(profile[5].count[2], 0);
    SEQAN_ASSERT_EQ(profile[5].count[3], 2);
    SEQAN_ASSERT_EQ(profile[5].count[4], 0);
}

// Call addToProfile() three times.
SEQAN_DEFINE_TEST(test_align_profile_add_to_profile_multiple_banded)
{
    typedef seqan::ProfileChar<seqan::Dna, int> TDnaProfile;
    typedef seqan::String<TDnaProfile> TProfileString;

    TProfileString profile = "CGAT";
    seqan::DnaString seq1 = "CGGAAT";
    seqan::DnaString seq2 = "CGGAT";
    seqan::DnaString seq3 = "AGAAT";

    addToProfile(profile, seq1, -3, 3);
    addToProfile(profile, seq2, -3, 3);
    addToProfile(profile, seq3, -3, 3);

    SEQAN_ASSERT_EQ(length(profile), 6u);

    SEQAN_ASSERT_EQ(profile[0].count[0], 1);
    SEQAN_ASSERT_EQ(profile[0].count[1], 3);
    SEQAN_ASSERT_EQ(profile[0].count[2], 0);
    SEQAN_ASSERT_EQ(profile[0].count[3], 0);
    SEQAN_ASSERT_EQ(profile[0].count[4], 0);

    SEQAN_ASSERT_EQ(profile[1].count[0], 0);
    SEQAN_ASSERT_EQ(profile[1].count[1], 0);
    SEQAN_ASSERT_EQ(profile[1].count[2], 4);
    SEQAN_ASSERT_EQ(profile[1].count[3], 0);
    SEQAN_ASSERT_EQ(profile[1].count[4], 0);

    SEQAN_ASSERT_EQ(profile[2].count[0], 1);
    SEQAN_ASSERT_EQ(profile[2].count[1], 0);
    SEQAN_ASSERT_EQ(profile[2].count[2], 2);
    SEQAN_ASSERT_EQ(profile[2].count[3], 0);
    SEQAN_ASSERT_EQ(profile[2].count[4], 1);

    SEQAN_ASSERT_EQ(profile[3].count[0], 1);
    SEQAN_ASSERT_EQ(profile[3].count[1], 0);
    SEQAN_ASSERT_EQ(profile[3].count[2], 0);
    SEQAN_ASSERT_EQ(profile[3].count[3], 0);
    SEQAN_ASSERT_EQ(profile[3].count[4], 3);

    SEQAN_ASSERT_EQ(profile[4].count[0], 4);
    SEQAN_ASSERT_EQ(profile[4].count[1], 0);
    SEQAN_ASSERT_EQ(profile[4].count[2], 0);
    SEQAN_ASSERT_EQ(profile[4].count[3], 0);
    SEQAN_ASSERT_EQ(profile[4].count[4], 0);

    SEQAN_ASSERT_EQ(profile[5].count[0], 0);
    SEQAN_ASSERT_EQ(profile[5].count[1], 0);
    SEQAN_ASSERT_EQ(profile[5].count[2], 0);
    SEQAN_ASSERT_EQ(profile[5].count[3], 4);
    SEQAN_ASSERT_EQ(profile[5].count[4], 0);
}

SEQAN_BEGIN_TESTSUITE(test_align_profile)
{
    SEQAN_CALL_TEST(test_align_profile_profile_test);

    SEQAN_CALL_TEST(test_align_profile_align_profile_sequence);
    SEQAN_CALL_TEST(test_align_profile_align_profile_sequence_frac);

    SEQAN_CALL_TEST(test_align_profile_add_to_profile);
    SEQAN_CALL_TEST(test_align_profile_add_to_profile_multiple);
    SEQAN_CALL_TEST(test_align_profile_add_to_profile_banded);
    SEQAN_CALL_TEST(test_align_profile_add_to_profile_multiple_banded);
}
SEQAN_END_TESTSUITE
