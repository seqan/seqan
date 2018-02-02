// ==========================================================================
//                        test_alignment_dp_profile.h
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

#ifndef SANDBOX_RMAERKER_TESTS_ALIGN2_TEST_ALIGNMENT_DP_PROFILE_H_
#define SANDBOX_RMAERKER_TESTS_ALIGN2_TEST_ALIGNMENT_DP_PROFILE_H_

#include <seqan/basic.h>

#include <seqan/align.h>

void testAlignmentDPProfileIsFreeEndGap()
{

    using namespace seqan;

    {
        typedef DPProfile_<GlobalAlignment_<>, LinearGaps, TracebackOn<> > TDPProfile;

        bool result0 = IsFreeEndGap_<TDPProfile, DPFirstRow>::VALUE;
        bool result1 = IsFreeEndGap_<TDPProfile, DPFirstColumn>::VALUE;
        bool result2 = IsFreeEndGap_<TDPProfile, DPLastRow>::VALUE;
        bool result3 = IsFreeEndGap_<TDPProfile, DPLastColumn>::VALUE;

        SEQAN_ASSERT_EQ(result0, false);
        SEQAN_ASSERT_EQ(result1, false);
        SEQAN_ASSERT_EQ(result2, false);
        SEQAN_ASSERT_EQ(result3, false);
    }

    {
        typedef DPProfile_<GlobalAlignment_<>, LinearGaps, TracebackOn<> > const TDPProfile;

        bool result0 = IsFreeEndGap_<TDPProfile, DPFirstRow>::VALUE;
        bool result1 = IsFreeEndGap_<TDPProfile, DPFirstColumn>::VALUE;
        bool result2 = IsFreeEndGap_<TDPProfile, DPLastRow>::VALUE;
        bool result3 = IsFreeEndGap_<TDPProfile, DPLastColumn>::VALUE;

        SEQAN_ASSERT_EQ(result0, false);
        SEQAN_ASSERT_EQ(result1, false);
        SEQAN_ASSERT_EQ(result2, false);
        SEQAN_ASSERT_EQ(result3, false);
    }

    {
        typedef DPProfile_<GlobalAlignment_<FreeEndGaps_<True, False, True, False> >, LinearGaps, TracebackOn<> > TDPProfile;

        bool result0 = IsFreeEndGap_<TDPProfile, DPFirstRow>::VALUE;
        bool result1 = IsFreeEndGap_<TDPProfile, DPFirstColumn>::VALUE;
        bool result2 = IsFreeEndGap_<TDPProfile, DPLastRow>::VALUE;
        bool result3 = IsFreeEndGap_<TDPProfile, DPLastColumn>::VALUE;

        SEQAN_ASSERT_EQ(result0, true);
        SEQAN_ASSERT_EQ(result1, false);
        SEQAN_ASSERT_EQ(result2, true);
        SEQAN_ASSERT_EQ(result3, false);
    }

    {
        typedef DPProfile_<GlobalAlignment_<FreeEndGaps_<True, False, True, False> >, LinearGaps, TracebackOn<> > const TDPProfile;

        bool result0 = IsFreeEndGap_<TDPProfile, DPFirstRow>::VALUE;
        bool result1 = IsFreeEndGap_<TDPProfile, DPFirstColumn>::VALUE;
        bool result2 = IsFreeEndGap_<TDPProfile, DPLastRow>::VALUE;
        bool result3 = IsFreeEndGap_<TDPProfile, DPLastColumn>::VALUE;

        SEQAN_ASSERT_EQ(result0, true);
        SEQAN_ASSERT_EQ(result1, false);
        SEQAN_ASSERT_EQ(result2, true);
        SEQAN_ASSERT_EQ(result3, false);
    }

    {
        typedef DPProfile_<GlobalAlignment_<FreeEndGaps_<True, True, True, True> >, LinearGaps, TracebackOn<> > TDPProfile;

        bool result0 = IsFreeEndGap_<TDPProfile, DPFirstRow>::VALUE;
        bool result1 = IsFreeEndGap_<TDPProfile, DPFirstColumn>::VALUE;
        bool result2 = IsFreeEndGap_<TDPProfile, DPLastRow>::VALUE;
        bool result3 = IsFreeEndGap_<TDPProfile, DPLastColumn>::VALUE;

        SEQAN_ASSERT_EQ(result0, true);
        SEQAN_ASSERT_EQ(result1, true);
        SEQAN_ASSERT_EQ(result2, true);
        SEQAN_ASSERT_EQ(result3, true);
    }

    {
        typedef DPProfile_<GlobalAlignment_<FreeEndGaps_<True, True, True, True> >, LinearGaps, TracebackOn<> > const TDPProfile;

        bool result0 = IsFreeEndGap_<TDPProfile, DPFirstRow>::VALUE;
        bool result1 = IsFreeEndGap_<TDPProfile, DPFirstColumn>::VALUE;
        bool result2 = IsFreeEndGap_<TDPProfile, DPLastRow>::VALUE;
        bool result3 = IsFreeEndGap_<TDPProfile, DPLastColumn>::VALUE;

        SEQAN_ASSERT_EQ(result0, true);
        SEQAN_ASSERT_EQ(result1, true);
        SEQAN_ASSERT_EQ(result2, true);
        SEQAN_ASSERT_EQ(result3, true);
    }

    {
        typedef DPProfile_<LocalAlignment_<>, LinearGaps, TracebackOn<> > TDPProfile;

        bool result0 = IsFreeEndGap_<TDPProfile, DPFirstRow>::VALUE;
        bool result1 = IsFreeEndGap_<TDPProfile, DPFirstColumn>::VALUE;
        bool result2 = IsFreeEndGap_<TDPProfile, DPLastRow>::VALUE;
        bool result3 = IsFreeEndGap_<TDPProfile, DPLastColumn>::VALUE;

        SEQAN_ASSERT_EQ(result0, true);
        SEQAN_ASSERT_EQ(result1, true);
        SEQAN_ASSERT_EQ(result2, true);
        SEQAN_ASSERT_EQ(result3, true);
    }

    {
        typedef DPProfile_<LocalAlignment_<>, LinearGaps, TracebackOn<> > const TDPProfile;

        bool result0 = IsFreeEndGap_<TDPProfile, DPFirstRow>::VALUE;
        bool result1 = IsFreeEndGap_<TDPProfile, DPFirstColumn>::VALUE;
        bool result2 = IsFreeEndGap_<TDPProfile, DPLastRow>::VALUE;
        bool result3 = IsFreeEndGap_<TDPProfile, DPLastColumn>::VALUE;

        SEQAN_ASSERT_EQ(result0, true);
        SEQAN_ASSERT_EQ(result1, true);
        SEQAN_ASSERT_EQ(result2, true);
        SEQAN_ASSERT_EQ(result3, true);
    }

    {
        bool result0 = IsFreeEndGap_<Nothing, DPFirstRow>::VALUE;
        bool result1 = IsFreeEndGap_<Nothing, DPFirstColumn>::VALUE;
        bool result2 = IsFreeEndGap_<Nothing, DPLastRow>::VALUE;
        bool result3 = IsFreeEndGap_<Nothing, DPLastColumn>::VALUE;

        SEQAN_ASSERT_EQ(result0, false);
        SEQAN_ASSERT_EQ(result1, false);
        SEQAN_ASSERT_EQ(result2, false);
        SEQAN_ASSERT_EQ(result3, false);
    }

}

void testAlignmentDPProfileIsGlobal()
{
    using namespace seqan;

    {
        typedef DPProfile_<GlobalAlignment_<>, LinearGaps, TracebackOn<> > TDPProfile;

        bool result0 = IsGlobalAlignment_<TDPProfile>::VALUE;
        bool result1 = IsGlobalAlignment_<GlobalAlignment_<> >::VALUE;

        SEQAN_ASSERT_EQ(result0, true);
        SEQAN_ASSERT_EQ(result1, true);
    }

    {
        typedef DPProfile_<LocalAlignment_<>, AffineGaps, TracebackOff> TDPProfile;

        bool result0 = IsGlobalAlignment_<TDPProfile>::VALUE;
        bool result1 = IsGlobalAlignment_<LocalAlignment_<> >::VALUE;

        SEQAN_ASSERT_EQ(result0, false);
        SEQAN_ASSERT_EQ(result1, false);
    }

    {
        typedef DPProfile_<GlobalAlignment_<>, LinearGaps, TracebackOn<> > const TDPProfile;

        bool result0 = IsGlobalAlignment_<TDPProfile>::VALUE;
        bool result1 = IsGlobalAlignment_<GlobalAlignment_<> >::VALUE;

        SEQAN_ASSERT_EQ(result0, true);
        SEQAN_ASSERT_EQ(result1, true);
    }

    {
        typedef DPProfile_<LocalAlignment_<>, AffineGaps, TracebackOff> const TDPProfile;

        bool result0 = IsGlobalAlignment_<TDPProfile>::VALUE;
        bool result1 = IsGlobalAlignment_<LocalAlignment_<> >::VALUE;

        SEQAN_ASSERT_EQ(result0, false);
        SEQAN_ASSERT_EQ(result1, false);
    }
}

void testAlignmentDPProfileIsLocal()
{
    using namespace seqan;

    {
        typedef DPProfile_<GlobalAlignment_<>, LinearGaps, TracebackOn<> > TDPProfile;

        bool result0 = IsLocalAlignment_<TDPProfile>::VALUE;
        bool result1 = IsLocalAlignment_<GlobalAlignment_<> >::VALUE;

        SEQAN_ASSERT_EQ(result0, false);
        SEQAN_ASSERT_EQ(result1, false);
    }

    {
        typedef DPProfile_<LocalAlignment_<>, AffineGaps, TracebackOff> TDPProfile;

        bool result0 = IsLocalAlignment_<TDPProfile>::VALUE;
        bool result1 = IsLocalAlignment_<LocalAlignment_<> >::VALUE;

        SEQAN_ASSERT_EQ(result0, true);
        SEQAN_ASSERT_EQ(result1, true);
    }

    {
        typedef DPProfile_<GlobalAlignment_<>, LinearGaps, TracebackOn<> > const TDPProfile;

        bool result0 = IsLocalAlignment_<TDPProfile>::VALUE;
        bool result1 = IsLocalAlignment_<GlobalAlignment_<> >::VALUE;

        SEQAN_ASSERT_EQ(result0, false);
        SEQAN_ASSERT_EQ(result1, false);
    }

    {
        typedef DPProfile_<LocalAlignment_<>, AffineGaps, TracebackOff> const TDPProfile;

        bool result0 = IsLocalAlignment_<TDPProfile>::VALUE;
        bool result1 = IsLocalAlignment_<LocalAlignment_<> >::VALUE;

        SEQAN_ASSERT_EQ(result0, true);
        SEQAN_ASSERT_EQ(result1, true);
    }
}

void testAlignmetnDPProfileIsTracebackEnabled()
{
    using namespace seqan;

    {
        typedef DPProfile_<GlobalAlignment_<>, LinearGaps, TracebackOn<> > TDPProfile;

        bool result0 = IsTracebackEnabled_<TDPProfile>::VALUE;
        bool result1 = IsTracebackEnabled_<TracebackOn<> >::VALUE;

        SEQAN_ASSERT_EQ(result0, true);
        SEQAN_ASSERT_EQ(result1, true);
    }

    {
        typedef DPProfile_<LocalAlignment_<>, AffineGaps, TracebackOff> TDPProfile;

        bool result0 = IsTracebackEnabled_<TDPProfile>::VALUE;
        bool result1 = IsTracebackEnabled_<TracebackOff>::VALUE;

        SEQAN_ASSERT_EQ(result0, false);
        SEQAN_ASSERT_EQ(result1, false);
    }

    {
        typedef DPProfile_<GlobalAlignment_<>, LinearGaps, TracebackOn<> > const TDPProfile;

        bool result0 = IsTracebackEnabled_<TDPProfile>::VALUE;
        bool result1 = IsTracebackEnabled_<TracebackOn<> >::VALUE;

        SEQAN_ASSERT_EQ(result0, true);
        SEQAN_ASSERT_EQ(result1, true);
    }

    {
        typedef DPProfile_<LocalAlignment_<>, AffineGaps, TracebackOff> const TDPProfile;

        bool result0 = IsTracebackEnabled_<TDPProfile>::VALUE;
        bool result1 = IsTracebackEnabled_<TracebackOff>::VALUE;

        SEQAN_ASSERT_EQ(result0, false);
        SEQAN_ASSERT_EQ(result1, false);
    }
}

SEQAN_DEFINE_TEST(test_alignment_dp_profile_is_free_end_gaps)
{
    testAlignmentDPProfileIsFreeEndGap();
}

SEQAN_DEFINE_TEST(test_alignment_dp_profile_is_global_alignment)
{
    testAlignmentDPProfileIsGlobal();
}

SEQAN_DEFINE_TEST(test_alignment_dp_profile_is_local_alignment)
{
    testAlignmentDPProfileIsLocal();
}

SEQAN_DEFINE_TEST(test_alignment_dp_profile_is_traceback_enabled)
{
    testAlignmetnDPProfileIsTracebackEnabled();
}

#endif  // #ifndef SANDBOX_RMAERKER_TESTS_ALIGN2_TEST_ALIGNMENT_DP_PROFILE_H_
