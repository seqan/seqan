// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
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
// Tests for the paramChooser.h header.
//
// Currently, we only test the I/O routines from this file.
// ==========================================================================

#undef SEQAN_ENABLE_TESTING
#define SEQAN_ENABLE_TESTING 1

#include <sstream>

#include <seqan/basic.h>
#include <seqan/file.h>

#include "../../../core/apps/splazers/paramChooser.h"

SEQAN_DEFINE_TEST(test_param_chooser_quality_distribution_from_prb_file)
{
    seqan::ParamChooserOptions pmOptions;
    pmOptions.totalN = 4;         // Read length.
    pmOptions.qualityCutoff = 10;  // Ignore reads with smaller average quality than threshold.
    
    seqan::String<double> qualDist;
    
    std::stringstream ss;
    ss << "40 -40 -40 -40    10 20 30 40    30 30 20 10    10 20 10 10\n"
       << "10 -40 -40 -40    10 10 0 10      9  0  0  0     5  5  5 10\n"
       << "40 -40 -40 -40    10 20 30 40    30 30 20 10    10 20 10 10";
    ss.seekg(0);
    
    SEQAN_ASSERT_EQ(qualityDistributionFromPrbFile(ss, qualDist, pmOptions), 0);
    
    SEQAN_ASSERT_EQ(length(qualDist), 4u);
    SEQAN_ASSERT_IN_DELTA(qualDist[0], 9.999e-05, 1.0e-07);
    SEQAN_ASSERT_IN_DELTA(qualDist[1], 9.999e-05, 1.0e-07);
    SEQAN_ASSERT_IN_DELTA(qualDist[2], 0.000999001, 1.0e-07);
    SEQAN_ASSERT_IN_DELTA(qualDist[3], 0.00990099, 1.0e-07);
}

SEQAN_BEGIN_TESTSUITE(test_param_chooser)
{
    SEQAN_CALL_TEST(test_param_chooser_quality_distribution_from_prb_file);
}
SEQAN_END_TESTSUITE
