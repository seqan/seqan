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
// Tests for the paramChooser.h header.
//
// Currently, we only test the I/O routines from this file.
// ==========================================================================

#undef SEQAN_ENABLE_TESTING
#define SEQAN_ENABLE_TESTING 1

#include <sstream>

#include <seqan/basic.h>
#include <seqan/file.h>

#include "../splazers/paramChooser.h"

SEQAN_DEFINE_TEST(test_param_chooser_quality_distribution_from_prb_file)
{
    seqan::ParamChooserOptions pmOptions;
    pmOptions.totalN = 4;         // Read length.
    pmOptions.verbose = false;
    pmOptions.qualityCutoff = 10;  // Ignore reads with smaller average quality than threshold.
    
    seqan::String<double> qualDist;
    
    std::stringstream ss;
    ss << "40 -40 -40 -40    10 20 30 40    30 30 20 10    10 20 10 10\n"
       << "10 -40 -40 -40    10 10 0 10      9  0  0  0     5  5  5 10\n"
       << "40 -40 -40 -40    10 20 30 40    30 30 20 10    10 20 10 10";
    ss.seekg(0);
    
    qualityDistributionFromPrbFile(ss, qualDist, pmOptions);
    
    SEQAN_ASSERT_EQ(length(qualDist), 4u);
    SEQAN_ASSERT_IN_DELTA(qualDist[0], 9.999e-05, 1.0e-07);
    SEQAN_ASSERT_IN_DELTA(qualDist[1], 9.999e-05, 1.0e-07);
    SEQAN_ASSERT_IN_DELTA(qualDist[2], 0.000999001, 1.0e-07);
    SEQAN_ASSERT_IN_DELTA(qualDist[3], 0.00990099, 1.0e-07);
}

SEQAN_DEFINE_TEST(test_param_chooser_quality_distribution_from_fastq_file)
{
    seqan::ParamChooserOptions pmOptions;
    pmOptions.totalN = 4;         // Read length.
    pmOptions.verbose = false;
    pmOptions.qualityCutoff = 10;  // Ignore reads with smaller average quality than threshold.
    
    seqan::String<double> qualDist;
    
    std::stringstream ss;
    ss << "@1\n"
       << "CGAT\n"
       << "+\n"
       << "II?5\n"
       // << "@2\n"  // see prb test, but from FASTQ does not kick out bad reads
       // << "CGAT\n"
       // << "+\n"
       // << "+!*+\n"
       << "@3\n"
       << "CGAT\n"
       << "+\n"
       << "II?5";
    ss.seekg(0);
    
    qualityDistributionFromFastQFile(ss, qualDist, pmOptions);
    
    SEQAN_ASSERT_EQ(length(qualDist), 4u);
    SEQAN_ASSERT_IN_DELTA(qualDist[0], 9.999e-05, 1.0e-07);
    SEQAN_ASSERT_IN_DELTA(qualDist[1], 9.999e-05, 1.0e-07);
    SEQAN_ASSERT_IN_DELTA(qualDist[2], 0.000999001, 1.0e-07);
    SEQAN_ASSERT_IN_DELTA(qualDist[3], 0.00990099, 1.0e-07);
}


SEQAN_DEFINE_TEST(test_param_chooser_parse_gapped_params)
{
    seqan::ParamChooserOptions pmOptions;
    pmOptions.totalN = 4;         // Read length.
    pmOptions.verbose = false;
    pmOptions.chooseOneGappedOnly = false;
    pmOptions.chooseUngappedOnly = false;
    pmOptions.minThreshold = 3;
    pmOptions.optionLossRate = 0.01;
    
    seqan::RazerSOptions<> rOptions;

    std::stringstream ss;
    ss << "errors    shape               t        lossrate   PM\n"
       << "\n"
       << "0    11000000100100101        3        0          40895\n"
       << "1    11000000100100101        2        0          492286\n"
       << "2    11000000100100101        1        0.0293446  48417798\n"
       << "3    11000000100100101        1        0.195237   48417798\n"
       << "0    1111100001              10        0          40484\n"
       << "1    1111100001               7        0.163535   48622\n"
       << "1    1111100001               6        0.0497708  51290\n"
       << "1    1111100001               5        0          61068\n"
       << "1    1111101                  4        0          61068";

    ss.seekg(0);
    SEQAN_ASSERT(parseGappedParams(rOptions, ss, pmOptions));
    SEQAN_ASSERT_EQ(rOptions.shape, "1111100001");
    SEQAN_ASSERT_EQ(rOptions.threshold, 10);

    pmOptions.chooseOneGappedOnly = true;
    ss.clear();
    ss.seekg(0);
    SEQAN_ASSERT(parseGappedParams(rOptions, ss, pmOptions));
    SEQAN_ASSERT_EQ(rOptions.shape, "1111100001");
    SEQAN_ASSERT_EQ(rOptions.threshold, 10);

    pmOptions.chooseUngappedOnly = true;
    ss.clear();
    ss.seekg(0);
    SEQAN_ASSERT_NOT(parseGappedParams(rOptions, ss, pmOptions));
}

SEQAN_DEFINE_TEST(test_param_parse_shapes_from_file)
{
    seqan::ParamChooserOptions pmOptions;

    std::stringstream ss;
    ss << "    11000000100100101        \n"
       << "1111100001\n"
       << "1111100001";
    ss.seekg(0);

    seqan::String<seqan::CharString> shapes;
    SEQAN_ASSERT_EQ(parseShapesFromFile(shapes, ss, pmOptions), 3);

    SEQAN_ASSERT_EQ(length(shapes), 3u);
    SEQAN_ASSERT_EQ(shapes[0], "11000000100100101");
    SEQAN_ASSERT_EQ(shapes[1], "1111100001");
    SEQAN_ASSERT_EQ(shapes[2], "1111100001");
}

SEQAN_BEGIN_TESTSUITE(test_param_chooser)
{
    SEQAN_CALL_TEST(test_param_chooser_quality_distribution_from_prb_file);
    SEQAN_CALL_TEST(test_param_chooser_quality_distribution_from_fastq_file);

    SEQAN_CALL_TEST(test_param_chooser_parse_gapped_params);
    SEQAN_CALL_TEST(test_param_parse_shapes_from_file);
}
SEQAN_END_TESTSUITE
