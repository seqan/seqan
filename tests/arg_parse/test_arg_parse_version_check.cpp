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
// Author: Svenja Mehringer <svenja.mehringer@fu-berlin.de>
// ==========================================================================

// Locally enable version check for this test.
#if defined(SEQAN_DISABLE_VERSION_CHECK)
#undef SEQAN_DISABLE_VERSION_CHECK
#endif
#if defined(SEQAN_VERSION_CHECK_OPT_IN)
#undef SEQAN_VERSION_CHECK_OPT_IN
#endif

#define SEQAN_DEBUG
#define SEQAN_TEST_VERSION_CHECK_

#include "test_arg_parse_version_check.h"

using namespace seqan;

SEQAN_BEGIN_TESTSUITE(test_arg_parse)
{
    // tests to ensure that version check tests can run properly
    SEQAN_CALL_TEST(test_path_availability);
    SEQAN_CALL_TEST(test_delete_version_files);
    SEQAN_CALL_TEST(test_create_files);

    // version check tests
    // IMPORTANT: there always needs to be the test 'test_delete_version_files'
    //            in between tests to ensure that no former files interfere with
    //            the test results.
    SEQAN_CALL_TEST(test_option_on);
    SEQAN_CALL_TEST(test_delete_version_files);
    SEQAN_CALL_TEST(test_option_off);
    SEQAN_CALL_TEST(test_delete_version_files);
    SEQAN_CALL_TEST(test_smaller_seqan_version);
    SEQAN_CALL_TEST(test_delete_version_files);
    SEQAN_CALL_TEST(test_delete_version_files);
    SEQAN_CALL_TEST(test_smaller_app_version);
    SEQAN_CALL_TEST(test_delete_version_files);
    SEQAN_CALL_TEST(test_greater_app_version);
    SEQAN_CALL_TEST(test_delete_version_files);
    SEQAN_CALL_TEST(test_time_out);

    // clean up 
    _removeFilesFromPath();
}
SEQAN_END_TESTSUITE
