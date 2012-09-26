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
// Author: @@Your Name@@ <@@Your Email@@>
// ==========================================================================
// @@Description of what is tested here@@
// ==========================================================================

// @@Lines and sections with two at-chars are comments that should either
// be removed or replaced by tests defined by you.@@

#include <seqan/basic.h>  // @@Includes testing infrastructure.@@
#include <seqan/file.h>   // @@Required to print strings in tests.@@

// @@Create one header with tests for each of your headers under test.@@
#include "test_template_strings.h"
#include "test_template_others.h"


SEQAN_BEGIN_TESTSUITE(test_template) {
    // Call tests.
    SEQAN_CALL_TEST(test_template_strings_example1);
    SEQAN_CALL_TEST(test_template_others_example1);

    // Verify checkpoints.
    // @@
    // Remove this line and accordingly add a call for each of your headers
    // under test.  When the SeqAn testing mode is enabled, each call to
    // SEQAN_CHECKPOINTS will be registered with the testing system.
    // Verification of checkpoints means that we will check for each registered
    // checkpoint to be hit.
    // @@
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/sequence/lexical.h");
}
SEQAN_END_TESTSUITE
