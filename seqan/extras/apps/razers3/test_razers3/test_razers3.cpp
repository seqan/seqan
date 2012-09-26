/*==========================================================================
  SeqAn - The Library for Sequence Analysis
  http://www.seqan.de
  ===========================================================================
  Copyright (C) 2010

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

  ===========================================================================
  Author: @@Your Name@@ <@@Your Email@@>
  ===========================================================================
  @@Description of what is tested here@@
  ===========================================================================
*/

// @@Lines and sections with two at-chars are comments that should either
// be removed or replaced by tests defined by you.@@

#include <seqan/basic.h>  // @@Includes testing infrastructure.@@
#include "test_framework.h"
#include <seqan/file.h>   // @@Required to print strings in tests.@@

#include <iostream>
#include <fstream>
#include <seqan/find.h>
#include <seqan/index.h>
#include <seqan/store.h>

#include "../razers.h"

// @@Create one header with tests for each of your headers under test.@@
#include "test_razers3_others.h"


SEQAN_BEGIN_TESTSUITE(test_razers3)
{
    // Call tests.
    SEQAN_CALL_TEST(test_split_algorithm);

    SEQAN_CALL_TEST(test_to_bucket);

    SEQAN_CALL_TEST(test_radix_pass);

    SEQAN_CALL_TEST(test_radix_sort);

    // Verify checkpoints.
    // @@
    // Remove this line and accordingly add a call for each of your headers
    // under test.  When the SeqAn testing mode is enabled, each call to
    // SEQAN_CHECKPOINTS will be registered with the testing system.
    // Verification of checkpoints means that we will check for each registered
    // checkpoint to be hit.
    // @@
    SEQAN_VERIFY_CHECKPOINTS("projects/library/apps/razers3/razers_parallel_reads.h");
}
SEQAN_END_TESTSUITE
