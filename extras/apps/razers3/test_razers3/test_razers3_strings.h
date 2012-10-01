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
  @@By its name, this would test the header template/strings.h.@@
  ===========================================================================
*/

// @@Lines and sections with two at-chars are comments that should either
// be removed or replaced by tests defined by you.@@

// @@Replace symbol accordingly to your file name!@@
#ifndef TEST_TEMPLATE_TEST_TEMPLATE_STRINGS_H_
#define TEST_TEMPLATE_TEST_TEMPLATE_STRINGS_H_

#include <seqan/basic.h>     // @@For the testing infrastructure.@@
#include <seqan/sequence.h>  // @@Replace with your header under test.@@


// @@Simple test on string operations and demonstration of assertions.
// Replace this with your own test.@@
SEQAN_DEFINE_TEST(test_template_strings_example1)
{
    using namespace seqan;

    // Define some constant test data for comparison...
    CharString const kStr1 = "test 1";
    CharString const kStr2 = "test 2";

    // Append to a string and make equality assertion on the result.
    CharString myStr = "test ";
    append(myStr, "1");
    SEQAN_ASSERT_EQ(kStr1, myStr);

    // Demonstration of other assertions.
    SEQAN_ASSERT_GT(kStr2, myStr);
    SEQAN_ASSERT_GEQ(kStr2, myStr);
    SEQAN_ASSERT_LT(myStr, kStr2);
    SEQAN_ASSERT_LEQ(myStr, kStr2);
}

// @@Replace symbol accordingly to your file name!@@
#endif  // TEST_TEMPLATE_TEST_TEMPLATE_STRINGS_H_
