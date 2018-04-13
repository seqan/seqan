// ==========================================================================
//                                 arg_parse
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
// Author: Stephan Aiche <stephan.aiche@fu-berlin.de>
// ==========================================================================

#define SEQAN_DEBUG


#include <seqan/basic.h>
#include <seqan/stream.h>

#include "test_arg_parse.h"
#include "test_arg_parse_argument.h"
#include "test_arg_parse_option.h"
#include "test_argument_parser.h"
#include "test_arg_parse_ctd_support.h"

using namespace seqan;

SEQAN_BEGIN_TESTSUITE(test_arg_parse)
{
    // Test an unset value
    SEQAN_CALL_TEST(test_unset_value);
    SEQAN_CALL_TEST(test_unset_values);

    // Call the int option tests
    SEQAN_CALL_TEST(test_int_short_argument);
    SEQAN_CALL_TEST(test_int_long_argument);
    SEQAN_CALL_TEST(test_non_int_argument);

    // Call the int64 option tests
    SEQAN_CALL_TEST(test_int64_short_argument);
    SEQAN_CALL_TEST(test_int64_long_argument);
    SEQAN_CALL_TEST(test_non_int64_argument);

    // Call the double option tests
    SEQAN_CALL_TEST(test_double_short_argument);
    SEQAN_CALL_TEST(test_double_long_argument);
    SEQAN_CALL_TEST(test_non_double_argument);
    SEQAN_CALL_TEST(test_double_scientific_notation);

    // Call the string option tests
    SEQAN_CALL_TEST(test_string_short_argument);
    SEQAN_CALL_TEST(test_string_long_argument);

    // Test missing argument with string option
    SEQAN_CALL_TEST(test_string_missing_argument);

    // Test list of arguments
    SEQAN_CALL_TEST(test_string_list);

    // Test min/max restrictions
    SEQAN_CALL_TEST(test_min_max_double_values_in_range);
    SEQAN_CALL_TEST(test_min_max_double_values_to_small);
    SEQAN_CALL_TEST(test_min_max_double_values_to_big);

    SEQAN_CALL_TEST(test_min_max_int_values_in_range);
    SEQAN_CALL_TEST(test_min_max_int_values_to_small);
    SEQAN_CALL_TEST(test_min_max_int_values_to_big);

    // Test allowed values.
    SEQAN_CALL_TEST(test_allowed_values_contained);
    SEQAN_CALL_TEST(test_allowed_values_not_contained);

    // Test file types
    SEQAN_CALL_TEST(test_input_file_short);
    SEQAN_CALL_TEST(test_input_file_long);
    SEQAN_CALL_TEST(test_input_file_missing);
    SEQAN_CALL_TEST(test_input_file_invalid_type);
    SEQAN_CALL_TEST(test_input_file_valid_type);
    SEQAN_CALL_TEST(test_input_file_extension);

    SEQAN_CALL_TEST(test_input_file_auto_file_ext_option);  // automatic --${name}-file-ext option.

    SEQAN_CALL_TEST(test_output_file_short);
    SEQAN_CALL_TEST(test_output_file_long);
    SEQAN_CALL_TEST(test_output_file_missing);
    SEQAN_CALL_TEST(test_output_file_invalid_type);
    SEQAN_CALL_TEST(test_output_file_valid_type);
    SEQAN_CALL_TEST(test_output_file_extension);
    SEQAN_CALL_TEST(test_output_file_extension_targz);
    SEQAN_CALL_TEST(test_output_file_explicit_extension_valid);
    SEQAN_CALL_TEST(test_output_file_explicit_extension_invalid);

    // Test for arguments.
    SEQAN_CALL_TEST(test_argument_string);
    SEQAN_CALL_TEST(test_argument_not_all_set);
    SEQAN_CALL_TEST(test_argument_double);
    SEQAN_CALL_TEST(test_argument_not_a_double);

    SEQAN_CALL_TEST(test_argument_auto_file_ext_option);  // automatic --${name}-file-ext option.

    // Test list of n-tuples arguments
    SEQAN_CALL_TEST(test_int_list_option);
    SEQAN_CALL_TEST(test_double_list_option);
    SEQAN_CALL_TEST(test_double_list_option_not_enough_arguments);

    // Test bools
    SEQAN_CALL_TEST(test_boolean_argument_on);
    SEQAN_CALL_TEST(test_boolean_argument_off);
    SEQAN_CALL_TEST(test_boolean_flags);
    SEQAN_CALL_TEST(test_combined_boolean_flags);
    SEQAN_CALL_TEST(test_long_short_flag_name);

    // cmd argument tests
    SEQAN_CALL_TEST(test_argument_string_type);
    SEQAN_CALL_TEST(test_argument_int_type);
    SEQAN_CALL_TEST(test_argument_int64_type);
    SEQAN_CALL_TEST(test_argument_double_type);
    SEQAN_CALL_TEST(test_argument_inputfile_type);
    SEQAN_CALL_TEST(test_argument_outputfile_type);
    SEQAN_CALL_TEST(test_argument_inputprefix_type);
    SEQAN_CALL_TEST(test_argument_outputprefix_type);
    SEQAN_CALL_TEST(test_argument_label);
    SEQAN_CALL_TEST(test_argument_invalid_cast);
    SEQAN_CALL_TEST(test_argument_min_max_boundaries);
    SEQAN_CALL_TEST(test_argument_valid_values);
    SEQAN_CALL_TEST(test_argument_valid_values_directories);

    SEQAN_CALL_TEST(test_argument_parser);
    SEQAN_CALL_TEST(test_parse_non_const_cstring);

    // default value test
    SEQAN_CALL_TEST(test_default_value);

    // conversion tests
    SEQAN_CALL_TEST(test_isDouble);
    SEQAN_CALL_TEST(test_isInt);

    // ctd tests
    SEQAN_CALL_TEST(test_arg_parse_ctd_support);
}
SEQAN_END_TESTSUITE
