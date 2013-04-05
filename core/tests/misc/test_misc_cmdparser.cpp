// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
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
// Tests for misc/misc_cmdparser.h.
// ==========================================================================

#include <seqan/misc/misc_cmdparser.h>  // Header under test.

#define SEQAN_DEBUG


#include <seqan/basic.h>
#include <iostream>

using namespace std;
using namespace seqan;

const char * A_INT_0 = "test";
const char * A_INT_1 = "-i";
const char * A_INT_2 = "--integer";
const char * A_INT_3 = "1";
const char * A_INT_4 = "-i1";
const char * A_INT_5 = "not-an-int";

const char * A_DOUBLE_0 = "test";
const char * A_DOUBLE_1 = "-d";
const char * A_DOUBLE_2 = "--double";
const char * A_DOUBLE_3 = "1.56";
const char * A_DOUBLE_4 = "-d1.56";
const char * A_DOUBLE_5 = "not-a-double";
const char * A_DOUBLE_6 = "6.0221418e23";

const char * A_STRING_0 = "test";
const char * A_STRING_1 = "-s";
const char * A_STRING_2 = "--string";
const char * A_STRING_3 = "this-is-a-string-value";
const char * A_STRING_4 = "-sthis-is-a-string-value";

const char * A_IN_FILE_0 = "test";
const char * A_IN_FILE_1 = "-i";
const char * A_IN_FILE_2 = "--in";
const char * A_IN_FILE_3 = "input.fasta";
const char * A_IN_FILE_4 = "-iinput.fasta";

const char * A_OUT_FILE_0 = "test";
const char * A_OUT_FILE_1 = "-o";
const char * A_OUT_FILE_2 = "--out";
const char * A_OUT_FILE_3 = "output.fasta";
const char * A_OUT_FILE_4 = "-ooutput.fasta";

const char * A_ARGUMENT_0 = "test";
const char * A_ARGUMENT_1 = "argument1";
const char * A_ARGUMENT_2 = "argument2";
const char * A_ARGUMENT_3 = "argument3";


// moved initialization of cmd parser out of the test functions
// to have single place to change in case of interface changes 
// or test extensions
void testInitDoubleParser(CommandLineParser & parser)
{
    addOption(parser, CommandLineOption("d", "double", "set a double option", OptionType::Double | OptionType::Label));
}

void testInitIntegerParser(CommandLineParser & parser)
{
    addOption(parser, CommandLineOption("i", "integer", "set an integer option", OptionType::Integer | OptionType::Label));
}

void testInitStringParser(CommandLineParser & parser)
{
    addOption(parser, CommandLineOption("s", "string", "set a string option", OptionType::String | OptionType::Label | OptionType::List));
}

void testInFileTypeParser(CommandLineParser & parser)
{
    addOption(parser, CommandLineOption("i", "in", "set a input file", OptionType::Label | OptionType::INPUTFILE));
}

void testOutFileTypeParser(CommandLineParser & parser)
{
    addOption(parser, CommandLineOption("o", "out", "set a output file", OptionType::Label | OptionType::OUTPUTFILE));
}

SEQAN_DEFINE_TEST(test_int_short_argument)
{

    CommandLineParser parser;
    testInitIntegerParser(parser);

    int argc = 3;
    const char * argv[3] = {A_INT_0, A_INT_1, A_INT_3};

    std::stringstream error_stream;

    SEQAN_ASSERT(parse(parser, argc, argv, error_stream));

    int integerValue = 0;
    SEQAN_ASSERT(getOptionValue(parser, "integer", integerValue));
    SEQAN_ASSERT_EQ(integerValue, 1);
    SEQAN_ASSERT_EQ(error_stream.str(), "");
}

SEQAN_DEFINE_TEST(test_int_long_argument)
{

    CommandLineParser parser;
    testInitIntegerParser(parser);

    int argc = 3;
    const char * argv[3] = {A_INT_0, A_INT_2, A_INT_3};

    std::stringstream error_stream;

    SEQAN_ASSERT(parse(parser, argc, argv, error_stream));

    int integerValue = 0;
    getOptionValue(parser, "integer", integerValue);
    SEQAN_ASSERT_EQ(integerValue, 1);
    SEQAN_ASSERT_EQ(error_stream.str(), "");
}

SEQAN_DEFINE_TEST(test_int_merged_short_argument)
{

    CommandLineParser parser;
    testInitIntegerParser(parser);

    int argc = 2;
    const char * argv[2] = {A_INT_0, A_INT_4};

    std::stringstream error_stream;

    SEQAN_ASSERT(parse(parser, argc, argv, error_stream));

    int integerValue = 0;
    getOptionValue(parser, "integer", integerValue);
    SEQAN_ASSERT_EQ(integerValue, 1);
    SEQAN_ASSERT_EQ(error_stream.str(), "");
}

SEQAN_DEFINE_TEST(test_non_int_argument)
{

    CommandLineParser parser;
    testInitIntegerParser(parser);

    int argc = 3;
    const char * argv[3] = {A_INT_0, A_INT_1, A_INT_5};

    std::stringstream error_stream;

    SEQAN_ASSERT(!parse(parser, argc, argv, error_stream));
    SEQAN_ASSERT_EQ(error_stream.str(), "test: \"not-an-int\" is not a valid integer value for '-i, --integer'\n");
}

SEQAN_DEFINE_TEST(test_double_short_argument)
{

    CommandLineParser parser;
    testInitDoubleParser(parser);

    int argc = 3;
    const char * argv[3] = {A_DOUBLE_0, A_DOUBLE_1, A_DOUBLE_3};

    std::stringstream error_stream;

    SEQAN_ASSERT(parse(parser, argc, argv, error_stream));

    double doubleValue = 0.0;
    getOptionValue(parser, "double", doubleValue);
    SEQAN_ASSERT_EQ(doubleValue, 1.56);
    SEQAN_ASSERT_EQ(error_stream.str(), "");

}

SEQAN_DEFINE_TEST(test_double_long_argument)
{

    CommandLineParser parser;
    testInitDoubleParser(parser);

    int argc = 3;
    const char * argv[3] = {A_DOUBLE_0, A_DOUBLE_2, A_DOUBLE_3};

    std::stringstream error_stream;

    SEQAN_ASSERT(parse(parser, argc, argv, error_stream));

    double doubleValue = 0.0;
    getOptionValue(parser, "double", doubleValue);
    SEQAN_ASSERT_EQ(doubleValue, 1.56);
    SEQAN_ASSERT_EQ(error_stream.str(), "");
}

SEQAN_DEFINE_TEST(test_double_merged_short_argument)
{
    CommandLineParser parser;
    testInitDoubleParser(parser);

    int argc = 2;
    const char * argv[2] = {A_DOUBLE_0, A_DOUBLE_4};

    std::stringstream error_stream;

    SEQAN_ASSERT(parse(parser, argc, argv, error_stream));

    double doubleValue = 0.0;
    getOptionValue(parser, "double", doubleValue);
    SEQAN_ASSERT_EQ(doubleValue, 1.56);
    SEQAN_ASSERT_EQ(error_stream.str(), "");
}

SEQAN_DEFINE_TEST(test_non_double_argument)
{

    CommandLineParser parser;
    testInitDoubleParser(parser);

    int argc = 3;
    const char * argv[3] = {A_DOUBLE_0, A_DOUBLE_1, A_DOUBLE_5};

    std::stringstream error_stream;

    SEQAN_ASSERT(!parse(parser, argc, argv, error_stream));
    SEQAN_ASSERT_EQ(error_stream.str(), "test: \"not-a-double\" is not a valid double value for '-d, --double'\n");
}

SEQAN_DEFINE_TEST(test_double_scientific_notation)
{

    CommandLineParser parser;
    testInitDoubleParser(parser);

    int argc = 3;
    const char * argv[3] = {A_DOUBLE_0, A_DOUBLE_1, A_DOUBLE_6};

    std::stringstream error_stream;

    double doubleValue = 0.0;
    SEQAN_ASSERT(parse(parser, argc, argv, error_stream));
    getOptionValue(parser, "double", doubleValue);
    SEQAN_ASSERT_EQ(doubleValue, 6.0221418e23);
    SEQAN_ASSERT_EQ(error_stream.str(), "");
}

SEQAN_DEFINE_TEST(test_string_short_argument)
{

    CommandLineParser parser;
    testInitStringParser(parser);

    int argc = 3;
    const char * argv[3] = {A_STRING_0, A_STRING_1, A_STRING_3};

    std::stringstream error_stream;

    SEQAN_ASSERT(parse(parser, argc, argv, error_stream));

    CharString value;
    getOptionValue(parser, "string", value);
    SEQAN_ASSERT_EQ(value, "this-is-a-string-value");
    SEQAN_ASSERT_EQ(error_stream.str(), "");
}

SEQAN_DEFINE_TEST(test_string_long_argument)
{

    CommandLineParser parser;
    testInitStringParser(parser);

    int argc = 3;
    const char * argv[3] = {A_STRING_0, A_STRING_2, A_STRING_3};

    std::stringstream error_stream;
    SEQAN_ASSERT(parse(parser, argc, argv, error_stream));

    CharString value;
    getOptionValue(parser, "string", value);
    SEQAN_ASSERT_EQ(value, "this-is-a-string-value");
    SEQAN_ASSERT_EQ(error_stream.str(), "");
}

SEQAN_DEFINE_TEST(test_string_merged_short_argument)
{

    CommandLineParser parser;
    testInitStringParser(parser);

    int argc = 2;
    const char * argv[2] = {A_STRING_0, A_STRING_4};

    std::stringstream error_stream;

    SEQAN_ASSERT(parse(parser, argc, argv, error_stream));

    CharString value;
    getOptionValue(parser, "string", value);
    SEQAN_ASSERT_EQ(value, "this-is-a-string-value");
    SEQAN_ASSERT_EQ(error_stream.str(), "");

}

SEQAN_DEFINE_TEST(test_string_missing_argument)
{

    CommandLineParser parser;
    testInitStringParser(parser);

    int argc = 2;
    const char * argv[2] = {A_STRING_0, A_STRING_2};

    std::stringstream error_stream;

    SEQAN_ASSERT(!parse(parser, argc, argv, error_stream));
    SEQAN_ASSERT_EQ(error_stream.str(), "test: '-s, --string' requires 1 value(s)\n");
}

SEQAN_DEFINE_TEST(test_string_list)
{
    CommandLineParser parser;
    testInitStringParser(parser);

    int argc = 7;
    const char * argv[7] = {A_STRING_0, A_STRING_1, A_STRING_3, A_STRING_2, A_STRING_3, A_STRING_1, A_STRING_3};

    std::stringstream error_stream;

    SEQAN_ASSERT(parse(parser, argc, argv, error_stream));

    String<CharString> const & values = getOptionValues(parser, "string");

    SEQAN_ASSERT_EQ(length(values), 3u);

    for (typename Size<String<CharString> >::Type i = 0; i < length(values); ++i)
    {
        SEQAN_ASSERT_EQ(value(values, i), "this-is-a-string-value");
    }

    SEQAN_ASSERT_EQ(error_stream.str(), "");
}

SEQAN_DEFINE_TEST(test_min_max_double_values_in_range)
{
    CommandLineParser parser;
    testInitDoubleParser(parser);

    setMinValue(parser, "double", "1.0");
    setMaxValue(parser, "double", "2.0");

    int argc = 3;
    const char * argv[3] = {A_DOUBLE_0, A_DOUBLE_2, A_DOUBLE_3};

    std::stringstream error_stream;

    SEQAN_ASSERT(parse(parser, argc, argv, error_stream));

    double doubleValue = 0.0;
    getOptionValue(parser, "double", doubleValue);
    SEQAN_ASSERT_EQ(doubleValue, 1.56);
    SEQAN_ASSERT_EQ(error_stream.str(), "");
}

SEQAN_DEFINE_TEST(test_min_max_double_values_to_small)
{
    CommandLineParser parser;
    testInitDoubleParser(parser);

    setMinValue(parser, "double", "1.6");

    int argc = 3;
    const char * argv[3] = {A_DOUBLE_0, A_DOUBLE_2, A_DOUBLE_3};

    std::stringstream error_stream;

    SEQAN_ASSERT(!parse(parser, argc, argv, error_stream));
    SEQAN_ASSERT_EQ(error_stream.str(), "-d, --double: given argument \"1.56\" is not in the required range [1.6:+inf]\n");
}

SEQAN_DEFINE_TEST(test_min_max_double_values_to_big)
{
    CommandLineParser parser;
    testInitDoubleParser(parser);

    setMaxValue(parser, "double", "1.5");

    int argc = 3;
    const char * argv[3] = {A_DOUBLE_0, A_DOUBLE_2, A_DOUBLE_3};

    std::stringstream error_stream;

    SEQAN_ASSERT(!parse(parser, argc, argv, error_stream));
    SEQAN_ASSERT_EQ(error_stream.str(), "-d, --double: given argument \"1.56\" is not in the required range [-inf:1.5]\n");
}

SEQAN_DEFINE_TEST(test_min_max_int_values_in_range)
{
    CommandLineParser parser;
    testInitIntegerParser(parser);

    setMinValue(parser, "integer", "-10");
    setMaxValue(parser, "integer", "2");

    int argc = 3;
    const char * argv[3] = {A_INT_0, A_INT_2, A_INT_3};

    std::stringstream error_stream;

    SEQAN_ASSERT(parse(parser, argc, argv, error_stream));

    int integerValue = 0;
    getOptionValue(parser, "integer", integerValue);
    SEQAN_ASSERT_EQ(integerValue, 1);
    SEQAN_ASSERT_EQ(error_stream.str(), "");
}

SEQAN_DEFINE_TEST(test_min_max_int_values_to_small)
{
    CommandLineParser parser;
    testInitIntegerParser(parser);

    setMinValue(parser, "integer", "3");

    int argc = 3;
    const char * argv[3] = {A_INT_0, A_INT_2, A_INT_3};

    std::stringstream error_stream;

    SEQAN_ASSERT(!parse(parser, argc, argv, error_stream));
    SEQAN_ASSERT_EQ(error_stream.str(), "-i, --integer: given argument \"1\" is not in the required range [3:+inf]\n");
}

SEQAN_DEFINE_TEST(test_min_max_int_values_to_big)
{
    CommandLineParser parser;
    testInitIntegerParser(parser);

    setMaxValue(parser, "integer", "-3");

    int argc = 3;
    const char * argv[3] = {A_INT_0, A_INT_2, A_INT_3};

    std::stringstream error_stream;

    SEQAN_ASSERT(!parse(parser, argc, argv, error_stream));
    SEQAN_ASSERT_EQ(error_stream.str(), "-i, --integer: given argument \"1\" is not in the required range [-inf:-3]\n");
}

SEQAN_DEFINE_TEST(test_allowed_values_contained)
{
    CommandLineParser parser;
    testInitStringParser(parser);

    setValidValues(parser, "string", CharString("a b c this-is-a-string-value"));

    int argc = 2;
    const char * argv[2] = {A_STRING_0, A_STRING_4};

    std::stringstream error_stream;

    SEQAN_ASSERT(parse(parser, argc, argv, error_stream));

    CharString value;
    getOptionValue(parser, "string", value);
    SEQAN_ASSERT_EQ(value, "this-is-a-string-value");
    SEQAN_ASSERT_EQ(error_stream.str(), "");
}

SEQAN_DEFINE_TEST(test_allowed_values_not_contained)
{
    CommandLineParser parser;
    testInitStringParser(parser);

    setValidValues(parser, "string", CharString("a b c"));

    int argc = 2;
    const char * argv[2] = {A_STRING_0, A_STRING_4};

    std::stringstream error_stream;

    SEQAN_ASSERT(!parse(parser, argc, argv, error_stream));
    SEQAN_ASSERT_EQ(error_stream.str(), "-s, --string: given argument \"this-is-a-string-value\" is not a valid value [a, b, c]\n");
}

SEQAN_DEFINE_TEST(test_input_file_short)
{
    CommandLineParser parser;
    testInFileTypeParser(parser);

    int argc = 3;
    const char * argv[3] = {A_IN_FILE_0, A_IN_FILE_1, A_IN_FILE_3};

    std::stringstream error_stream;

    SEQAN_ASSERT(parse(parser, argc, argv, error_stream));

    CharString value;
    SEQAN_ASSERT(getOptionValue(parser, "in", value));
    SEQAN_ASSERT_EQ(value, "input.fasta");
    SEQAN_ASSERT_EQ(error_stream.str(), "");
}

SEQAN_DEFINE_TEST(test_input_file_long)
{
    CommandLineParser parser;
    testInFileTypeParser(parser);

    int argc = 3;
    const char * argv[3] = {A_IN_FILE_0, A_IN_FILE_2, A_IN_FILE_3};

    std::stringstream error_stream;

    SEQAN_ASSERT(parse(parser, argc, argv, error_stream));

    CharString value;
    SEQAN_ASSERT(getOptionValue(parser, "in", value));
    SEQAN_ASSERT_EQ(value, "input.fasta");
    SEQAN_ASSERT_EQ(error_stream.str(), "");
}

SEQAN_DEFINE_TEST(test_input_file_merged)
{
    CommandLineParser parser;
    testInFileTypeParser(parser);

    int argc = 2;
    const char * argv[2] = {A_IN_FILE_0, A_IN_FILE_4};

    std::stringstream error_stream;

    SEQAN_ASSERT(parse(parser, argc, argv, error_stream));

    CharString value;
    SEQAN_ASSERT(getOptionValue(parser, "in", value));
    SEQAN_ASSERT_EQ(value, "input.fasta");
    SEQAN_ASSERT_EQ(error_stream.str(), "");
}

SEQAN_DEFINE_TEST(test_input_file_missing)
{
    CommandLineParser parser;
    testInFileTypeParser(parser);

    int argc = 2;
    const char * argv[2] = {A_IN_FILE_0, A_IN_FILE_1};

    std::stringstream error_stream;

    SEQAN_ASSERT(!parse(parser, argc, argv, error_stream));
    SEQAN_ASSERT_EQ(error_stream.str(), "test: '-i, --in' requires 1 value(s)\n");
}

SEQAN_DEFINE_TEST(test_input_file_invalid_type)
{
    CommandLineParser parser;
    testInFileTypeParser(parser);

    setValidValues(parser, "in", CharString("FASTA fa"));

    int argc = 2;
    const char * argv[2] = {A_IN_FILE_0, A_IN_FILE_4};

    std::stringstream error_stream;

    SEQAN_ASSERT(!parse(parser, argc, argv, error_stream));
    SEQAN_ASSERT_EQ(error_stream.str(), "-i, --in: given argument \"input.fasta\" is not a valid file type [FASTA, fa]\n");
}

SEQAN_DEFINE_TEST(test_input_file_valid_type)
{
    CommandLineParser parser;
    testInFileTypeParser(parser);

    setValidValues(parser, "in", CharString("fasta FASTA fa"));

    int argc = 2;
    const char * argv[2] = {A_IN_FILE_0, A_IN_FILE_4};

    std::stringstream error_stream;

    SEQAN_ASSERT(parse(parser, argc, argv, error_stream));

    CharString value;
    SEQAN_ASSERT(getOptionValue(parser, "in", value));
    SEQAN_ASSERT_EQ(value, "input.fasta");
    SEQAN_ASSERT_EQ(error_stream.str(), "");
}

SEQAN_DEFINE_TEST(test_output_file_short)
{
    CommandLineParser parser;
    testOutFileTypeParser(parser);

    int argc = 3;
    const char * argv[3] = {A_OUT_FILE_0, A_OUT_FILE_1, A_OUT_FILE_3};

    std::stringstream error_stream;

    SEQAN_ASSERT(parse(parser, argc, argv, error_stream));

    CharString value;
    SEQAN_ASSERT(getOptionValue(parser, "out", value));
    SEQAN_ASSERT_EQ(value, "output.fasta");
    SEQAN_ASSERT_EQ(error_stream.str(), "");
}

SEQAN_DEFINE_TEST(test_output_file_long)
{
    CommandLineParser parser;
    testOutFileTypeParser(parser);

    int argc = 3;
    const char * argv[3] = {A_OUT_FILE_0, A_OUT_FILE_2, A_OUT_FILE_3};

    std::stringstream error_stream;

    SEQAN_ASSERT(parse(parser, argc, argv, error_stream));

    CharString value;
    SEQAN_ASSERT(getOptionValue(parser, "out", value));
    SEQAN_ASSERT_EQ(value, "output.fasta");
    SEQAN_ASSERT_EQ(error_stream.str(), "");
}

SEQAN_DEFINE_TEST(test_output_file_merged)
{
    CommandLineParser parser;
    testOutFileTypeParser(parser);

    int argc = 2;
    const char * argv[2] = {A_OUT_FILE_0, A_OUT_FILE_4};

    std::stringstream error_stream;

    SEQAN_ASSERT(parse(parser, argc, argv, error_stream));

    CharString value;
    SEQAN_ASSERT(getOptionValue(parser, "out", value));
    SEQAN_ASSERT_EQ(value, "output.fasta");
    SEQAN_ASSERT_EQ(error_stream.str(), "");
}

SEQAN_DEFINE_TEST(test_output_file_missing)
{
    CommandLineParser parser;
    testOutFileTypeParser(parser);

    int argc = 2;
    const char * argv[2] = {A_OUT_FILE_0, A_OUT_FILE_1};

    std::stringstream error_stream;

    SEQAN_ASSERT(!parse(parser, argc, argv, error_stream));
    SEQAN_ASSERT_EQ(error_stream.str(), "test: '-o, --out' requires 1 value(s)\n");
}

SEQAN_DEFINE_TEST(test_output_file_invalid_type)
{
    CommandLineParser parser;
    testOutFileTypeParser(parser);

    setValidValues(parser, "out", CharString("FASTA fa"));

    int argc = 2;
    const char * argv[2] = {A_OUT_FILE_0, A_OUT_FILE_4};

    std::stringstream error_stream;

    SEQAN_ASSERT(!parse(parser, argc, argv, error_stream));
    SEQAN_ASSERT_EQ(error_stream.str(), "-o, --out: given argument \"output.fasta\" is not a valid file type [FASTA, fa]\n");
}

SEQAN_DEFINE_TEST(test_output_file_valid_type)
{
    CommandLineParser parser;
    testOutFileTypeParser(parser);

    setValidValues(parser, "out", CharString("fasta FASTA fa"));

    int argc = 2;
    const char * argv[2] = {A_OUT_FILE_0, A_OUT_FILE_4};

    std::stringstream error_stream;

    SEQAN_ASSERT(parse(parser, argc, argv, error_stream));

    CharString value;
    SEQAN_ASSERT(getOptionValue(parser, "out", value));
    SEQAN_ASSERT_EQ(value, "output.fasta");
    SEQAN_ASSERT_EQ(error_stream.str(), "");
}

SEQAN_DEFINE_TEST(test_commandline_argument_not_enough_arguments)
{
    CommandLineParser parser;
    requiredArguments(parser, 2);

    int argc = 2;
    const char * argv[2] = {A_ARGUMENT_0, A_ARGUMENT_1};

    std::stringstream error_stream;

    SEQAN_ASSERT(!parse(parser, argc, argv, error_stream));
    SEQAN_ASSERT_EQ(error_stream.str(), "");
}

SEQAN_DEFINE_TEST(test_commandline_argument_required_arguments_set)
{
    CommandLineParser parser;
    requiredArguments(parser, 2);

    int argc = 3;
    const char * argv[3] = {A_ARGUMENT_0, A_ARGUMENT_1, A_ARGUMENT_2};

    std::stringstream error_stream;

    SEQAN_ASSERT(parse(parser, argc, argv, error_stream));
    SEQAN_ASSERT_EQ(error_stream.str(), "");

    String<CharString> arguments = getArgumentValues(parser);
    SEQAN_ASSERT_EQ(length(arguments), 2u);
    SEQAN_ASSERT_EQ(value(arguments, 0), A_ARGUMENT_1);
    SEQAN_ASSERT_EQ(value(arguments, 1), A_ARGUMENT_2);
}


SEQAN_DEFINE_TEST(test_isDouble)
{
    CharString a = "this is not a double";
    CharString b = "2.5";
    CharString c = "-45.5245";
    CharString d = "6.0221418e23";
    CharString e = "-45.5245aeeeb";

    SEQAN_ASSERT(!_isDouble(a));
    SEQAN_ASSERT(_isDouble(b));
    SEQAN_ASSERT(_isDouble(c));
    SEQAN_ASSERT(_isDouble(d));
    SEQAN_ASSERT(!_isDouble(e));
}

SEQAN_DEFINE_TEST(test_isInt)
{
    CharString a = "this is not an int";
    CharString b = "2";
    CharString c = "-4253252";
    CharString d = "6aaefgeag";

    SEQAN_ASSERT(!_isInt(a));
    SEQAN_ASSERT(_isInt(b));
    SEQAN_ASSERT(_isInt(c));
    SEQAN_ASSERT(!_isInt(d));
}

SEQAN_BEGIN_TESTSUITE(test_misc_cmdparser)
{
    // Call the int option tests
    SEQAN_CALL_TEST(test_int_short_argument);
    SEQAN_CALL_TEST(test_int_long_argument);
    SEQAN_CALL_TEST(test_int_merged_short_argument);
    SEQAN_CALL_TEST(test_non_int_argument);

    // Call the double option tests
    SEQAN_CALL_TEST(test_double_short_argument);
    SEQAN_CALL_TEST(test_double_long_argument);
    SEQAN_CALL_TEST(test_double_merged_short_argument);
    SEQAN_CALL_TEST(test_non_double_argument);
    SEQAN_CALL_TEST(test_double_scientific_notation);

    // Call the string option tests
    SEQAN_CALL_TEST(test_string_short_argument);
    SEQAN_CALL_TEST(test_string_long_argument);
    SEQAN_CALL_TEST(test_string_merged_short_argument);

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
    SEQAN_CALL_TEST(test_input_file_merged);
    SEQAN_CALL_TEST(test_input_file_missing);
    SEQAN_CALL_TEST(test_input_file_invalid_type);
    SEQAN_CALL_TEST(test_input_file_valid_type);

    SEQAN_CALL_TEST(test_output_file_short);
    SEQAN_CALL_TEST(test_output_file_long);
    SEQAN_CALL_TEST(test_output_file_merged);
    SEQAN_CALL_TEST(test_output_file_missing);
    SEQAN_CALL_TEST(test_output_file_invalid_type);
    SEQAN_CALL_TEST(test_output_file_valid_type);

    // test command line arguments
    SEQAN_CALL_TEST(test_commandline_argument_not_enough_arguments);
    SEQAN_CALL_TEST(test_commandline_argument_required_arguments_set);

    // conversion tests
    SEQAN_CALL_TEST(test_isDouble);
    SEQAN_CALL_TEST(test_isInt);
}
SEQAN_END_TESTSUITE
