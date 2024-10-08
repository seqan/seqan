// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2024, Knut Reinert, FU Berlin
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
// Tests for arg_parse/arg_parse_ctd_support.h.
// ==========================================================================

#ifndef SEQAN_TESTS_ARG_PARSE_TEST_ARG_PARSE_CTD_SUPPORT_H_
#define SEQAN_TESTS_ARG_PARSE_TEST_ARG_PARSE_CTD_SUPPORT_H_

#include <seqan/basic.h>

#include <seqan/arg_parse/argument_parser.h>
#include <seqan/arg_parse/arg_parse_ctd_support.h>
#include <seqan/arg_parse/arg_parse_parse.h>

SEQAN_DEFINE_TEST(test_arg_parse_ctd_support)
{
    // define a minimal argument parser
    seqan2::ArgumentParser parser;

    setAppName(parser, "test_app");
    setVersion(parser, "0.1.0");
    setShortDescription(parser, "This is a test-app.");
    setCategory(parser, "SeqAn/Testing");
    addDescription(parser, "This is the first line of our test description.");
    addDescription(parser, "The second one contains formating <\\fIbla\\fP>.");

    addOption(parser, seqan2::ArgParseOption("b", "bool", "set a bool option", seqan2::ArgParseArgument::BOOL));
    addOption(parser, seqan2::ArgParseOption("d", "double", "set a double option", seqan2::ArgParseArgument::DOUBLE));
    addOption(parser, seqan2::ArgParseOption("i", "integer", "set an integer option", seqan2::ArgParseArgument::INTEGER));
    setMinValue(parser, "i", "1");
    setMaxValue(parser, "i", "10");
    addOption(parser, seqan2::ArgParseOption("j", "int64", "set a 64 bit integer option", seqan2::ArgParseArgument::INT64));
    addOption(parser, seqan2::ArgParseOption("s", "string", "set a string option", seqan2::ArgParseArgument::STRING, "", true));
    setValidValues(parser, "s", "a b c");
    addOption(parser, seqan2::ArgParseOption("f", "in", "set an input file", seqan2::ArgParseArgument::INPUT_FILE));
    setValidValues(parser, "f", "fasta");
    addOption(parser, seqan2::ArgParseOption("o", "out", "set an output file", seqan2::ArgParseArgument::OUTPUT_FILE));
    setValidValues(parser, "o", "sam");
    addOption(parser, seqan2::ArgParseOption("ip", "input-prefix-option", "set an input prefix", seqan2::ArgParseArgument::INPUT_PREFIX));
    setValidValues(parser, "ip", "btx");
    addOption(parser, seqan2::ArgParseOption("op", "output-prefix-option", "set an output prefix", seqan2::ArgParseArgument::OUTPUT_PREFIX));
    setValidValues(parser, "output-prefix-option", "blub");

    addOption(parser, seqan2::ArgParseOption("hi", "hidden", "a hidden option - will not appear in the ctd", seqan2::ArgParseArgument::STRING));
    hideOption(parser, "hi");

    addOption(parser, seqan2::ArgParseOption("ad", "advanced", "an advanced option - will appear as advanced in the ctd, too", seqan2::ArgParseArgument::STRING));
    setAdvanced(parser, "advanced");

    addArgument(parser, seqan2::ArgParseArgument(ArgParseArgument::DOUBLE, "DOUBLE"));
    setHelpText(parser, 0, "Double Argument");
    addArgument(parser, seqan2::ArgParseArgument(ArgParseArgument::STRING, "STRING"));
    setHelpText(parser, 1, "String Argument");
    addArgument(parser, seqan2::ArgParseArgument(ArgParseArgument::STRING, "DOC"));
    setHelpText(parser, 2, "Documentated Argument with \\fBformating\\fP");
    addArgument(parser, seqan2::ArgParseArgument(ArgParseArgument::OUTPUT_FILE, "OUTPUT-FILE"));
    setHelpText(parser, 3, "Testing output file arguments");

    // export ctd
    seqan2::CharString outPath = SEQAN_TEMP_FILENAME();
    append(outPath, ".ctd");

    std::ofstream ofstream(toCString(outPath), std::ofstream::out | std::ofstream::binary);
    writeCTD(parser, ofstream);
    ofstream.close();

    // compare ctd to expected
    seqan2::CharString goldPath = getAbsolutePath("/tests/arg_parse/test_app.ctd");

    SEQAN_ASSERT(seqan2::_compareTextFilesAlt(toCString(outPath), toCString(goldPath)));
}

#endif // SEQAN_TESTS_ARG_PARSE_TEST_ARG_PARSE_CTD_SUPPORT_H_
