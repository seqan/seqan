// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2012, Knut Reinert, FU Berlin
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
// Tests for arg_parse/arg_parse_option.h.
// ==========================================================================

#ifndef SEQAN_CORE_TESTS_ARG_PARSE_TEST_ARG_PARSE_CTD_SUPPORT_H_
#define SEQAN_CORE_TESTS_ARG_PARSE_TEST_ARG_PARSE_CTD_SUPPORT_H_

#include <seqan/basic.h>

#include <seqan/arg_parse/argument_parser.h>
#include <seqan/arg_parse/arg_parse_ctd_support.h>
#include <seqan/arg_parse/arg_parse_parse.h>

SEQAN_DEFINE_TEST(test_arg_parse_ctd_support)
{
    // define a minimal argument parser
    seqan::ArgumentParser parser;

    setAppName(parser, "test_app");
    setVersion(parser, "0.1.0");
    setShortDescription(parser, "This is a test-app.");
    setCategory(parser, "SeqAn/Testing");

    addOption(parser, seqan::ArgParseOption("d", "double", "set a double option", seqan::ArgParseArgument::DOUBLE));
    addOption(parser, seqan::ArgParseOption("i", "integer", "set an integer option", seqan::ArgParseArgument::INTEGER));
    addOption(parser, seqan::ArgParseOption("s", "string", "set a string option", seqan::ArgParseArgument::STRING, "", true));
    addOption(parser, seqan::ArgParseOption("f", "in", "set an input file", seqan::ArgParseArgument::INPUTFILE));
    addOption(parser, seqan::ArgParseOption("o", "out", "set an output file", seqan::ArgParseArgument::OUTPUTFILE));

    // export ctd
    seqan::CharString outPath = SEQAN_TEMP_FILENAME();
    append(outPath, ".ctd");

    std::ofstream ofstream(toCString(outPath));
    writeCTD(parser, ofstream);
    ofstream.close();

    // compare ctd to expected
    seqan::CharString goldPath = SEQAN_PATH_TO_ROOT();
    append(goldPath, "/core/tests/arg_parse/test_app.ctd");

    SEQAN_ASSERT(seqan::_compareTextFiles(toCString(outPath), toCString(goldPath)));
}

#endif // SEQAN_CORE_TESTS_ARG_PARSE_TEST_ARG_PARSE_CTD_SUPPORT_H_
