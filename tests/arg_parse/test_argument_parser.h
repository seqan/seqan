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
// Author: Stephan Aiche <stephan.aiche@fu-berlin.de>
// ==========================================================================
// Tests for arg_parse/arg_parse_option.h.
// ==========================================================================

#ifndef SEQAN_TESTS_ARG_PARSE_TEST_ARGUMENT_PARSER_H_
#define SEQAN_TESTS_ARG_PARSE_TEST_ARGUMENT_PARSER_H_

#include <seqan/basic.h>

#include <seqan/arg_parse/argument_parser.h>
#include <seqan/arg_parse/arg_parse_ctd_support.h>
#include <seqan/arg_parse/arg_parse_parse.h>

using namespace seqan;

// TODO(aiche): write tests for basic functions in ArgumentParser
SEQAN_DEFINE_TEST(test_argument_parser)
{
    ArgumentParser parser;
    addOption(parser, ArgParseOption("i", "integer", "help of an integer option", ArgParseArgument::INTEGER, "", true));
    _assignArgumentValue(getOption(parser, "i"), "10");
    _assignArgumentValue(getOption(parser, "i"), "12");
}


#endif // SEQAN_TESTS_ARG_PARSE_TEST_ARGUMENT_PARSER_H_
