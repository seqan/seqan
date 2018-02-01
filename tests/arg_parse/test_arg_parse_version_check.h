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
// Author: Svenja Mehringer <ssvenja.mehringer@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_TESTS_ARG_PARSE_TEST_ARG_PARSE_VERSION_CHECK_H_
#define SEQAN_TESTS_ARG_PARSE_TEST_ARG_PARSE_VERSION_CHECK_H_

#include <seqan/arg_parse.h>

#include <iostream>
#include <fstream>
#include <chrono>

struct TestVersionCheck_
{
    static const std::chrono::duration<long int>::rep TIME_NOW;

    static const std::string PATH; // See arg_parse_version_check.h for _getPath()
    static const std::string APP_NAME;
    static const std::string APP_VERSION_FILENAME;
    static const std::string APP_TIMESTAMP_FILENAME;

    static constexpr const char * const OPTION_VERSION_CHECK = "--version-check";
    static constexpr const char * const OPTION_OFF = "OFF";
    static constexpr const char * const OPTION_ON = "ON";
};

const std::chrono::duration<long int>::rep TestVersionCheck_::TIME_NOW = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now().time_since_epoch()).count();
const std::string TestVersionCheck_::PATH = seqan::_getPath();
const std::string TestVersionCheck_::APP_NAME = "test_version_check_" + std::to_string(TIME_NOW); // avoid name conflicts
const std::string TestVersionCheck_::APP_VERSION_FILENAME     = TestVersionCheck_::PATH + "/" +
                                                                TestVersionCheck_::APP_NAME +
                                                                static_cast<std::string>(".version");
const std::string TestVersionCheck_::APP_TIMESTAMP_FILENAME   = TestVersionCheck_::PATH + "/" +
                                                                TestVersionCheck_::APP_NAME +
                                                                static_cast<std::string>("_usr.timestamp");

namespace seqan {

inline ArgumentParser::ParseResult _simulateArgumentParser(String<CharString> & stream_result,
                                                           bool & app_call_succeeded,
                                                           int argc,
                                                           const char ** argv)
{
    ArgumentParser parser;
    setAppName(parser, TestVersionCheck_::APP_NAME);
    setVersion(parser, "2.3.4");

    std::stringstream err_stream;
    std::stringstream out_stream;

    ArgumentParser::ParseResult res = parse(parser, argc, argv, out_stream, err_stream);

    appendValue(stream_result, out_stream.str());
    appendValue(stream_result, err_stream.str());

    // call future.get() to artificially wait for the thread to finish and avoid
    // any interference with following tests
    if(parser.appVersionCheckFuture.valid())
            app_call_succeeded = parser.appVersionCheckFuture.get();

    return res;
}

inline void _removeFilesFromPath()
{
    std::remove(TestVersionCheck_::APP_VERSION_FILENAME.c_str());
    std::remove(TestVersionCheck_::APP_TIMESTAMP_FILENAME.c_str());
}

template <typename TMessage>
inline void _createFile(std::string const & filename, TMessage const & message)
{
    std::ofstream out_file(filename.c_str());
    if (out_file.is_open())
    {
        out_file << message;
        out_file.close();
    }
}

// even if the homedir might not be writable at least the tmp dir should be
SEQAN_DEFINE_TEST(test_path_availability)
{
    SEQAN_ASSERT(!TestVersionCheck_::PATH.empty()); // TODO:: correct when function is working!
}

// other tests might fail if the removal of intermediate results is not possible
SEQAN_DEFINE_TEST(test_delete_version_files)
{
    _removeFilesFromPath();
    SEQAN_ASSERT(!fileExists(TestVersionCheck_::APP_VERSION_FILENAME.c_str()));
    SEQAN_ASSERT(!fileExists(TestVersionCheck_::APP_TIMESTAMP_FILENAME.c_str()));
}

SEQAN_DEFINE_TEST(test_create_files)
{
    _createFile(TestVersionCheck_::APP_VERSION_FILENAME.c_str(), "20.5.9");
    _createFile(TestVersionCheck_::APP_TIMESTAMP_FILENAME.c_str(), TestVersionCheck_::TIME_NOW);

    SEQAN_ASSERT(fileExists(TestVersionCheck_::APP_VERSION_FILENAME.c_str()));
    SEQAN_ASSERT(fileExists(TestVersionCheck_::APP_TIMESTAMP_FILENAME.c_str()));

    _removeFilesFromPath(); // clear files again
    SEQAN_ASSERT(!fileExists(TestVersionCheck_::APP_VERSION_FILENAME.c_str()));
    SEQAN_ASSERT(!fileExists(TestVersionCheck_::APP_TIMESTAMP_FILENAME.c_str()));
}

SEQAN_DEFINE_TEST(test_option_on)
{
    int argc(3);
    const char * argv[3] = {TestVersionCheck_::APP_NAME.c_str(), TestVersionCheck_::OPTION_VERSION_CHECK,
                            TestVersionCheck_::OPTION_ON};
    String<CharString> stream_result;
    bool app_call_succeeded(false);

    ArgumentParser::ParseResult res = _simulateArgumentParser(stream_result, app_call_succeeded, argc, argv);

    SEQAN_ASSERT_EQ(res, ArgumentParser::PARSE_OK);
    SEQAN_ASSERT_EQ(stream_result[0], "");
    SEQAN_ASSERT_EQ(stream_result[1], "");

    // make sure that all files now exist
    SEQAN_ASSERT(fileExists(TestVersionCheck_::APP_TIMESTAMP_FILENAME.c_str()));
    if (app_call_succeeded)
        SEQAN_ASSERT(fileExists(TestVersionCheck_::APP_VERSION_FILENAME.c_str()));
}

SEQAN_DEFINE_TEST(test_option_off)
{
    int argc(3);
    const char * argv[3] = {TestVersionCheck_::APP_NAME.c_str(), TestVersionCheck_::OPTION_VERSION_CHECK,
                            TestVersionCheck_::OPTION_OFF};
    String<CharString> stream_result;
    bool app_call_succeeded(false);

    ArgumentParser::ParseResult res = _simulateArgumentParser(stream_result, app_call_succeeded, argc, argv);

    SEQAN_ASSERT_EQ(res, ArgumentParser::PARSE_OK);
    SEQAN_ASSERT_EQ(stream_result[0], "");
    SEQAN_ASSERT_EQ(stream_result[1], "");

    // make sure that no files exist
    SEQAN_ASSERT(!fileExists(TestVersionCheck_::APP_TIMESTAMP_FILENAME.c_str()));
    SEQAN_ASSERT(!fileExists(TestVersionCheck_::APP_VERSION_FILENAME.c_str()));
}

// case: the current argument parser has a smaller seqan version than is present in the version file
SEQAN_DEFINE_TEST(test_smaller_seqan_version)
{
    int argc(3);
    const char * argv[3] = {TestVersionCheck_::APP_NAME.c_str(), TestVersionCheck_::OPTION_VERSION_CHECK,
                            TestVersionCheck_::OPTION_ON};
    String<CharString> stream_result;
    bool app_call_succeeded(false);

    // create version file with euqal app version and a greater seqan version than the current (2.3.4)
    std::stringstream version_info;
    version_info << "2.3.4" << std::endl << "20.5.9";
    _createFile(TestVersionCheck_::APP_VERSION_FILENAME.c_str(), version_info.str());

    // create timestamp file that dates one day before current to trigger a message
    _createFile(TestVersionCheck_::APP_TIMESTAMP_FILENAME.c_str(), TestVersionCheck_::TIME_NOW - 86401); // one day = 86400 seconds

    ArgumentParser::ParseResult res = _simulateArgumentParser(stream_result, app_call_succeeded, argc, argv);

    SEQAN_ASSERT_EQ(res, ArgumentParser::PARSE_OK);
    SEQAN_ASSERT_EQ(stream_result[0], "");
    SEQAN_ASSERT_EQ(stream_result[1], VersionControlTags_<>::MESSAGE_SEQAN_UPDATE);
}

// case: the current argument parser has a smaller app version than is present in the version file
SEQAN_DEFINE_TEST(test_smaller_app_version)
{
    int argc(3);
    const char * argv[3] = {TestVersionCheck_::APP_NAME.c_str(), TestVersionCheck_::OPTION_VERSION_CHECK,
                            TestVersionCheck_::OPTION_ON};
    String<CharString> stream_result;
    bool app_call_succeeded(false);

    // create version file with equal seqan version and a greater app version than the current (2.3.4)
    std::string seqan_version = std::to_string(SEQAN_VERSION_MAJOR) + "." +
                                std::to_string(SEQAN_VERSION_MINOR) + "." +
                                std::to_string(SEQAN_VERSION_PATCH);
    std::stringstream version_info;
    version_info << "20.5.9" << std::endl << seqan_version;
    _createFile(TestVersionCheck_::APP_VERSION_FILENAME.c_str(), version_info.str());

    // create timestamp file that dates one day before current to trigger a message
    _createFile(TestVersionCheck_::APP_TIMESTAMP_FILENAME.c_str(), TestVersionCheck_::TIME_NOW - 86401); // one day = 86400 seconds

    ArgumentParser::ParseResult res = _simulateArgumentParser(stream_result, app_call_succeeded, argc, argv);

    SEQAN_ASSERT_EQ(res, ArgumentParser::PARSE_OK);
    SEQAN_ASSERT_EQ(stream_result[0], "");
    SEQAN_ASSERT_EQ(stream_result[1], VersionControlTags_<>::MESSAGE_APP_UPDATE);
}

// case: the current argument parser has a greater app version than is present in the version file
SEQAN_DEFINE_TEST(test_greater_app_version)
{
    int argc(3);
    const char * argv[3] = {TestVersionCheck_::APP_NAME.c_str(), TestVersionCheck_::OPTION_VERSION_CHECK,
                            TestVersionCheck_::OPTION_ON};
    String<CharString> stream_result;
    bool app_call_succeeded(false);

    // create version file with equal seqan version and a smaller app version than the current (2.3.4)
    std::string seqan_version = std::to_string(SEQAN_VERSION_MAJOR) + "." +
                                std::to_string(SEQAN_VERSION_MINOR) + "." +
                                std::to_string(SEQAN_VERSION_PATCH);
    std::stringstream version_info;
    version_info << "1.5.9" << std::endl << seqan_version;
    _createFile(TestVersionCheck_::APP_VERSION_FILENAME.c_str(), version_info.str());

    // create timestamp file that dates one day before current to trigger a message
    _createFile(TestVersionCheck_::APP_TIMESTAMP_FILENAME.c_str(), TestVersionCheck_::TIME_NOW - 86401); // one day = 86400 seconds
 
    ArgumentParser::ParseResult res = _simulateArgumentParser(stream_result, app_call_succeeded, argc, argv);

    SEQAN_ASSERT_EQ(res, ArgumentParser::PARSE_OK);
    SEQAN_ASSERT_EQ(stream_result[0], "");
    SEQAN_ASSERT_EQ(stream_result[1], VersionControlTags_<>::MESSAGE_REGISTERED_APP_UPDATE);
}

SEQAN_DEFINE_TEST(test_time_out)
{
    int argc(3);
    const char * argv[3] = {TestVersionCheck_::APP_NAME.c_str(), TestVersionCheck_::OPTION_VERSION_CHECK,
                            TestVersionCheck_::OPTION_ON};
    String<CharString> stream_result;
    bool app_call_succeeded(false);

    // create timestamp files
    _createFile(TestVersionCheck_::APP_TIMESTAMP_FILENAME.c_str(), TestVersionCheck_::TIME_NOW);

    ArgumentParser::ParseResult res = _simulateArgumentParser(stream_result, app_call_succeeded, argc, argv);

    SEQAN_ASSERT_EQ(res, ArgumentParser::PARSE_OK);
    SEQAN_ASSERT_EQ(stream_result[0], "");
    SEQAN_ASSERT_EQ(stream_result[1], "");

    SEQAN_ASSERT(!fileExists(TestVersionCheck_::APP_VERSION_FILENAME.c_str()));
}

} // namespace seqan

#endif  // SEQAN_TESTS_ARG_PARSE_TEST_ARG_PARSE_VERSION_CHECK_H_
