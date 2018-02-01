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

#ifndef SEQAN_INCLUDE_ARG_PARSE_VERSION_CHECK_H_
#define SEQAN_INCLUDE_ARG_PARSE_VERSION_CHECK_H_

#include <sys/stat.h>

#include <iostream>
#include <fstream>
#include <regex>
#include <future>
#include <chrono>

namespace seqan
{

// ==========================================================================
// Forwards
// ==========================================================================

// NOTE(rrahn): In-file forward for function call operator.
struct VersionCheck;
inline void _checkForNewerVersion(VersionCheck &, std::promise<bool>);
inline std::string _getPath();
constexpr const char * _getOS();
constexpr const char * _getBitSys();

// ==========================================================================
// Tags, Classes, Enums
// ==========================================================================

template <typename TVoidSpec = void>
struct VersionControlTags_
{
    static constexpr char const * const SEQAN_NAME         = "seqan";
    static constexpr char const * const UNREGISTERED_APP   = "UNREGISTERED_APP";

    static constexpr char const * const MESSAGE_SEQAN_UPDATE =
        "[SEQAN INFO] :: A new SeqAn version is available online.\n"
        "[SEQAN INFO] :: Please visit www.seqan.de for an update or inform the developer of this app.\n"
        "[SEQAN INFO] :: If you don't wish to receive further notifications, set --version-check OFF.\n\n";
    static constexpr char const * const MESSAGE_APP_UPDATE =
        "[APP INFO] :: A new version of this application is now available.\n"
        "[APP INFO] :: Visit www.seqan.de for updates of official SeqAn applications.\n"
        "[APP INFO] :: If you don't wish to receive further notifications, set --version-check OFF.\n\n";
    static constexpr char const * const MESSAGE_UNREGISTERED_APP =
        "[SEQAN INFO] :: Thank you for using SeqAn!\n"
        "[SEQAN INFO] :: Do you wish to register your app for update notifications?\n"
        "[SEQAN INFO] :: Just send an email to support@seqan.de with your app name and version number.\n"
        "[SEQAN INFO] :: If you don't wish to receive further notifications, set --version-check OFF.\n\n";
    static constexpr char const * const MESSAGE_REGISTERED_APP_UPDATE =
        "[APP INFO] :: We noticed the app version you use is newer than the one registered with us.\n"
        "[APP INFO] :: Please send us an email with the new version so we can correct it (support@seqan.de)\n\n";
};

template <typename TVoidSpec>
constexpr char const * const VersionControlTags_<TVoidSpec>::SEQAN_NAME;
template <typename TVoidSpec>
constexpr char const * const VersionControlTags_<TVoidSpec>::UNREGISTERED_APP;
template <typename TVoidSpec>
constexpr char const * const VersionControlTags_<TVoidSpec>::MESSAGE_SEQAN_UPDATE;
template <typename TVoidSpec>
constexpr char const * const VersionControlTags_<TVoidSpec>::MESSAGE_APP_UPDATE;
template <typename TVoidSpec>
constexpr char const * const VersionControlTags_<TVoidSpec>::MESSAGE_UNREGISTERED_APP;
template <typename TVoidSpec>
constexpr char const * const VersionControlTags_<TVoidSpec>::MESSAGE_REGISTERED_APP_UPDATE;

struct VersionCheck
{
    // ----------------------------------------------------------------------------
    // Member Variables
    // ----------------------------------------------------------------------------
    std::string _url;
    std::string _name;
    std::string _version = "0.0.0";
    std::string _program;
    std::string _command;
    std::string _path = _getPath();
    std::string _timestamp_filename;
    std::ostream & errorStream;

    // ----------------------------------------------------------------------------
    // Constructors
    // ----------------------------------------------------------------------------

    VersionCheck(std::string name,
                 std::string const & version,
                 std::ostream & errorStream) :
        _name{std::move(name)},
        errorStream(errorStream)
    {
        std::smatch versionMatch;
#if defined(NDEBUG) || defined(SEQAN_TEST_VERSION_CHECK_)
        _timestamp_filename = _path + "/" + _name + "_usr.timestamp";
#else
        _timestamp_filename = _path + "/" + _name + "_dev.timestamp";
#endif
        if (!version.empty() &&
            std::regex_search(version, versionMatch, std::regex("^([[:digit:]]+\\.[[:digit:]]+\\.[[:digit:]]+).*")))
        {
            _version = versionMatch.str(1); // in case the git revision number is given take only version number
        }
        _url = static_cast<std::string>("http://seqan-update.informatik.uni-tuebingen.de/check/SeqAn_") + _getOS() + _getBitSys() + _name + "_" + _version;
        _getProgram();
        _updateCommand();
    }

    // ----------------------------------------------------------------------------
    // Member Functions
    // ----------------------------------------------------------------------------

#if defined(STDLIB_VS)
    void _getProgram()
    {
        _program = "powershell.exe -NoLogo -NonInteractive -Command \"& {Invoke-WebRequest -erroraction 'silentlycontinue' -OutFile";
    }
#else  // Unix based platforms.
    void _getProgram()
    {
        // ask if system call for version or help is successfull
        if (!system("wget --version > /dev/null 2>&1"))
            _program = "wget -q -O";
        else if (!system("curl --version > /dev/null 2>&1"))
            _program =  "curl -o";
#ifndef __linux  // ftp call does not work on linux
        else if (!system("which ftp > /dev/null 2>&1"))
            _program =  "ftp -Vo";
#endif
        else
            _program.clear();
    }
#endif  // defined(STDLIB_VS)

    void _updateCommand()
    {
        if (!_program.empty())
        {
            _command = _program + " " + _path + "/" + _name + ".version " + _url;
#if defined(STDLIB_VS)
            _command = _command + "; exit  [int] -not $?}\" > nul 2>&1";
#else
            _command = _command + " > /dev/null 2>&1";
#endif
        }
    }

    inline void operator()(std::promise<bool> versionCheckProm)
    {
        _checkForNewerVersion(*this, std::move(versionCheckProm));
    }
};

// ==========================================================================
// Metafunctions
// ==========================================================================

// ==========================================================================
// Functions
// ==========================================================================

// ----------------------------------------------------------------------------
// Function setURL()
// ----------------------------------------------------------------------------

inline void setURL(VersionCheck & me, std::string url)
{
    std::swap(me._url, url);
    me._updateCommand();
}

// ----------------------------------------------------------------------------
// Function _getOS()
// ----------------------------------------------------------------------------

constexpr const char * _getOS()
{
    //get system information
#ifdef __linux
    return "Linux";
#elif __APPLE__
    return "MacOS";
#elif defined(STDLIB_VS)
    return "Windows";
#elif __FreeBSD__
    return "FreeBSD";
#elif __OpenBSD__
    return "OpenBSD";
#else
    return "unknown";
#endif
}

// ----------------------------------------------------------------------------
// Function _getBitSys()
// ----------------------------------------------------------------------------

constexpr const char * _getBitSys()
{
#if SEQAN_IS_32_BIT
    return "_32_";
#else
     return "_64_";
#endif
}

// ----------------------------------------------------------------------------
// Function _checkWritability()
// ----------------------------------------------------------------------------

#if defined(STDLIB_VS)
inline bool _checkWritability(std::string const & path)
{
    DWORD ftyp = GetFileAttributesA(path.c_str());
    if (ftyp == INVALID_FILE_ATTRIBUTES || !(ftyp & FILE_ATTRIBUTE_DIRECTORY))
    {
        if (!CreateDirectory(path.c_str(), NULL))
            return false;
    }

    HANDLE dummyFile; // check writablity by trying to create a file in GENERIC_WRITE mode
    std::string fileName(path + "/dummy.txt");
    dummyFile = CreateFile(fileName.c_str(), GENERIC_WRITE, 0, NULL, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);

    if (dummyFile == INVALID_HANDLE_VALUE)
        return false;

    CloseHandle(dummyFile);
    bool successful_deletion = DeleteFile(fileName.c_str());
    SEQAN_ASSERT(successful_deletion);
    if (!successful_deletion)
        return false;

    return true;
}
#else
inline bool _checkWritability(std::string const & path)
{
    struct stat d_stat;
    if (stat(path.c_str(), &d_stat) < 0)
    {
        // try to make dir
        std::string makeDir("mkdir -p " + path);
        if (system(makeDir.c_str()))
            return false; // could not create home dir

        if (stat(path.c_str(), &d_stat) < 0) // repeat stat
            return false;
    }

    if (!(d_stat.st_mode & S_IWUSR))
        return false; // dir not writable

    return true;
}
#endif

// ----------------------------------------------------------------------------
// Function _getPath()
// ----------------------------------------------------------------------------

inline std::string _getPath()
{
    std::string path;
#if defined(STDLIB_VS)
    path = std::string(getenv("UserProfile")) + "/.config/seqan";
#else
    path = std::string(getenv("HOME")) + "/.config/seqan";
#endif

    // check if user has permission to write to home path
    if (!_checkWritability(path))
    {
#if defined(STDLIB_VS)
        TCHAR tmp_path [MAX_PATH];
        if (GetTempPath(MAX_PATH, tmp_path) != 0)
        {
            path = tmp_path;
        }
        else
        { //GetTempPath() returns 0 on failure
            path.clear();
        }
# else // unix
        path = "/tmp";
#endif
    }
    return path;
}

// ----------------------------------------------------------------------------
// Function _getFileTimeDiff()
// ----------------------------------------------------------------------------

inline double _getFileTimeDiff(VersionCheck const & me)
{
    double curr = static_cast<double>(std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    std::ifstream timestamp_file;
    timestamp_file.open(me._timestamp_filename.c_str());

    if (timestamp_file.is_open())
    {
        std::string str_time;
        std::getline(timestamp_file, str_time);
        timestamp_file.close();
        double d_time;
        lexicalCast(d_time, str_time);
        return curr - d_time;
    }
    return curr;
}

// ----------------------------------------------------------------------------
// Function _getNumbersFromString()
// ----------------------------------------------------------------------------

inline String<int> _getNumbersFromString(std::string const & str)
{
    String<int> numbers;
    StringSet<std::string> set;
    strSplit(set, str, EqualsChar<'.'>(), false);
    for (auto & num : set)
    {
        appendValue(numbers, lexicalCast<int>(num));
    }
    return numbers;
}

// ----------------------------------------------------------------------------
// Function _readVersionString()
// ----------------------------------------------------------------------------

inline void _readVersionStrings(std::vector<std::string> & versions, std::string const & version_file)
{
    std::ifstream myfile;
    myfile.open(version_file.c_str());
    std::string app_version;
    std::string seqan_version;
    if (myfile.is_open())
    {
        std::getline(myfile, app_version); // get first line which should only contain the version number of the app

#if !defined(NDEBUG) || defined(SEQAN_TEST_VERSION_CHECK_)
        if (app_version == VersionControlTags_<>::UNREGISTERED_APP)
            versions[0] = app_version;
#endif // !defined(NDEBUG) || defined(SEQAN_TEST_VERSION_CHECK_)

        if (std::regex_match(app_version, std::regex("^[[:digit:]]+\\.[[:digit:]]+\\.[[:digit:]]+$")))
            versions[0] = app_version;

        std::getline(myfile, seqan_version); // get second line which should only contain the version number of seqan

        if (std::regex_match(seqan_version, std::regex("^[[:digit:]]+\\.[[:digit:]]+\\.[[:digit:]]+$")))
            versions[1] = seqan_version;

        myfile.close();
    }
}

// ----------------------------------------------------------------------------
// Function _callServer()
// ----------------------------------------------------------------------------

inline void _callServer(VersionCheck const me, std::promise<bool> prom)
{
    // update timestamp
    std::ofstream timestamp_file(me._timestamp_filename.c_str());
    if (timestamp_file.is_open())
    {
        timestamp_file << std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now().time_since_epoch()).count();
        timestamp_file.close();
    }

    // system call
    // http response is stored in a file '.config/seqan/{app_name}_version'
    if (system(me._command.c_str()))
        prom.set_value(false);
    else
        prom.set_value(true);
}

// ----------------------------------------------------------------------------
// Function checkForNewerVersion()
// ----------------------------------------------------------------------------

inline void _checkForNewerVersion(VersionCheck & me, std::promise<bool> prom)
{
    if (me._path.empty()) // neither home dir nor temp dir are writable
    {
        prom.set_value(false);
        return;
    }

    std::string version_filename(me._path + "/" + me._name + ".version");
    double min_time_diff(86400);                                 // one day = 86400 seonds
    double file_time_diff(_getFileTimeDiff(me)); // time difference in seconds

    if (file_time_diff < min_time_diff)
    {
        prom.set_value(false); // only check for newer version once a day
        return;
    }

    std::vector<std::string> str_server_versions{"", ""};
    _readVersionStrings(str_server_versions, version_filename);

#if !defined(NDEBUG) || defined(SEQAN_TEST_VERSION_CHECK_) // only check seqan version in debug or testing mode
    if (!str_server_versions[1].empty()) // seqan version
    {
        std::string seqan_version = std::to_string(SEQAN_VERSION_MAJOR) + "." +
                                    std::to_string(SEQAN_VERSION_MINOR) + "." +
                                    std::to_string(SEQAN_VERSION_PATCH);
        Lexical<> version_comp(_getNumbersFromString(seqan_version), _getNumbersFromString(str_server_versions[1]));

        if (isLess(version_comp))
            me.errorStream << VersionControlTags_<>::MESSAGE_SEQAN_UPDATE;
    }
    if (str_server_versions[0] == VersionControlTags_<>::UNREGISTERED_APP)
        me.errorStream << VersionControlTags_<>::MESSAGE_UNREGISTERED_APP;
#endif

    if (!str_server_versions[0].empty() & !(str_server_versions[0] == VersionControlTags_<>::UNREGISTERED_APP)) // app version
    {
        Lexical<> version_comp(_getNumbersFromString(me._version), _getNumbersFromString(str_server_versions[0]));

#if defined(NDEBUG) || defined(SEQAN_TEST_VERSION_CHECK_) // only check app version in release or testing mode
        if (isLess(version_comp))
            me.errorStream << VersionControlTags_<>::MESSAGE_APP_UPDATE;
#endif // defined(NDEBUG) || defined(SEQAN_TEST_VERSION_CHECK_)

#if !defined(NDEBUG) || defined(SEQAN_TEST_VERSION_CHECK_) // only notify developer that app version should be updated on server
        if (isGreater(version_comp))
            me.errorStream << VersionControlTags_<>::MESSAGE_REGISTERED_APP_UPDATE;
#endif // !defined(NDEBUG) || defined(SEQAN_TEST_VERSION_CHECK_)
    }

    if (me._program.empty())
    {
        prom.set_value(false);
        return;
    }

    // launch a seperate thread to not defer runtime.
    std::thread(_callServer, me, std::move(prom)).detach();
}

} // namespace seqan
#endif //SEQAN_INCLUDE_ARG_PARSE_VERSION_CHECK_H_
