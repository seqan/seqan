// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2016, Knut Reinert, FU Berlin
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

#include <iostream>
#include <fstream>
#include <string>
#include <sys/stat.h>
#include <time.h>
#include <stdio.h>
#include <regex>
#include <future>
#include <sys/types.h> 
#include <utime.h>
#include <chrono>

namespace seqan
{

// ==========================================================================
// Forwards
// ==========================================================================

// NOTE(rrahn): In-file forward for function call operator.
struct VersionCheck;
inline bool _checkForNewerVersion(VersionCheck &);
inline std::string _getOS();
inline std::string _getPath();
inline std::string _getBitSys();

// ==========================================================================
// Tags, Classes, Enums
// ==========================================================================

struct VersionControlTags
{
    constexpr static const char * SEQAN_NAME = "seqan";
    constexpr static const char * UNREGISTERED_APP = "UNREGISTERED APP";
};

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
    std::string _website = "https://github.com/seqan/seqan/tree/master/apps";
    std::string _path = _getPath();

    // ----------------------------------------------------------------------------
    // Constructors
    // ----------------------------------------------------------------------------

    VersionCheck(std::string const & name,
                 std::string const & version,
                 std::string const & website)
    {
        _name = name;
        if (!version.empty() &&
            std::regex_match(version, std::regex("^[[:digit:]]+\\.[[:digit:]]+\\.[[:digit:]]+.*")))
            _version = version.substr(0,5); // in case the git revision number is given take only version number
        if(!website.empty())
            _website = website;
        _url = "http://www.seqan.de/version_check/SeqAn_" + _getOS() + _getBitSys() + _name + "_" + _version;
        _getProgram();
        _updateCommand();
    }

    // ----------------------------------------------------------------------------
    // Member Functions
    // ----------------------------------------------------------------------------

#if defined(PLATFORM_WINDOWS)
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
#endif  // defined(PLATFORM_WINDOWS)

    void _updateCommand()
    {
        if (!_program.empty())
        {
            _command = _program + " " + _path + "/" + _name + ".version " + _url;
#if defined(PLATFORM_WINDOWS)
            _command = _command + "; exit  [int] -not $?}\" > nul 2>&1";
#endif
        }
    }

    inline void operator()(std::promise<bool> versionCheckProm)
    {
        versionCheckProm.set_value(_checkForNewerVersion(*this));
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
inline std::string _getOS()
{
    //get system information
    std::string os;
#ifdef __linux
    os = "Linux";
#elif __APPLE__
    os = "MacOS";
#elif defined(PLATFORM_WINDOWS)
    os = "Windows";
#elif __FreeBSD__
    os = "FreeBSD";
#elif __OpenBSD__
    os = "OpenBSD";
#else
    os = "unknown";
#endif
    return os;
}

// ----------------------------------------------------------------------------
// Function _getPath()
// ----------------------------------------------------------------------------
inline std::string _getPath()
{
    std::string path;
#if defined(PLATFORM_WINDOWS)
    path = std::string(getenv("UserProfile")) + "/.config/seqan";
#else
    path = std::string(getenv("HOME")) + "/.config/seqan";
#endif
    return path;
}

// ----------------------------------------------------------------------------
// Function _getPath()
// ----------------------------------------------------------------------------
inline std::string _getBitSys()
{
    std::string bitSys;

#if SEQAN_IS_32_BIT
    bitSys = "_32_";
#else
    bitSys = "_64_";
#endif
    return bitSys;
}

// ----------------------------------------------------------------------------
// Function _checkWritability()
// ----------------------------------------------------------------------------
#if defined(PLATFORM_WINDOWS)
inline bool _checkWritability(std::string const & path)
{
    DWORD ftyp = GetFileAttributesA(path.c_str());
    if (ftyp == INVALID_FILE_ATTRIBUTES || !(ftyp & FILE_ATTRIBUTE_DIRECTORY))
    {
        if (!CreateDirectory(path.c_str(), NULL))
            return false;
    }

    HANDLE dummyFile; // check writablity by trying to create a file in GENERIC_WRITE mode
    std::string fileName = path + "/dummy.txt";
    dummyFile = CreateFile(fileName.c_str(), GENERIC_WRITE, 0, NULL, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);

    if (dummyFile == INVALID_HANDLE_VALUE)
        return false;
    
    CloseHandle(dummyFile);
    DeleteFile(fileName.c_str());
    return true;
}
#else
inline bool _checkWritability(std::string const & path)
{
    struct stat d_stat;
    if (stat(path.c_str(), &d_stat) < 0)
    {
        // try to make dir
        std::string makeDir = "mkdir -p " + path;
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
// Function _getFileTimeDiff()
// ----------------------------------------------------------------------------
inline  double _getFileTimeDiff(std::string const & timestamp_filename)
{
    double curr = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    std::ifstream timestamp_file;
    timestamp_file.open(timestamp_filename.c_str());

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
inline void _getNumbersFromString(String<int> & numbers, std::string const & str)
{
    std::string number;
    std::istringstream iss(str);
    while (std::getline(iss, number, '.'))
    {
        if (!number.empty())
        {
            appendValue(numbers, atoi(toCString(number)));
        }
    }
}

// ----------------------------------------------------------------------------
// Function _getVersionNumbers()
// ----------------------------------------------------------------------------
inline bool _readVersionNumbers(String<int> & version_numbers, std::string const & version_file)
{
    std::ifstream myfile;
    myfile.open(version_file.c_str());
    std::string line;
    if (myfile.is_open())
    {
        std::getline(myfile,line); // get first line which should only contain the version number
        if (std::regex_match(line, std::regex("^[[:digit:]]+\\.[[:digit:]]+\\.[[:digit:]]+$")))
        {
            _getNumbersFromString(version_numbers, line);
            myfile.close();
        }
        else if (line == VersionControlTags::UNREGISTERED_APP)
        {
            std::cerr << "[SEQAN INFO] :: Thank you for using SeqAn!\n"
                      << "[SEQAN INFO] :: You might want to regsiter you app for support and version check features?!\n"
                      << "[SEQAN INFO] :: Just send us an email to seqan@team.fu-berlin.de with your app name and version number.\n"
                      << "[SEQAN INFO] :: If you don't want to recieve this message anymore set --version_check OFF\n\n";
            myfile.close();
            return false;
        }
        else
        {
            myfile.close();
            return false;
        }
    }
    else
    {
        return false;
    }
    return true;
}

// ----------------------------------------------------------------------------
// Function _isSmaller()
// ----------------------------------------------------------------------------
inline bool _isSmaller(String<int> & left, String<int> & right)
{
    for (unsigned i = 0; i < length(left); ++i)
    {
        if (left[i] < right[i])
        {
            return true;
        }
    }
    return false;
}

// ----------------------------------------------------------------------------
// Function _callServer()
// ----------------------------------------------------------------------------

inline bool _callServer(VersionCheck const & me)
{
    // update timestamp
    std::string timestamp_filename = me._path + "/" + me._name + ".timestamp";
    std::ofstream timestamp_file(timestamp_filename.c_str());
    if (timestamp_file.is_open())
    {
        timestamp_file << std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now().time_since_epoch()).count();
        timestamp_file.close();
    }

    // system call
    // http response is stored in a file '.config/seqan/{app_name}_version'
    if (system(me._command.c_str()))
        return false;

    return true;
}

// ----------------------------------------------------------------------------
// Function checkForNewerVersion()
// ----------------------------------------------------------------------------

inline bool _checkForNewerVersion(VersionCheck & me)
{
    if (!_checkWritability(me._path))
    {
#if defined(PLATFORM_WINDOWS)
        TCHAR tmp_path [MAX_PATH];
        if (GetTempPath (MAX_PATH, tmp_path) != 0)
            me._path = tmp_path;
        else //GetTempPath() returns 0 on failure
            return false;
# else // unix
        me._path = "/tmp";
#endif
        me._updateCommand();
    }

    std::string version_filename   = me._path + "/" + me._name + ".version";
    std::string timestamp_filename = me._path + "/" + me._name + ".timestamp";
    double min_time_diff = 86400;                                 // one day = 86400 seonds
    double file_time_diff = _getFileTimeDiff(timestamp_filename); // time difference in seconds

    if (file_time_diff < min_time_diff)
        return false; // only check for newer version once a day

    String<int> new_ver;
    if (_readVersionNumbers(new_ver, version_filename))
    {
        String<int> old_ver;
        _getNumbersFromString(old_ver, me._version);

        if (_isSmaller(old_ver, new_ver))
        {
            if(me._name == VersionControlTags::SEQAN_NAME)
            {
                std::cerr << "[SEQAN INFO] :: There is a newer SeqAn version available : SeqAn "
                          << new_ver[0] << "." << new_ver[1] << "." << new_ver[2] << " Go to "
                          << me._website << "\n"
                          << "[SEQAN INFO] :: If you don't want to recieve this message again set --version-check APP_ONLY" << "\n\n";
            }
            else
            {
                std::cerr << "[APP INFO] :: There is a newer version available: " << me._name << " " 
                          << new_ver[0] << "." << new_ver[1] << "." << new_ver[2] << "\n"
                          << "[APP INFO] :: Check out " << me._website << "\n"
                          << "[APP INFO] :: If you don't want to recieve this message again set --version_check OFF" << "\n\n";
            }
        }
        else if (_isSmaller(new_ver, old_ver))
        {
            std::cerr << "[APP INFO] :: We noticed your app version (" << me._version << ") is newer than the one registered.\n"
                      << "[APP INFO] :: If you are the developer of this app, please send us an email to update your version info (support@seqan.de)" 
                      << "[APP INFO] :: If not, you might want to contact the developer."
                      << "\n\n";
        }
    }

    if (me._program.empty())
        return false;

    // launch a seperate thread to not defer runtime
    // std::cout << me._command << std::endl;
    return _callServer(me);
}

} // namespace seqan
#endif //SEQAN_INCLUDE_ARG_PARSE_VERSION_CHECK_H_
