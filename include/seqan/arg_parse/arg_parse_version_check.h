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

// ==========================================================================
// Tags, Classes, Enums
// ==========================================================================

struct VersionCheck
{
    // ----------------------------------------------------------------------------
    // Member Variables
    // ----------------------------------------------------------------------------

    std::string _url = "http://openms-update.informatik.uni-tuebingen.de/check/OpenMS_KNIME_";
    std::string _path = std::string(getenv("HOME")) + "/.config/seqan";
    std::string _name;
    std::string _version = "0.0.0";
    std::string _program;
    std::string _command;
    std::string _website = "www.seqan.de/applications";
    std::future<bool> & _fut;

    //get system information
#ifdef __linux
    std::string _os = "Linux";
#elif __APPLE__
    _std::string os = "MacOS"
#elif __WINDOWS__
    _std::string os = "Windows"
#elif __FreeBSD__
    _std::string os = "FreeBSD";
#elif __OpenBSD__
    std::string _os = "OpenBSD";
#else
    std::string _os = "unknown";
#endif

    // ----------------------------------------------------------------------------
    // Constructors
    // ----------------------------------------------------------------------------

    VersionCheck(std::future<bool> & fut,
                 std::string const & name,
                 std::string const & version,
                 std::string const & website):
    _fut(fut)
    {
        //_fut = &fut;
        _name = name;
        if (!version.empty())
            _version = version;
        if(!website.empty())
            _website = website;
        _getProgram();
        _updateCommand();
    }

    VersionCheck(VersionCheck const & rhs):
    _fut(rhs._fut)
    {
        _name    = rhs._name;
        _version = rhs._version;
        _website = rhs._website;
        _program = rhs._program;
        _command = rhs._command;
    }

    VersionCheck(VersionCheck && rhs) :
        _fut(rhs._fut)
    {
        swap(rhs);
    }

    VersionCheck & operator=(VersionCheck rhs)
    {
        swap(rhs);
        return *this;
    }

    // ----------------------------------------------------------------------------
    // Member Functions
    // ----------------------------------------------------------------------------

    inline void swap(VersionCheck & rhs)
    {
        std::swap(_name,    rhs._name);
        std::swap(_version, rhs._version);
        std::swap(_website, rhs._website);
        std::swap(_program, rhs._program);
        std::swap(_command, rhs._command);
        std::swap(_fut,     rhs._fut);
    }

#ifdef __unix
    void _getProgram()
    {
        // ask if system call for version or help is successfull
        if (!system("wget --version > /dev/null 2>&1"))
            _program = "wget -q -O";
        else if (!system("curl --version > /dev/null 2>&1"))
            _program =  "curl -o";
#ifndef __linux
        // ftp call does not work on linux
        else if (!system("which ftp > /dev/null 2>&1"))
            _program =  "ftp -Vo";
#endif
        else
            _program.clear();
    }
#else // windows
    void _getProgramm()
    {
        _program.clear();
    }
#endif

    void _updateCommand()
    {
        if (!_program.empty())
            _command = _program + " " + _path + "/" + _name + "_version " + _url + _os + "_64_"+ _name + "_" + _version;
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
// Function _checkWritability()
// ----------------------------------------------------------------------------
inline bool _checkWritability(std::string path)
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

// ----------------------------------------------------------------------------
// Function _getFileTimeDiff()
// ----------------------------------------------------------------------------
inline double _getFileTimeDiff(std::string version_file)
{
    time_t curr;
    time(&curr); // get current time

    struct stat t_stat;
    if (stat(version_file.c_str(), &t_stat) < 0)  // file does not exist
        return -1;
    else
        return (difftime(curr, t_stat.st_mtime)); // returns time in seconds
}

// ----------------------------------------------------------------------------
// Function _getNumbersFromString()
// ----------------------------------------------------------------------------
inline void _getNumbersFromString(seqan::String<int> & numbers, std::string const & str)
{
    std::string number;
    std::istringstream iss(str);
    while (std::getline(iss, number, '.'))
    {
        if (!number.empty())
        {
            seqan::appendValue(numbers, atoi(seqan::toCString(number)));
        }
    }
}

// ----------------------------------------------------------------------------
// Function _getVersionNumbers()
// ----------------------------------------------------------------------------
inline bool _readVersionNumbers(seqan::String<int> & version_numbers, std::string const & version_file)
{
    std::ifstream myfile;
    myfile.open(version_file.c_str());
    std::string line;
    if (myfile.is_open())
    {
        std::getline(myfile,line); // get first line which should only contain the version number
        if (std::regex_match(line, std::regex("^[[:digit:]]\\.[[:digit:]]\\.[[:digit:]]$")))
        {
            _getNumbersFromString(version_numbers, line);
            myfile.close();
        }
        else if (line == "UNREGISTERED APP")
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
inline bool _isSmaller(seqan::String<int> & left, seqan::String<int> & right)
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
inline bool _callServer(VersionCheck me)
{
    // system call
    // http response is stored in a file '.config/seqan/{app_name}_version'
    //std::cout << "I will perform the following command:\n" << me._command << std::endl;
    if (system(seqan::toCString(me._command)))
    {
        std::ofstream version_file (seqan::toCString(me._path + "/" + me._name + "_version"));
        if (version_file.is_open())
        {
            version_file << "UNABLE TO CALL HOME.\n";
            version_file.close();
        }
        return false;
    }
    return true;
}

// ----------------------------------------------------------------------------
// Function checkForNewerVersion()
// ----------------------------------------------------------------------------
inline bool checkForNewerVersion(VersionCheck & me)
{
    std::string version_file = me._path + "/" + me._name + "_version";
    double min_time_diff = 86400;                           // one day = 86400 seonds
    double file_time_diff = _getFileTimeDiff(version_file); // check for time the version file was last updated

    if (file_time_diff > -1) // file exists. TODO:: ask if there should be a time limit here too
    {
        seqan::String<int> new_ver;
        if (_readVersionNumbers(new_ver, version_file))
        {
            seqan::String<int> old_ver;
            _getNumbersFromString(old_ver, me._version);

            if (_isSmaller(old_ver, new_ver))
            {
                if(me._name == "seqan")
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
            else if (_isSmaller(new_ver, old_ver)) // can only happen for registered app that is developed further
            {
                std::cerr << "[APP INFO] :: Thank you for registering " << me._name << ".\n" 
                          << "[APP INFO] :: We noticed you developed a newer version of your app (" << me._version << ")\n"
                          << "[APP INFO] :: Please send us an email to update your version info (seqan@team.fu-berlin.de)" << "\n\n";
            }
        }
    }

    if (me._program.empty())
        return false;

    if (file_time_diff < min_time_diff && file_time_diff > -1)
        return false;

    if (!_checkWritability(me._path))
        me._path = "/tmp";

    // launch a seperate thread to not defer runtime
    me._fut = std::async(std::launch::async, _callServer, me);
    return true;
}

#endif //SEQAN_INCLUDE_ARG_PARSE_VERSION_CHECK_H_
