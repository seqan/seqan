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

#ifndef SEQAN_INCLUDE_ARG_PARSE_CALLING_HOME_H_
#define SEQAN_INCLUDE_ARG_PARSE_CALLING_HOME_H_

#include <iostream>
#include <fstream>
#include <string>
#include <sys/stat.h>
#include <time.h>
#include <stdio.h>

// ==========================================================================
// Tags, Classes, Enums
// ==========================================================================

struct CallingHome
{
    std::string _url = "http://openms-update.informatik.uni-tuebingen.de/check/OpenMS_KNIME_";
    std::string _homeDir = getenv("HOME");
    std::string _path = _homeDir + "/.config/seqan";
    std::string _appName;
    std::string _version;
    std::string _program;
    std::string _command;

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


#ifdef __unix
    void _getProgram()
    {
        // ask if system call for version or help is successfull
        if (!system("wget --version > /dev/null 2>&1"))
        {
            _program = "wget -q -O";
        }
        else if (!system("curl --version > /dev/null 2>&1"))
        {
            _program =  "curl -o";
        }
#ifndef __linux
        // ftp call does not work on linux
        else if (!system("which ftp > /dev/null 2>&1"))
        {
            _program =  "ftp -Vo";
        }
#endif
        else
        {
            _program.clear();
        }
    }
#else // windows
    void _getProgramm()
    {
        _program.clear();
    }
#endif

    //TODO::(smehringer) make void
    void _updateCommand()
    {
        if (!_program.empty())
            _command = _program + " " + _path + "/" + _appName + "_version " + _url + _os + "_64_"+ _appName + "_" + _version;
    }

    CallingHome(std::string const & name, std::string const & version)
    {
        _appName = name;
        if (version.empty())
            _version = "0.0.0";
        else
            _version = version;

        _getProgram();
        _updateCommand();
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
inline void setURL(CallingHome & me, std::string url)
{
    std::swap(me._url, url);
    me._updateCommand();
}

// ----------------------------------------------------------------------------
// Function _checkWritability()
// ----------------------------------------------------------------------------
inline bool _checkWritability(CallingHome const & me)
{
    struct stat d_stat;
    if (stat(me._path.c_str(), &d_stat) < 0)
    {
        // try to make dir
        std::string makeDir = "mkdir -p " + me._path;
        if (system(makeDir.c_str()))
            return false; // could not create home dir

        if (stat(me._path.c_str(), &d_stat) < 0) // repeat stat
            return false;
    }

    if (!(d_stat.st_mode & S_IWUSR))
        return false; // dir not writable

    return true;
}

// ----------------------------------------------------------------------------
// Function _checkDate()
// ----------------------------------------------------------------------------
inline bool _checkDate(CallingHome const & me)
{
    unsigned min_time_diff = 86400; // (seconds) minimum time difference before next call home

    time_t curr;
    time(&curr); // get current time

    struct stat t_stat; // get last modified time of version file (if existing)
    std::string version_file = me._path + "/" + me._appName + "_version";
    if (stat(version_file.c_str(), &t_stat) < 0)
        return true; // no file found so calling home shall be performed

    if (difftime(curr, t_stat.st_mtime) < min_time_diff)
    {
        std::cout << "Version file " + me._appName + "_version is up to date.";
        return false;
    }

    return true;
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
        if (line != "UNABLE TO CALL HOME.")
        {
            // TODO:: try catch just in case ?
            _getNumbersFromString(version_numbers, line);
            myfile.close();
        }
        else
        {
            myfile.close();
            //std::cout << "found UNABLE TO CALL HOME" << std::endl;
            return false;
        }
    }
    else
    {
        return false;
    }
    //std::cout << "length " << length(version_numbers) << std::endl;
    if (length(version_numbers) != 3)
        return false;

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
// Function checkForNewerVersion()
// ----------------------------------------------------------------------------
inline bool checkForNewerVersion(CallingHome const & me)
{
    struct stat t_stat;
    std::string version_file = me._path + "/" + me._appName + "_version";
    if (stat(version_file.c_str(), &t_stat) < 0) // check if file exists
        return false;

    seqan::String<int> new_ver;
    if (!_readVersionNumbers(new_ver, version_file))
        return false;

    seqan::String<int> old_ver;
    _getNumbersFromString(old_ver, me._version);

    if (_isSmaller(old_ver, new_ver))
        std::cerr << "There is a newer version of " << me._appName << " available!\n";
    // TODO:: supply webiste adress for newer version
    // default : www.seqan.de/applications

    return true;
}

// ----------------------------------------------------------------------------
// Function callHome()
// ----------------------------------------------------------------------------
inline bool callHome(CallingHome const & me)
{
    if (!_checkWritability(me))
        return false;

    if (!_checkDate(me))
        return false;

    if (me._program.empty())
        return false;

    if (me._command.empty())
        return false;

    // system call
    // http response is stored in a file '.{app_name}_version'
    std::cout << "I will perform the following command:\n" << me._command << std::endl;
    if (system(seqan::toCString(me._command)))
    {
        std::ofstream version_file (seqan::toCString(me._path + "/" + me._appName + "_version"));
        if (version_file.is_open())
        {
            version_file << "UNABLE TO CALL HOME.\n";
            version_file.close();
        }
        return false;
    }
    return true;
}

#endif //SEQAN_INCLUDE_ARG_PARSE_CALLING_HOME_H_
