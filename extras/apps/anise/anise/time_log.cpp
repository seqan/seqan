// ==========================================================================
//                                   ANISE
// ==========================================================================
// Copyright (c) 2006-2014, Knut Reinert, FU Berlin
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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

#include "time_log.h"

#include <fstream>
#include <chrono>
#include <mutex>
#include <ctime>
#include <iomanip>

#include "anise/temporary_file_manager.h"
#include "anise/file_name_tokens.h"
#include "anise/time_log.h"

// --------------------------------------------------------------------------
// Class TimeLogImpl
// --------------------------------------------------------------------------

class TimeLogImpl
{
public:
    // Log for a given token, step no, and site id.  Optionally, specify a duration in seconds and a comment text.  The
    // current time stamp will be written to the file as well.
    void log(char const * token, int stepNo = -1, int siteID = -1, double seconds = -1,
             char const * text = "")
    {
        std::lock_guard<std::mutex> lock(mutex);

        std::time_t rawtime;
        std::time(&rawtime);
        std::tm *local_time = std::localtime(&rawtime);
        char buffer[100];
        snprintf(buffer, 100, "%4d/%02d/%02d %02d:%02d:%02d", 1900 + local_time->tm_year,
                 local_time->tm_mon + 1, local_time->tm_mday, local_time->tm_hour,
                 local_time->tm_min, local_time->tm_sec);
        f << buffer << "\t" << token << "\t";
        if (stepNo >= 0)
            f << stepNo << "\t";
        else
            f << "<no step>\t";
        if (siteID >= 0)
            f << siteID << "\t";
        else
            f << "<no site>\t";
        if (seconds != -1)
            f << seconds << "\t";
        else
            f << "<no time>\t";
        if (*text)
            f << "text\t";
        else
            f << ".";
        f << "\n" << std::flush;
    }

    void open(TemporaryFileManager & tmpMgr)
    {
        tmpMgr.open(f, std::ios::binary | std::ios::out | std::ios::app,
                    TIME_LOG_TOKEN, TIME_LOG_EXT, -1, -1);
        f << "#TIME STAMP\tTOKEN\tSTEP\tSITE\tSECONDS\tTEXT\n";
    }

private:

    // Mutex to use for serializing access.
    std::mutex mutex;
    // The file to append to.
    std::fstream f;
};

// --------------------------------------------------------------------------
// Class TimeLog
// --------------------------------------------------------------------------

TimeLog::TimeLog() : impl(new TimeLogImpl)
{}

TimeLog::~TimeLog()
{}

void TimeLog::open(TemporaryFileManager & tmpMgr)
{
    impl->open(tmpMgr);
}

void TimeLog::log(char const * token, int stepNo, int siteID, double seconds,
                  char const * text)
{
    impl->log(token, stepNo, siteID, seconds, text);
}

// ----------------------------------------------------------------------------
// Function withTimeLog()
// ----------------------------------------------------------------------------

void withTimeLog(char const * tokenPrefix, int stepNo, int siteID, std::function<void()> fun)
{
    using namespace std::chrono;

    // Log START event.
    std::string tokenBegin = tokenPrefix;
    tokenBegin.append(" START");
    system_clock::time_point startTime = system_clock::now();
    TimeLog::instance().log(tokenBegin.c_str(), stepNo, siteID);

    // Call function.
    fun();

    // Log END event.
    system_clock::time_point endTime = system_clock::now();
    double secs = duration_cast<seconds>(endTime - startTime).count();
    std::string tokenEnd = tokenPrefix;
    tokenEnd.append(" END");
    TimeLog::instance().log(tokenEnd.c_str(), stepNo, siteID, secs);
}

void withTimeLog(char const * tokenPrefix, int stepNo, std::function<void()> fun)
{
    withTimeLog(tokenPrefix, stepNo, -1, fun);
}

void withTimeLog(char const * tokenPrefix, std::function<void()> fun)
{
    withTimeLog(tokenPrefix, -1, -1, fun);
}

