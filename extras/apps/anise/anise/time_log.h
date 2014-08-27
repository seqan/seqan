// ==========================================================================
//                                   ANISE
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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

#ifndef SANDBOX_HERBARIUM_APPS_ANISE_ANISE_TIME_LOG_H_
#define SANDBOX_HERBARIUM_APPS_ANISE_ANISE_TIME_LOG_H_

#include <memory>
#include <functional>  // for std::function<> in withTimeLog()

// ============================================================================
// Forwards
// ============================================================================

class TimeLogImpl;
class TemporaryFileManager;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

class TimeLog
{
public:

    TimeLog();
    ~TimeLog();  // for pimpl

    // Open the time log.
    void open(TemporaryFileManager & tmpMgr);

    // Log for a given token, step no, and site id.  Optionally, specify a duration in seconds and a comment text.  The
    // current time stamp will be written to the file as well.
    void log(char const * token, int stepNo = -1, int siteID = -1, double seconds = -1,
             char const * text = "");

    // Return singleton instance.
    static TimeLog & instance()
    {
        static TimeLog i;
        return i;
    }

private:

    std::unique_ptr<TimeLogImpl> impl;
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function withTimeLog()
// ----------------------------------------------------------------------------

// Run the given callable with a wrapper.
//
// This wrapper will perform logging.

void withTimeLog(char const * tokenPrefix, std::function<void()> fun);
void withTimeLog(char const * tokenPrefix, int stepNo, std::function<void()> fun);
void withTimeLog(char const * tokenPrefix, int stepNo, int siteID,
                 std::function<void()> fun);

#endif  // #ifndef SANDBOX_HERBARIUM_APPS_ANISE_ANISE_TIME_LOG_H_
