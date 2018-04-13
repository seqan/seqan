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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Simple helper class for creating timelines.
// ==========================================================================

#ifndef APP_RAZERS_PROFILE_TIMELINE_H_
#define APP_RAZERS_PROFILE_TIMELINE_H_

#include <iostream>

#ifdef STDLIB_VS
#include <process.h>
#else  // #ifdef STDLIB_VS
#include <unistd.h>
#endif  // #ifdef STDLIB_VS

#ifdef _OPENMP
#include <omp.h>
#endif

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Typedes, Classes, Enums
// ============================================================================

struct TimelineEntry_
{
    unsigned entryType;
    bool isBegin;
    double timestamp;

    TimelineEntry_() :
        entryType(0), isBegin(false), timestamp(-1) {}

    TimelineEntry_(bool isBegin_, unsigned entryType_, double timestamp_) :
        entryType(entryType_), isBegin(isBegin_), timestamp(timestamp_) {}
};

class Timeline
{
public:
    enum
    {
        INITIAL_SIZE = 1024
    };

    double initTimestamp;
    String<Pair<CharString, CharString> > _taskTypeNames;  // (short name, long name)
    String<String<TimelineEntry_> > _entries;

    // TODO(holtgrew): Move to global function?
    static Timeline & instance()
    {
        static Timeline tl;
        return tl;
    }

private:
    // Can only construct in instance()!
    Timeline() :
        initTimestamp(sysTime()) {}

    // No copy-construction, no assignment, we are a Singleton!
    Timeline(Timeline const &);
    void operator=(Timeline const &);
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

inline
unsigned
timelineAddTaskType(CharString const & shortName, CharString const & longName)
{
    appendValue(Timeline::instance()._taskTypeNames, Pair<CharString, CharString>(shortName, longName));
    return length(Timeline::instance()._taskTypeNames) - 1;
}

inline
unsigned
timelineAddTaskType(CharString const & shortName)
{
    return timelineAddTaskType(shortName, shortName);
}

inline
void
initTimeline(int threadCount)
{
    resize(Timeline::instance()._entries, threadCount);
    for (int i = 0; i < threadCount; ++i)
        reserve(Timeline::instance()._entries[i], static_cast<int>(Timeline::INITIAL_SIZE), Generous());
    timelineAddTaskType("WAIT");
}

inline
void
initTimeline()
{
    initTimeline(omp_get_max_threads());
}

inline
double
timelineBeginTask(unsigned taskTypeNo)
{
    SEQAN_ASSERT_LT_MSG(taskTypeNo, length(Timeline::instance()._taskTypeNames), "Too large task type no!");
    double timestamp = sysTime();
    appendValue(Timeline::instance()._entries[omp_get_thread_num()], TimelineEntry_(true, taskTypeNo, timestamp));
    return timestamp;
}

inline
double
timelineEndTask(unsigned taskTypeNo)
{
    SEQAN_ASSERT_LT_MSG(taskTypeNo, length(Timeline::instance()._taskTypeNames), "Too large task type no!");
    double timestamp = sysTime();
    appendValue(Timeline::instance()._entries[omp_get_thread_num()], TimelineEntry_(false, taskTypeNo, timestamp));
    return timestamp;
}

inline
void
dumpTimeline(char const * path, bool appendPid)
{
    const double gapIgnore = 0.0001;
    char * pathBuffer = new char[strlen(path) + 30];
    strcpy(pathBuffer, path);
    if (appendPid)
    {
#ifdef STDLIB_VS
        int pid = _getpid();
#else // #ifdef STDLIB_VS
        int pid = getpid();
#endif // #ifdef STDLIB_VS
        char buffer[30];
        snprintf(buffer, 30, "%d", pid);
        strcat(pathBuffer, ".");
        strcat(pathBuffer, buffer);
    }
    Timeline const & timeline = Timeline::instance();
    FILE * fp = fopen(pathBuffer, "w");
    fprintf(fp, "@SQN:PROFILE\n");
    fprintf(fp, "@TIME\t%f\t%f\n", timeline.initTimestamp, sysTime());

    // Dump event types;
    for (unsigned i = 0; i < length(timeline._taskTypeNames); ++i)
        fprintf(fp, "@EVENT\t%u\t%s\t%s\n", i, toCString(timeline._taskTypeNames[i].i1), toCString(timeline._taskTypeNames[i].i2));

    // Dump events.
    const char * arr[] = {"END", "BEGIN"};
    for (unsigned threadId = 0; threadId < length(timeline._entries); ++threadId)
    {
        for (unsigned i = 0; i < length(timeline._entries[threadId]); ++i)
        {
            while (true)
            {
                if (i + 1 < length(timeline._entries[threadId]) &&
                    !timeline._entries[threadId][i].isBegin &&
                    timeline._entries[threadId][i + 1].isBegin &&
                    timeline._entries[threadId][i].entryType == timeline._entries[threadId][i + 1].entryType &&
                    fabs(timeline._entries[threadId][i + 1].timestamp - timeline._entries[threadId][i].timestamp) < gapIgnore)
                {
                    i += 2;
                }
                else
                {
                    break;
                }
            }
            fprintf(fp, "%u\t%s\t%d\t%f\n", threadId, arr[timeline._entries[threadId][i].isBegin], timeline._entries[threadId][i].entryType, timeline._entries[threadId][i].timestamp);
        }
    }

    fclose(fp);
    delete[] pathBuffer;
}

inline
void
dumpTimeline(char const * path)
{
    dumpTimeline(path, false);
}

}  // namespace seqan

#endif  // #ifndef APP_RAZERS_PROFILE_TIMELINE_H_
