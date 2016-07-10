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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_HEADER_SYSTEM_EVENT_WIN_H
#define SEQAN_HEADER_SYSTEM_EVENT_WIN_H

namespace seqan {

// NOTE(rrahn): This calss implements an synchronization event for asynchronous file I/O on windows platforms.
// This file is required by file_async.
#ifdef STDLIB_VS

static SECURITY_ATTRIBUTES EventDefaultAttributes =
{
    sizeof(SECURITY_ATTRIBUTES),
    NULL,
    true
};

struct Event        // this class mustn't exceed the size of HANDLE (needed by waitForAll/Any)
{
    typedef HANDLE Handle;
    enum : unsigned {Infinite = INFINITE};
    Handle hEvent;

    Event() :
        hEvent(NULL) {}

    Event(BOOL initial)
    {
        SEQAN_DO_SYS2(open(initial), "Could not create Event");
    }

    // Move constructors
    Event(Event & other, Move) :
        hEvent(other.hEvent)
    {
        other.hEvent = NULL;
    }

    Event(Event && other) :
        hEvent(other.hEvent)
    {
        other.hEvent = NULL;
    }

    ~Event()
    {
        if (*this) SEQAN_DO_SYS2(close(), "Could not destroy Event");
    }

    inline Event & operator=(Event const & origin)
    {
        // resource sharing is not yet supported (performance reason)
        // it needs a reference counting technique
        if (origin)
        {
            hEvent = origin.hEvent;
            const_cast<Event &>(origin).hEvent = NULL;
        }
        else
            hEvent = NULL;
        return *this;
    }

    inline bool open(BOOL initial = FALSE)
    {
        return (hEvent = CreateEvent(&EventDefaultAttributes, TRUE, initial, NULL)) != NULL;
    }

    inline bool close()
    {
        bool success = CloseHandle(hEvent);
        hEvent = NULL;
        return success;
    }

    inline bool wait(DWORD timeoutMilliSec, bool & inProgress)
    {
        if (!hEvent)
        {
            inProgress = false;
            return true;
        }
        DWORD result = WaitForSingleObject(hEvent, timeoutMilliSec);
        inProgress = (result == WAIT_TIMEOUT);
        return result == WAIT_OBJECT_0 || inProgress;
    }

    inline bool wait()
    {
        bool dummy;
        return wait(Infinite, dummy);
    }

    inline bool signal()
    {
        return SetEvent(hEvent) != 0;
    }

    inline bool reset()
    {
        return ResetEvent(hEvent) != 0;
    }

    inline operator bool() const
    {
        return hEvent != NULL;
    }

private:

    Event(Event const &) :
        hEvent(NULL)
    {
        // we only support move construction (no copy-construction)
    }
};


//////////////////////////////////////////////////////////////////////////////
// global event functions

inline void reset(Event & e)
{
    e.reset();
}

inline bool waitForAll(Event eventList[], DWORD count, DWORD timeoutMilliSec)
{
    return WaitForMultipleObjects(count, &eventList[0].hEvent, true, timeoutMilliSec) != WAIT_TIMEOUT;
}

inline bool waitForAll(Event eventList[], DWORD count)
{
    return waitForAll(eventList, count, Event::Infinite);
}

inline int waitForAny(Event eventList[], DWORD count, DWORD timeoutMilliSec)
{
    DWORD result = WaitForMultipleObjects(count, &eventList[0].hEvent, false, timeoutMilliSec);

    if (/*result >= WAIT_OBJECT_0 && */ result < WAIT_OBJECT_0 + count)
        return result - WAIT_OBJECT_0;

    return -1;
}

inline int waitForAny(Event eventList[], DWORD count)
{
    return waitForAny(eventList, count, Event::Infinite);
}

template <>
struct HasMoveConstructor<Event> : True {};

//////////////////////////////////////////////////////////////////////////////
// global event functions

inline bool open(Event & e, bool initial)
{
    return e.open(initial);
}

inline bool open(Event & e)
{
    return open(e, false);
}

inline bool close(Event & e)
{
    return e.close();
}

inline bool waitFor(Event & e)
{
    return e.wait();
}

template <typename TTime>
inline bool waitFor(Event & e, TTime timeoutMilliSec, bool & inProgress)
{
        #ifdef disabledSEQAN_PROFILE
    double begin = sysTime();
    bool b = e.wait(timeoutMilliSec, inProgress);
    double end = sysTime();
    if (begin != end)
        std::cerr << "waitTime: " << end - begin << std::endl;
    return b;

        #else
    return e.wait(timeoutMilliSec, inProgress);

        #endif
}

inline bool signal(Event & e)
{
    return e.signal();
}

#endif  // #ifdef STDLIB_VS
}

#endif
