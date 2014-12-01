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

#ifndef SANDBOX_HERBARIUM_APPS_ANISE_SHARED_SLIDING_WINDOW_H_
#define SANDBOX_HERBARIUM_APPS_ANISE_SHARED_SLIDING_WINDOW_H_

#include <algorithm>
#include <set>
#include <vector>

#include <seqan/basic.h>

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class SlidingWindowInterval
// ----------------------------------------------------------------------------

// Interval in SlidingWindowAlgo.

struct SlidingWindowInterval
{
    int rID;
    int beginPos;
    int endPos;
    unsigned itemID;

    explicit SlidingWindowInterval(int rID = -1, int beginPos = -1, int endPos = -1, unsigned itemID = -1) :
            rID(rID), beginPos(beginPos), endPos(endPos), itemID(itemID)
    {}

    bool operator<(SlidingWindowInterval const & other) const
    {
        return mkTuple() < other.mkTuple();
    }

private:
    std::tuple<int, int, int, unsigned> mkTuple() const
    {
        return std::make_tuple(rID, beginPos, endPos, itemID);
    }
};

// ----------------------------------------------------------------------------
// Class SlidingWindowAlgo
// ----------------------------------------------------------------------------

// TODO(holtgrew): Line-sweep is the better name.

class SlidingWindowAlgo
{
public:
    // (rID, pos, isOpen)
    typedef std::tuple<int, int, bool> TLocation;
    typedef std::set<unsigned>::const_iterator TIterator;

    SlidingWindowAlgo() : currentLocation(-1, -1, false), it(events.begin())
    {}

    SlidingWindowAlgo(std::vector<SlidingWindowInterval> const & intervals) :
            events(makeEvents(intervals)), currentLocation(-1, -1, false), it(events.begin())
    {
        if (!events.empty())
            currentLocation = std::make_tuple(events.front().rID, events.front().pos, false);
    }

    // Return current location.
    TLocation location() const
    {
        return currentLocation;
    }

    // Returns true if at end.
    bool atEnd() const
    {
        return (it == events.end());
    }

    // Advance to next location and return it.
    TLocation advance()
    {
        if (atEnd())
            return currentLocation;

        TLocation loc = currentLocation;
        if (!std::get<2>(loc))
        {
            // Is after the end events and before the begin events.  Advance only in the isBegin flag.
            std::get<2>(loc) = true;
            for (; it != events.end() && it->location() <= loc; ++it)
                processEvent(*it);
        }
        else
        {
            // Is after both end and begin events, advance to next available.
            std::get<1>(loc) += 1;
            std::get<2>(loc) = false;
            for (; it != events.end() && it->location() < loc; ++it)
                processEvent(*it);
        }

        if (it == events.end())
        {
            currentLocation = events.back().location();
            std::get<2>(currentLocation) = true;
        }
        else if (!std::get<2>(loc))
        {
            std::get<0>(currentLocation) = it->rID;
            std::get<1>(currentLocation) = it->pos;
            std::get<2>(currentLocation) = false;
        }
        else
        {
            currentLocation = loc;
        }
        return currentLocation;
    }

    // Advance to the given location, at least as far as this location.
    TLocation advanceTo(int rID, int pos, bool isBegin)
    {
        while (!atEnd() && location() < std::make_tuple(rID, pos, isBegin))
            advance();
        return location();
    }

    // Begin / end iterator to the currently active/open item ids.
    TIterator activeBegin() const
    {
        return openEvents.begin();
    }

    TIterator activeEnd() const
    {
        return openEvents.end();
    }

    // Count currently active/open events.
    size_t activeCount() const
    {
        return openEvents.size();
    }

    // Return item ID for event ID.
    unsigned getItemID(unsigned eventID)
    {
        return events[eventID].itemID;
    }

private:

    // Forward declaration.
    struct Event;

    void processEvent(Event const & event)
    {
        if (event.isBegin)
        {
            SEQAN_ASSERT_NOT(openEvents.count(event.itemID));
            openEvents.insert(event.itemID);
        }
        else
        {
            openEvents.erase(event.itemID);
        }
    }

    // Make events from sliding window intervals.
    std::vector<Event> makeEvents(std::vector<SlidingWindowInterval> const & intervals) const
    {
        std::vector<Event> result;

        for (auto const & i : intervals)
        {
            result.push_back(Event(i.rID, i.beginPos, true, i.itemID));
            result.push_back(Event(i.rID, i.endPos, false, i.itemID));
        }

        std::sort(result.begin(), result.end());
        return result;
    }

    // The current events.
    std::vector<Event> events;
    // The currently open event ids.
    std::set<unsigned> openEvents;
    // The current location.
    TLocation currentLocation;
    // The current item in events.
    std::vector<Event>::const_iterator it;

    struct Event
    {
        int rID;
        int pos;
        bool isBegin;
        unsigned itemID;  // ID of item

        explicit
        Event(int rID = -1, int pos = -1, int isBegin = false, unsigned itemID = -1) :
            rID(rID), pos(pos), isBegin(isBegin), itemID(itemID)
        {}

        bool operator<(Event const & other) const
        {
            return location() < other.location();
        }

        std::tuple<int, int, bool> location() const
        {
            return std::make_tuple(rID, pos, isBegin);
        }
    };
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

#endif  // #ifndef SANDBOX_HERBARIUM_APPS_ANISE_SHARED_SLIDING_WINDOW_H_
