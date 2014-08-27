// ==========================================================================
//                                  ANISE
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

#include "progress_indicator.h"

// TODO(holtgrew): Better fallback for non-terminal mode.
#include <seqan/misc/misc_terminal.h>

#include <atomic>
#include <mutex>
#include <iostream>
#include <chrono>
#include <stdexcept>
#include <thread>
#include <sstream>
#include <cmath>

namespace  // anonymous namespace
{

// ----------------------------------------------------------------------------
// Class ProgressState
// ----------------------------------------------------------------------------

// Stores the progress state.

class ProgressState
{
    // Advance progress in [min, max], current position is current.
    std::atomic_llong min, max, current;
    // Mutex for updating.
    std::recursive_mutex mutex;

public:

    explicit ProgressState(long long int min = 0, long long int max = 100) : min(min), max(max), current(min)
    {}

    ProgressState(long long int min, long long int max, long long int current) : min(min), max(max), current(current)
    {}

    // Advance to the given position, must be in [min, max].

    virtual void advanceTo(long long int dest)
    {
        std::lock_guard<std::recursive_mutex> lock(mutex);

        if (dest < min.load() || dest > max.load())
            return;  // invalid progress target
            //throw std::runtime_error("Invalid progress target.");
        current = dest;
    }

    // Advance forwards (or backwards) by delta steps.

    virtual void advanceBy(long long int delta)
    {
        std::lock_guard<std::recursive_mutex> lock(mutex);

        long long int dest = current.load() + delta;
        advanceTo(dest);
    }

    // Force advancing to completion.

    virtual void complete()
    {
        std::lock_guard<std::recursive_mutex> lock(mutex);

        current.store(max.load());
    }

    // Return the min/max/current step count.

    long long int minSteps() const
    {
        return min;
    }

    long long int maxSteps() const
    {
        return max;
    }

    long long int currentSteps() const
    {
        return current;
    }

    // Return completion ratio.

    double completionRatio() const
    {
        if (max == min)
            return 0;
        return 1.0 * (current - min) / (max - min);
    }
};

// ----------------------------------------------------------------------------
// Class EstimationImpossible
// ----------------------------------------------------------------------------

// Handles the time estimation.

class EstimationImpossible : public std::runtime_error
{
public:
    explicit EstimationImpossible(std::string const & whatArg) : std::runtime_error(whatArg)
    {}
};

// ----------------------------------------------------------------------------
// Class EstimatingProgressState
// ----------------------------------------------------------------------------

class EstimatingProgressState : public ProgressState
{
    // The start time of the process.
    std::chrono::time_point<std::chrono::system_clock> startTime;

public:

    explicit EstimatingProgressState(long long int min = 0, long long int max = 100) :
        ProgressState(min, max), startTime(std::chrono::system_clock::now())
    {}

    EstimatingProgressState(long long int min, long long int max, long long int current) :
        ProgressState(min, max, current), startTime(std::chrono::system_clock::now())
    {}

    // Get elapsed time in seconds.

    std::chrono::duration<double> elapsed() const
    {
        return std::chrono::system_clock::now() - startTime;
    }

    // Get estimate based on assumed linear progress.
    //
    // Throws EstimationImpossible in case that there would be a division by zero.

    std::chrono::duration<double> estimateTotal() const
    {
        long long int stepsDone = currentSteps() - minSteps();
        if (!stepsDone)
            return std::chrono::duration<double>();  // throw EstimationImpossible("No step completed yet.");
        long long int stepsToDo = maxSteps() - minSteps();
        if (!stepsToDo)
            return startTime - startTime;  // nothing to do
        double doneRatio = 1.0 * stepsDone / stepsToDo;
        return elapsed() / doneRatio;
    }
};

// Display the bar part for a progress bar.

class ProgressBarDisplay
{
    long long int width;  // number of characters to use
    char barChar;  // Bar character to use
    char emptyChar;  // Empty character to use.
    char first, last;  // start/end characters, '\0' if none

public:

    explicit
    ProgressBarDisplay(long long int width = 50, char barChar = '#', char first = '|',
                       char last = '|') :
        width(width), barChar(barChar), emptyChar(' '), first(first), last(last)
    {}

    // Prlong long int bar for completion ratio.

    void print(std::ostream & out, double ratio)
    {
        // Number of bar characters to use.
        long long int count = width - (first != '\0') - (first != '\0');

        if (first != '\0')
            out << first;

        out << std::string(ceil(ratio * count), barChar);
        out << std::string(count - ceil(ratio * count), emptyChar);

        if (last != '\0')
            out << last;
    }
};

}  // anonymous namespace

// ----------------------------------------------------------------------------
// Class ProgressBarImpl
// ----------------------------------------------------------------------------

// Displays a progress bar.
//
// Requires a terminal to write to.

class ProgressBarImpl : public EstimatingProgressState
{
    std::ostream * out;  // out stream to write to
    long long int lineWidth;  // current line width

    std::atomic<bool> done;  // display newline if not done on destruction

    // Mutex for the display.
    std::recursive_mutex mutex;

public:

    bool active;
    std::string label;
    ProgressBarDisplay display;
    // Thread used for updating the display.
    std::thread updateThread;
    // Does output go to a terminal (i.e. '\r' works?)
    bool isTerminal;

    explicit ProgressBarImpl(std::ostream & out, long long int min = 0, long long int max = 100, bool active = true) :
            EstimatingProgressState(min, max), out(&out), lineWidth(0), done(false), active(active),
            updateThread([this,active]() { if (active) this->runUpdate(); }), isTerminal(seqan::isTerminal())
    {}

    ProgressBarImpl(std::ostream & out, long long int min, long long int max, long long int current, bool active = true) :
            EstimatingProgressState(min, max, current), out(&out), lineWidth(0), done(false), active(active),
            updateThread([this,active]() { if (active) this->runUpdate(); }), isTerminal(seqan::isTerminal())
    {}

    ~ProgressBarImpl()
    {
        done = true;
        updateThread.join();
    }

    void advanceTo(long long int dest)
    {
        ProgressState::advanceTo(dest);
        if (!isTerminal)
            updateDisplay();
    }

    void advanceBy(long long int delta)
    {
        ProgressState::advanceBy(delta);
        if (!isTerminal)
            updateDisplay();
    }

    // Force advancing to completion.

    void complete()
    {
        ProgressState::complete();
        if (!isTerminal)
            updateDisplay();
    }

    // Update display.
    void updateDisplay()
    {
        if (!active)
            return;

        std::lock_guard<std::recursive_mutex> lock(mutex);

        std::stringstream ss;

        // Estimate time.
        double estimated = 0;
        try
        {
            estimated = estimateTotal().count();
        }
        catch (EstimationImpossible e)
        {
            estimated = -1;
        }

        ss << label << " ";
        display.print(ss, completionRatio());
        if (estimated < 0)
            ss << " " << currentSteps() - minSteps() << "/" << maxSteps() - minSteps() << " [" << fmtDuration(elapsed()) << "]";
        else
            ss << " " << currentSteps() - minSteps() << "/" << maxSteps() - minSteps() << " [" << fmtDuration(elapsed())
           << "/" << fmtDuration(estimateTotal()) << "]";

        // clearOutLine();
        writeLine(ss.str());
    }

    // Convert duration to HH:MM:SS display.
    std::string fmtDuration(std::chrono::duration<double> const & duration)
    {
        double dbl = duration.count();
        long long int hours = dbl / 60 / 60;
        dbl -= hours * 60 * 60;
        long long int minutes = dbl / 60;
        dbl -= minutes * 60;
        long long int seconds = dbl;

        char buffer[100];
        snprintf(buffer, 100, "%02lld:%02lld:%02lld", hours, minutes, seconds);
        return buffer;
    }

    // Finish display.
    void finish()
    {
        if (done || !active)
        {
            done = true;
            return;
        }
        // Complete progress and update display one last time.
        complete();
        updateDisplay();
        // Write newline.
        if (isTerminal)
            *out << "\n";
        lineWidth = 0;
    }

private:

    // Clear the output line.
    void clearOutLine()
    {
        if (isTerminal)
        {
            *out << '\r';
            for (long long int i = 0; i < lineWidth; ++i)
                *out << ' ';
            *out << '\r';
        }
        else
        {
            *out << '\n';
        }
        *out << std::flush;
        lineWidth = 0;
    }

    // Write to out, updating lineWidth.
    void writeLine(std::string line)
    {
        std::string filler((line.size() < (size_t)lineWidth) ? lineWidth - line.size() : 0, ' ');
        if (isTerminal)
            *out << '\r';
        *out << line << filler;
        if (!isTerminal)
            *out << '\n';
        *out << std::flush;
        lineWidth = std::max((long long int)line.size(), lineWidth);
    }

    // Called for updating the display.
    void runUpdate()
    {
        if (!isTerminal)
            return;  // no update in case of no terminal
        while (!done && maxSteps() != currentSteps())
        {
            this->updateDisplay();
            std::chrono::milliseconds dura(500);
            std::this_thread::sleep_for(dura);
        }
    }
};

// ----------------------------------------------------------------------------
// Class ProgressBar
// ----------------------------------------------------------------------------

ProgressBar::ProgressBar(std::ostream & out, long long int min, long long int max, bool active) :
        impl(new ProgressBarImpl(out, min, max, active))
{}

ProgressBar::ProgressBar(std::ostream & out, long long int min, long long int max, long long int current, bool active) :
        impl(new ProgressBarImpl(out, min, max, current, active))
{}

ProgressBar::~ProgressBar()
{}

void ProgressBar::setLabel(char const * label)
{
    impl->label = label;
}

void ProgressBar::setActive(bool state)
{
    impl->active = state;
}

void ProgressBar::advanceBy(long long int delta)
{
    impl->advanceBy(delta);
}

void ProgressBar::advanceTo(long long int dest)
{
    impl->advanceTo(dest);
}

void ProgressBar::updateDisplay()
{
    impl->updateDisplay();
}

void ProgressBar::finish()
{
    impl->finish();
}
