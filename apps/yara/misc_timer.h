// ==========================================================================
//                      Yara - Yet Another Read Aligner
// ==========================================================================
// Copyright (c) 2011-2018, Enrico Siragusa, FU Berlin
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
//     * Neither the name of Enrico Siragusa or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ENRICO SIRAGUSA OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================

#ifndef APP_YARA_MISC_TIMER_H_
#define APP_YARA_MISC_TIMER_H_

using namespace seqan;

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class Timer
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec = void>
struct Timer
{
    TValue _begin;
    TValue _end;

    Timer() :
        _begin(0),
        _end(0)
    {};
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function start()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
inline void start(Timer<TValue, TSpec> & timer)
{
    timer._begin = sysTime();
}

// ----------------------------------------------------------------------------
// Function stop()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
inline void stop(Timer<TValue, TSpec> & timer)
{
    timer._end = sysTime();
}

// ----------------------------------------------------------------------------
// Function getValue()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
inline TValue getValue(Timer<TValue, TSpec> const & timer)
{
    return timer._end - timer._begin;
}

template <typename TValue, typename TSpec>
inline TValue getValue(Timer<TValue, TSpec> & timer)
{
    return getValue(static_cast<Timer<TValue, TSpec> const &>(timer));
}

// ----------------------------------------------------------------------------
// Function operator<<
// ----------------------------------------------------------------------------

template <typename TStream, typename TValue, typename TSpec>
inline TStream & operator<<(TStream & os, Timer<TValue, TSpec> const & timer)
{
    os << getValue(timer) << " sec";
    return os;
}

// ----------------------------------------------------------------------------
// Function printRuler()
// ----------------------------------------------------------------------------

template <typename TStream>
inline void printRuler(TStream & os)
{
    os << std::endl
       << "================================================================================"
       << std::endl << std::endl;
}

#endif // APP_YARA_MISC_TIMER_H_
