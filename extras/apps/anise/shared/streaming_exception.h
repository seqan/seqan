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

// TODO(holtgrew): Can we move more stuff from header to .cpp file?

#ifndef SANDBOX_HERBARIUM_APPS_ANISE_STREAMING_EXCEPTION_H_
#define SANDBOX_HERBARIUM_APPS_ANISE_STREAMING_EXCEPTION_H_

#include <stdexcept>
#include <memory>
#include <sstream>

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class StreamingException
// ----------------------------------------------------------------------------

// Handy Exception class that allow to use the following usage:
//
// throw StreamingException() << "Text " << 10 << "!";

class StreamingException : public std::runtime_error
{
public:

    StreamingException() : std::runtime_error(""), ss(new std::stringstream)
    {}

    StreamingException(StreamingException const & other) : std::runtime_error(""), ss(std::move(other.ss))
    {}

    template <typename T>
    StreamingException & operator<<(T const & x)
    {
        *ss << x;
        return *this;
    }

    virtual const char * what() const throw()
    {
        s = ss->str();
        return s.c_str();
    }

private:

    mutable std::unique_ptr<std::stringstream> ss;
    mutable std::string s;
};

// ----------------------------------------------------------------------------
// Class AniseIOException
// ----------------------------------------------------------------------------

class AniseIOException : public StreamingException
{
public:

    AniseIOException() : StreamingException()
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

#endif  // #ifndef SANDBOX_HERBARIUM_APPS_ANISE_STREAMING_EXCEPTION_H_
