// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2016, Knut Reinert, FU Berlin
// Copyright (c) 2013 NVIDIA Corporation
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

#ifndef PLATFORM_WINDOWS
#define PLATFORM_WINDOWS

#define PLATFORM_WINDOWS_VS

// ==========================================================================
// Compiler Defines For MSVC.
// ==========================================================================

// Make <windows.h> not define min() and max() as macros.
#ifndef NOMINMAX
#define NOMINMAX
#endif  // #ifndef NOMINMAX

// ==========================================================================
// Disable Warnings
// ==========================================================================

// Disable warning for identifer name truncation.  There is not much we can
// do about this.  Boost also has this problem and they chose to suppress
// it globally.  So did we.
//
// Documentation of C4504 from Microsoft:
//   http://msdn.microsoft.com/en-us/library/074af4b6%28v=vs.80%29.aspx
// Boost Warnings Guidelines:
//   https://svn.boost.org/trac/boost/wiki/Guidelines/WarningsGuidelines
#pragma warning( disable : 4503 )

// Disable warning for C++ compliant behaviour for default-initializing
// arrays in classes.
//
// Documentation of C4345 from Microsoft:
//   http://msdn.microsoft.com/de-de/library/1ywe7hcy(v=vs.80).aspx
// Documentation of C4351 from Microsoft:
//   http://msdn.microsoft.com/en-us/library/wewb47ee(v=vs.80).aspx
#pragma warning( disable : 4345 )
#pragma warning( disable : 4351 )

// Disable warning for "this" used in derived c'tor
// Documentation of C4355 from Microsoft:
//   https://msdn.microsoft.com/en-us/library/3c594ae3(v=vs.100).aspx
#pragma warning( disable : 4355 )

// ==========================================================================
// Visual Studio Specific Workarounds.
// ==========================================================================

// Workaround for missing round() from C99 in Visual Studio.
template <typename T>
inline T round(T const & x)
{
    return static_cast<T>(floor(x + 0.5));
}

#endif  // #ifndef PLATFORM_WINDOWS
