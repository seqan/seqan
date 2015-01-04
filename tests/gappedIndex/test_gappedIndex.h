// ==========================================================================
//                             test_gappedIndex.h
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
// Author: Your Name <your.email@example.net>
// ==========================================================================

#ifndef CORE_TESTS_GAPPEDINDEX_TEST_GAPPEDINDEX_H_
#define CORE_TESTS_GAPPEDINDEX_TEST_GAPPEDINDEX_H_

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct _ShapeDefs
{
    CyclicShape<FixedShape<0,GappedShape<HardwiredShape<> >, 1> >       S_10;
    CyclicShape<FixedShape<0,GappedShape<HardwiredShape<1,2> >, 1> >    S_11010;
    CyclicShape<FixedShape<0,GappedShape<HardwiredShape<1,1,1> >, 2> >  S_111100;
    CyclicShape<FixedShape<0,GappedShape<HardwiredShape<4> >, 0> >      S_10001;
    CyclicShape<FixedShape<1,GappedShape<HardwiredShape<> >, 0> >       S_01;
    CyclicShape<FixedShape<2,GappedShape<HardwiredShape<1> >, 0> >      S_0011;
    
    CyclicShape<GenericShape> s_10;
    CyclicShape<GenericShape> s_11010;
    CyclicShape<GenericShape> s_111100;
    CyclicShape<GenericShape> s_10001;
    CyclicShape<GenericShape> s_01;
    CyclicShape<GenericShape> s_0011;
    
    _ShapeDefs()
    {
        stringToCyclicShape(s_10,       "10");
        stringToCyclicShape(s_11010,    "11010");
        stringToCyclicShape(s_111100,   "111100");
        stringToCyclicShape(s_10001,    "10001");
        stringToCyclicShape(s_01,       "01");
        stringToCyclicShape(s_0011,     "0011");
    }
};


// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

#endif  // #ifndef CORE_TESTS_GAPPEDINDEX_TEST_GAPPEDINDEX_H_
