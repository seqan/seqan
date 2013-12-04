// ==========================================================================
//                             fixed_cylcic_shape
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
// Author: Sascha Meiers <meiers@inf.fu-berlin.de>
// ==========================================================================
// Demo on how to use a Fixed-at-compile-time CyclicShape
// ==========================================================================

#include <iostream>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>      // For printing SeqAn Strings.
#include <seqan/modifier.h>

using namespace seqan;
int main(int argc, char const ** argv)
{
    //![Define FixedCyclicShape]
    typedef GappedShape<HardwiredShape<1, 1, 3> >       TInnerShape; // 111001
    // You can also use predefied Shapes such as Patternhunter
    
    typedef CyclicShape<FixedShape<0, TInnerShape, 1> > TShape;      // 1110010
    TShape shape;
    //![Define FixedCyclicShape]
    
    // print cyclic Shape
    CharString out;
    cyclicShapeToString(out, shape);
    std::cout << "shape: " << out << std::endl;
    std::cout << "weight: " << (unsigned)weight(shape)  // alternative: WEIGHT<TShape>::VALUE
    << ", span: " << shape.span << std::endl;           // alternative: TShape::span
    
    //![Define FixedCyclicShape Modified String]
    CharString str = "das ist das Haus vom Nikolaus";
    ModifiedString<CharString, ModCyclicShape<TShape> > modStr (str);
    std::cout << str << " => " << modStr << std::endl;
    //![Define FixedCyclicShape Modified String]
    return 0;
}
