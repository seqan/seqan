// ==========================================================================
//                            generic_cyclic_shape
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
// Demo on how to use a generic CyclicShape
// ==========================================================================

#include <iostream>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>      // For printing SeqAn Strings.
#include <seqan/modifier.h>

using namespace seqan;
int main(int argc, char const ** argv)
{
    //![Define GenericCyclicShape]
    typedef CyclicShape<GenericShape> TShape;
    TShape shape;
    stringToCyclicShape(shape, "1110010");
    //![Define GenericCyclicShape]
    
    // print cyclic Shape
    CharString out;
    cyclicShapeToString(out, shape);
    std::cout << "shape: " << out << std::endl;
    std::cout << "weight: " << weight(shape)
              << ", span: " << shape.span << std::endl;
    
    //![Define GenericCyclicShape Modified String]
    CharString str = "das ist das Haus vom Nikolaus";
    ModifiedString<CharString, ModCyclicShape<TShape> > modStr (str, shape);
    std::cout << str << " => " << modStr << std::endl;
    //![Define GenericCyclicShape Modified String]
    
    
    //![CyclicShape Care Positions]
    std::cout << std::endl << "relative care positions: ";
    for(unsigned i=0; i<weight(shape); ++i)
        std::cout << shape.diffs[i] << ",";

    std::cout << std::endl << "absolute care positions: ";
    String<int> carePos;
    carePositions(carePos, shape);

    for(unsigned i=0; i<weight(shape); ++i)
        std::cout << carePos[i] << ",";
    std::cout << std::endl;
    //![CyclicShape Care Positions]
    
    return 0;
}
