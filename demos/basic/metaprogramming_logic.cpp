// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Example program used in the documentation of the metaprogramming functions.
// ==========================================================================

#include <iostream>

#include <seqan/basic.h>

// Demonstration of using True and False as base classes.
//![inheriting from true false]
template <typename T>
struct IsInt32 : seqan::False {};

template <>
struct IsInt32<int> : seqan::True {};
//![inheriting from true false]

// Helper functions for print tags True/False.
//![true false print helpers]
void printBoolType(seqan::True const & /*tag*/)
{
    std::cout << "true" << std::endl;
}

void printBoolType(seqan::False const & /*tag*/)
{
    std::cout << "false" << std::endl;
}
//![true false print helpers]

int main()
{
    // Example for using the tags True and False.
//![tags true false]
    std::cout << seqan::False::VALUE << "\n"                                         // => "0"
              << seqan::True::VALUE << "\n"                                          // => "1"
              << seqan::IsSameType<seqan::False, seqan::False::Type>::VALUE << "\n"; // => "1"
//![tags true false]

    // And an example for using IsInt32.
//![using isint32]
    std::cout << IsInt32<bool>::VALUE << "\n"  // => "0"
              << IsInt32<int>::VALUE << "\n";  // => "1"
//![using isint32]

    // Demonstration of using Eval.
//![print bool type eval]
    printBoolType(seqan::Eval<true>::Type());  // => "true"
    printBoolType(seqan::Eval<false>::Type()); // => "false"
//![print bool type eval]

    // Demonstration of using Not.
//![print bool type not]
    printBoolType(seqan::Not<seqan::False>::Type());  // => "true"
    printBoolType(seqan::Not<seqan::True>::Type());   // => "false"
//![print bool type not]

    // Demonstration of using NotC.
//![print bool type notc]
    printBoolType(seqan::NotC<false>::Type());  // => "true"
    printBoolType(seqan::NotC<true>::Type());   // => "false"
//![print bool type notc]

    // Demonstration of using Or.
//![print bool type or]
    printBoolType(seqan::Or<seqan::False, seqan::False>::Type());  // => "false"
    printBoolType(seqan::Or<seqan::False, seqan::True >::Type());  // => "true"
    printBoolType(seqan::Or<seqan::True,  seqan::False>::Type());  // => "true"
    printBoolType(seqan::Or<seqan::True,  seqan::True >::Type());  // => "true"
//![print bool type or]

    // Demonstration of using OrC.
//![print bool type orc]
    printBoolType(seqan::OrC<false, false>::Type());  // => "false"
    printBoolType(seqan::OrC<false, true >::Type());  // => "true"
    printBoolType(seqan::OrC<true,  false>::Type());  // => "true"
    printBoolType(seqan::OrC<true,  true >::Type());  // => "true"
//![print bool type orc]

    // Demonstration of using And.
//![print bool type and]
    printBoolType(seqan::And<seqan::False, seqan::False>::Type());  // => "false"
    printBoolType(seqan::And<seqan::False, seqan::True >::Type());  // => "true"
    printBoolType(seqan::And<seqan::True,  seqan::False>::Type());  // => "true"
    printBoolType(seqan::And<seqan::True,  seqan::True >::Type());  // => "true"
//![print bool type and]

    // Demonstration of using AndC.
//![print bool type andc]
    printBoolType(seqan::AndC<false, false>::Type());  // => "false"
    printBoolType(seqan::AndC<false, true >::Type());  // => "true"
    printBoolType(seqan::AndC<true,  false>::Type());  // => "true"
    printBoolType(seqan::AndC<true,  true >::Type());  // => "true"
//![print bool type andc]

    // Demonstration of using If.
//![print bool type if]
    printBoolType(seqan::If<seqan::True, seqan::True,  seqan::False>::Type());  // => "true"
    printBoolType(seqan::If<seqan::True, seqan::False, seqan::True >::Type());  // => "false"
//![print bool type if]

    // Demonstration of using IfC.
//![print bool type ifc]
    printBoolType(seqan::If<seqan::True, seqan::True,  seqan::False>::Type());  // => "true"
    printBoolType(seqan::If<seqan::True, seqan::False, seqan::True >::Type());  // => "false"
//![print bool type ifc]

    // Demonstration of the shortcut-to-::Type feature.
//![shortcut to type feature]
    typedef seqan::And<seqan::Or<seqan::True, seqan::False>, seqan::True> TResult;
    printBoolType(TResult());  // => "true"
//![shortcut to type feature]

    return 0;
}
