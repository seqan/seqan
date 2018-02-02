// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
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
// Author: Rene Rahn <rene.rahn@fu-berlin.de>
// ==========================================================================

#include <vector>
#include <array>

#include <seqan/sequence.h>
#include <seqan/stream.h>

using namespace seqan;

int main()
{
    std::vector<int> vec = {10, 12, 14};
    DnaString str = "AGT";
    std::array<double, 3> arr = { {3.14, 2.71, 1.41} };

    auto zipCont = makeZipView(vec, str, infix(str, 0, 3), arr, reverseString(str));

    // Range based for-loop.
    std::cout << "Using range based for-loop!" << std::endl;
    for (auto elem : zipCont)
    {
        std::cout << (std::get<0>(elem) += 3) << ",  " << std::get<1>(elem) << ", " << std::get<2>(elem) << ", " << std::get<3>(elem) << ", " << std::get<4>(elem) << std::endl;
    }

    // Using iterator.
    std::cout << "\nUsing iterator!" << std::endl;
    for (auto it = begin(zipCont, Standard()); it != end(zipCont, Standard()); ++it)
    {
        std::cout << (std::get<0>(*it) + 3) << ",  " << std::get<1>(*it) << ", " << std::get<2>(*it) << ", " << std::get<3>(*it) << ", " << std::get<4>(*it) << std::endl;
    }

    // Using value function and position.
    std::cout << "\nUsing value and position!" << std::endl;
    for (unsigned it = 0; it < length(zipCont); ++it)
    {
        std::cout << (std::get<0>(value(zipCont,it)) -= 3) << ",  " << std::get<1>(value(zipCont,it)) << ", " << std::get<2>(value(zipCont,it)) << ", " << std::get<3>(value(zipCont,it)) << ", " << std::get<4>(value(zipCont, it)) << std::endl;
    }
    return 0;
}
