// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
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

#include <iostream>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>      // For printing SeqAn Strings.
#include <seqan/misc/edit_environment.h>

int main()
{
    using namespace seqan;

    Dna5String original = "CGAT";

    typedef StringEnumerator<Dna5String, EditEnvironment<HammingDistance, 2> > THammingEnumerator;
    typedef Iterator<THammingEnumerator>::Type THammingIterator;
    std::cerr << "Enumerating Hamming distance environment of " << original << " of distance 2\n";
    THammingEnumerator hammingEnumerator(original);
    for (THammingIterator itH = begin(hammingEnumerator); !atEnd(itH); goNext(itH))
        std::cout << *itH << '\n';

    typedef StringEnumerator<Dna5String, EditEnvironment<LevenshteinDistance, 2> > TEditEnumerator;
    typedef Iterator<TEditEnumerator>::Type TEditIterator;
    std::cerr << "\nEnumerating edit distance environment of " << original << " of distance 1-2\n";
    TEditEnumerator editEnumerator(original);
    for (TEditIterator itE = begin(editEnumerator); !atEnd(itE); goNext(itE))
        std::cout << *itE << '\n';

    return 0;
}
