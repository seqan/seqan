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
// Demo for using the repeat finding algorithm (period 1 case).  The first
// sequence is read from a FASTA file and all repeats are printed.
// ==========================================================================

// Comment out/uncomment for using sequential/parallel repeat finding code.
// Parallel code is enabled automatically when OpenMP is available, use this
// preprocessor define to disable it in any case.
// #define SEQAN_ENABLE_PARALLELISM 0

#include <fstream>

#include <seqan/basic.h>
#include <seqan/index.h>
#include <seqan/seq_io.h>

int main(int argc, char ** argv)
{
    using namespace seqan;

    // Check parameters.
    if (argc != 2)
    {
        std::cerr << "USAGE: find_repeats IN.fasta";
        return 1;
    }

    // Load first sequence from fasta file.
    std::ifstream inFile(argv[1]);
    if (!inFile.good())
    {
        std::cerr << "Could not open " << argv[1] << std::endl;
        return 1;
    }
    RecordReader<std::ifstream, SinglePass<> > reader(inFile);
    CharString id;
    Dna5String seq;
    if (readRecord(id, seq, reader, Fasta()))
    {
        std::cerr << "Could not read sequence!" << std::endl;
        return 1;
    }

    // Find repeats and print them.
    String<Repeat<unsigned, unsigned> > repeats;
    findRepeats(repeats, seq, 1000);
    std::cerr << "# of repeats: " << length(repeats) << std::endl;
    for (unsigned i = 0; i < length(repeats); ++i)
    {
        std::cerr << "i == " << i << ", beginPosition = " << repeats[i].beginPosition
                  << ", endPosition = " << repeats[i].endPosition
                  << ", period = " << repeats[i].period << std::endl;
    }

    return 0;
}
