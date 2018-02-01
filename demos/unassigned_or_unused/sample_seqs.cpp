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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de.>
// ==========================================================================
// Sample NUM reads from input file FILE.
// ==========================================================================

#include <iostream>
#include <set>
#include <random>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/seq_io.h>

using namespace seqan;

const unsigned SEED = 42;

int main(int argc, char const ** argv)
{
    // Check command line count.
    if (argc != 3)
    {
        std::cerr << "Invalid number of arguments.\n"
                  << "USAGE: sample_seqs IN.{fasta,fastq} NUM\n";
        return 1;
    }

    // Get number of sequences to sample.
    unsigned num = atoi(argv[2]);

    // Open file.
    std::cerr << "Opening file..." << std::endl;
    SeqFileIn seqFileIn;
    if (!open(seqFileIn, argv[1]))
    {
        std::cerr << "Could not open file " << argv[1] << "\n";
        return 1;
    }

    String<CharString> ids;
    String<String<Dna5Q> > seqs;
    readRecords(ids, seqs, seqFileIn);

    // Sanity check on number to sample.
    if (num > length(seqs))
    {
        std::cerr << "Request to sample more reads than there actually are!" << std::endl;
        return 0;
    }
    if (2 * num > length(seqs))
    {
        std::cerr << "WARNING There are " << length(seqs)
                  << " reads, we want to sample " << num
                  << ", it make quite some time!" << std::endl;
    }

    // Now, sample reads to pick.
    std::mt19937 rng(SEED);
    std::uniform_int_distribution<unsigned> pdf(0, length(seqs) - 1);
    std::cerr << "Sampling ids..." << std::endl;
    std::set<unsigned> sampledIds;
    while (sampledIds.size() < num)
    {
        unsigned x = pdf(rng);
        sampledIds.insert(x);
    }

    // Finally, sample reads.
    SeqFileOut seqFileOut(std::cout, Fastq());
    std::set<unsigned>::iterator it = sampledIds.begin();
    for (; it != sampledIds.end(); ++it)
        writeRecord(seqFileOut, ids[*it], seqs[*it]);

    return 0;
}
