// ==========================================================================
//                         Mason - A Read Simulator
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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

#include "simulate_genome.h"
#include <random>

// ----------------------------------------------------------------------------
// Function simulateGenome()
// ----------------------------------------------------------------------------

// Simulate a genome given the simulation options.
//
// The resulting sequence is written to stream.

int simulateGenome(seqan::SeqFileOut & stream, MasonSimulateGenomeOptions const & options)
{
    // Initialize std generator and distribution
    std::mt19937 generator(42); // seed 100
    std::uniform_real_distribution<double> distribution(0, 1);
    auto randomNumber = std::bind ( distribution, generator );

    seqan::CharString id;
    seqan::Dna5String contig;

    for (unsigned i = 0; i < length(options.contigLengths); ++i)
    {
        clear(id);
        clear(contig);

        std::stringstream ss;
        ss << (i + 1);
        id = ss.str();

        std::cerr << "contig " << id << " ...";

        for (int j = 0; j < options.contigLengths[i];)
        {
            double x = randomNumber();
            if (x < 0.25)
                appendValue(contig, 'A');
            else if (x < 0.5)
                appendValue(contig, 'C');
            else if (x < 0.75)
                appendValue(contig, 'G');
            else if (x < 1.0)
                appendValue(contig, 'T');
            else
                continue;  // Redraw.
            ++j;
        }

        try
        {
            writeRecord(stream, id, contig);
        }
        catch (seqan::IOError const & ioErr)
        {
            std::cerr << "\nERROR: Could not write contig " << id << " to output file.\n";
            return 1;
        }

        std::cerr << " DONE\n";
    }

    return 0;
}

// ----------------------------------------------------------------------------
// Function simulateGenome()
// ----------------------------------------------------------------------------

int simulateGenome(char const * filename, MasonSimulateGenomeOptions const & options)
{
    seqan::SeqFileOut stream;
    if (!open(stream, filename))
    {
        std::cerr << "ERROR: Could not open " << filename << "for writing!\n";
        return 1;
    }

    return simulateGenome(stream, options);
}
