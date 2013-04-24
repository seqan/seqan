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
// Demo of using seq_io module.
//
// We will read FASTA and FASTQ sequence files and write out the id and the
// sequence in TSV format.  Works with any combination of FASTA, FASTQ and
// gzip/bzip2 compression.  Plain text files and gzip files can also be read
// from stdin.
//
// USAGE: example_seq_io FILE
//        example_seq_io - <FILE
// ==========================================================================

#include <iostream>

#include <seqan/seq_io.h>
#include <seqan/sequence.h>

int main(int argc, char const ** argv)
{
    if (argc != 2)
    {
        std::cerr << "ERROR: Invalid arguments.\n"
                  << "USAGE: example_seq_io FILE\n";
        return 1;
    }

    seqan::SequenceStream seqStream(argv[1]);
    if (!isGood(seqStream))
    {
        std::cerr << "ERROR: Could not open file " << argv[1] << " for reading.\n";
        return 1;
    }
    seqan::CharString id, seq;

    int res = 0;
    while (!atEnd(seqStream))
    {
        if ((res = readRecord(id, seq, seqStream)) != 0)
        {
            std::cerr << "ERROR: There was an error reading from sequence file.\n";
            return 1;
        }

        std::cout << id << "\t" << seq << "\n";
    }
    
    return 0;
}
