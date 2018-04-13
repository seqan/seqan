// ==========================================================================
//                              parse_lm_stellar
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

#include <iostream>
#include <fstream>

#include <seqan/basic.h>
#include <seqan/stream.h>
#include <seqan/parse_lm.h>

using namespace seqan;

int main(int argc, char const ** argv)
{
    // Check arguments.
    if (argc != 2)
    {
        std::cerr << "Invalid argument count.\n"
                  << "USAGE parse_lm_stellar FILE.gff\n";
        return 1;
    }

    // Open input file.
    std::fstream inStream(argv[1], std::ios::in | std::ios::binary);
    if (!inStream.good())
    {
        std::cerr << "Could not open " << argv[1] << "!\n";
        return 1;
    }

    // Read local matches in GFF Stellar format.
    DirectionIterator<std::fstream, Input>::Type reader(inStream);
    LocalMatchStore<> lmStore;
    int i = 0;
    try
    {
        while (!atEnd(reader))
        {
            readRecord(lmStore, reader, StellarGff());
            i++;
        }
    }
    catch (std::runtime_error & e)
    {
        std::cerr << "Invalid Stellar GFF record #" << i << ": " << e.what() << '\n';
        return 1;
    }

    // Dump records to stdout again.
    std::cout << "# Reverse strandness is indicated by end < begin\n"
              << "#db\tdb_beg\tdb_end\t+/-\tquery\tq_beg\tq_end\tcigar\n";
    for (unsigned i = 0; i < length(lmStore.matchStore); ++i)
    {

        std::cout << lmStore.sequenceNameStore[lmStore.matchStore[i].subjectId] << "\t"
                  << lmStore.matchStore[i].subjectBeginPos << "\t"
                  << lmStore.matchStore[i].subjectEndPos << "\t";
        if (lmStore.matchStore[i].subjectBeginPos < lmStore.matchStore[i].subjectEndPos)
            std::cout << "+\t";
        else
            std::cout << "-\t";
        std::cout << lmStore.sequenceNameStore[lmStore.matchStore[i].queryId] << "\t"
                  << lmStore.matchStore[i].queryBeginPos << "\t"
                  << lmStore.matchStore[i].queryEndPos << "\t";

        if (length(lmStore.cigarStore) > lmStore.matchStore[i].id)
        {
            String<CigarElement<> > const & cigar = lmStore.cigarStore[lmStore.matchStore[i].id];
            for (unsigned i = 0; i < length(cigar); ++i)
                std::cout << cigar[i].count << cigar[i].operation;
        }
        std::cout << '\n';
    }

    return 0;
}
