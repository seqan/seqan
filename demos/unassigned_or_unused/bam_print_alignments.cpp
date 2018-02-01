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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Iterate over BAM alignment, dump lines in SAM format, followed by the read
// alignment.
//
// USAGE: bam_print_alignments REF.fasta ALIGN.bam
// ==========================================================================

#include <iostream>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>      // For printing SeqAn Strings.
#include <seqan/align.h>
#include <seqan/bam_io.h>
#include <seqan/seq_io.h>

using namespace seqan;

#if SEQAN_HAS_ZLIB

void trimSeqHeaderToId(CharString & header)
{
    unsigned i = 0;
    for (; i < length(header); ++i)
        if (isspace(header[i]))
            break;
    resize(header, i);
}

int main(int argc, char const ** argv)
{
    // Check command line arguments.
    if (argc != 3)
    {
        std::cerr << "USAGE: bam_print_alignments REF.fasta FILE.bam" << std::endl;
        return 1;
    }

    // Read FASTA file.
    std::cerr << "Reading FASTA " << argv[1] << std::endl;
    StringSet<CharString> refNameStore;
    StringSet<Dna5String> seqs;
    SeqFileIn inSeq;
    if (!open(inSeq, argv[1]))
    {
        std::cerr << "Could not open FASTA file " << argv[1] << std::endl;
        return 1;
    }
    readRecords(refNameStore, seqs, inSeq);
    for (unsigned i = 0; i < length(refNameStore); ++i)
        trimSeqHeaderToId(refNameStore[i]);

    // Open BGZF stream.
    std::cerr << "Opening BAM " << argv[2] << std::endl;
    BamFileIn bamFileIn;
    if (!open(bamFileIn, argv[2]))
    {
        std::cerr << "[ERROR] Could not open SAM/BAM file" << argv[2] << std::endl;
        return 1;
    }

    // Read Header.
    std::cerr << "Reading Header" << std::endl;
    NameStoreCache<StringSet<CharString> > refNameStoreCache(refNameStore);
    BamIOContext<StringSet<CharString> > context(refNameStore, refNameStoreCache);
    BamHeader header;
    readHeader(header, bamFileIn);

    // Stream through file, getting alignment and dumping it.
    std::cerr << "Reading Alignments..." << std::endl;
    Align<Dna5String> align;
    BamAlignmentRecord record;
    BamFileOut samOut(std::cout, Sam());
    while (!atEnd(bamFileIn))
    {
        readRecord(record, bamFileIn);

        if (record.rID == BamAlignmentRecord::INVALID_REFID)
            continue;  // Skip * reference.

        // Convert BAM record to alignment.
        bamRecordToAlignment(align, seqs[record.rID], record);
        // Dump record as SAM and the alignment.
        writeRecord(samOut, record);
        std::cout << align << std::endl;
    }

    return 0;
}

#else  // #if SEQAN_HAS_ZLIB

int main()
{
    std::cerr << "zlib is required for bam_print_alignment demo." << std::endl;

    return 1;
}

#endif  // #if SEQAN_HAS_ZLIB
