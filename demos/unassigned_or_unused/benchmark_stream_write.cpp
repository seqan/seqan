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
// Author: Hannes Hauswedell <hauswedell@mi.fu-berlin.de>
// ==========================================================================
// Very simple benchmarking tool for stream writing
// ==========================================================================

#include <cstdio>
#include <fstream>
#if SEQAN_HAS_ZLIB
#include <zlib.h>
#endif  // #if SEQAN_HAS_ZLIB
#if SEQAN_HAS_BZIP2
#include <bzlib.h>
#endif  // #if SEQAN_HAS_BZIP2

#include <seqan/basic.h>
#include <seqan/stream.h>
#include <seqan/stream.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>

#if defined(STDLIB_VS)
#define UNLINK _unlink
#else
#define UNLINK unlink
#endif

const int MB = 1024 * 1024;

using namespace seqan;

template <typename TMetas, typename TSeqs>
void constructFastaStrings(TMetas & metas, TSeqs & seqs)
{
    resize(metas, 4);
    resize(seqs, 4);

    reserve(seqs[0], MB);
    metas[0] = "1MB of as";
    for (int i = 0; i < MB; ++i)
        appendValue(seqs[0], 'a');

    reserve(seqs[1], 32 * MB);
    metas[1] = "32MB of as";
    for (int i = 0; i < 32 * MB; ++i)
        appendValue(seqs[1], 'c');

    reserve(seqs[2], 256 * MB);
    metas[2] = "256MB of gs";
    for (int i = 0; i < 256 * MB; ++i)
        appendValue(seqs[2], 'g');

    reserve(seqs[3], 512 * MB);
    metas[3] = "512MB of ts";
    for (int i = 0; i < 512 * MB; ++i)
        appendValue(seqs[3], 't');
}

template <typename TMetas, typename TSeqs>
void doIt(TMetas const & metas, TSeqs const & seqs)
{
    double before = 0;
    double after = 0;
    {
        CharString tempFilename = SEQAN_TEMP_FILENAME();
        char filenameBuffer[1000];
        strncpy(filenameBuffer, toCString(tempFilename), 999);

        std::fstream fileStream(filenameBuffer, std::ios_base::in | std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
        SeqFileOut file(fileStream, Fasta());
        context(file).options = SequenceOutputOptions(0);

        before = sysTime();
        std::cerr << "Writing with new IO (no linebreaks, fstream) .... " << std::flush;
        for (int i = 0; i < 4; ++i)
            writeRecord(file, metas[i], seqs[i]);
        close(file);
        after = sysTime();
        std::cerr << "completed in " << after - before << "s\n" << std::flush;
        UNLINK(filenameBuffer);
    }
    {
        CharString tempFilename = SEQAN_TEMP_FILENAME();
        char filenameBuffer[1000];
        strncpy(filenameBuffer, toCString(tempFilename), 999);

        String<char, MMap<> > mmapString;
        open(mmapString, filenameBuffer);

        before = sysTime();
        std::cerr << "Writing with new IO (no linebreaks, MMAP) .... " << std::flush;
        for (int i = 0; i < 4; ++i)
            writeRecord(mmapString, metas[i], seqs[i], Fasta(), SequenceOutputOptions(0));
        close(mmapString);
        after = sysTime();
        std::cerr << "completed in " << after - before << "s\n" << std::flush;
        UNLINK(filenameBuffer);
    }
//    {
//        CharString tempFilename = SEQAN_TEMP_FILENAME();
//        char filenameBuffer[1000];
//        strncpy(filenameBuffer, toCString(tempFilename), 999);
//        FILE * file = fopen(filenameBuffer, "wb+");
//
//        before = sysTime();
//        std::cerr << "Writing with new IO (no linebreaks,  cstdio) .... " << std::flush;
//        for (int i = 0; i < 4; ++i)
//            writeRecord(file, metas[i], seqs[i], Fasta(), SequenceOutputOptions(0));
//        fclose(file);
//        after = sysTime();
//        std::cerr << "completed in " << after - before << "s\n"<< std::flush;
//        UNLINK(filenameBuffer);
//    }
    {
        CharString tempFilename = SEQAN_TEMP_FILENAME();
        char filenameBuffer[1000];
        strncpy(filenameBuffer, toCString(tempFilename), 999);

        std::fstream file(filenameBuffer, std::ios_base::in | std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
        SEQAN_ASSERT(file.is_open());
        DirectionIterator<std::fstream, Output>::Type iter(file);

        before = sysTime();
        std::cerr << "Writing with new IO (with linebreaks, fstream) .... " << std::flush;
        for (int i = 0; i < 4; ++i)
            writeRecord(iter, metas[i], seqs[i], Fasta());
        close(file);
        after = sysTime();
        std::cerr << "completed in " << after - before << "s\n" << std::flush;
        UNLINK(filenameBuffer);
    }
    {
        CharString tempFilename = SEQAN_TEMP_FILENAME();
        char filenameBuffer[1000];
        strncpy(filenameBuffer, toCString(tempFilename), 999);

        String<char, MMap<> > mmapString;
        open(mmapString, filenameBuffer);

        before = sysTime();
        std::cerr << "Writing with new IO (with linebreaks, MMAP) .... to " << std::flush;
        for (int i = 0; i < 4; ++i)
            writeRecord(mmapString, metas[i], seqs[i], Fasta());
        close(mmapString);
        after = sysTime();
        std::cerr << "completed in " << after - before << "s\n" << std::flush;
        UNLINK(filenameBuffer);
    }
//    {
//        CharString tempFilename = SEQAN_TEMP_FILENAME();
//        char filenameBuffer[1000];
//        strncpy(filenameBuffer, toCString(tempFilename), 999);
//        FILE * file = fopen(filenameBuffer, "wb+");
//
//        before = sysTime();
//        std::cerr << "Writing with new IO (with linebreaks, cstdio) .... " << std::flush;
//        for (int i = 0; i < 4; ++i)
//            writeRecord(file, metas[i], seqs[i], Fasta());
//        fclose(file);
//        after = sysTime();
//        std::cerr << "completed in " << after - before << "s\n"<< std::flush;
//        UNLINK(filenameBuffer);
//    }
    {
        CharString tempFilename = SEQAN_TEMP_FILENAME();
        char filenameBuffer[1000];
        strncpy(filenameBuffer, toCString(tempFilename), 999);

        std::fstream file(filenameBuffer, std::ios_base::in | std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
        SEQAN_ASSERT(file.is_open());
        DirectionIterator<std::fstream, Output>::Type iter(file);

        before = sysTime();
        std::cerr << "Writing with old IO.... to " << std::flush;
        for (int i = 0; i < 4; ++i)
            writeRecord(iter, seqs[i], metas[i], Fasta());
        close(file);
        after = sysTime();
        std::cerr << "completed in " << after - before << "s\n" << std::flush;
        UNLINK(filenameBuffer);
    }
}

int main(int argc, char const ** argv)
{
    (void)argc;
    (void)argv;

    {
        std::cerr << "START\n -- using CharStrings\n" << std::flush;
        StringSet<CharString> metas;
        StringSet<CharString> seqs;

        double before = sysTime();
        std::cerr << "Constructing Strings ...." << std::flush;
        constructFastaStrings(metas, seqs);
        double after = sysTime();
        std::cerr << "completed in " << after - before << "s\n" << std::flush;

        for (int i = 0; i < 5; ++i)
        {
            std::cerr << "RUN No" << i << "\n" << std::flush;
            doIt(metas, seqs);
        }
    }
    {
        std::cerr << "START\n -- using DnaStrings\n" << std::flush;
        StringSet<CharString> metas;
        StringSet<DnaString> seqs;

        double before = sysTime();
        std::cerr << "Constructing Strings ...." << std::flush;
        constructFastaStrings(metas, seqs);
        double after = sysTime();
        std::cerr << "completed in " << after - before << "s\n" << std::flush;

        for (int i = 0; i < 5; ++i)
        {
            std::cerr << "RUN No" << i << "\n" << std::flush;
            doIt(metas, seqs);
        }
    }
    return 0;
}
