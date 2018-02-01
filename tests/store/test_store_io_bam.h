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
// Tests for the SeqAn model store, I/O functionality.
// ==========================================================================

#include <seqan/basic.h>  // For test functionality.
#include <seqan/store.h>  // Header under test.

using namespace seqan;

SEQAN_DEFINE_TEST(test_store_io_bam_read)
{
  /*
    // Construct name to reference FASTA files.
    char fastaBuffer[1023];
    strcpy(fastaBuffer, getAbsolutePath(""));
    strcat(fastaBuffer, "/projects/tests/store/toy.fa");
    // Construct file name to SAM file.
    char samBuffer[1023];
    strcpy(samBuffer, getAbsolutePath(""));
    strcat(samBuffer, "/projects/tests/store/toy.sam");
    // Construct file name to BAM file.
    char bamBuffer[1023];
    strcpy(bamBuffer, getAbsolutePath(""));
    strcat(bamBuffer, "/projects/tests/store/toy.bam");

    // Load FragmentStore from SAM and BAM file.
    FragmentStore<> samStore;
    loadContigs(samStore, fastaBuffer);
    FILE * samFp = fopen(samBuffer, "rb");
    read(samFp, samStore, Sam());
    fclose(samFp);

    FragmentStore<> bamStore;
    loadContigs(bamStore, fastaBuffer);
    samfile_t * bamFp = samopen(bamBuffer, "rb", 0);
    read(bamFp, bamStore, Bam());
    samclose(bamFp);

    // Check that the stores are the same.
    SEQAN_ASSERT_EQ(length(samStore.alignedReadStore), length(bamStore.alignedReadStore));
    for (unsigned i = 0; i < length(samStore.alignedReadStore); ++i)
        SEQAN_ASSERT(samStore.alignedReadStore[i] == bamStore.alignedReadStore[i]);
    SEQAN_ASSERT_EQ(length(samStore.alignedReadTagStore), length(bamStore.alignedReadTagStore));
    for (unsigned i = 0; i < length(samStore.alignedReadTagStore); ++i)
        SEQAN_ASSERT(samStore.alignedReadTagStore[i] == bamStore.alignedReadTagStore[i]);
    SEQAN_ASSERT_EQ(length(samStore.alignQualityStore), length(bamStore.alignQualityStore));
    for (unsigned i = 0; i < length(samStore.alignQualityStore); ++i)
        SEQAN_ASSERT(samStore.alignQualityStore[i] == bamStore.alignQualityStore[i]);
    SEQAN_ASSERT_EQ(length(samStore.contigFileStore), length(bamStore.contigFileStore));
    for (unsigned i = 0; i < length(samStore.contigFileStore); ++i)
        SEQAN_ASSERT(samStore.contigFileStore[i] == bamStore.contigFileStore[i]);
    SEQAN_ASSERT_EQ(length(samStore.contigNameStore), length(bamStore.contigNameStore));
    for (unsigned i = 0; i < length(samStore.contigNameStore); ++i)
        SEQAN_ASSERT(samStore.contigNameStore[i] == bamStore.contigNameStore[i]);
    SEQAN_ASSERT_EQ(length(samStore.contigStore), length(bamStore.contigStore));
    for (unsigned i = 0; i < length(samStore.contigStore); ++i)
        SEQAN_ASSERT(samStore.contigStore[i] == bamStore.contigStore[i]);
    SEQAN_ASSERT_EQ(length(samStore.libraryNameStore), length(bamStore.libraryNameStore));
    for (unsigned i = 0; i < length(samStore.libraryNameStore); ++i)
        SEQAN_ASSERT(samStore.libraryNameStore[i] == bamStore.libraryNameStore[i]);
    SEQAN_ASSERT_EQ(length(samStore.libraryStore), length(bamStore.libraryStore));
    for (unsigned i = 0; i < length(samStore.libraryStore); ++i)
        SEQAN_ASSERT(samStore.libraryStore[i] == bamStore.libraryStore[i]);
    SEQAN_ASSERT_EQ(length(samStore.matePairNameStore), length(bamStore.matePairNameStore));
    for (unsigned i = 0; i < length(samStore.matePairNameStore); ++i)
        SEQAN_ASSERT(samStore.matePairNameStore[i] == bamStore.matePairNameStore[i]);
    SEQAN_ASSERT_EQ(length(samStore.matePairStore), length(bamStore.matePairStore));
    for (unsigned i = 0; i < length(samStore.matePairStore); ++i)
        SEQAN_ASSERT(samStore.matePairStore[i] == bamStore.matePairStore[i]);
    SEQAN_ASSERT_EQ(length(samStore.readNameStore), length(bamStore.readNameStore));
    for (unsigned i = 0; i < length(samStore.readNameStore); ++i)
        SEQAN_ASSERT(samStore.readNameStore[i] == bamStore.readNameStore[i]);
    SEQAN_ASSERT_EQ(length(samStore.readSeqStore), length(bamStore.readSeqStore));
    for (unsigned i = 0; i < length(samStore.readSeqStore); ++i)
        SEQAN_ASSERT(samStore.readSeqStore[i] == bamStore.readSeqStore[i]);
    SEQAN_ASSERT_EQ(length(samStore.readStore), length(bamStore.readStore));
    for (unsigned i = 0; i < length(samStore.readStore); ++i)
        SEQAN_ASSERT(samStore.readStore[i] ==  bamStore.readStore[i]);

    // TODO(holtgrew): Actually check contents.
    // */
}

SEQAN_DEFINE_TEST(test_store_io_bam_write)
{
    // Construct name to reference FASTA files.
    char fastaBuffer[1023];
    strcpy(fastaBuffer, getAbsolutePath(""));
    strcat(fastaBuffer, "/projects/tests/store/toy.fa");
    // strcat(fastaBuffer, "/projects/tests/store/ex1.fa");
    // Construct file name to SAM file.
    char samBuffer[1023];
    strcpy(samBuffer, getAbsolutePath(""));
    strcat(samBuffer, "/projects/tests/store/toy.sam");
    // strcat(samBuffer, "/projects/tests/store/ex1.sam");
    // Construct path to a temporary output file.
    char tmpBuffer[1023];
    strcpy(tmpBuffer, SEQAN_TEMP_FILENAME());

    // Load FragmentStore from SAM and BAM file.
    FragmentStore<> samStore;
    loadContigs(samStore, fastaBuffer);
    FILE * samFp = fopen(samBuffer, "rb");
    read(samFp, samStore, Sam());
    fclose(samFp);

    // Write out BAM file.
    // write(tmpBuffer, samStore, Bam());
    FILE * outFb = fopen(tmpBuffer, "wb");
    write(outFb, samStore, Sam());
    fclose(outFb);
}
