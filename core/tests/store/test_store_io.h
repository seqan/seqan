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
// Tests for the SeqAn model store, I/O functionality.
// ==========================================================================

#include <seqan/basic.h>  // For test functionality.
#include <seqan/store.h>  // Header under test.

#include <seqan/misc/misc_svg.h>

using namespace seqan;

SEQAN_DEFINE_TEST(test_store_io_read_ucsc_known_genes)
{
    // The file contains 13 annotations in total which will be checked line
    // after line.
    seqan::CharString ucscPath = SEQAN_PATH_TO_ROOT();
    append(ucscPath, "/core/tests/store/example_known_genes.tsv");

    std::fstream f(toCString(ucscPath), std::ios::binary | std::ios::in);
    SEQAN_ASSERT(f.good());

    seqan::FragmentStore<> store;

    read(f, store, seqan::Ucsc());
    f.close();

    Iterator<FragmentStore<>, AnnotationTree<> >::Type it;
    it = begin(store, AnnotationTree<>());

    SEQAN_ASSERT_EQ(getType(it), "<root>");
    SEQAN_ASSERT_EQ(getAnnotation(it).beginPos, maxValue(getAnnotation(it).beginPos));
    SEQAN_ASSERT_EQ(getAnnotation(it).endPos, maxValue(getAnnotation(it).endPos));
    SEQAN_ASSERT_EQ(value(it), 0u);
    SEQAN_ASSERT_EQ(getAnnotation(it).parentId, maxValue(getAnnotation(it).parentId));
    SEQAN_ASSERT_EQ(getParentName(it), "<root>");
    goNext(it);

    SEQAN_ASSERT_EQ(getType(it), "mRNA");
    SEQAN_ASSERT_EQ(getAnnotation(it).beginPos, 33031813);
    SEQAN_ASSERT_EQ(getAnnotation(it).endPos, 33026870);
    SEQAN_ASSERT_EQ(value(it), 1u);
    SEQAN_ASSERT_EQ(getAnnotation(it).parentId, 0u);
    SEQAN_ASSERT_EQ(getParentName(it), "<root>");
    goNext(it);

    SEQAN_ASSERT_EQ(getType(it), "CDS");
    SEQAN_ASSERT_EQ(getAnnotation(it).beginPos, 33026870);
    SEQAN_ASSERT_EQ(getAnnotation(it).endPos, 33026870);
    SEQAN_ASSERT_EQ(value(it), 2u);
    SEQAN_ASSERT_EQ(getAnnotation(it).parentId, 1u);
    SEQAN_ASSERT_EQ(getParentName(it), "uc002yoz.1");
    goNext(it);

    SEQAN_ASSERT_EQ(getType(it), "exon");
    SEQAN_ASSERT_EQ(getAnnotation(it).beginPos, 33027740);
    SEQAN_ASSERT_EQ(getAnnotation(it).endPos, 33026870);
    SEQAN_ASSERT_EQ(value(it), 3u);
    SEQAN_ASSERT_EQ(getAnnotation(it).parentId, 1u);
    SEQAN_ASSERT_EQ(getParentName(it), "uc002yoz.1");
    goNext(it);

    SEQAN_ASSERT_EQ(getType(it), "exon");
    SEQAN_ASSERT_EQ(getAnnotation(it).beginPos, 33030540);
    SEQAN_ASSERT_EQ(getAnnotation(it).endPos, 33030246);
    SEQAN_ASSERT_EQ(value(it), 4u);
    SEQAN_ASSERT_EQ(getAnnotation(it).parentId, 1u);
    SEQAN_ASSERT_EQ(getParentName(it), "uc002yoz.1");
    goNext(it);

    SEQAN_ASSERT_EQ(getType(it), "exon");
    SEQAN_ASSERT_EQ(getAnnotation(it).beginPos, 33031813);
    SEQAN_ASSERT_EQ(getAnnotation(it).endPos, 33031709);
    SEQAN_ASSERT_EQ(value(it), 5u);
    SEQAN_ASSERT_EQ(getAnnotation(it).parentId, 1u);
    SEQAN_ASSERT_EQ(getParentName(it), "uc002yoz.1");
    goNext(it);

    SEQAN_ASSERT_EQ(getType(it), "mRNA");
    SEQAN_ASSERT_EQ(getAnnotation(it).beginPos, 33031934);
    SEQAN_ASSERT_EQ(getAnnotation(it).endPos, 33041243);
    SEQAN_ASSERT_EQ(value(it), 6u);
    SEQAN_ASSERT_EQ(getAnnotation(it).parentId, 0u);
    SEQAN_ASSERT_EQ(getParentName(it), "<root>");
    goNext(it);

    SEQAN_ASSERT_EQ(getType(it), "CDS");
    SEQAN_ASSERT_EQ(getAnnotation(it).beginPos, 33032082);
    SEQAN_ASSERT_EQ(getAnnotation(it).endPos, 33040891);
    SEQAN_ASSERT_EQ(value(it), 7u);
    SEQAN_ASSERT_EQ(getAnnotation(it).parentId, 6u);
    SEQAN_ASSERT_EQ(getParentName(it), "uc002ypa.3");
    goNext(it);

    SEQAN_ASSERT_EQ(getType(it), "exon");
    SEQAN_ASSERT_EQ(getAnnotation(it).beginPos, 33031934);
    SEQAN_ASSERT_EQ(getAnnotation(it).endPos, 33032154);
    SEQAN_ASSERT_EQ(value(it), 8u);
    SEQAN_ASSERT_EQ(getAnnotation(it).parentId, 6u);
    SEQAN_ASSERT_EQ(getParentName(it), "uc002ypa.3");
    goNext(it);

    SEQAN_ASSERT_EQ(getType(it), "exon");
    SEQAN_ASSERT_EQ(getAnnotation(it).beginPos, 33036102);
    SEQAN_ASSERT_EQ(getAnnotation(it).endPos, 33036199);
    SEQAN_ASSERT_EQ(value(it), 9u);
    SEQAN_ASSERT_EQ(getAnnotation(it).parentId, 6u);
    SEQAN_ASSERT_EQ(getParentName(it), "uc002ypa.3");
    goNext(it);

    SEQAN_ASSERT_EQ(getType(it), "exon");
    SEQAN_ASSERT_EQ(getAnnotation(it).beginPos, 33038761);
    SEQAN_ASSERT_EQ(getAnnotation(it).endPos, 33038831);
    SEQAN_ASSERT_EQ(value(it), 10u);
    SEQAN_ASSERT_EQ(getAnnotation(it).parentId, 6u);
    SEQAN_ASSERT_EQ(getParentName(it), "uc002ypa.3");
    goNext(it);

    SEQAN_ASSERT_EQ(getType(it), "exon");
    SEQAN_ASSERT_EQ(getAnnotation(it).beginPos, 33039570);
    SEQAN_ASSERT_EQ(getAnnotation(it).endPos, 33039688);
    SEQAN_ASSERT_EQ(value(it), 11u);
    SEQAN_ASSERT_EQ(getAnnotation(it).parentId, 6u);
    SEQAN_ASSERT_EQ(getParentName(it), "uc002ypa.3");
    goNext(it);

    SEQAN_ASSERT_EQ(getType(it), "exon");
    SEQAN_ASSERT_EQ(getAnnotation(it).beginPos, 33040783);
    SEQAN_ASSERT_EQ(getAnnotation(it).endPos, 33041243);
    SEQAN_ASSERT_EQ(value(it), 12u);
    SEQAN_ASSERT_EQ(getAnnotation(it).parentId, 6u);
    SEQAN_ASSERT_EQ(getParentName(it), "uc002ypa.3");
    goNext(it);

    SEQAN_ASSERT(atEnd(it));
}

SEQAN_DEFINE_TEST(test_store_io_write_ucsc_known_genes)
{
    seqan::CharString ucscPath = SEQAN_PATH_TO_ROOT();
    append(ucscPath, "/core/tests/store/example_known_genes.tsv");

    std::fstream fin(toCString(ucscPath), std::ios::binary | std::ios::in);
    SEQAN_ASSERT(fin.good());

    seqan::FragmentStore<> store;
    read(fin, store, seqan::Ucsc());

    seqan::CharString outPath  = SEQAN_TEMP_FILENAME();
    append(outPath, ".tsv");
    std::fstream fout(toCString(outPath), std::ios::binary | std::ios::out);
    write(fout, store, seqan::Ucsc());
    fout.close();

    seqan::CharString goldPath = SEQAN_PATH_TO_ROOT();
    append(goldPath, "/core/tests/store/example_known_genes.tsv");

    SEQAN_ASSERT(seqan::_compareTextFiles(toCString(outPath), toCString(goldPath)));
}

SEQAN_DEFINE_TEST(test_store_io_read_gff)
{
    seqan::CharString gffPath = SEQAN_PATH_TO_ROOT();
    append(gffPath, "/core/tests/store/example_gff.tsv");

    std::fstream f(toCString(gffPath), std::ios::binary | std::ios::in);
    SEQAN_ASSERT(f.good());

    seqan::FragmentStore<> store;

    read(f, store, seqan::Gff());
    f.close();

    Iterator<FragmentStore<>, AnnotationTree<> >::Type it;
    it = begin(store, AnnotationTree<>());

    SEQAN_ASSERT_EQ(getType(it), "<root>");
    SEQAN_ASSERT_EQ(getAnnotation(it).beginPos, maxValue(getAnnotation(it).beginPos));
    SEQAN_ASSERT_EQ(getAnnotation(it).endPos, maxValue(getAnnotation(it).endPos));
    SEQAN_ASSERT_EQ(value(it), 0u);
    SEQAN_ASSERT_EQ(getAnnotation(it).parentId, 4294967295);
    SEQAN_ASSERT_EQ(getParentName(it), "<root>");
    goNext(it);

    SEQAN_ASSERT_EQ(getType(it), "mRNA");
    SEQAN_ASSERT_EQ(getAnnotation(it).beginPos, 1299u);
    SEQAN_ASSERT_EQ(getAnnotation(it).endPos, 9000u);
    SEQAN_ASSERT_EQ(value(it), 1u);
    SEQAN_ASSERT_EQ(getAnnotation(it).parentId, 0u);
    SEQAN_ASSERT_EQ(getParentName(it), "<root>");
    goNext(it);

    SEQAN_ASSERT_EQ(getType(it), "exon");
    SEQAN_ASSERT_EQ(getAnnotation(it).beginPos, 1299u);
    SEQAN_ASSERT_EQ(getAnnotation(it).endPos, 1500);
    SEQAN_ASSERT_EQ(value(it), 2u);
    SEQAN_ASSERT_EQ(getAnnotation(it).parentId, 1);
    SEQAN_ASSERT_EQ(getParentName(it), "mrna0001");
    goNext(it);

    SEQAN_ASSERT_EQ(getType(it), "exon");
    SEQAN_ASSERT_EQ(getAnnotation(it).beginPos, 1049u);
    SEQAN_ASSERT_EQ(getAnnotation(it).endPos, 1500u);
    SEQAN_ASSERT_EQ(value(it), 3u);
    SEQAN_ASSERT_EQ(getAnnotation(it).parentId, 1);
    SEQAN_ASSERT_EQ(getParentName(it), "mrna0001");
}

SEQAN_DEFINE_TEST(test_store_io_write_gff)
{
    seqan::CharString ucscPath = SEQAN_PATH_TO_ROOT();
    append(ucscPath, "/core/tests/store/example_gff.tsv");

    std::fstream fin(toCString(ucscPath), std::ios::binary | std::ios::in);
    SEQAN_ASSERT(fin.good());

    seqan::FragmentStore<> store;
    read(fin, store, seqan::Gff());

    seqan::CharString outPath  = SEQAN_TEMP_FILENAME();
    append(outPath, ".tsv");
    std::cerr << outPath << std::endl;
    std::fstream fout(toCString(outPath), std::ios::binary | std::ios::out);
    write(fout, store, seqan::Gff());
    fout.close();

    seqan::CharString goldPath = SEQAN_PATH_TO_ROOT();
    append(goldPath, "/core/tests/store/example_gff.tsv");

    SEQAN_ASSERT(seqan::_compareTextFiles(toCString(outPath), toCString(goldPath)));
}

SEQAN_DEFINE_TEST(test_store_io_read_gtf)
{
    seqan::CharString gffPath = SEQAN_PATH_TO_ROOT();
    append(gffPath, "/core/tests/store/example_gtf.tsv");

    std::fstream f(toCString(gffPath), std::ios::binary | std::ios::in);
    SEQAN_ASSERT(f.good());

    seqan::FragmentStore<> store;

    read(f, store, seqan::Gtf());
    f.close();

    Iterator<FragmentStore<>, AnnotationTree<> >::Type it;
    it = begin(store, AnnotationTree<>());

    SEQAN_ASSERT_EQ(getType(it), "<root>");
    SEQAN_ASSERT_EQ(getAnnotation(it).beginPos, maxValue(getAnnotation(it).beginPos));
    SEQAN_ASSERT_EQ(getAnnotation(it).endPos, maxValue(getAnnotation(it).endPos));
    SEQAN_ASSERT_EQ(value(it), 0u);
    SEQAN_ASSERT_EQ(getAnnotation(it).parentId, 4294967295);
    SEQAN_ASSERT_EQ(getParentName(it), "<root>");
    goNext(it);

    SEQAN_ASSERT_EQ(getType(it), "gene");
    SEQAN_ASSERT_EQ(getAnnotation(it).beginPos, 13182);
    SEQAN_ASSERT_EQ(getAnnotation(it).endPos, 5140);
    SEQAN_ASSERT_EQ(value(it), 1u);
    SEQAN_ASSERT_EQ(getAnnotation(it).parentId, 0);
    SEQAN_ASSERT_EQ(getParentName(it), "<root>");
    goNext(it);

    SEQAN_ASSERT_EQ(getType(it), "mRNA");
    SEQAN_ASSERT_EQ(getAnnotation(it).beginPos, 13182);
    SEQAN_ASSERT_EQ(getAnnotation(it).endPos, 5140);
    SEQAN_ASSERT_EQ(value(it), 2u);
    SEQAN_ASSERT_EQ(getAnnotation(it).parentId, 1);
    SEQAN_ASSERT_EQ(getParentName(it), "1");
    goNext(it);

    SEQAN_ASSERT_EQ(getType(it), "inter");
    SEQAN_ASSERT_EQ(getAnnotation(it).beginPos, 8522);
    SEQAN_ASSERT_EQ(getAnnotation(it).endPos, 5140);
    SEQAN_ASSERT_EQ(value(it), 3u);
    SEQAN_ASSERT_EQ(getAnnotation(it).parentId, 2);
    SEQAN_ASSERT_EQ(getParentName(it), "2");
    goNext(it);

    SEQAN_ASSERT_EQ(getType(it), "inter_CNS");
    SEQAN_ASSERT_EQ(getAnnotation(it).beginPos, 9711);
    SEQAN_ASSERT_EQ(getAnnotation(it).endPos, 8522);
    SEQAN_ASSERT_EQ(value(it), 4u);
    SEQAN_ASSERT_EQ(getAnnotation(it).parentId, 2);
    SEQAN_ASSERT_EQ(getParentName(it), "2");
    goNext(it);

    SEQAN_ASSERT_EQ(getType(it), "inter");
    SEQAN_ASSERT_EQ(getAnnotation(it).beginPos, 13182);
    SEQAN_ASSERT_EQ(getAnnotation(it).endPos, 9711);
    SEQAN_ASSERT_EQ(value(it), 5u);
    SEQAN_ASSERT_EQ(getAnnotation(it).parentId, 2);
    SEQAN_ASSERT_EQ(getParentName(it), "2");
    goNext(it);

    SEQAN_ASSERT_EQ(getType(it), "gene");
    SEQAN_ASSERT_EQ(getAnnotation(it).beginPos, 70151);
    SEQAN_ASSERT_EQ(getAnnotation(it).endPos, 65148);
    SEQAN_ASSERT_EQ(value(it), 6u);
    SEQAN_ASSERT_EQ(getAnnotation(it).parentId, 0);
    SEQAN_ASSERT_EQ(getParentName(it), "<root>");
    goNext(it);

    SEQAN_ASSERT_EQ(getType(it), "mRNA");
    SEQAN_ASSERT_EQ(getAnnotation(it).beginPos, 70151);
    SEQAN_ASSERT_EQ(getAnnotation(it).endPos, 65148);
    SEQAN_ASSERT_EQ(value(it), 7u);
    SEQAN_ASSERT_EQ(getAnnotation(it).parentId, 6);
    SEQAN_ASSERT_EQ(getParentName(it), "140.000");
    goNext(it);

    SEQAN_ASSERT_EQ(getType(it), "3UTR");
    SEQAN_ASSERT_EQ(getAnnotation(it).beginPos, 65487);
    SEQAN_ASSERT_EQ(getAnnotation(it).endPos, 65148);
    SEQAN_ASSERT_EQ(value(it), 8u);
    SEQAN_ASSERT_EQ(getAnnotation(it).parentId, 7);
    SEQAN_ASSERT_EQ(getParentName(it), "140.000.1");
    goNext(it);

    SEQAN_ASSERT_EQ(getType(it), "3UTR");
    SEQAN_ASSERT_EQ(getAnnotation(it).beginPos, 66992);
    SEQAN_ASSERT_EQ(getAnnotation(it).endPos, 66822);
    SEQAN_ASSERT_EQ(value(it), 9u);
    SEQAN_ASSERT_EQ(getAnnotation(it).parentId, 7);
    SEQAN_ASSERT_EQ(getParentName(it), "140.000.1");
    goNext(it);

    SEQAN_ASSERT_EQ(getType(it), "stop_codon");
    SEQAN_ASSERT_EQ(getAnnotation(it).beginPos, 66995);
    SEQAN_ASSERT_EQ(getAnnotation(it).endPos, 66992);
    SEQAN_ASSERT_EQ(value(it), 10u);
    SEQAN_ASSERT_EQ(getAnnotation(it).parentId, 7);
    SEQAN_ASSERT_EQ(getParentName(it), "140.000.1");
    goNext(it);

    SEQAN_ASSERT_EQ(getType(it), "CDS");
    SEQAN_ASSERT_EQ(getAnnotation(it).beginPos, 66999);
    SEQAN_ASSERT_EQ(getAnnotation(it).endPos, 66995);
    SEQAN_ASSERT_EQ(value(it), 11u);
    SEQAN_ASSERT_EQ(getAnnotation(it).parentId, 7);
    SEQAN_ASSERT_EQ(getParentName(it), "140.000.1");
    goNext(it);

    SEQAN_ASSERT_EQ(getType(it), "intron_CNS");
    SEQAN_ASSERT_EQ(getAnnotation(it).beginPos, 70151);
    SEQAN_ASSERT_EQ(getAnnotation(it).endPos, 70102);
    SEQAN_ASSERT_EQ(value(it), 12u);
    SEQAN_ASSERT_EQ(getAnnotation(it).parentId, 7);
    SEQAN_ASSERT_EQ(getParentName(it), "140.000.1");
}

SEQAN_DEFINE_TEST(test_store_io_write_gtf)
{
    seqan::CharString ucscPath = SEQAN_PATH_TO_ROOT();
    append(ucscPath, "/core/tests/store/example_gtf.tsv");

    std::fstream fin(toCString(ucscPath), std::ios::binary | std::ios::in);
    SEQAN_ASSERT(fin.good());

    seqan::FragmentStore<> store;
    read(fin, store, seqan::Gtf());

    seqan::CharString outPath  = SEQAN_TEMP_FILENAME();
    append(outPath, ".tsv");
    std::cerr << outPath << std::endl;
    std::fstream fout(toCString(outPath), std::ios::binary | std::ios::out);
    write(fout, store, seqan::Gtf());
    fout.close();

    seqan::CharString goldPath = SEQAN_PATH_TO_ROOT();
    append(goldPath, "/core/tests/store/example_gtf.tsv");

    SEQAN_ASSERT(seqan::_compareTextFiles(toCString(outPath), toCString(goldPath)));
}



SEQAN_DEFINE_TEST(test_store_io_sam)
{
    FragmentStore<> store;
    char buffer[1023];

    // 1. LOAD CONTIGS
    strcpy(buffer, SEQAN_PATH_TO_ROOT());
    strcat(buffer, "/core/tests/store/ex1.fa");

    loadContigs(store, buffer);

    // 2. LOAD SAM ALIGNMENTS
    strcpy(buffer, SEQAN_PATH_TO_ROOT());
    strcat(buffer, "/core/tests/store/ex1.sam.copy");
    MultiSeqFile sam1;
    open(sam1.concat, buffer);
    split(sam1, Raw());

    {
        // read reference Sam from file
        std::ifstream samFile(buffer);
        SEQAN_ASSERT(samFile);
        read(samFile, store, Sam());
    }

    //AlignedReadLayout layout;
    //layoutAlignment(layout, store);
    //printAlignment(std::cout, Raw(), layout, store, 0, 0, 1000, 0, 1000);

    // 3. WRITE SAM ALIGNMENTS
    strcpy(buffer, SEQAN_TEMP_FILENAME());
    {
        // write Sam to temp file
        std::ofstream samFileOut(buffer);
        SEQAN_ASSERT(samFileOut);
        write(samFileOut, store, Sam());
    }

    // 4. COMPARE BOTH SAM FILES
    MultiSeqFile sam2;
    open(sam2.concat, buffer);
    split(sam2, Raw());

    SEQAN_ASSERT(!empty(sam1));
    SEQAN_ASSERT(!empty(sam2));
    for (unsigned i = 0; i < length(sam1); ++i)
    {
        if (sam1[i] != sam2[i])
        {
            std::cout << "    \t" << sam1[i] << std::endl;
            std::cout << " != \t" << sam2[i] << std::endl;
            SEQAN_ASSERT_FAIL("Files differ in line %d.", i);
        }
    }
}
