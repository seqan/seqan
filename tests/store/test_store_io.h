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

#include <seqan/misc/svg.h>

using namespace seqan;

SEQAN_DEFINE_TEST(test_store_io_read_ucsc_known_genes)
{
    // The file contains 13 annotations in total which will be checked line
    // after line.
    seqan::CharString ucscPath = seqan::getAbsolutePath("/tests/store/example_knownGene.txt");

    UcscFileIn fin(toCString(ucscPath));
    seqan::FragmentStore<> store;
    readRecords(store, fin);

    Iterator<FragmentStore<>, AnnotationTree<> >::Type it;
    it = begin(store, AnnotationTree<>());

    SEQAN_ASSERT_EQ(getType(it), "<root>");
    SEQAN_ASSERT_EQ(getAnnotation(it).beginPos, std::numeric_limits<decltype(getAnnotation(it).beginPos)>::max());
    SEQAN_ASSERT_EQ(getAnnotation(it).endPos, std::numeric_limits<decltype(getAnnotation(it).endPos)>::max());
    SEQAN_ASSERT_EQ(value(it), 0u);
    SEQAN_ASSERT_EQ(getAnnotation(it).parentId, std::numeric_limits<decltype(getAnnotation(it).parentId)>::max());
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

SEQAN_DEFINE_TEST(test_store_io_read_ucsc_known_genes_and_isoforms)
{
    // The file contains 13 annotations in total which will be checked line
    // after line.
    seqan::CharString ucscGenesPath = getAbsolutePath("/tests/store/example_knownGene.txt");
    seqan::CharString ucscIsoformsPath = getAbsolutePath("/tests/store/example_knownIsoforms.txt");

    UcscFileIn finGenes(toCString(ucscGenesPath));
    UcscFileIn finIsoforms(toCString(ucscIsoformsPath));
    seqan::FragmentStore<> store;
    readRecords(store, finGenes);
    readRecords(store, finIsoforms);

    Iterator<FragmentStore<>, AnnotationTree<> >::Type it;
    it = begin(store, AnnotationTree<>());

    SEQAN_ASSERT_EQ(getType(it), "<root>");
    SEQAN_ASSERT_EQ(getAnnotation(it).beginPos, std::numeric_limits<decltype(getAnnotation(it).beginPos)>::max());
    SEQAN_ASSERT_EQ(getAnnotation(it).endPos, std::numeric_limits<decltype(getAnnotation(it).endPos)>::max());
    SEQAN_ASSERT_EQ(value(it), 0u);
    SEQAN_ASSERT_EQ(getAnnotation(it).parentId, std::numeric_limits<decltype(getAnnotation(it).parentId)>::max());
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

    SEQAN_ASSERT_EQ(getType(it), "gene");
    SEQAN_ASSERT_EQ(value(it), 13u);
    SEQAN_ASSERT_EQ(getAnnotation(it).parentId, 0u);
    SEQAN_ASSERT_EQ(getParentName(it), "<root>");
    goNext(it);

    SEQAN_ASSERT_EQ(getType(it), "mRNA");
    SEQAN_ASSERT_EQ(value(it), 14u);
    SEQAN_ASSERT_EQ(getAnnotation(it).parentId, 13u);
    SEQAN_ASSERT_EQ(getParentName(it), "GENE1");
    goNext(it);

    SEQAN_ASSERT_EQ(getType(it), "gene");
    SEQAN_ASSERT_EQ(value(it), 15u);
    SEQAN_ASSERT_EQ(getAnnotation(it).parentId, 0u);
    SEQAN_ASSERT_EQ(getParentName(it), "<root>");
    goNext(it);

    SEQAN_ASSERT_EQ(getType(it), "mRNA");
    SEQAN_ASSERT_EQ(value(it), 16u);
    SEQAN_ASSERT_EQ(getAnnotation(it).parentId, 15u);
    SEQAN_ASSERT_EQ(getParentName(it), "GENE2");
    goNext(it);

    SEQAN_ASSERT_EQ(getType(it), "gene");
    SEQAN_ASSERT_EQ(value(it), 17u);
    SEQAN_ASSERT_EQ(getAnnotation(it).parentId, 0u);
    SEQAN_ASSERT_EQ(getParentName(it), "<root>");
    goNext(it);

    SEQAN_ASSERT_EQ(getType(it), "mRNA");
    SEQAN_ASSERT_EQ(value(it), 18u);
    SEQAN_ASSERT_EQ(getAnnotation(it).parentId, 17u);
    SEQAN_ASSERT_EQ(getParentName(it), "GENE3");
    goNext(it);

    SEQAN_ASSERT(atEnd(it));
}

SEQAN_DEFINE_TEST(test_store_io_write_ucsc_known_genes)
{
    seqan::CharString ucscPath = getAbsolutePath("/tests/store/example_knownGene.txt");

    UcscFileIn fin(toCString(ucscPath));
    seqan::FragmentStore<> store;
    readRecords(store, fin);

    seqan::CharString outPath  = SEQAN_TEMP_FILENAME();
    append(outPath, ".knownGene.txt");
    UcscFileOut fout(toCString(outPath));
    writeRecords(fout, store);
    close(fout);

    seqan::CharString goldPath = getAbsolutePath("/tests/store/example_knownGene.txt");

    SEQAN_ASSERT(seqan::_compareTextFilesAlt(toCString(outPath), toCString(goldPath)));
}

SEQAN_DEFINE_TEST(test_store_io_read_gff)
{
    seqan::CharString gffPath = getAbsolutePath("/tests/store/example.gff");

    GffFileIn f(toCString(gffPath));
    typedef typename seqan::FragmentStore<>::TAnnotationStoreElement::TId TId;

    seqan::FragmentStore<> store;

    readRecords(store, f);

    Iterator<FragmentStore<>, AnnotationTree<> >::Type it;
    it = begin(store, AnnotationTree<>());

    SEQAN_ASSERT_EQ(getType(it), "<root>");
    SEQAN_ASSERT_EQ(getAnnotation(it).beginPos, std::numeric_limits<decltype(getAnnotation(it).beginPos)>::max());
    SEQAN_ASSERT_EQ(getAnnotation(it).endPos, std::numeric_limits<decltype(getAnnotation(it).endPos)>::max());
    SEQAN_ASSERT_EQ(value(it), 0u);
    SEQAN_ASSERT_EQ(getAnnotation(it).parentId, std::numeric_limits<TId>::max());
    SEQAN_ASSERT_EQ(getParentName(it), "<root>");
    goNext(it);

    SEQAN_ASSERT_EQ(getType(it), "mRNA");
    SEQAN_ASSERT_EQ(getAnnotation(it).beginPos, 1299);
    SEQAN_ASSERT_EQ(getAnnotation(it).endPos, 9000);
    SEQAN_ASSERT_EQ(value(it), 1u);
    SEQAN_ASSERT_EQ(getAnnotation(it).parentId, 0u);
    SEQAN_ASSERT_EQ(getParentName(it), "<root>");
    goNext(it);

    SEQAN_ASSERT_EQ(getType(it), "exon");
    SEQAN_ASSERT_EQ(getAnnotation(it).beginPos, 1299);
    SEQAN_ASSERT_EQ(getAnnotation(it).endPos, 1500);
    SEQAN_ASSERT_EQ(value(it), 2u);
    SEQAN_ASSERT_EQ(getAnnotation(it).parentId, 1u);
    SEQAN_ASSERT_EQ(getParentName(it), "mrna0001");
    goNext(it);

    SEQAN_ASSERT_EQ(getType(it), "exon");
    SEQAN_ASSERT_EQ(getAnnotation(it).beginPos, 1049);
    SEQAN_ASSERT_EQ(getAnnotation(it).endPos, 1500    );
    SEQAN_ASSERT_EQ(value(it), 3u);
    SEQAN_ASSERT_EQ(getAnnotation(it).parentId, 1u);
    SEQAN_ASSERT_EQ(getParentName(it), "mrna0001");
}

SEQAN_DEFINE_TEST(test_store_io_write_gff)
{
    seqan::CharString goldPath = getAbsolutePath("/tests/store/example.gff");

    GffFileIn fin(toCString(goldPath));
    seqan::FragmentStore<> store;
    readRecords(store, fin);

    seqan::CharString outPath = SEQAN_TEMP_FILENAME();
    append(outPath, ".gff");
    GffFileOut fout(toCString(outPath));
    writeRecords(fout, store);
    close(fout);

    SEQAN_ASSERT(seqan::_compareTextFilesAlt(toCString(outPath), toCString(goldPath)));
}

SEQAN_DEFINE_TEST(test_store_io_read_gtf)
{
    typedef typename seqan::FragmentStore<>::TAnnotationStoreElement::TId TId;

    seqan::CharString gtfPath = getAbsolutePath("/tests/store/example.gtf");

    GffFileIn fin(toCString(gtfPath));
    seqan::FragmentStore<> store;
    readRecords(store, fin);

    Iterator<FragmentStore<>, AnnotationTree<> >::Type it;
    it = begin(store, AnnotationTree<>());

    SEQAN_ASSERT_EQ(getType(it), "<root>");
    SEQAN_ASSERT_EQ(getAnnotation(it).beginPos, std::numeric_limits<decltype(getAnnotation(it).beginPos)>::max());
    SEQAN_ASSERT_EQ(getAnnotation(it).endPos, std::numeric_limits<decltype(getAnnotation(it).endPos)>::max());
    SEQAN_ASSERT_EQ(value(it), 0u);
    SEQAN_ASSERT_EQ(getAnnotation(it).parentId, std::numeric_limits<TId>::max());
    SEQAN_ASSERT_EQ(getParentName(it), "<root>");
    goNext(it);

    SEQAN_ASSERT_EQ(getType(it), "gene");
    SEQAN_ASSERT_EQ(getAnnotation(it).beginPos, 13182);
    SEQAN_ASSERT_EQ(getAnnotation(it).endPos, 5140);
    SEQAN_ASSERT_EQ(value(it), 1u);
    SEQAN_ASSERT_EQ(getAnnotation(it).parentId, 0u);
    SEQAN_ASSERT_EQ(getParentName(it), "<root>");
    goNext(it);

    SEQAN_ASSERT_EQ(getType(it), "mRNA");
    SEQAN_ASSERT_EQ(getAnnotation(it).beginPos, 13182);
    SEQAN_ASSERT_EQ(getAnnotation(it).endPos, 5140);
    SEQAN_ASSERT_EQ(value(it), 2u);
    SEQAN_ASSERT_EQ(getAnnotation(it).parentId, 1u);
    SEQAN_ASSERT_EQ(getParentName(it), "gene1");
    goNext(it);

    SEQAN_ASSERT_EQ(getType(it), "inter");
    SEQAN_ASSERT_EQ(getAnnotation(it).beginPos, 8522);
    SEQAN_ASSERT_EQ(getAnnotation(it).endPos, 5140);
    SEQAN_ASSERT_EQ(value(it), 3u);
    SEQAN_ASSERT_EQ(getAnnotation(it).parentId, 2u);
    SEQAN_ASSERT_EQ(getParentName(it), "trans2");
    goNext(it);

    SEQAN_ASSERT_EQ(getType(it), "inter_CNS");
    SEQAN_ASSERT_EQ(getAnnotation(it).beginPos, 9711);
    SEQAN_ASSERT_EQ(getAnnotation(it).endPos, 8522);
    SEQAN_ASSERT_EQ(value(it), 4u);
    SEQAN_ASSERT_EQ(getAnnotation(it).parentId, 2u);
    SEQAN_ASSERT_EQ(getParentName(it), "trans2");
    goNext(it);

    SEQAN_ASSERT_EQ(getType(it), "inter");
    SEQAN_ASSERT_EQ(getAnnotation(it).beginPos, 13182);
    SEQAN_ASSERT_EQ(getAnnotation(it).endPos, 9711);
    SEQAN_ASSERT_EQ(value(it), 5u);
    SEQAN_ASSERT_EQ(getAnnotation(it).parentId, 2u);
    SEQAN_ASSERT_EQ(getParentName(it), "trans2");
    goNext(it);

    SEQAN_ASSERT_EQ(getType(it), "gene");
    SEQAN_ASSERT_EQ(getAnnotation(it).beginPos, 70151);
    SEQAN_ASSERT_EQ(getAnnotation(it).endPos, 65148);
    SEQAN_ASSERT_EQ(value(it), 6u);
    SEQAN_ASSERT_EQ(getAnnotation(it).parentId, 0u);
    SEQAN_ASSERT_EQ(getParentName(it), "<root>");
    goNext(it);

    SEQAN_ASSERT_EQ(getType(it), "mRNA");
    SEQAN_ASSERT_EQ(getAnnotation(it).beginPos, 70151);
    SEQAN_ASSERT_EQ(getAnnotation(it).endPos, 65148);
    SEQAN_ASSERT_EQ(value(it), 7u);
    SEQAN_ASSERT_EQ(getAnnotation(it).parentId, 6u);
    SEQAN_ASSERT_EQ(getParentName(it), "140.000");
    goNext(it);

    SEQAN_ASSERT_EQ(getType(it), "3UTR");
    SEQAN_ASSERT_EQ(getAnnotation(it).beginPos, 65487);
    SEQAN_ASSERT_EQ(getAnnotation(it).endPos, 65148);
    SEQAN_ASSERT_EQ(value(it), 8u);
    SEQAN_ASSERT_EQ(getAnnotation(it).parentId, 7u);
    SEQAN_ASSERT_EQ(getParentName(it), "140.000.1");
    goNext(it);

    SEQAN_ASSERT_EQ(getType(it), "3UTR");
    SEQAN_ASSERT_EQ(getAnnotation(it).beginPos, 66992);
    SEQAN_ASSERT_EQ(getAnnotation(it).endPos, 66822);
    SEQAN_ASSERT_EQ(value(it), 9u);
    SEQAN_ASSERT_EQ(getAnnotation(it).parentId, 7u);
    SEQAN_ASSERT_EQ(getParentName(it), "140.000.1");
    goNext(it);

    SEQAN_ASSERT_EQ(getType(it), "stop_codon");
    SEQAN_ASSERT_EQ(getAnnotation(it).beginPos, 66995);
    SEQAN_ASSERT_EQ(getAnnotation(it).endPos, 66992);
    SEQAN_ASSERT_EQ(value(it), 10u);
    SEQAN_ASSERT_EQ(getAnnotation(it).parentId, 7u);
    SEQAN_ASSERT_EQ(getParentName(it), "140.000.1");
    goNext(it);

    SEQAN_ASSERT_EQ(getType(it), "CDS");
    SEQAN_ASSERT_EQ(getAnnotation(it).beginPos, 66999);
    SEQAN_ASSERT_EQ(getAnnotation(it).endPos, 66995);
    SEQAN_ASSERT_EQ(value(it), 11u);
    SEQAN_ASSERT_EQ(getAnnotation(it).parentId, 7u);
    SEQAN_ASSERT_EQ(getParentName(it), "140.000.1");
    goNext(it);

    SEQAN_ASSERT_EQ(getType(it), "intron_CNS");
    SEQAN_ASSERT_EQ(getAnnotation(it).beginPos, 70151);
    SEQAN_ASSERT_EQ(getAnnotation(it).endPos, 70102);
    SEQAN_ASSERT_EQ(value(it), 12u);
    SEQAN_ASSERT_EQ(getAnnotation(it).parentId, 7u);
    SEQAN_ASSERT_EQ(getParentName(it), "140.000.1");
}

SEQAN_DEFINE_TEST(test_store_io_write_gtf)
{
    seqan::CharString goldPath = getAbsolutePath("/tests/store/example.gtf");

    GffFileIn fin(toCString(goldPath));
    seqan::FragmentStore<> store;
    readRecords(store, fin);

    seqan::CharString outPath = SEQAN_TEMP_FILENAME();
    append(outPath, ".gtf");
    GffFileOut fout(toCString(outPath));
    writeRecords(fout, store);
    close(fout);

    SEQAN_ASSERT(seqan::_compareTextFilesAlt(toCString(outPath), toCString(goldPath)));
}

// Read in SAM file, write out SAM file.
SEQAN_DEFINE_TEST(test_store_io_sam)
{
    FragmentStore<> store;

    // 1. LOAD CONTIGS
    std::string goldPathRef = getAbsolutePath("/tests/store/ex1.fa");
    loadContigs(store, toCString(goldPathRef));

    // 2. LOAD SAM ALIGNMENTS
    std::string goldPathSam = getAbsolutePath("/tests/store/ex1.copy.sam");
    BamFileIn inFile(toCString(goldPathSam));
    readRecords(store, inFile);

    // 3. WRITE SAM ALIGNMENTS
    std::string testPathSam = (std::string)SEQAN_TEMP_FILENAME() + ".sam";
    BamFileOut outFile(toCString(testPathSam));
    writeRecords(outFile, store);
    close(outFile);

    SEQAN_ASSERT(seqan::_compareTextFilesAlt(toCString(goldPathSam), toCString(testPathSam)));
}

SEQAN_DEFINE_TEST(test_store_io_sam2)
{
    FragmentStore<> store;

    // 1. LOAD CONTIGS
    std::string goldPathRef = getAbsolutePath("/tests/store/ex1.fa");
    loadContigs(store, toCString(goldPathRef));

    // 2. LOAD SAM ALIGNMENTS
    std::string goldPathSam = getAbsolutePath("/tests/store/ex1.copy.sam");
    BamFileIn inFile(toCString(goldPathSam));
    readRecords(store, inFile);

    // 3. WRITE SAM ALIGNMENTS
    std::string testPathSam = (std::string)SEQAN_TEMP_FILENAME() + ".sam";
    BamFileOut outFile(toCString(testPathSam));
    writeRecords(outFile, store);
    close(outFile);

    SEQAN_ASSERT(seqan::_compareTextFilesAlt(toCString(goldPathSam), toCString(testPathSam)));
}

template <typename TFragStore>
void _appendReadAlignments(TFragStore &store, char const *fileName)
{
    using namespace seqan;

    std::string str = getAbsolutePath(fileName);
    BamFileIn inFile(str.c_str());
    readRecords(store, inFile);
  }

template <typename TFragStore>
void _writeStore(TFragStore &store, std::string const &outPath, char const *suffix)
{
    AlignedReadLayout layout;
    layoutAlignment(layout, store);

    std::string outPathTxt = outPath + suffix;
    std::ofstream file(toCString(outPathTxt));
    printAlignment(file, layout, store, 0, 0, 1030, 0, 36);
    printAlignment(file, layout, store, 1, 0, 1030, 0, 36);
    file.close();

    std::string goldFileName = (std::string)"/tests/store/ex1.splitmerge" + (std::string)suffix;
    std::string goldPathTxt = seqan::getAbsolutePath(goldFileName.c_str());
    SEQAN_ASSERT(seqan::_compareTextFilesAlt(toCString(outPathTxt), toCString(goldPathTxt)));
}

SEQAN_DEFINE_TEST(test_store_io_split_sam)
{
    using namespace seqan;

    FragmentStore<> store;

    // 1. LOAD CONTIGS
    std::string fastaFileName = seqan::getAbsolutePath("/tests/store/ex1.fa");
    loadContigs(store, toCString(fastaFileName));

    std::string outPath = (std::string)SEQAN_TEMP_FILENAME();

    _appendReadAlignments(store, "/tests/store/ex1_a1.sam");
    _writeStore(store, outPath, ".1.txt");
    _appendReadAlignments(store, "/tests/store/ex1_a2.sam");
    _writeStore(store, outPath, ".2.txt");
    _appendReadAlignments(store, "/tests/store/ex1_a3.sam");
    _writeStore(store, outPath, ".3.txt");
    _appendReadAlignments(store, "/tests/store/ex1_b.sam");
    _writeStore(store, outPath, ".4.txt");

    std::string outPathSam = outPath + ".sam";
    BamFileOut outFile(outPathSam.c_str());
    writeRecords(outFile, store);
    close(outFile);

    std::string goldPathSam = seqan::getAbsolutePath("/tests/store/ex1.splitmerge.sam");
    SEQAN_ASSERT(seqan::_compareTextFilesAlt(toCString(outPathSam), toCString(goldPathSam)));
}

#if SEQAN_HAS_ZLIB

// Read in BAM file, write out SAM file.
SEQAN_DEFINE_TEST(test_store_io_read_bam)
{
    FragmentStore<> store;

    // 1. LOAD CONTIGS
    std::string fastaFileName = seqan::getAbsolutePath("/tests/store/ex1.fa");
    loadContigs(store, toCString(fastaFileName));

    // 2. LOAD BAM ALIGNMENTS
    std::string bamFileName = seqan::getAbsolutePath("/tests/store/ex1.bam");

    // Read reference Sam from file.
    {
        BamFileIn inFile(toCString(bamFileName));
        readRecords(store, inFile);
    }

//    AlignedReadLayout layout;
//    layoutAlignment(layout, store);
//    printAlignment(std::cout, layout, store, 0, 0, 1000, 0, 1000);

    // 3. WRITE SAM ALIGNMENTS
    std::string outFileName = (std::string)SEQAN_TEMP_FILENAME() + ".sam";
    // Write Sam to temp file.
    BamFileOut outFile(toCString(outFileName));
    writeRecords(outFile, store);
    close(outFile);

    // 4. COMPARE BOTH SAM FILES
    CharString samFileName = seqan::getAbsolutePath("/tests/store/ex1.copy.sam");
    SEQAN_ASSERT(seqan::_compareTextFilesAlt(toCString(samFileName), toCString(outFileName)));
}

#endif  // #if SEQAN_HAS_ZLIB


// Read AMOS and check for some basic properties.  Then, write out as SAM and verify with the expected result.
SEQAN_DEFINE_TEST(test_store_io_read_amos)
{
    // Get path to input file.
    std::string inPath = seqan::getAbsolutePath("/tests/store/toy.amos");
    // Get path to temporary file.
    std::string outPathSam = (std::string)SEQAN_TEMP_FILENAME() + ".sam";
    std::string outPathFasta = (std::string)SEQAN_TEMP_FILENAME() + ".fa";

    // Read in AMOS.
    seqan::FragmentStore<> store;
    std::fstream fAmosIn(toCString(inPath), std::ios::binary | std::ios::in);
    SEQAN_ASSERT(fAmosIn.good());
    int res = read(store, fAmosIn, seqan::Amos());
    SEQAN_ASSERT_EQ(res, 0);

    // Write out contigs and SAM file.
    BamFileOut fSamOut(toCString(outPathSam));
    writeRecords(fSamOut, store);
    close(fSamOut);

    SeqFileOut fFastaOut(toCString(outPathFasta));
    writeContigs(fFastaOut, store);
    close(fFastaOut);

    // Compare result.
    SEQAN_ASSERT_EQ(length(store.contigNameStore), 2u);
    SEQAN_ASSERT_EQ(length(store.readNameStore), 11u);
    SEQAN_ASSERT_EQ(length(store.matePairNameStore), 1u);
    SEQAN_ASSERT_EQ(length(store.matePairStore), 1u);
    SEQAN_ASSERT_EQ(length(store.alignedReadStore), 12u);

    seqan::CharString goldPathSam = seqan::getAbsolutePath("/tests/store/amos_to_sam_result.sam");
    SEQAN_ASSERT(seqan::_compareTextFilesAlt(toCString(outPathSam), toCString(goldPathSam)));
    seqan::CharString goldPathFasta = seqan::getAbsolutePath("/tests/store/amos_to_sam_result.fasta");


    SEQAN_ASSERT(seqan::_compareTextFilesAlt(toCString(outPathFasta), toCString(goldPathFasta)));
}

// Read SAM and write out as AMOS.  The resulting AMOS file is compared to a gold standard file.
SEQAN_DEFINE_TEST(test_store_io_write_amos)
{
    // Get path to input files.
    std::string inPathSam = seqan::getAbsolutePath("/tests/store/ex1.copy.sam");
    std::string inPathFasta = seqan::getAbsolutePath("/tests/store/ex1.fa");
    // Get path to temporary file.
    std::string outPathAmos = SEQAN_TEMP_FILENAME();

    // Read in SAM and FASTA.
    seqan::FragmentStore<> store;
    loadContigs(store, toCString(inPathFasta));
    BamFileIn fSamIn(toCString(inPathSam));
    readRecords(store, fSamIn);

    // Write out AMOS file.
    std::fstream fAmosOut(toCString(outPathAmos), std::ios::binary | std::ios::out);
    write(fAmosOut, store, seqan::Amos());
    fAmosOut.close();

    // Compare result.
    seqan::CharString goldPathAmos = seqan::getAbsolutePath("/tests/store/sam_to_amos_result.amos");
    SEQAN_ASSERT(seqan::_compareTextFilesAlt(toCString(outPathAmos), toCString(goldPathAmos)));
}

// Read SAM and write out as AMOS.  The resulting AMOS file is compared to a gold standard file.
SEQAN_DEFINE_TEST(test_store_io_readwrite_amos)
{
    // Get path to input files.
    std::string goldPathAmos = seqan::getAbsolutePath("/tests/store/toy.amos");
    // Get path to temporary file.
    std::string outPathAmos = SEQAN_TEMP_FILENAME();

    // Read in AMOS.
    seqan::FragmentStore<> store;
    std::fstream fAmosIn(toCString(goldPathAmos), std::ios::binary | std::ios::in);
    SEQAN_ASSERT(fAmosIn.good());
    int res = read(store, fAmosIn, seqan::Amos());
    SEQAN_ASSERT_EQ(res, 0);

    // Write out AMOS file.
    std::fstream fAmosOut(toCString(outPathAmos), std::ios::binary | std::ios::out);
    write(fAmosOut, store, seqan::Amos());
    fAmosOut.close();

    // Compare result.
    SEQAN_ASSERT(seqan::_compareTextFilesAlt(toCString(outPathAmos), toCString(goldPathAmos)));
}
