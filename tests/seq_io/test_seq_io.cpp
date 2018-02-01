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

#include <seqan/basic.h>
#include <seqan/stream.h>

#include <seqan/seq_io.h>

#include "test_seq_io_generic.h"
#include "test_stream_write_fasta.h"

#include "test_genomic_region.h"
#include "test_fai_index.h"

#include "test_sequence_file.h"
#include "test_stream_read_embl.h"
#include "test_stream_read_genbank.h"

#include "test_read_bam.h"
#include "test_write_bam.h"

#include "test_tag_select_intersect.h"

SEQAN_BEGIN_TESTSUITE(test_seq_io)
{
    // Test recognition of supported file types.
    SEQAN_CALL_TEST(test_seq_io_sequence_file_recognize_file_type_gz_fasta);
    SEQAN_CALL_TEST(test_seq_io_sequence_file_recognize_file_type_bz2_fasta);

    // Test recognition of supported file formats.
    SEQAN_CALL_TEST(test_seq_io_sequence_file_recognize_file_format_text_fasta);
    SEQAN_CALL_TEST(test_seq_io_sequence_file_recognize_file_format_text_fastq);

    // Test transefer of format from input file to output file (tag_select_intersect).
    SEQAN_CALL_TEST(test_tag_select_intersect);

    // Test reading with different interfaces.
    SEQAN_CALL_TEST(test_seq_io_sequence_file_read_record_text_fasta);
    SEQAN_CALL_TEST(test_seq_io_sequence_file_read_all_text_fasta);

    // Test writing with different interfaces.
    SEQAN_CALL_TEST(test_seq_io_sequence_file_write_record_text_fasta);
    SEQAN_CALL_TEST(test_seq_io_sequence_file_write_all_text_fasta);
    SEQAN_CALL_TEST(test_seq_io_sequence_file_write_record_text_fastq_no_qual);
    SEQAN_CALL_TEST(test_seq_io_sequence_file_write_record_text_fastq_with_qual);
    SEQAN_CALL_TEST(test_seq_io_sequence_file_write_all_text_fastq_no_qual);
    SEQAN_CALL_TEST(test_seq_io_sequence_file_write_all_text_fastq_with_qual);

    // Test isOpen functionality
    SEQAN_CALL_TEST(test_seq_io_sequence_file_isOpen_fileIn);
    SEQAN_CALL_TEST(test_seq_io_sequence_file_isOpen_fileOut);

    // Test parsing for GenomicRegion.
    SEQAN_CALL_TEST(test_seq_io_genomic_region_default_constructed);
    SEQAN_CALL_TEST(test_seq_io_genomic_region_from_string);
    SEQAN_CALL_TEST(test_seq_io_genomic_region_clear);
    SEQAN_CALL_TEST(test_seq_io_genomic_region_parse_chrom);
    SEQAN_CALL_TEST(test_seq_io_genomic_region_parse_chrom_begin);
    SEQAN_CALL_TEST(test_seq_io_genomic_region_parse_chrom_begin_end);
    SEQAN_CALL_TEST(test_seq_io_genomic_region_to_string_interval);
    SEQAN_CALL_TEST(test_seq_io_genomic_region_to_string_point);

    // Test FaiIndex.
    SEQAN_CALL_TEST(test_seq_io_genomic_fai_index_build);
    SEQAN_CALL_TEST(test_seq_io_genomic_fai_index_write);
    SEQAN_CALL_TEST(test_seq_io_genomic_fai_index_read);
    SEQAN_CALL_TEST(test_seq_io_genomic_fai_index_read_sequence);
    SEQAN_CALL_TEST(test_seq_io_genomic_fai_index_read_region);

    // Tests for EMBL
    SEQAN_CALL_TEST(test_stream_read_embl_single_char_array_stream);
    SEQAN_CALL_TEST(test_stream_read_embl_record_char_array_stream);
    SEQAN_CALL_TEST(test_stream_read_embl_single_mmap);
    SEQAN_CALL_TEST(test_stream_read_embl_single_batch_mmap);

    // Tests for GenBank
    SEQAN_CALL_TEST(test_stream_read_genbank_single_char_array_stream);
    SEQAN_CALL_TEST(test_stream_read_genbank_record_char_array_stream);
    SEQAN_CALL_TEST(test_stream_read_genbank_single_mmap);
    SEQAN_CALL_TEST(test_stream_read_genbank_single_batch_mmap);

    // Tests for BAM-File
    SEQAN_CALL_TEST(test_seq_io_bam_file_sam_read_sequences);
    SEQAN_CALL_TEST(test_seq_io_bam_file_sam_read_sequences_and_qualities);
    SEQAN_CALL_TEST(test_seq_io_bam_file_sam_write_sequences);
    SEQAN_CALL_TEST(test_seq_io_bam_file_sam_write_sequences_and_qualities);

#if SEQAN_HAS_ZLIB
    SEQAN_CALL_TEST(test_seq_io_bam_file_bam_read_sequences);
    SEQAN_CALL_TEST(test_seq_io_bam_file_bam_read_sequences_and_qualities);
    SEQAN_CALL_TEST(test_seq_io_bam_file_bam_write_sequences);
    SEQAN_CALL_TEST(test_seq_io_bam_file_bam_write_sequences_and_qualities);
#endif
}

SEQAN_END_TESTSUITE
