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

#include <seqan/basic.h>
#include <seqan/file.h>

#include <seqan/seq_io.h>

//#include "test_seq_io_generic.h"
//
//#include "test_sequence_stream.h"
//#include "test_genomic_region.h"
//#include "test_fai_index.h"
//#include "test_stream_record_reader_fasta.h"
//#include "test_stream_record_reader_fastq.h"
//#include "test_stream_guess_stream_format.h"
#include "test_stream_write_fasta.h"
//#include "test_stream_read_embl.h"
//#include "test_stream_read_genbank.h"
//#include "test_stream_read_auto_format.h"

SEQAN_BEGIN_TESTSUITE(test_seq_io)
{
#if 0
    // Test simple readFasta() function.
    SEQAN_CALL_TEST(test_seq_io_read_fasta);

    // Test recognition of supported file types.
    SEQAN_CALL_TEST(test_seq_io_sequence_stream_recognize_file_type_text_fasta);
    SEQAN_CALL_TEST(test_seq_io_sequence_stream_recognize_file_type_gz_fasta);
    SEQAN_CALL_TEST(test_seq_io_sequence_stream_recognize_file_type_bz2_fasta);

    // Test recognition of supported file formats.
    SEQAN_CALL_TEST(test_seq_io_sequence_stream_recognize_file_format_text_fasta);
    SEQAN_CALL_TEST(test_seq_io_sequence_stream_recognize_file_format_text_fastq);

    // Test reading with different interfaces.
    SEQAN_CALL_TEST(test_seq_io_sequence_stream_read_record_text_fasta);
    SEQAN_CALL_TEST(test_seq_io_sequence_stream_read_batch_text_fasta);
    SEQAN_CALL_TEST(test_seq_io_sequence_stream_read_all_text_fasta);

    // Test writing with different interfaces.
    SEQAN_CALL_TEST(test_seq_io_sequence_stream_write_record_text_fasta);
    SEQAN_CALL_TEST(test_seq_io_sequence_stream_write_all_text_fasta);
    SEQAN_CALL_TEST(test_seq_io_sequence_stream_write_record_text_fastq_no_qual);
    SEQAN_CALL_TEST(test_seq_io_sequence_stream_write_record_text_fastq_with_qual);
    SEQAN_CALL_TEST(test_seq_io_sequence_stream_write_all_text_fastq_no_qual);
    SEQAN_CALL_TEST(test_seq_io_sequence_stream_write_all_text_fastq_with_qual);

    // Test parsing for GenomicRegion.
    SEQAN_CALL_TEST(test_seq_io_genomic_region_default_constructed);
    SEQAN_CALL_TEST(test_seq_io_genomic_region_from_string);
    SEQAN_CALL_TEST(test_seq_io_genomic_region_clear);
    SEQAN_CALL_TEST(test_seq_io_genomic_region_parse_chrom);
    SEQAN_CALL_TEST(test_seq_io_genomic_region_parse_chrom_begin);
    SEQAN_CALL_TEST(test_seq_io_genomic_region_parse_chrom_begin_end);

    // Test FaiIndex.
    SEQAN_CALL_TEST(test_seq_io_genomic_fai_index_build);
    SEQAN_CALL_TEST(test_seq_io_genomic_fai_index_write);
    SEQAN_CALL_TEST(test_seq_io_genomic_fai_index_read);
    SEQAN_CALL_TEST(test_seq_io_genomic_fai_index_read_sequence);
    SEQAN_CALL_TEST(test_seq_io_genomic_fai_index_read_region);

    // -------------- File format specific code ------------------

    /* RecordReader is tested for each Format with fstream and mmap-backend.
     * Every backend is tested with single-pass record-reading, those that have
     * double-pass implementations are also tested for double-pass
     * record-reading and batch reading (double-pass whole document). For the
     * last case an additional test for concat-direct-strings is performed. */

    // Tests for FASTA
    SEQAN_CALL_TEST(test_stream_record_reader_fasta_single_fstream);
    SEQAN_CALL_TEST(test_stream_record_reader_fasta_double_fstream);
    SEQAN_CALL_TEST(test_stream_record_reader_fasta_batch_single_fstream);
    SEQAN_CALL_TEST(test_stream_record_reader_fasta_batch_double_fstream);
    SEQAN_CALL_TEST(test_stream_record_reader_fasta_single_mmap);
    SEQAN_CALL_TEST(test_stream_record_reader_fasta_double_mmap);
    SEQAN_CALL_TEST(test_stream_record_reader_fasta_batch_single_mmap);
    SEQAN_CALL_TEST(test_stream_record_reader_fasta_batch_single_concat_mmap);
    SEQAN_CALL_TEST(test_stream_record_reader_fasta_batch_double_mmap);
    SEQAN_CALL_TEST(test_stream_record_reader_fasta_batch_double_concat_mmap);

    // Tests for FASTQ
    SEQAN_CALL_TEST(test_stream_record_reader_fastq_single_fstream);
    SEQAN_CALL_TEST(test_stream_record_reader_fastq_double_fstream);
    SEQAN_CALL_TEST(test_stream_record_reader_fastq_batch_single_fstream);
    SEQAN_CALL_TEST(test_stream_record_reader_fastq_batch_double_fstream);
    SEQAN_CALL_TEST(test_stream_record_reader_fastq_single_mmap);
    SEQAN_CALL_TEST(test_stream_record_reader_fastq_double_mmap);
    SEQAN_CALL_TEST(test_stream_record_reader_fastq_batch_single_mmap);
    SEQAN_CALL_TEST(test_stream_record_reader_fastq_batch_single_concat_mmap);
    SEQAN_CALL_TEST(test_stream_record_reader_fastq_batch_double_mmap);
    SEQAN_CALL_TEST(test_stream_record_reader_fastq_batch_double_concat_mmap);
    SEQAN_CALL_TEST(test_stream_record_reader_fastq_check_stream_format);

    // Tests for EMBL
    SEQAN_CALL_TEST(test_stream_read_embl_single_char_array_stream);
    SEQAN_CALL_TEST(test_stream_read_embl_record_char_array_stream);
    SEQAN_CALL_TEST(test_stream_read_embl_batch_char_array_stream);
    SEQAN_CALL_TEST(test_stream_read_embl_single_mmap);
    SEQAN_CALL_TEST(test_stream_read_embl_single_batch_mmap);
    SEQAN_CALL_TEST(test_stream_read_embl_single_batch_concat_mmap);

    // Tests for GenBank
    SEQAN_CALL_TEST(test_stream_read_genbank_single_char_array_stream);
    SEQAN_CALL_TEST(test_stream_read_genbank_record_char_array_stream);
    SEQAN_CALL_TEST(test_stream_read_genbank_batch_char_array_stream);
    SEQAN_CALL_TEST(test_stream_read_genbank_single_mmap);
    SEQAN_CALL_TEST(test_stream_read_genbank_single_batch_mmap);
    SEQAN_CALL_TEST(test_stream_read_genbank_single_batch_concat_mmap);

    // TODO Tests for other formats once they are supported

    // Tests for file format auto-detection
    SEQAN_CALL_TEST(test_stream_guess_stream_format_auto_fasta);
    SEQAN_CALL_TEST(test_stream_guess_stream_format_auto_fastq);
    SEQAN_CALL_TEST(test_stream_guess_stream_format_auto_bogus);

    // Tests for reading with automatic file format detection.
    SEQAN_CALL_TEST(test_stream_read_record_auto_format_quals_fasta);
    SEQAN_CALL_TEST(test_stream_read_record_auto_format_quals_fastq);
    SEQAN_CALL_TEST(test_stream_read_record_auto_format_no_quals_fasta);
    SEQAN_CALL_TEST(test_stream_read_record_auto_format_no_quals_fastq);

    SEQAN_CALL_TEST(test_stream_read_auto_format_quals_fasta);
    SEQAN_CALL_TEST(test_stream_read_auto_format_quals_fastq);
    SEQAN_CALL_TEST(test_stream_read_auto_format_no_quals_fasta);
    SEQAN_CALL_TEST(test_stream_read_auto_format_no_quals_fastq);

    // Tests for FASTA with Amino Acid alphabet.
    SEQAN_CALL_TEST(test_stream_record_reader_fasta_protein_single_fstream);
    SEQAN_CALL_TEST(test_stream_record_reader_fasta_annotated_protein_single_fstream);
#endif
    /* Tests for writing file formats.
     * Only fstream used as a representative */

    // Tests for FASTA
    SEQAN_CALL_TEST(test_stream_write_record_fasta_default);
    SEQAN_CALL_TEST(test_stream_write_record_fasta_nolinebreaks);
    SEQAN_CALL_TEST(test_stream_write_record_fastq_default_separate_qual);
    SEQAN_CALL_TEST(test_stream_write_record_fastq_default_qual_in_seq);
    SEQAN_CALL_TEST(test_stream_write_record_fastq_linebreaks_qualmeta);
//    SEQAN_CALL_TEST(test_stream_write2_fasta_default);
//    SEQAN_CALL_TEST(test_stream_write2_fastq_default_separate_qual);
//    SEQAN_CALL_TEST(test_stream_write2_fastq_default_qual_in_seq);
}
SEQAN_END_TESTSUITE
