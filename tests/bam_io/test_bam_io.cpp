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
#include <seqan/sequence.h>
#include <seqan/stream.h>

#include "test_bam_alignment_record.h"
#include "test_bam_file.h"
#include "test_bam_header_record.h"
#include "test_bam_io_context.h"
#include "test_bam_sam_conversion.h"
#include "test_bam_tags_dict.h"
#include "test_read_sam.h"
#include "test_write_sam.h"
#include "test_read_bam.h"
#include "test_write_bam.h"

#if SEQAN_HAS_ZLIB
#include "test_bam_index.h"
#endif

SEQAN_BEGIN_TESTSUITE(test_bam_io)
{
    // Test BamAlignmentRecord class.
    SEQAN_CALL_TEST(test_bam_io_bam_alignment_record_class);
    SEQAN_CALL_TEST(test_bam_io_bam_alignment_record_has_flag_multiple);
    SEQAN_CALL_TEST(test_bam_io_bam_alignment_record_has_flag_unmapped);
    SEQAN_CALL_TEST(test_bam_io_bam_alignment_record_has_flag_next_unmapped);
    SEQAN_CALL_TEST(test_bam_io_bam_alignment_record_has_flag_rc);
    SEQAN_CALL_TEST(test_bam_io_bam_alignment_record_has_flag_next_rc);
    SEQAN_CALL_TEST(test_bam_io_bam_alignment_record_has_flag_first);
    SEQAN_CALL_TEST(test_bam_io_bam_alignment_record_has_flag_last);
    SEQAN_CALL_TEST(test_bam_io_bam_alignment_record_has_flag_secondary);
    SEQAN_CALL_TEST(test_bam_io_bam_alignment_record_has_flag_qc_no_pass);
    SEQAN_CALL_TEST(test_bam_io_bam_alignment_record_has_flag_duplicate);
    SEQAN_CALL_TEST(test_bam_io_bam_alignment_record_has_flag_supplementary);

    // Test BamHeader and BamHeaderRecord classes.
    SEQAN_CALL_TEST(test_bam_io_bam_header_class);
    SEQAN_CALL_TEST(test_bam_io_bam_header_typedefs);
    SEQAN_CALL_TEST(test_bam_io_bam_header_record_class);
    SEQAN_CALL_TEST(test_bam_io_bam_header_record_typedefs);
    SEQAN_CALL_TEST(test_bam_io_bam_header_record_find_tag_key);
    SEQAN_CALL_TEST(test_bam_io_bam_header_record_get_tag_value);

    // Test BamIoContext.
    SEQAN_CALL_TEST(test_bam_io_bam_io_context_standalone);
    SEQAN_CALL_TEST(test_bam_io_bam_io_context_fragment_store);

    // Test BAM<->SAM tag conversion.
    SEQAN_CALL_TEST(test_assign_tags_bam_to_sam_two_tags);
    SEQAN_CALL_TEST(test_assign_tags_bam_to_sam_type_A);
    SEQAN_CALL_TEST(test_assign_tags_bam_to_sam_type_c);
    SEQAN_CALL_TEST(test_assign_tags_bam_to_sam_type_C);
    SEQAN_CALL_TEST(test_assign_tags_bam_to_sam_type_s);
    SEQAN_CALL_TEST(test_assign_tags_bam_to_sam_type_S);
    SEQAN_CALL_TEST(test_assign_tags_bam_to_sam_type_i);
    SEQAN_CALL_TEST(test_assign_tags_bam_to_sam_type_I);
    SEQAN_CALL_TEST(test_assign_tags_bam_to_sam_type_f);
    SEQAN_CALL_TEST(test_assign_tags_bam_to_sam_type_Z);
    SEQAN_CALL_TEST(test_assign_tags_bam_to_sam_type_H);
    SEQAN_CALL_TEST(test_assign_tags_bam_to_sam_type_Bc);
    SEQAN_CALL_TEST(test_assign_tags_bam_to_sam_type_BC);
    SEQAN_CALL_TEST(test_assign_tags_bam_to_sam_type_Bs);
    SEQAN_CALL_TEST(test_assign_tags_bam_to_sam_type_BS);
    SEQAN_CALL_TEST(test_assign_tags_bam_to_sam_type_Bi);
    SEQAN_CALL_TEST(test_assign_tags_bam_to_sam_type_BI);
    SEQAN_CALL_TEST(test_assign_tags_bam_to_sam_type_Bf);
    SEQAN_CALL_TEST(test_assign_tags_sam_to_bam_two_tags);
    SEQAN_CALL_TEST(test_assign_tags_sam_to_bam_type_A);
    SEQAN_CALL_TEST(test_assign_tags_sam_to_bam_type_i);
    SEQAN_CALL_TEST(test_assign_tags_sam_to_bam_type_f);
    SEQAN_CALL_TEST(test_assign_tags_sam_to_bam_type_Z);
    SEQAN_CALL_TEST(test_assign_tags_sam_to_bam_type_H);
    SEQAN_CALL_TEST(test_assign_tags_sam_to_bam_type_Bc);
    SEQAN_CALL_TEST(test_assign_tags_sam_to_bam_type_BC);
    SEQAN_CALL_TEST(test_assign_tags_sam_to_bam_type_Bs);
    SEQAN_CALL_TEST(test_assign_tags_sam_to_bam_type_BS);
    SEQAN_CALL_TEST(test_assign_tags_sam_to_bam_type_Bi);
    SEQAN_CALL_TEST(test_assign_tags_sam_to_bam_type_BI);
    SEQAN_CALL_TEST(test_assign_tags_sam_to_bam_type_Bf);

    // Test BamTagsDict.
    SEQAN_CALL_TEST(test_bam_tags_dict_get_type_size);
    SEQAN_CALL_TEST(test_bam_tags_dict_get_type);
    SEQAN_CALL_TEST(test_bam_tags_dict_length);
    SEQAN_CALL_TEST(test_bam_tags_dict_extract_value_type_A);
    SEQAN_CALL_TEST(test_bam_tags_dict_extract_value_type_c);
    SEQAN_CALL_TEST(test_bam_tags_dict_extract_value_type_C);
    SEQAN_CALL_TEST(test_bam_tags_dict_extract_value_type_s);
    SEQAN_CALL_TEST(test_bam_tags_dict_extract_value_type_S);
    SEQAN_CALL_TEST(test_bam_tags_dict_extract_value_type_i);
    SEQAN_CALL_TEST(test_bam_tags_dict_extract_value_type_I);
    SEQAN_CALL_TEST(test_bam_tags_dict_get_value_type_A);
    SEQAN_CALL_TEST(test_bam_tags_dict_get_value_type_c);
    SEQAN_CALL_TEST(test_bam_tags_dict_get_value_type_C);
    SEQAN_CALL_TEST(test_bam_tags_dict_get_value_type_s);
    SEQAN_CALL_TEST(test_bam_tags_dict_get_value_type_S);
    SEQAN_CALL_TEST(test_bam_tags_dict_get_value_type_i);
    SEQAN_CALL_TEST(test_bam_tags_dict_get_value_type_I);
    SEQAN_CALL_TEST(test_bam_tags_dict_get_value_type_f);
    SEQAN_CALL_TEST(test_bam_tags_dict_get_value_type_Z);
    SEQAN_CALL_TEST(test_bam_tags_dict_get_value_type_H);
    SEQAN_CALL_TEST(test_bam_tags_dict_get_value_type_Bc);
    SEQAN_CALL_TEST(test_bam_tags_dict_get_value_type_BC);
    SEQAN_CALL_TEST(test_bam_tags_dict_get_value_type_Bs);
    SEQAN_CALL_TEST(test_bam_tags_dict_get_value_type_BS);
    SEQAN_CALL_TEST(test_bam_tags_dict_get_value_type_Bi);
    SEQAN_CALL_TEST(test_bam_tags_dict_get_value_type_BI);
    SEQAN_CALL_TEST(test_bam_tags_dict_get_value_type_Bf);
    SEQAN_CALL_TEST(test_bam_tags_dict_erase_tag);
    SEQAN_CALL_TEST(test_bam_tags_dict_set_tag_value);
    SEQAN_CALL_TEST(test_bam_tags_dict_append_tag_value);
    SEQAN_CALL_TEST(test_bam_tags_dict_const_bam_tags_sequence);
    
    // Test SAM I/O.
    SEQAN_CALL_TEST(test_bam_io_sam_read_header);
    SEQAN_CALL_TEST(test_bam_io_sam_read_alignment);
    SEQAN_CALL_TEST(test_bam_io_sam_write_header);
    SEQAN_CALL_TEST(test_bam_io_sam_write_alignment);

    // Test BAM I/O.
    SEQAN_CALL_TEST(test_bam_io_bam_read_header);
    SEQAN_CALL_TEST(test_bam_io_bam_read_alignment);
    SEQAN_CALL_TEST(test_bam_io_bam_write_header);
    SEQAN_CALL_TEST(test_bam_io_bam_write_alignment);

    // Test isOpen.
    SEQAN_CALL_TEST(test_bam_io_bam_file_isOpen_fileIn);
    SEQAN_CALL_TEST(test_bam_io_bam_file_isOpen_fileOut);

#if SEQAN_HAS_ZLIB
    // Test BamStream class.
    SEQAN_CALL_TEST(test_bam_io_bam_file_sam_file_size);
    SEQAN_CALL_TEST(test_bam_io_bam_file_sam_read_header);
    SEQAN_CALL_TEST(test_bam_io_bam_file_sam_read_records);
    SEQAN_CALL_TEST(test_bam_io_bam_file_sam_write_header);
    SEQAN_CALL_TEST(test_bam_io_bam_file_sam_write_records);

    SEQAN_CALL_TEST(test_bam_io_bam_file_bam_file_size);
    SEQAN_CALL_TEST(test_bam_io_bam_file_bam_read_header);
    SEQAN_CALL_TEST(test_bam_io_bam_file_bam_read_records);
    SEQAN_CALL_TEST(test_bam_io_bam_file_bam_read_ex1);
    SEQAN_CALL_TEST(test_bam_io_bam_file_bam_write_header);
    SEQAN_CALL_TEST(test_bam_io_bam_file_bam_write_records);
    SEQAN_CALL_TEST(test_bam_io_bam_file_bam_file_seek);

    // Issue 489
    SEQAN_CALL_TEST(test_bam_io_sam_file_issue_489);

    // Test BAM indices.
    SEQAN_CALL_TEST(test_bam_io_bam_index_save);
//  TODO(dadi): uncomment when BamIndex.build index is fixed
//    SEQAN_CALL_TEST(test_bam_io_bam_index_build);
    SEQAN_CALL_TEST(test_bam_io_bam_index_open);
#endif
}
SEQAN_END_TESTSUITE
