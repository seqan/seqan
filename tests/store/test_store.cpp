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
// Tests for the SeqAn module store.
// ==========================================================================

#include <seqan/basic.h>
#include "test_store_io.h"

SEQAN_BEGIN_TESTSUITE(test_store)
{
    // the UCSC knownGene format
    SEQAN_CALL_TEST(test_store_io_read_ucsc_known_genes);
    SEQAN_CALL_TEST(test_store_io_read_ucsc_known_genes_and_isoforms);
    SEQAN_CALL_TEST(test_store_io_write_ucsc_known_genes);

    // the gff format
    SEQAN_CALL_TEST(test_store_io_read_gff);
    SEQAN_CALL_TEST(test_store_io_write_gff);

    // the gtf format
    SEQAN_CALL_TEST(test_store_io_read_gtf);
    SEQAN_CALL_TEST(test_store_io_write_gtf);

    // Tests for the AMOS format.
    SEQAN_CALL_TEST(test_store_io_readwrite_amos);
    SEQAN_CALL_TEST(test_store_io_read_amos);
    SEQAN_CALL_TEST(test_store_io_write_amos);

    // Tests for the SAM/BAM format.
    SEQAN_CALL_TEST(test_store_io_sam);
    SEQAN_CALL_TEST(test_store_io_sam2);
    SEQAN_CALL_TEST(test_store_io_split_sam);
#if SEQAN_HAS_ZLIB
    SEQAN_CALL_TEST(test_store_io_read_bam);
#endif  // #if SEQAN_HAS_ZLIB
}
SEQAN_END_TESTSUITE
