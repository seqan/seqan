// ==========================================================================
//                                  fm_index
// ==========================================================================
// Copyright (c) 2006-2011, Knut Reinert, FU Berlin
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
// Author: Jochen Singer <jochen.singer@fu-berlin.de>
// ==========================================================================
// Facade header for module sequence.
// ==========================================================================

#ifndef EXTRAS_INCLUDE_SEQAN_INDEX_H_
#define EXTRAS_INCLUDE_SEQAN_INDEX_H_

// ===========================================================================
// Prerequisites.
// ===========================================================================

#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/index.h>

// ===========================================================================
// First Header Group.
// ===========================================================================

#include <seqan/index/index_fm_rank_support_bit_string.h>
#include <seqan/index/index_fm_rank_support_bit_string_iterator.h>
#include <seqan/index/index_fm_sparse_string.h>
#include <seqan/index/index_fm_compressed_sa.h>
#include <seqan/index/index_fm_compressed_sa_iterator.h>
#include <seqan/index/index_fm_prefix_sum_table.h>
#include <seqan/index/index_fm_lf_table.h>
#include <seqan/index/index_fm_right_array_binary_tree.h>
#include <seqan/index/index_fm_right_array_binary_tree_iterator.h>
#include <seqan/index/index_fm_wavelet_tree.h>
#include <seqan/index/index_fm.h>
#include <seqan/index/index_fm_stree.h>

#endif  // SANDBOX_SINGER_INCLUDE_SEQAN_FM_INDEX_H_
