// ==========================================================================
//                                  fm_sequence
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
// Facade header for module fm_sequence.
// ==========================================================================

#ifndef SANDBOX_SINGER_INCLUDE_SEQAN_FM_SEQUENCE_H_
#define SANDBOX_SINGER_INCLUDE_SEQAN_FM_SEQUENCE_H_

// ===========================================================================
// Prerequisites.
// ===========================================================================

#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/sequence.h>
#include <seqan/index.h>

// ===========================================================================
// First Header Group.
// ===========================================================================

//#include "seqan/fm_sequence/prefix_sum_table.h"
//#include "seqan/fm_sequence/bit_string.h"
//#include "seqan/fm_sequence/bit_string_beta.h"
//#include "seqan/fm_sequence/sparse_string.h"
//#include "seqan/fm_sequence/compressed_sa.h"
//#include "seqan/fm_sequence/compressed_sa_impl.h"
//#include "seqan/fm_sequence/wavelet_tree.h"
//#include "seqan/fm_sequence/wavelet_tree_structure.h"

#include "seqan/fm_sequence/bit_string_beta.h"
#include "seqan/fm_sequence/sparse_string_beta.h"
#include "seqan/fm_sequence/compressed_sa_beta.h"
#include "seqan/fm_sequence/compressed_sa_iterator_beta.h"
#include "seqan/fm_sequence/lf_table_beta.h"
#include "seqan/fm_sequence/prefix_sum_table_beta.h"
#include "seqan/fm_sequence/wavelet_tree_structure_beta.h"
#include "seqan/fm_sequence/wavelet_tree_structure_iterator_beta.h"
#include "seqan/fm_sequence/wavelet_tree_beta.h"

#endif  // SANDBOX_SINGER_INCLUDE_SEQAN_FM_SEQUENCE_H_
