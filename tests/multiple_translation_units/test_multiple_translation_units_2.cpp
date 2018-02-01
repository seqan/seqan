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
// Author: Jochen Singer <jochen.singer@fu-berlin.de>
// ==========================================================================

// This test simply checks whether all functions are inline or templates.

#include <seqan/align.h>
#include <seqan/align_extend.h>
#include <seqan/align_parallel.h>
#include <seqan/align_profile.h>
#include <seqan/align_split.h>
#include <seqan/alignment_free.h>
#include <seqan/arg_parse.h>
#include <seqan/bam_io.h>
#include <seqan/basic.h>
#include <seqan/bed_io.h>
#include <seqan/consensus.h>
#include <seqan/file.h>
#include <seqan/find.h>
#include <seqan/gff_io.h>
#include <seqan/graph_algorithms.h>
#include <seqan/graph_align.h>
#include <seqan/graph_msa.h>
#include <seqan/graph_types.h>
#include <seqan/index.h>
#include <seqan/journaled_set.h>
#include <seqan/map.h>
#include <seqan/math.h>
#include <seqan/modifier.h>
#include <seqan/parallel.h>
#include <seqan/parse_lm.h>
#include <seqan/pipe.h>
#include <seqan/platform.h>
#include <seqan/random.h>
#include <seqan/realign.h>
#include <seqan/reduced_aminoacid.h>
#include <seqan/roi_io.h>
#include <seqan/score.h>
#include <seqan/seeds.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <seqan/sequence_journaled.h>
#include <seqan/simple_intervals_io.h>
#include <seqan/statistics.h>
#include <seqan/store.h>
#include <seqan/stream.h>
#include <seqan/system.h>
#include <seqan/translation.h>
#include <seqan/ucsc_io.h>
#include <seqan/vcf_io.h>
#include <seqan/version.h>
