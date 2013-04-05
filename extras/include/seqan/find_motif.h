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

#ifndef SEQAN_HEADER_FIND_MOTIF_H
#define SEQAN_HEADER_FIND_MOTIF_H

//____________________________________________________________________________
// prerequisites

#include <algorithm>
#include <iterator>
#include <numeric>
#include <vector>
#include <set>
#include <cfloat>

#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/index/shape_base.h>
#include <seqan/index/shape_gapped.h>
#include <seqan/random.h>

//____________________________________________________________________________

#include <seqan/find_motif/sequence_model_types.h>

#include <seqan/find_motif/pseudocount_base.h>
#include <seqan/find_motif/pseudocount_mode_c.h>
#include <seqan/find_motif/pseudocount_mode_p.h>

#include <seqan/find_motif/frequency_distribution.h>
#include <seqan/find_motif/profile.h>

#include <seqan/find_motif/find_motif_base.h>
#include <seqan/find_motif/find_motif_pms1.h>
#include <seqan/find_motif/find_motif_pmsp.h>
#include <seqan/find_motif/find_motif_projection.h>
#include <seqan/find_motif/find_motif_epatternbranching.h>
#include <seqan/find_motif/em_algorithm.h>


#endif //#ifndef SEQAN_HEADER_...
