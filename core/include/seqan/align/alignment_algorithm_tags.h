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
// We put the tag definition into its own header so we can include them
// independently from the algorithms.
// ==========================================================================

#ifndef SEQAN_CORE_INCLUDE_SEQAN_ALIGN_ALIGNMENT_ALGORITHM_TAGS_H_
#define SEQAN_CORE_INCLUDE_SEQAN_ALIGN_ALIGNMENT_ALGORITHM_TAGS_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Global Alignment Algorithm Tags
// ----------------------------------------------------------------------------

// TODO(holtgrew): Rename MyersBitVector to Myers? Clashes with find module at the moment.

/**
.Tag.Pairwise Global Alignment Algorithms
..cat:Alignments
..summary:Tags used for selecting pairwise global alignment algorithms.
..tag
...Gotoh:Gotoh's for affine gap costs.
...NeedlemanWunsch:The Needleman-Wunsch algorithm for linear gap costs.
...Hirschberg:Hirschberg's algorithm using linear space.
...MyersBitVector:Myer's bit-vector algorithm.
...MyersHirschberg:Combination of Myer's and Hirschberg's algorithm.
..see:Function.globalAlignment
..see:Function.globalAlignmentScore
..include:seqan/align.h
*/

struct Gotoh_;
typedef Tag<Gotoh_> Gotoh;

struct NeedlemanWunsch_;
typedef Tag<NeedlemanWunsch_> NeedlemanWunsch;

struct Hirschberg_;
typedef Tag<Hirschberg_> Hirschberg;

struct MyersBitVector_;
typedef Tag<MyersBitVector_> MyersBitVector;

struct MyersHirschberg_;
typedef Tag<MyersHirschberg_> MyersHirschberg;

// ----------------------------------------------------------------------------
// Local Alignment Algorithm Tags
// ----------------------------------------------------------------------------

/**
.Tag.Pairwise Local Alignment Algorithms
..cat:Alignments
..summary:Tags used for selecting pairwise global alignment algorithms.
..tag
...SmithWaterman:Smith-Waterman algorithm for local alignments.
...WatermanEggert:Smith-Waterman algorithm with declumping to identify suboptimal local alignments.
..see:Function.localAlignment
..see:Class.LocalAlignmentEnumerator
..include:seqan/align.h
*/

struct SmithWaterman_;
typedef Tag<SmithWaterman_> SmithWaterman;

struct WatermanEggert_;
typedef Tag<WatermanEggert_> WatermanEggert;

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

}  // namespace seqan

#endif  // #ifndef SEQAN_CORE_INCLUDE_SEQAN_ALIGN_ALIGNMENT_ALGORITHM_TAGS_H_
