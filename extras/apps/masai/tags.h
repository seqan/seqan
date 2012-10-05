// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
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
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================
// This file contains tags.
// ==========================================================================

#ifndef SEQAN_EXTRAS_MASAI_TAGS_H_
#define SEQAN_EXTRAS_MASAI_TAGS_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

using namespace seqan;

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct All_;
struct AllBest_;
struct KBest_;
struct AnyBest_;

typedef Tag<All_>       All;
typedef Tag<AllBest_>   AllBest;
typedef Tag<KBest_>     KBest;
typedef Tag<AnyBest_>   AnyBest;

// ============================================================================

struct SingleEnd_;
struct PairedEnd_;

typedef Tag<SingleEnd_>     SingleEnd;
typedef Tag<PairedEnd_>     PairedEnd;

struct LeftMate_;
struct RightMate_;

typedef Tag<LeftMate_>      LeftMate;
typedef Tag<RightMate_>     RightMate;

struct LeftFile_;
struct RightFile_;

typedef Tag<LeftFile_>      LeftFile;
typedef Tag<RightFile_>     RightFile;

// ============================================================================

struct SingleBacktracking_;
struct MultipleBacktracking_;

typedef Tag<SingleBacktracking_>    SingleBacktracking;
typedef Tag<MultipleBacktracking_>  MultipleBacktracking;


// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

#endif  // #ifndef SEQAN_EXTRAS_MASAI_TAGS_H_
