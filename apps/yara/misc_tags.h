// ==========================================================================
//                      Yara - Yet Another Read Aligner
// ==========================================================================
// Copyright (c) 2011-2018, Enrico Siragusa, FU Berlin
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
//     * Neither the name of Enrico Siragusa or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ENRICO SIRAGUSA OR THE FU BERLIN BE LIABLE
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

#ifndef APP_YARA_MISC_TAGS_H_
#define APP_YARA_MISC_TAGS_H_

#include <seqan/basic.h>

using namespace seqan;

// ============================================================================
// Enums
// ============================================================================

enum MappingMode
{
    STRATA, ALL
};

enum LibraryOrientation
{
    FWD_REV, FWD_FWD, REV_REV, ANY
};

enum SecondaryAlignments
{
    TAG, RECORD, OMIT
};

enum Sensitivity
{
    LOW, HIGH, FULL
};

// ============================================================================
// Tags
// ============================================================================

// ----------------------------------------------------------------------------
// Mapping Strategy Tags
// ----------------------------------------------------------------------------

struct Strata_;
struct All_;

typedef Tag<Strata_>    Strata;
typedef Tag<All_>       All;

// ----------------------------------------------------------------------------
// Pairing Strategy Tags
// ----------------------------------------------------------------------------

struct AnchorOne_;
struct AnchorBoth_;

typedef Tag<AnchorBoth_>    AnchorOne;
typedef Tag<AnchorOne_>     AnchorBoth;

// ----------------------------------------------------------------------------
// Sequencing Technologies Tags
// ----------------------------------------------------------------------------

struct SingleEnd_;
struct PairedEnd_;

typedef Tag<SingleEnd_>     SingleEnd;
typedef Tag<PairedEnd_>     PairedEnd;

// ----------------------------------------------------------------------------
// Paired-End / Mate-Pairs Tags
// ----------------------------------------------------------------------------

struct LeftMate_;
struct RightMate_;

typedef Tag<LeftMate_>      LeftMate;
typedef Tag<RightMate_>     RightMate;

struct FwdRev_;
struct RevFwd_;
struct FwdFwd_;
struct RevRev_;

typedef Tag<FwdRev_>        FwdRev;
typedef Tag<RevFwd_>        RevFwd;
typedef Tag<FwdFwd_>        FwdFwd;
typedef Tag<RevRev_>        RevRev;

#endif  // #ifndef APP_YARA_MISC_TAGS_H_
