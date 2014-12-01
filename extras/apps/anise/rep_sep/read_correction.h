// ==========================================================================
//                                   ANISE
// ==========================================================================
// Copyright (c) 2006-2014, Knut Reinert, FU Berlin
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

#ifndef SANDBOX_HERBARIUM_APPS_ANISE_REP_SEP_READ_CORRECTION_H_
#define SANDBOX_HERBARIUM_APPS_ANISE_REP_SEP_READ_CORRECTION_H_

#include <functional>
#include <vector>

#include "asm/frag_store.h"

// ============================================================================
// External Forwards
// ============================================================================

namespace seqan { class BamAlignmentRecord; }

namespace rep_sep {

// ============================================================================
// Forwards
// ============================================================================

class FeatureMap;
struct FeatureReadSet;
struct ReadSeparatorOptions;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

typedef std::function<bool(seqan::BamAlignmentRecord const &)> TIsReadCorrectable;

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function performReadCorrection()
// ----------------------------------------------------------------------------

// Correct sequences of reads in records using the consensus and alignment from store.
//
// Only reads with "mS" (KEY_BIRTH_STEP) tag value <= maxStepNo are corrected.  Columns stored in FeatureMap are
// corrected using the features in readSet.  Note that this can actually assign the same read base multiple times but we
// allow for such few errors.

void performReadCorrection(std::vector<seqan::BamAlignmentRecord> & records,
                           TFragmentStore & store,
                           FeatureMap const & featureMap,
                           TIsReadCorrectable isReadCorrectable);

}  // namespace rep_sep

#endif  // #ifndef SANDBOX_HERBARIUM_APPS_ANISE_REP_SEP_READ_CORRECTION_H_
