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

#ifndef SANDBOX_HERBARIUM_APPS_ANISE_ASM_BOG_UNITIGGER_H_
#define SANDBOX_HERBARIUM_APPS_ANISE_ASM_BOG_UNITIGGER_H_

#include "asm/unitigger_impl.h"

#include "asm/best_overlap_graph.h"

namespace assembler {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class BogUnitigger
// ----------------------------------------------------------------------------

// Unitigging using a CABOG-like best overlap approach.

class BogUnitigger : public UnitiggerImpl
{
public:
    BogUnitigger(bool logging = false) :
            withTimer([](char const *, std::function<void()> fn) {
                    return fn(); }),
            logging(logging)
    {}

    BogUnitigger(std::function<void(char const *, std::function<void()>)> withTimer,
                  bool logging = false) : withTimer(withTimer), logging(logging)
    {}

    std::unique_ptr<ContigGraph> run(seqan::StringSet<seqan::Dna5String> const & seqs,
                                     std::vector<Overlap> const & overlaps) const;

protected:

    // Return vector with read lengths.
    std::vector<unsigned> computeReadLengths(
            seqan::StringSet<seqan::Dna5String> const & seqs) const;

    // Compute ContigGraph from BogPathSet and set of initial overlaps.
    std::unique_ptr<ContigGraph> computeContigGraph(BogPathSet const & bogPathSet,
                                                    ContainmentFilter const & filter,
                                                    std::vector<Overlap> const & overlaps,
                                                    seqan::StringSet<seqan::Dna5String> const & seqs) const;

    // Augment ContigGraph with the contained alignments in filter.
    std::unique_ptr<ContigGraph> augmentContigGraph(std::unique_ptr<ContigGraph> in,
                                                    ContainmentFilter const & filter,
                                                    seqan::StringSet<seqan::Dna5String> const & seqs) const;

    // Function to use for timer.
    std::function<void(char const *, std::function<void()>)> withTimer;
    // Flag that enables/disables logging.
    bool logging;
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

}  // namespace assembler

#endif  // #ifndef SANDBOX_HERBARIUM_APPS_ANISE_ASM_BOG_UNITIGGER_H_
