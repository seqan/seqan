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

#ifndef SANDBOX_HERBARIUM_APPS_ANISE_REP_SEP_CLIQUE_SET_H_
#define SANDBOX_HERBARIUM_APPS_ANISE_REP_SEP_CLIQUE_SET_H_

#include <vector>

#include <boost/dynamic_bitset.hpp>

#include "rep_sep/clique.h"

namespace rep_sep {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class CliqueSet
// ----------------------------------------------------------------------------

// The set of cliques.

class CliqueSet
{
public:
    explicit CliqueSet(unsigned numReads) : numReads_(numReads) {}
    

    // Process the given read, adding to cliques or creating new one.
    //
    // nbh is the neighbourhood of read (overlaps, EXCLUDING the read itself, sorted ascendingly), readSet is the read
    // set.  The entries of nbh are ids/indices into readSet.
    void processRead(Read const & read,
                     boost::dynamic_bitset<> const & nbh,
                     FeatureReadSet const & readSet,
                     bool logging = false);

    std::vector<Clique> const & cliques() const { return cliques_; }
    std::vector<unsigned> const & activeCliques() const { return activeCliques_; }

private:
    // Deactivate cliques that end left of (contigID, beginPos).
    void deactivateCliques(unsigned contigID, int beginPos);

    // Add a new clique and mark it as active.
    void addClique(Clique const & clique);

    // Set of cliques and the ids of the active cliques.
    std::vector<Clique> cliques_;
    std::vector<unsigned> activeCliques_;

    // Number of reads managed.
    unsigned numReads_;
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

}  // namespace rep_sep

#endif  // #ifndef SANDBOX_HERBARIUM_APPS_ANISE_REP_SEP_CLIQUE_SET_H_
