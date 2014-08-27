// ==========================================================================
//                                   ANISE
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

#include "clique_set.h"

#include <algorithm>
#include <iostream>

#include <boost/dynamic_bitset.hpp>

namespace rep_sep {

namespace {  // anonymous namespace

// ----------------------------------------------------------------------------
// Function eraseIf()
// ----------------------------------------------------------------------------

template <typename TValue, typename TPred>
void eraseIf(std::vector<TValue> & v, TPred pred)
{
    v.erase(std::remove_if(v.begin(), v.end(), pred), v.end());
}

}  // anonymous namespace

// --------------------------------------------------------------------------
// Class CliqueSet
// --------------------------------------------------------------------------

void CliqueSet::deactivateCliques(unsigned contigID, int beginPos)
{
    auto deactivate = [&](unsigned cliqueID) {
        return (cliques_[cliqueID].contigID != contigID ||
                cliques_[cliqueID].endPos <= beginPos);
    };
    eraseIf(activeCliques_, deactivate);
}

void CliqueSet::addClique(Clique const & clique)
{
    activeCliques_.push_back(cliques_.size());
    cliques_.push_back(clique);
}

void CliqueSet::processRead(Read const & read,
                            boost::dynamic_bitset<> const & nbh,
                            FeatureReadSet const & readSet,
                            bool logging)
{
    if (logging)
        std::cerr << "Processing read " << read << "\n\t#active cliques: " << activeCliques_.size() << "\n";
    // Deactivate cliques that cannot become active any more.
    deactivateCliques(read.contigID, read.beginPos);

    bool intersectedWithAny = false;  // set to true below in case of intersection
    std::vector<Clique> newCliques;  // tentative new cliques to be added below

    // Try to combine read and nbh with any current clique.
    std::vector<unsigned> modifiedCliqueIDs;  // collect modified clique IDs
    boost::dynamic_bitset<> intersectionResult;  // buffer used below
    for (auto cliqueID : activeCliques_)
    {
        auto & clique = cliques_[cliqueID];

        if (logging)
            std::cerr << "Considering clique " << clique << "\n";

        intersectionResult.clear();
        clique.intersection(intersectionResult, nbh);
        auto num = intersectionResult.count();
        if (!num)
            continue;  // no intersection yet
        intersectedWithAny = (intersectedWithAny || (num > 0));

        if (logging)
            std::cerr << "num == " << num << "\n";
        if (num == clique.size())
        {
            clique.addRead(read);
            modifiedCliqueIDs.push_back(cliqueID);
            if (logging)
                std::cerr << "clique.addRead(read)\t(" << __LINE__ << ")\n"
                          << " => clique == " << clique << "\n";
        }
        else  // (num < clique.size())
        {
            Clique newClique(read.contigID, intersectionResult, readSet);
            newClique.addRead(read);
            if (logging)
                std::cerr << "newClique == " << newClique << "\n";
            newCliques.push_back(newClique);
        }
    }

    std::sort(modifiedCliqueIDs.begin(), modifiedCliqueIDs.end());

    // Remove duplicates (by readIDs list) first.
    auto ltReadIDs = [](Clique const & lhs, Clique const & rhs) { return (lhs.bitSet < rhs.bitSet); };
    std::sort(newCliques.begin(), newCliques.end(), ltReadIDs);
    auto eqReadIDs = [](Clique const & lhs, Clique const & rhs) { return (lhs.bitSet == rhs.bitSet); };
    newCliques.erase(std::unique(newCliques.begin(), newCliques.end(), eqReadIDs), newCliques.end());

    // Purge duplicates and subsumed from new cliques.
    for (auto it = newCliques.begin(); it != newCliques.end(); ++it)
    {
        auto includes2 = [&](unsigned cliqueID) { return cliques_[cliqueID].includes(*it); };
        auto includes = [&](Clique const & other) { return other.includes(*it); };
        if (logging)
        {
            std::cerr << "Checking clique " << *it << "\n"
                      << "\tincluded modified?\t" << std::any_of(modifiedCliqueIDs.begin(), modifiedCliqueIDs.end(), includes2) << "\n"
                      << "\tincluded new 1/2?\t" << std::any_of(newCliques.begin(), it, includes) << "\n"
                      << "\tincluded new 2/2?\t" << std::any_of(std::next(it), newCliques.end(), includes) << "\n";
        }
        if (!std::any_of(modifiedCliqueIDs.begin(), modifiedCliqueIDs.end(), includes2) &&
            !std::any_of(newCliques.begin(), it, includes) &&
            !std::any_of(std::next(it), newCliques.end(), includes))
        {
            if (logging)
                std::cerr << "addClique(" << *it << ")\t(" << __LINE__ << ")\n";
            addClique(*it);
        }
    }

    // Add new singleton clique in case it did not intersect with any.
    if (!intersectedWithAny)
    {
        Clique clique(read, numReads_);
        if (logging)
            std::cerr << "addClique(" << clique << ")\t(" << __LINE__ << ")\n";
        addClique(clique);
    }
}

}  // namespace rep_sep
