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

#include "feature_map.h"

#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <iterator>

namespace rep_sep {

// ----------------------------------------------------------------------------
// Function eraseIf()
// ----------------------------------------------------------------------------

template <typename TValue, typename TPred>
void eraseIf(std::vector<TValue> & v, TPred pred)
{
    v.erase(std::remove_if(v.begin(), v.end(), pred), v.end());
}

// ----------------------------------------------------------------------------
// Class FeatureDescription
// ----------------------------------------------------------------------------

std::ostream & operator<<(std::ostream & out, FeatureDescription const & desc)
{
    char const * label[2] = { "COLUMN", "LINK" };
    return out << "FeatureDescription(id=" << desc.id << ", kind=" << label[(int)desc.kind]
               << ", contigID=" << desc.contigID << ", pos=" << desc.pos
               << ", coverage=" << desc.coverage << ")";
}

// ----------------------------------------------------------------------------
// Class FeatureMap
// ----------------------------------------------------------------------------

namespace {  // anonymous namespace
}  // anonymous namespace

void FeatureMap::insert(FeatureDescription const & desc)
{
    needsRefresh = true;
    features.push_back(desc);
}

void FeatureMap::refresh()
{
    needsRefresh = false;

    std::sort(features.begin(), features.end());

    map.clear();
    for (unsigned i = 0; i < features.size(); ++i)
        map[features[i].id] = i;
}

std::pair<FeatureMap::TIterator, FeatureMap::TIterator> FeatureMap::contigFeatures(unsigned contigID) const
{
    if (needsRefresh)
        throw std::runtime_error("Requires refresh.");
    FeatureDescription val;
    val.kind = FeatureDescription::COLUMN;
    val.contigID = contigID;
    auto lt = [](FeatureDescription const & lhs, FeatureDescription const & rhs) {
        return (std::make_pair((int)lhs.kind, lhs.contigID) < std::make_pair((int)rhs.kind, rhs.contigID));
    };
    return std::equal_range(features.begin(), features.end(), val, lt);
}

std::pair<FeatureMap::TIterator, FeatureMap::TIterator> FeatureMap::contigFeatures() const
{
    if (needsRefresh)
        throw std::runtime_error("Requires refresh.");
    FeatureDescription val;
    val.kind = FeatureDescription::COLUMN;
    auto lt = [](FeatureDescription const & lhs, FeatureDescription const & rhs) { return ((int)lhs.kind < (int)rhs.kind); };
    return std::equal_range(features.begin(), features.end(), val, lt);
}

std::pair<FeatureMap::TIterator, FeatureMap::TIterator> FeatureMap::linkFeatures() const
{
    if (needsRefresh)
        throw std::runtime_error("Requires refresh.");
    FeatureDescription val;
    val.kind = FeatureDescription::LINK;
    auto lt = [](FeatureDescription const & lhs, FeatureDescription const & rhs) { return ((int)lhs.kind < (int)rhs.kind); };
    return std::equal_range(features.begin(), features.end(), val, lt);
}

void FeatureMap::applyContigIDMapping(std::vector<unsigned> const & oldToNew)
{
    eraseIf(features, [&oldToNew](FeatureDescription const & desc) {
            if (desc.kind != FeatureDescription::COLUMN)
                return false;
            else
                return (oldToNew.at(desc.contigID) == (unsigned)-1);
        });
    std::for_each(features.begin(), features.end(), [&oldToNew](FeatureDescription & desc) {
            if (desc.kind == FeatureDescription::COLUMN)
                desc.contigID = oldToNew.at(desc.contigID);
        });
    refresh();
}

FeatureMap::TIterator FeatureMap::find(unsigned featureID) const
{
    auto it = map.find(featureID);
    if (it == map.end())
        return features.end();
    else
        return features.begin() + it->second;
}

void FeatureMap::print(std::ostream & out) const
{
    out << "FeatureMap (#features=" << features.size() << ")\n";
    std::copy(features.begin(), features.end(), std::ostream_iterator<FeatureDescription>(out, "\n"));
}

}  // namespace rep_sep
