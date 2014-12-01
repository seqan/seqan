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

#ifndef SANDBOX_HERBARIUM_APPS_ANISE_REP_SEP_FEATURE_MAP_H_
#define SANDBOX_HERBARIUM_APPS_ANISE_REP_SEP_FEATURE_MAP_H_

#include <tuple>
#include <map>
#include <vector>
#include <iosfwd>

namespace rep_sep {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class FeatureDescription
// ----------------------------------------------------------------------------

// TODO(holtgrew): Add value distribution?

struct FeatureDescription
{
    static const unsigned INVALID = (unsigned)-1;
    // Describes the kind of the feature (1) separating column, (2) paired-end link between two contigs that is in
    // conflict with another link.
    enum Kind { COLUMN, LINK };

    FeatureDescription() = default;
    FeatureDescription(unsigned id, Kind kind) : id(id), kind(kind) {}
    FeatureDescription(unsigned id, Kind kind, unsigned contigID, int pos, unsigned coverage) :
            id(id), kind(kind), contigID(contigID), pos(pos), coverage(coverage)
    {}

    // The ID of the feature.
    unsigned id { INVALID };
    // The kind of feature.
    Kind kind { COLUMN };
    // The contig the feature is on.
    unsigned contigID { INVALID };
    // Position on the contig, if any (kind == COLUMN).
    int pos { 0 };
    // Number of reads covering this feature.
    unsigned coverage { 0 };

    // Lexicographic sorting.
    bool operator<(FeatureDescription const & other) const
    { return makeTuple() < other.makeTuple(); }

private:
    // Helper to create tuple for linear sorting/comparison.
    std::tuple<int, unsigned, int, unsigned> makeTuple() const
    { return std::make_tuple((int)kind, contigID, pos, coverage); }
};

std::ostream & operator<<(std::ostream & out, FeatureDescription const & desc);

// ----------------------------------------------------------------------------
// Class FeatureMap
// ----------------------------------------------------------------------------

class FeatureMap
{
public:
    typedef std::vector<FeatureDescription>::const_iterator TIterator;

    FeatureMap() = default;

    // Insert new feature into map, invalidates sorting in features and refresh() is required.
    void insert(FeatureDescription const & desc);
    // Refresh sorting in features.
    void refresh();

    // Returns features of kind COLUMN for the given contig / all contigs.
    std::pair<TIterator, TIterator> contigFeatures(unsigned contigID) const;
    std::pair<TIterator, TIterator> contigFeatures() const;
    // Returns features of kind LINK.
    std::pair<TIterator, TIterator> linkFeatures() const;

    // Returns number of features of the given kind.
    size_t numFeatures(FeatureDescription::Kind kind) const
    {
        auto range = (kind == FeatureDescription::COLUMN) ? contigFeatures() : linkFeatures();
        return (range.second - range.first);
    }

    // Returns total size.
    size_t size() const { return features.size(); }

    // The usual container iterators.
    TIterator begin() const { return features.begin(); }
    TIterator end() const { return features.end(); }

    // Apply old-to-new mapping.  In case of value being (unsigned)-1, entry is discarded.
    void applyContigIDMapping(std::vector<unsigned> const & oldToNew);

    // Find feature with the given ID and return end() in case of failure.
    TIterator find(unsigned featureID) const;

    void print(std::ostream & out) const;

    void clear() { features.clear(); refresh(); }

private:
    // Whether features requires sorting.
    bool needsRefresh { false };
    // Map from feature ID to feature value.
    std::map<unsigned, unsigned> map;
    // The features.
    std::vector<FeatureDescription> features;
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

}  // namespace rep_sep

#endif  // #ifndef SANDBOX_HERBARIUM_APPS_ANISE_REP_SEP_FEATURE_MAP_H_
