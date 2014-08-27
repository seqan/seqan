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

#ifndef SANDBOX_HERBARIUM_APPS_ANISE_REP_SEP_WEIGHTED_FEATURE_VECTOR_H_
#define SANDBOX_HERBARIUM_APPS_ANISE_REP_SEP_WEIGHTED_FEATURE_VECTOR_H_

#include <iosfwd>
#include <utility>
#include <vector>
#include <tuple>

#include "rep_sep/feature_vector.h"

namespace rep_sep {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class WeightedFeature
// ----------------------------------------------------------------------------

struct WeightedFeature : Feature
{
    int count { INVALID };

    WeightedFeature() = default;
    WeightedFeature(int id, int value, int count = 1) : Feature(id, value), count(count) {}
    WeightedFeature(Feature const & feature, int count = 1) : Feature(feature), count(count) {}

    bool operator==(WeightedFeature const & other) const
    {
        return (makeTuple() == other.makeTuple());
    }

    bool operator<(WeightedFeature const & other) const
    {
        return (makeTuple() < other.makeTuple());
    }

    void print(std::ostream & out) const;

private:

    std::tuple<int, int, int> makeTuple() const
    {
        return std::make_tuple(id, value, count);
    }
};

inline std::ostream & operator<<(std::ostream & out, WeightedFeature const & feature)
{
    feature.print(out);
    return out;
}

// ----------------------------------------------------------------------------
// Class WeightedFeatureVector
// ----------------------------------------------------------------------------

// A list of features, sorted by feature id.  Allows O(log n) query for feature value, conflict checks with other
// WeightedFeatureVector objects, adding of features, and merging with other vectors.

class WeightedFeatureVector
{
public:
    typedef std::vector<WeightedFeature>::const_iterator TIterator;

    WeightedFeatureVector() = default;

    TIterator find(int featureID) const;
    size_t count(int featureID) const { return (find(featureID) != end()); }

    void insert(WeightedFeature const & feature);
    TIterator begin() const;
    TIterator end() const;

    // Returns number of conflicts.
    int conflicts(WeightedFeatureVector const & other) const;

    size_t size() const { return features.size(); }
    bool empty() const { return features.empty(); }

    // Non-destructive and destructive merging.  Up to ignoreConflicts are ignored.
    WeightedFeatureVector & mergeWithThis(WeightedFeatureVector const & other,
                                          int ignoreConflicts = 0);

    bool operator==(WeightedFeatureVector const & other) const
    {
        return (features == other.features);
    }

    void print(std::ostream & out) const;

    void swap(WeightedFeatureVector & other)
    {
        using std::swap;
        swap(features, other.features);
    }

private:

    // The features.
    std::vector<WeightedFeature> features;
};

inline std::ostream & operator<<(std::ostream & out, WeightedFeatureVector const & vec)
{
    vec.print(out);
    return out;
}

inline void swap(WeightedFeatureVector & lhs, WeightedFeatureVector & rhs)
{
    lhs.swap(rhs);
}

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

}  // namespace rep_sep

#endif  // #ifndef SANDBOX_HERBARIUM_APPS_ANISE_REP_SEP_WEIGHTED_FEATURE_VECTOR_H_
