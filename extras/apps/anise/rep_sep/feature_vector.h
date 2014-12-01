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

#ifndef SANDBOX_HERBARIUM_APPS_ANISE_REP_SEP_FEATURE_VECTOR_H_
#define SANDBOX_HERBARIUM_APPS_ANISE_REP_SEP_FEATURE_VECTOR_H_

#include <iosfwd>
#include <utility>

#include <vector>

namespace rep_sep {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class Feature
// ----------------------------------------------------------------------------

// Represent one instance (value) of a feature with a given ID.

struct Feature
{
    static const int INVALID = -1;
    int id    { INVALID };
    int value { INVALID };

    Feature() = default;

    Feature(int id, int value) : id(id), value(value)
    {}

    bool operator==(Feature const & other) const
    {
        return (makeTuple() == other.makeTuple());
    }

    bool operator<(Feature const & other) const
    {
        return (makeTuple() < other.makeTuple());
    }

    void print(std::ostream & out) const;

private:

    std::pair<int, int> makeTuple() const
    {
        return std::make_pair(id, value);
    }
};

inline std::ostream & operator<<(std::ostream & out, Feature const & feature)
{
    feature.print(out);
    return out;
}

// ----------------------------------------------------------------------------
// Class FeatureVector
// ----------------------------------------------------------------------------

// A list of features, sorted by feature id.  Allows O(log n) query for feature value, conflict checks with other
// FeatureVector objects, adding of features, and merging with other vectors.

class FeatureVector
{
public:
    typedef std::vector<Feature>::const_iterator TIterator;

    FeatureVector() = default;

    TIterator find(int featureID) const;
    size_t count(int featureID) const { return (find(featureID) != end()); }

    void insert(Feature const & feature);
    TIterator begin() const;
    TIterator end() const;

    bool conflicts(FeatureVector const & other) const;

    size_t size() const { return features.size(); }
    bool empty() const { return features.empty(); }

    // Non-destructive and destructive merging.
    FeatureVector merge(FeatureVector const & other) const;
    FeatureVector & mergeWithThis(FeatureVector const & other);

    bool operator==(FeatureVector const & other) const
    {
        return (features == other.features);
    }

    void print(std::ostream & out) const;

    void swap(FeatureVector & other)
    {
        using std::swap;
        swap(features, other.features);
    }

private:

    // The features.
    std::vector<Feature> features;
};

inline std::ostream & operator<<(std::ostream & out, FeatureVector const & vec)
{
    vec.print(out);
    return out;
}

inline void swap(FeatureVector & lhs, FeatureVector & rhs)
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

#endif  // #ifndef SANDBOX_HERBARIUM_APPS_ANISE_REP_SEP_FEATURE_VECTOR_H_
