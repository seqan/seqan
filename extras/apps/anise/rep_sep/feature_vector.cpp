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

#include "feature_vector.h"

#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "rep_sep/utils.h"

namespace rep_sep {

namespace {    // anonymous namespace

// Compare feature based on ID only.

bool idsLessThan(Feature const & lhs, Feature const & rhs) { return (lhs.id < rhs.id); }

// Compare features based on value only.

bool valuesEqual(Feature const & lhs, Feature const & rhs) { return (lhs.value == rhs.value); }

}  // anonymous namespace

// --------------------------------------------------------------------------
// Class Feature
// --------------------------------------------------------------------------

void Feature::print(std::ostream & out) const
{
    out << "(id=" << id << ", value=" << value << ")";
}

// --------------------------------------------------------------------------
// Class FeatureVector
// --------------------------------------------------------------------------

FeatureVector::TIterator FeatureVector::find(int featureID) const
{
    auto lt = [](Feature const & lhs, int id) { return (lhs.id < id); };

    auto it = std::lower_bound(begin(), end(), featureID, lt);
    if (it != end() && (it->id == featureID))
        return it;
    else
        return end();
}

void FeatureVector::insert(Feature const & feature)
{
    auto it = std::lower_bound(features.begin(), features.end(), feature);
    if (it != end() && it->id == feature.id)
    {
        // Feature already exists, check that it has the same value (and ignore it) or throw an exception in case of
        // errors.
        if (it->value != feature.value)
            throw std::runtime_error("Feature already exists with different value!");
    }
    else
    {
        features.insert(it, feature);
    }
}

FeatureVector::TIterator FeatureVector::begin() const
{
    return features.begin();
}

FeatureVector::TIterator FeatureVector::end() const
{
    return features.end();
}

bool FeatureVector::conflicts(FeatureVector const & other) const
{
    return !predicateOverIntersection(begin(), end(), other.begin(), other.end(),
                                      idsLessThan, valuesEqual);
}

void FeatureVector::print(std::ostream & out) const
{
    out << "<";
    if (!features.empty())
        out << features[0];
    for (unsigned i = 1; i < features.size(); ++i)
    {
        out << ", ";
        out << features[i];
    }
    out << ">";
}

FeatureVector & FeatureVector::mergeWithThis(FeatureVector const & other)
{
    if (conflicts(other))
    {
        std::stringstream ss;
        ss << "Incompatible feature vectors: " << *this << " vs. " << other;
        throw std::runtime_error(ss.str());
    }

    decltype(features) tmp;
    std::set_union(features.begin(), features.end(), other.begin(), other.end(),
                   std::back_inserter(tmp), idsLessThan);
    using std::swap;
    swap(tmp, features);

    return *this;
}

FeatureVector FeatureVector::merge(FeatureVector const & other) const
{
    FeatureVector result(*this);
    result.mergeWithThis(other);
    return result;
}

}  // namespace rep_sep
