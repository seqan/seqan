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

#include "weighted_feature_vector.h"

#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include <seqan/basic.h>

#include "rep_sep/utils.h"

namespace rep_sep {

namespace {    // anonymous namespace

// Compare feature based on ID only.

bool idsLessThan(WeightedFeature const & lhs, WeightedFeature const & rhs) { return (lhs.id < rhs.id); }

// Compare features based on value only.

bool valuesEqual(WeightedFeature const & lhs, WeightedFeature const & rhs) { return (lhs.value == rhs.value); }

// Merge adjacent entries and sum up count if key/value equal.

void mergeAdjacent(std::vector<WeightedFeature> & features, int ignoreConflicts)
{
    if (features.empty())
        return;

    auto itL = features.begin();
    auto itR = std::next(itL);

    for (; itR != features.end(); ++itR)
    {
        SEQAN_ASSERT_LEQ(itL->id, itR->id);

        if (itL->id == itR->id)
        {
            if (itL->value != itR->value)
            {
                if (ignoreConflicts >= 0)
                    SEQAN_CHECK(std::min(itL->count, itR->count) <= ignoreConflicts,
                                "Too many conflicts min(...)=%d, ignoreConflicts=%d.",
                                std::min(itL->count, itR->count), ignoreConflicts);
                itL->count = std::max(itL->value, itR->value);
            }
            else
            {
                itL->count += itR->count;
            }
        }
        else
        {
            if (++itL != itR)
                *itL = *itR;
        }
    }

    features.resize(std::next(itL) - features.begin());
}

// Return number of conflicts, returns smaller count (e.g. non-consensus count).
//
// Based upon textbook std::set_union().

// TODO(holtgrew): This could use some love.

template <typename TIter1, typename TIter2, typename TComp, typename TPred>
int sumCountsOverConflicts(TIter1 first1, TIter2 last1,
                           TIter2 first2, TIter2 last2,
                           TComp comp, TPred pred)
{
    int result = 0;

    while (first1 != last1 && first2 != last2)
    {
        if (comp(*first1, *first2))
        {
            ++first1;
        }
        else
        {
            if (!comp(*first2, *first1))
            {
                if (!pred(*first1, *first2))
                    result += std::min(first1->count, first2->count);
                ++first1;
            }
            ++first2;
        }
    }

    return result;
}


}  // anonymous namespace

// --------------------------------------------------------------------------
// Class WeightedFeature
// --------------------------------------------------------------------------

void WeightedFeature::print(std::ostream & out) const
{
    out << "(id=" << id << ", value=" << value << ", count=" << count << ")";
}

// --------------------------------------------------------------------------
// Class WeightedFeatureVector
// --------------------------------------------------------------------------

WeightedFeatureVector::TIterator WeightedFeatureVector::find(int featureID) const
{
    auto lt = [](WeightedFeature const & lhs, int id) { return (lhs.id < id); };

    auto it = std::lower_bound(begin(), end(), featureID, lt);
    if (it != end() && (it->id == featureID))
        return it;
    else
        return end();
}

void WeightedFeatureVector::insert(WeightedFeature const & feature)
{
    auto it = std::lower_bound(features.begin(), features.end(), feature);
    if (it != end() && it->id == feature.id)
    {
        // WeightedFeature already exists, check that it has the same value (and ignore it) or throw an exception in case of
        // errors.
        if (it->value != feature.value)
            throw std::runtime_error("WeightedFeature already exists with different value!");
    }
    else
    {
        features.insert(it, feature);
    }
}

WeightedFeatureVector::TIterator WeightedFeatureVector::begin() const
{
    return features.begin();
}

WeightedFeatureVector::TIterator WeightedFeatureVector::end() const
{
    return features.end();
}

int WeightedFeatureVector::conflicts(WeightedFeatureVector const & other) const
{
    return sumCountsOverConflicts(begin(), end(), other.begin(), other.end(),
                                  idsLessThan, valuesEqual);
}

void WeightedFeatureVector::print(std::ostream & out) const
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

WeightedFeatureVector & WeightedFeatureVector::mergeWithThis(WeightedFeatureVector const & other,
                                                             int ignoreConflicts)
{
    if (ignoreConflicts >= 0 && conflicts(other) > ignoreConflicts)
    {
        std::stringstream ss;
        ss << "Incompatible weighted feature vectors: " << *this << " vs. " << other;
        throw std::runtime_error(ss.str());
    }

    decltype(features) tmp;
    std::merge(features.begin(), features.end(), other.begin(), other.end(),
               std::back_inserter(tmp), idsLessThan);
    mergeAdjacent(tmp, ignoreConflicts);
    using std::swap;
    swap(tmp, features);

    return *this;
}

}  // namespace rep_sep
