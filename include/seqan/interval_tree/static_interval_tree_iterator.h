// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
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
// Iterator for Static IntervalTree objects.
// ==========================================================================

#ifndef INCLUDE_SEQAN_INTERVAL_TREE_STATIC_INTERVAL_TREE_ITERATOR_H_
#define INCLUDE_SEQAN_INTERVAL_TREE_STATIC_INTERVAL_TREE_ITERATOR_H_

#include <vector>

#include <seqan/basic.h>
#include <seqan/sequence.h>

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

struct Static_;
typedef Tag<Static_> Static;

template <typename TValue, typename TSpec> class IntervalTree;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class Iter                                         [IntervalTree, non-const]
// ----------------------------------------------------------------------------

template <typename TValue>
class Iter<IntervalTree<TValue, Static>, void>
{
private:
    typedef typename std::vector<IntervalTreeEntry<TValue> >::iterator TIter_;

public:
    TIter_ it;

    Iter() : it() {}
    explicit Iter(TIter_ it) : it(it) {}

    bool operator==(Iter const & other)
    {
        return (it == other.it);
    }

    bool operator!=(Iter const & other)
    {
        return (it != other.it);
    }

    TValue & operator*()
    {
        return it->value;
    }

    TValue & operator*() const
    {
        return it->value;
    }

    TValue * operator->() const
    {
        return &it->value;
    }

    TValue * operator->()
    {
        return &it->value;
    }

    Iter operator++(int)
    {
        Iter copy(*this);
        operator++();
        return copy;
    }

    Iter & operator++()
    {
        ++it;
        return *this;
    }


    Iter operator--(int)
    {
        Iter copy(*this);
        operator--();
        return copy;
    }

    Iter & operator--()
    {
        --it;
        return *this;
    }
};

// ----------------------------------------------------------------------------
// Class Iter                                             [IntervalTree, const]
// ----------------------------------------------------------------------------

template <typename TValue>
class Iter<IntervalTree<TValue, Static> const, void>
{
private:
    typedef typename std::vector<IntervalTreeEntry<TValue> >::const_iterator TIter_;

public:
    TIter_ it;

    Iter() : it() {}
    explicit Iter(TIter_ it) : it(it) {}

    bool operator==(Iter const & other)
    {
        return (it == other.it);
    }

    bool operator!=(Iter const & other)
    {
        return (it != other.it);
    }

    TValue const & operator*()
    {
        return it->value;
    }

    TValue const & operator*() const
    {
        return it->value;
    }

    TValue const * operator->() const
    {
        return &it->value;
    }

    TValue const * operator->()
    {
        return &it->value;
    }

    Iter operator++(int)
    {
        Iter copy(*this);
        operator++();
        return copy;
    }

    Iter & operator++()
    {
        ++it;
        return *this;
    }


    Iter operator--(int)
    {
        Iter copy(*this);
        operator--();
        return copy;
    }

    Iter & operator--()
    {
        --it;
        return *this;
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

}  // namespace seqan

#endif  // INCLUDE_SEQAN_INTERVAL_TREE_STATIC_INTERVAL_TREE_ITERATOR_H_
