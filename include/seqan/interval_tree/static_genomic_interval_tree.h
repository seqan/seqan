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
// Static, immutable version of interval tree.
// ==========================================================================

// TODO(holtgrewe): Implement neighbor search in case of no overlapping intervals.

#ifndef INCLUDE_SEQAN_INTERVAL_TREE_STATIC_GENOMIC_INTERVAL_TREE_H_
#define INCLUDE_SEQAN_INTERVAL_TREE_STATIC_GENOMIC_INTERVAL_TREE_H_

#include <algorithm>
#include <vector>
#include <utility>

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include "base_genomic_interval_tree.h"
#include "interval_tree_entry.h"
#include "static_interval_tree.h"
#include "static_genomic_interval_tree_iterator.h"

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class StaticGenomicGenomicIntervalTree
// ----------------------------------------------------------------------------

// Tags for specializing GenomicIntervalTree<>.

struct Static_;
typedef Tag<Static_> Static;

/*!
 * @class StaticGenomicGenomicIntervalTree
 * @headerfile <seqan/interval_tree.h>
 * @brief Static @link GenomicIntervalTree @endlink implementation.
 *
 * @signature template <typename TValue>
 *            class GenomicIntervalTree<TValue, Static>;
 *
 * @tparam TValue Value to store in the tree.
 * @tparam TSpec  Specializing tag.
 *
 * The implementation is based on Cormen's augmented binary search trees as described in <a
 * href="http://en.wikipedia.org/wiki/Interval_tree#Augmented_tree">this Wikipedia article</a>.
 */
template <typename TValue>
class GenomicIntervalTree<TValue, Static>
{
public:
    GenomicIntervalTree() : totalSize(0) {}

    template <typename TContainer>
    GenomicIntervalTree(TContainer const & container)
    {
        initialize(container);
    }

    std::vector<IntervalTree<TValue, Static> > contigTrees;
    size_t totalSize;

    void updateTotalSize()
    {
        totalSize = 0;
        for (unsigned i = 0; i < contigTrees.size(); ++i)
            totalSize += length(contigTrees[i]);
    }

    template <typename TContainer> void initialize(TContainer const & container);
};

template <typename TValue>
template <typename TContainer>
void GenomicIntervalTree<TValue, Static>::initialize(TContainer const & container)
{
    contigTrees.resize(length(container));
    for (unsigned i = 0; i < length(container); ++i)
        contigTrees[i].initialize(container[i]);

    updateTotalSize();
}

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Size
// ----------------------------------------------------------------------------

/**
 * @mfn GenomicIntervalTree#Size
 * @brief The size type of the <tt>GenomicIntervalTree</tt>.
 */

template <typename TValue>
struct Size<GenomicIntervalTree<TValue, Static> >
{
    typedef typename std::vector<TValue>::size_type Type;
};

template <typename TValue>
struct Size<GenomicIntervalTree<TValue, Static> const>
{
    typedef typename std::vector<TValue>::size_type Type;
};

// ----------------------------------------------------------------------------
// Metafunction Iterator
// ----------------------------------------------------------------------------

/**
 * @mfn GenomicIntervalTree#Iterator
 * @brief The iterator type of the <tt>GenomicIntervalTree</tt>.
 */

template <typename TValue, typename TIterSpec>
struct Iterator<GenomicIntervalTree<TValue, Static>, TIterSpec>
{
    typedef Iter<GenomicIntervalTree<TValue, Static>, void> Type;
};

template <typename TValue, typename TIterSpec>
struct Iterator<GenomicIntervalTree<TValue, Static> const, TIterSpec>
{
    typedef Iter<GenomicIntervalTree<TValue, Static> const, void> Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function build()
// ----------------------------------------------------------------------------

template <typename TValue, typename TPos, typename TContainer>
void build(GenomicIntervalTree<TValue, Static> & tree,
           TContainer const & values)
{
    tree.initialize(values);
    tree.updateTotalSize();
}

template <typename TValue, typename TPos, typename TContainer>
void build(GenomicIntervalTree<TValue, Static> & tree,
           TContainer const & values,
           int rID)
{
    if (rID > (int)tree.contigTrees.size())
        tree.contigTrees[rID].initialize(values);
    tree.updateTotalSize();
}

// ----------------------------------------------------------------------------
// Function length()
// ----------------------------------------------------------------------------

// TODO(holtgrewe): Get rid of non-const length once all such functions are removed
template <typename TValue>
typename Size<GenomicIntervalTree<TValue, Static> >::Type
length(GenomicIntervalTree<TValue, Static> const & tree)
{
    return tree.totalSize;
}

template <typename TValue>
typename Size<GenomicIntervalTree<TValue, Static> >::Type
length(GenomicIntervalTree<TValue, Static> & tree)
{
    return tree.totalSize;
}

// ----------------------------------------------------------------------------
// Function begin()
// ----------------------------------------------------------------------------

// TODO(holtgrewe): For some reason, using TIterSpec template parameter instead of Standard/Rooted created ambiguous function calls and references to missing class specialization of IteratorDefaultImp_, same below for end().

template <typename TValue>
typename Iterator<GenomicIntervalTree<TValue, Static> const, Standard>::Type
begin(GenomicIntervalTree<TValue, Static> const & tree, Standard const & /*tag*/)
{
    return typename Iterator<GenomicIntervalTree<TValue, Static> const, Standard>::Type(tree, true);
}

template <typename TValue>
typename Iterator<GenomicIntervalTree<TValue, Static>, Standard>::Type
begin(GenomicIntervalTree<TValue, Static> & tree, Standard const & /*tag*/)
{
    return typename Iterator<GenomicIntervalTree<TValue, Static>, Standard>::Type(tree, true);
}

template <typename TValue>
typename Iterator<GenomicIntervalTree<TValue, Static> const, Rooted>::Type
begin(GenomicIntervalTree<TValue, Static> const & tree, Rooted const & /*tag*/)
{
    return typename Iterator<GenomicIntervalTree<TValue, Static> const, Standard>::Type(tree, true);
}

template <typename TValue>
typename Iterator<GenomicIntervalTree<TValue, Static>, Rooted>::Type
begin(GenomicIntervalTree<TValue, Static> & tree, Rooted const & /*tag*/)
{
    return typename Iterator<GenomicIntervalTree<TValue, Static>, Standard>::Type(tree, true);
}

// ----------------------------------------------------------------------------
// Function end()
// ----------------------------------------------------------------------------

template <typename TValue>
typename Iterator<GenomicIntervalTree<TValue, Static> const, Standard>::Type
end(GenomicIntervalTree<TValue, Static> const & tree, Standard const & /*tag*/)
{
    return typename Iterator<GenomicIntervalTree<TValue, Static> const, Standard>::Type(tree, false);
}

template <typename TValue>
typename Iterator<GenomicIntervalTree<TValue, Static>, Standard>::Type
end(GenomicIntervalTree<TValue, Static> & tree, Standard const & /*tag*/)
{
    return typename Iterator<GenomicIntervalTree<TValue, Static>, Standard>::Type(tree, false);
}

template <typename TValue>
typename Iterator<GenomicIntervalTree<TValue, Static> const, Rooted>::Type
end(GenomicIntervalTree<TValue, Static> const & tree, Rooted const & /*tag*/)
{
    return typename Iterator<GenomicIntervalTree<TValue, Static> const, Standard>::Type(tree, false);
}

template <typename TValue>
typename Iterator<GenomicIntervalTree<TValue, Static>, Rooted>::Type
end(GenomicIntervalTree<TValue, Static> & tree, Rooted const & /*tag*/)
{
    return typename Iterator<GenomicIntervalTree<TValue, Static>, Standard>::Type(tree, false);
}

// ----------------------------------------------------------------------------
// Function findOverlappingWithPoint()
// ----------------------------------------------------------------------------

// Guarantee: result is sorted by (beginPos, endPos), by search algorithm.

template <typename TValue, typename TPos, typename TResult>
void findOverlappingWithPoint(GenomicIntervalTree<TValue, Static> const & tree,
                              int rID,
                              TPos point,
                              TResult & result)
{
    typedef GenomicIntervalTree<TValue, Static> const                     TGenomicIntervalTree;
    typedef typename Position<TGenomicIntervalTree>::Type                 TPos2;
    typedef typename Size<TGenomicIntervalTree>::Type                     TSize;
    typedef typename Iterator<TGenomicIntervalTree const, Standard>::Type TIter;

    typedef typename Iterator<IntervalTree<TValue, Static> const, Standard>::Type TTreeIter;
    std::vector<TTreeIter> tmpResult;

    findOverlappingWithPoint(tree.contigTrees[rID], point, tmpResult);

    result.clear();
    result.reserve(tmpResult.size());
    for (unsigned i = 0; i < tmpResult.size(); ++i)
        result.push_back(TIter(tree, rID, tmpResult[i]));
}

// ----------------------------------------------------------------------------
// Function findOverlappingWithInterval()
// ----------------------------------------------------------------------------

// Guarantee: result is sorted by (beginPos, endPos), by search algorithm.

template <typename TValue, typename TPos, typename TResult>
void findOverlappingWithInterval(GenomicIntervalTree<TValue, Static> const & tree,
                                 int rID,
                                 TPos posBegin,
                                 TPos posEnd,
                                 TResult & result)
{
    typedef GenomicIntervalTree<TValue, Static> const                     TGenomicIntervalTree;
    typedef typename Position<TGenomicIntervalTree>::Type                 TPos2;
    typedef typename Size<TGenomicIntervalTree>::Type                     TSize;
    typedef typename Iterator<TGenomicIntervalTree const, Standard>::Type TIter;

    typedef typename Iterator<IntervalTree<TValue, Static> const, Standard>::Type TTreeIter;
    std::vector<TTreeIter> tmpResult;

    findOverlappingWithInterval(tree.contigTrees[rID], posBegin, posEnd, tmpResult);

    result.clear();
    result.reserve(tmpResult.size());
    for (unsigned i = 0; i < tmpResult.size(); ++i)
        result.push_back(TIter(tree, rID, tmpResult[i]));
}

}  // namespace seqan

#endif  // INCLUDE_SEQAN_INTERVAL_TREE_STATIC_GENOMIC_INTERVAL_TREE_H_
