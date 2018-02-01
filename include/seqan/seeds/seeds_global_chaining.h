// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
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
// Global chaining algorithms.
// ==========================================================================

#ifndef SEQAN_SEEDS_SEEDS_GLOBAL_CHAINING_H_
#define SEQAN_SEEDS_SEEDS_GLOBAL_CHAINING_H_

#include <map>

namespace seqan {

// ===========================================================================
// Enums, Tags, Classes, Specializations
// ===========================================================================

/*!
 * @defgroup GlobalChainingTags Global Chaining Tags
 * @brief Tags for specifying global chaining method.
 *
 * @tag GlobalChainingTags#SparseChaining
 * @brief Tag for global sparse chaining.
 */

struct SparseChaining_;
typedef Tag<SparseChaining_> SparseChaining;

// ===========================================================================
// Metafunctions
// ===========================================================================

// ===========================================================================
// Functions
// ===========================================================================

template <typename TIntermediate>
inline bool _checkScoreInvariant(TIntermediate const & list)
{
    if (list.empty())
        return true;
    auto it = list.begin();
    auto score = it->i2;
    ++it;
    for (; it != list.end(); ++it)
    {
        if (score > it->i2)
            return false;
        score = it->i2;
    }
    return true;
}

/*!
 * @fn chainSeedsGlobally
 * @headerfile <seqan/seeds.h>
 * @brief Global chaining of seeds.
 *
 * @signature void chainSeedsGlobally(target, seedSet, tag);
 *
 * @param[out] target  A @link ContainerConcept container @endlink to append the seeds to.
 * @param[in]  seedSet The @link SeedSet @endlink object to get the seeds from.
 * @param[in]  tag     The tag to select the algorithm with (currently only @link GlobalChainingTags#SparseChaining
 *                     SparseChaining @endlink is supported).
 *
 * Chaining of seeds between two sequences can be performed using sparse chaining as defined in (Gusfield, 1997).
 *
 * @section Example
 *
 * The following example demonstrates how to use the <tt>chainSeedsGlobally()</tt> function.  First, a @link SeedSet
 * @endlink is built and filled with @link SimpleSeed @endlink object.  Then, a @link String @endlink of @link
 * SimpleSeed @endlink objects is defined and filled using the <tt>chainSeedsGlobally()</tt> function.
 *
 * @include demos/dox/seeds/seeds_chaining.cpp
 *
 * The output is as follows.  Only the first and last seeds are written to <tt>std::cout</tt>.
 *
 * @include demos/dox/seeds/seeds_chaining.cpp.stdout
 *
 * @section References
 *
 * <ul>
 * <li>Dan Gusfield. Algorithms on Strings, Trees, and Sequences:  Computer Science and Computational Biology.
 *     Cambridge University Press, January 1997.</li>
 * </ul>
 */

// TODO(holtgrew): Implement scored!
// NOTE(rrahn): Some general notes regarding the sparse chaining algorithm by Gusfield.
// From the text it is hard to follow the correct algorithm, as one usually expects the y-coordinates to
// be increasing integers (sequence positions). However, the algorithm assumes the coordinates of the
// rectangles to be placed in the fourth quadrant of the cartesian coordinate system (negative y-coordinates).
// To adapt the algorithm for the positive sequence space we sort the y-coordinates in *L* in ascending order.
// We then can use lower_bound and upper_bound on the sorted set to find the corresponding seed j as described in
// the algorithm.
template <typename TTargetContainer, typename TSeed, typename TSeedSetSpec>
void
chainSeedsGlobally(
        TTargetContainer & target,
        SeedSet<TSeed, TSeedSetSpec> const & seedSet,
        SparseChaining const &)
{
    typedef typename Position<TSeed>::Type TPosition;
    typedef typename Size<TSeed>::Type TSize;

    // -----------------------------------------------------------------------
    // Preparation
    // -----------------------------------------------------------------------

    // We copy over the seeds from the seed set into an array of seeds.  We can then directly reference seed by their
    // index in this array which is simpler than handling iterators into the std::set<> of the seed set.
    String<TSeed> seeds;
    resize(seeds, length(seedSet));
    std::copy(seedSet._seeds.begin(), seedSet._seeds.end(), begin(seeds, Standard()));

    // -----------------------------------------------------------------------
    // Step 1: Generate the sorted list of interval points.
    // -----------------------------------------------------------------------

    // This list is I in Gusfield's description.  An interval point is a triple of (dimension 0 border value, is start
    // point, pointer to seed it belongs to).
    typedef Triple<TPosition, bool, unsigned> TIntervalPoint;
    typedef String<TIntervalPoint> TIntervalPoints;
    typedef typename Iterator<TIntervalPoints, Standard>::Type TIntervalPointsIterator;

    TIntervalPoints intervalPoints;
    // std::cout << ",--- high quality seeds" << std::endl;
    // TODO(holtgrew): Using something dense here for qualities and predecessors should be faster but requires seeds to be available by consecutive ids.
    std::map<unsigned, TSize> qualityOfChainEndingIn;
    std::map<unsigned, unsigned> predecessor;
    for (unsigned i = 0; i < length(seeds); ++i)
    {
        // Since we use gap space, we have to use "false" for end points so the lexical ordering gives us what the
        // sparse chaining algorithm expects.
        // std::cout << "| " << **it << std::endl;
        qualityOfChainEndingIn[i] = seedSize(seeds[i]);
        predecessor[i] = std::numeric_limits<unsigned>::max();
        appendValue(intervalPoints, TIntervalPoint(beginPositionH(seeds[i]), true, i));
        appendValue(intervalPoints, TIntervalPoint(endPositionH(seeds[i]), false, i));
    }
    // std::cout << "`--" << std::endl;
    std::sort(begin(intervalPoints, Standard()), end(intervalPoints, Standard()));
    // std::cout << ",--- interval points" << std::endl;
    // for (unsigned i = 0; i < length(intervalPoints); ++i) {
    //     std::cout << "| (" << intervalPoints[i].i1 << ", " << intervalPoints[i].i2 << ", " << intervalPoints[i].i3 << ")" << std::endl;
    // }
    // std::cout << "`---" << std::endl;

    // -----------------------------------------------------------------------
    // Step 2: Build the chain.
    // -----------------------------------------------------------------------
    // We build a list of "intermediate solutions" (referred to as *L* in the Gusfield book).
    // Each such solution is represented by the triple (end position in dim1,
    // value of best chain so far, last seed of the chain).
    typedef Triple<TPosition, TSize, unsigned> TIntermediateSolution;
    typedef std::multiset<TIntermediateSolution> TIntermediateSolutions;
    typedef typename TIntermediateSolutions::iterator TIntermediateSolutionsIterator;

    // For all interval points...
    TIntermediateSolutions intermediateSolutions;
    for (TIntervalPointsIterator it_k = begin(intervalPoints), itEnd = end(intervalPoints); it_k != itEnd; ++it_k) {
        // The seed belonging ot the interval point is seed k.
        TSeed const & seed_k = seeds[it_k->i3];

        // std::cout << "Processing interval point (" << it->i1 << ", " << it->i2 << ", " << it->i3 << ")" << std::endl;
        if (it_k->i2) {  // Is is begin point.
            // Find the closest seed j (in y-dimension) with an
            // entry in L whose end coordinate is less or equal the begin position of k.
            //
            // STL gives us upper_bound which returns a pointer to the
            // *first* one that compares greater than the reference
            // one.  Searching for the this one and decrementing the
            // result iterator gives the desired result.
            TIntermediateSolution referenceSolution(beginPositionV(seed_k), std::numeric_limits<TSize>::max(), std::numeric_limits<unsigned>::max());
            // std::cout << "    intermediateSolutions.upper_bound(" << beginPositionV(seed_k) << ")" << std::endl;
            TIntermediateSolutionsIterator it_j = intermediateSolutions.upper_bound(referenceSolution);

            // Special case not dealt with in the book: If L is empty or there is no chain
            // that ends before k begins, simply continue with the next point in I.
            if (intermediateSolutions.empty() || it_j == intermediateSolutions.begin())
            {
                continue;
            }
            // Go to the last value in L, i.e. l_j <= h_k.
            --it_j;
            // std::cout << "     --> " << seeds[it_j->i3] << std::endl;
            // Now, we have found such a seed j.
            SEQAN_ASSERT_LEQ(endPositionV(seeds[it_j->i3]), beginPositionV(seed_k));
            // Update the intermediate solution value for k and set predecessor.
            qualityOfChainEndingIn[it_k->i3] += it_j->i2;
            // std::cout << "  UPDATE qualityOfChainEndingIn[" << it->i3 << "] == " << qualityOfChainEndingIn[it->i3] << std::endl;
            predecessor[it_k->i3] = it_j->i3;
            // std::cout << "         predecessor[" << it->i3 << "] == " << itJ->i3 << std::endl;
        } else {  // Is end point.
            // Search for the first triple j in L with l_j >= l_k.
            // Or to put it in differently, find the first chain that ends
            // left and below the chain that ends in k. These are other possible solutions,
            // that either result in a better score or must be deleted as the score of the chain ending in k
            // is bigger. Hence, every other seed that could connect to both j and k, prefers the one with the higher
            // score.
            //
            // We can use the lower_bound, which gives the first triple j such that l_j >= l_k
            SEQAN_ASSERT_GT(endPositionV(seed_k), 0u);
            TIntermediateSolution referenceSolution(endPositionV(seed_k), 0, std::numeric_limits<unsigned>::max());
            TIntermediateSolutionsIterator it_j = intermediateSolutions.lower_bound(referenceSolution);

            // If there was a valid solution in L...
            if (it_j != intermediateSolutions.end())
            {
                // Found this intermediate solution.
                SEQAN_ASSERT_GEQ(it_j->i1, endPositionV(seed_k));
                TSeed const & seed_j = seeds[it_j->i3];
                // ... start a new chain at k if the vertical end of k is
                // above the vertical end of the chain ending in j or if
                // both k and j end at the same vertical position, while
                // the score of the chain ending in k is bigger than the
                // score of the chain ending in j.
                if (endPositionV(seed_j) > endPositionV(seed_k) ||
                    (endPositionV(seed_j) == endPositionV(seed_k) && qualityOfChainEndingIn[it_k->i3] > it_j->i2))
                {
                    TIntermediateSolution sol(endPositionV(seed_k), qualityOfChainEndingIn[it_k->i3], it_k->i3);
                    // std::cout << "  INSERT (" << sol.i1 << ", " << sol.i2 << ", " << sol.i3 << ")" << __LINE__  << std::endl;
                    intermediateSolutions.insert(sol);

                    // Delete all intermediate solutions where end1 >= end1 of k and have a lower score than k
                    // to ensure that the invariant of V(j) >= V(j'), with j' <= j holds.
                    // Roughly then, there is no chain ending in a seed below the seed_k, that has a lower score
                    // than the chain ending in seed_k. Thus the last value in `intermediateSolutions` will
                    // always point to the optimal chain.
                    TIntermediateSolutionsIterator itDel = intermediateSolutions.upper_bound(referenceSolution);
                    TIntermediateSolutionsIterator itDelEnd = intermediateSolutions.end();
                    while (itDel != itDelEnd)
                    {
                        TIntermediateSolutionsIterator ptr = itDel;
                        ++itDel;
                        if (qualityOfChainEndingIn[it_k->i3] > ptr->i2)
                        {
                            // std::cout << "  ERASE (" << ptr->i1 << ", " << ptr->i2 << ", " << ptr->i3 << ")" << std::endl;
                            intermediateSolutions.erase(ptr);
                        }
                    }
                }
            } // ... otherwise, add a triple for k in L if either L is empty or the last triple in
              // L has a lower score than the chain ending in k.
            else if (intermediateSolutions.empty() || (--it_j)->i2 < qualityOfChainEndingIn[it_k->i3])
            {
                // None found.  Insert a new triple for seed k.
                TIntermediateSolution sol(endPositionV(seed_k), qualityOfChainEndingIn[it_k->i3], it_k->i3);
                // std::cout << "  INSERT (" << sol.i1 << ", " << sol.i2 << ", " << sol.i3 << ") " << __LINE__ << std::endl;
                intermediateSolutions.insert(sol);
            }
            // Check if the invariant holds, that the scores in L are in a non-decreasing order.
            SEQAN_ASSERT(_checkScoreInvariant(intermediateSolutions));
        }
    }

    // std::cout << "Maximal quality: " << intermediateSolutions.rbegin()->i2 << std::endl;

    // -----------------------------------------------------------------------
    // Step 3: Write out the resulting chain.
    // -----------------------------------------------------------------------
    // TODO(holtgrew): We could use two different algorithms for target containers that are strings and those that are lists.
    clear(target);
    unsigned next = intermediateSolutions.rbegin()->i3;
    while (next != std::numeric_limits<unsigned>::max())
    {
        appendValue(target, seeds[next]);
        next = predecessor[next];
    }
    reverse(target);

    // Assert that the resulting chain is non-overlapping.
    #if SEQAN_ENABLE_DEBUG
    if (length(target) > 0u) {
        typedef typename Iterator<TTargetContainer, Standard>::Type TIterator;
        // std::cerr << ".-- Chain (" << __FILE__ << ":" << __LINE__ << "):" << std::endl;
        // for (TIterator it = begin(target, Standard()); it != end(target, Standard()); ++it)
        //     std::cerr << "| " << *it << std::endl;
        // std::cerr << "`--" << std::endl;
        TIterator itPrevious = begin(target, Standard());
        TIterator it = itPrevious;
        TIterator itEnd = end(target, Standard());
        // std::cout << *it << std::endl;
        ++it;
        for (; it != itEnd; ++it) {
            // std::cout << *it << std::endl;
            SEQAN_ASSERT_LEQ(endPositionH(*itPrevious), beginPositionH(*it));
            SEQAN_ASSERT_LEQ(endPositionV(*itPrevious), beginPositionV(*it));
            itPrevious = it;
        }
    }
    #endif  // #if SEQAN_ENABLE_DEBUG
}

}  // namespace seqan

#endif  // #ifndef SEQAN_SEEDS_SEEDS_GLOBAL_CHAINING_H_
