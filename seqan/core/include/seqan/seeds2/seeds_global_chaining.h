// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
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

/**
.Tag.Global Chaining
..summary:Tags for selecting the global chaining algorithm.
..cat:Seed Handling
..see:Function.chainSeedsGlobally
..tag:SparseChaining:
    Chaining as described in (Gusfield, 1997) section 13.3.
..include:seqan/seeds2.h
 */
struct SparseChaining_;
typedef Tag<SparseChaining_> SparseChaining;

// ===========================================================================
// Metafunctions
// ===========================================================================

// ===========================================================================
// Functions
// ===========================================================================

// TODO(holtgrew): Implement scored!
/**
.Function.chainSeedsGlobally
..summary:Global chaining of seeds.
..signature:chainSeedsGlobally(target, seedSet, tag)
..include:seqan/seeds2.h
*/
template <typename TTargetContainer, typename TSeedSpec, typename TSeedSetSpec, typename TSeedConfig>
void
chainSeedsGlobally(
        TTargetContainer & target,
        SeedSet<TSeedSpec, TSeedSetSpec, TSeedConfig> const & seedSet,
        SparseChaining const &)
{
    SEQAN_CHECKPOINT;

    typedef SeedSet<TSeedSpec, TSeedSetSpec, TSeedConfig> TSeedSet;
    typedef typename Value<TSeedSet>::Type TSeed;
    typedef typename Position<TSeed>::Type TPosition;
    typedef typename Size<TSeed>::Type TSize;
    typedef typename TSeedSet::THighQualitySeeds const THighQualitySeeds;
    typedef typename THighQualitySeeds::const_iterator THighQualitySeedsIterator;

    // -----------------------------------------------------------------------
    // Step 1: Generate the sorted list of interval points.
    // -----------------------------------------------------------------------
    // This list is I in Gusfield's description.  An interval point is
    // a triple of (dimension 0 border value, is start point, pointer
    // to seed it belongs to).
    typedef Triple<TPosition, bool, TSeed *> TIntervalPoint;
    typedef String<TIntervalPoint> TIntervalPoints;
    typedef typename Iterator<TIntervalPoints, Standard>::Type TIntervalPointsIterator;

    TIntervalPoints intervalPoints;
    // std::cout << ",--- high quality seeds" << std::endl;
    // TODO(holtgrew): Using something dense here for qualities and predecessors should be faster but requires seeds to be available by consecutive ids.
    std::map<TSeed *, TSize> qualityOfChainEndingIn;
    std::map<TSeed *, TSeed *> predecessor;
    for (THighQualitySeedsIterator it = seedSet._highQualitySeeds.begin(), itEnd = seedSet._highQualitySeeds.end(); it != itEnd; ++it) {
        // Since we use gap space, we have to use "false" for end
        // points so the lexical ordering gives us what the sparse
        // chaining algorithm expects.
        // std::cout << "| " << **it << std::endl;
        qualityOfChainEndingIn[*it] = getSeedSize(**it);
        appendValue(intervalPoints, TIntervalPoint(getBeginDim0(**it), true, *it));
        appendValue(intervalPoints, TIntervalPoint(getEndDim0(**it), false, *it));
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
    // We build a list of "intermediate solutions".  Each such
    // solution is represented by the triple (end position in dim1,
    // value of best chain so far, last seed of the chain).
    typedef Triple<TPosition, TSize, TSeed *> TIntermediateSolution;
    typedef std::multiset<TIntermediateSolution> TIntermediateSolutions;
    typedef typename TIntermediateSolutions::iterator TIntermediateSolutionsIterator;

    // For all interval points...
    TIntermediateSolutions intermediateSolutions;
    for (TIntervalPointsIterator it = begin(intervalPoints), itEnd = end(intervalPoints); it != itEnd; ++it) {
        // The seed belonging ot the interval point is seed k.
        TSeed const & seedK = value(value(it).i3);

        // std::cout << "Processing interval point (" << value(it).i1 << ", " << value(it).i2 << ", " << value(it).i3 << ")" << std::endl;
        if (value(it).i2) {  // Is is begin point.
            // Find the closest seed (in dimension 1) to seed k with an
            // entry in intermediateSolutions whose end coordinate in
            // dimension 1 is <= the begin coordinate in dimension 1
            // of seedK.
            //
            // STL gives us upper_bound which returns a pointer to the
            // *first* one that compares greater than the reference
            // one.  Searching for the this one and decrementing the
            // result iterator gives the desired result.
            TIntermediateSolution referenceSolution(getBeginDim1(seedK), 0, 0);
            // std::cout << "    intermediateSolutions.upper_bound(" << getBeginDim1(seedK) << ")" << std::endl;
            TIntermediateSolutionsIterator itJ = intermediateSolutions.upper_bound(referenceSolution);
            if (itJ == intermediateSolutions.begin()) {
                if (intermediateSolutions.size() > 0 &&
                    intermediateSolutions.rbegin()->i1 <= getBeginDim1(seedK)) {
                    itJ = intermediateSolutions.end();
                    --itJ;
                } else {
                    continue;
                }
            } else {
                SEQAN_ASSERT_GT(intermediateSolutions.size(), 0u);  // TODO(holtgrew): Remove this assertion?
                --itJ;
            }
            // std::cout << "     --> " << value(itJ->i3) << std::endl;
            // Now, we have found such a seed j.
            SEQAN_ASSERT_LEQ(getEndDim1(value(itJ->i3)), getEndDim1(seedK));
            // Update the intermediate solution value for k and set predecessor.
            qualityOfChainEndingIn[value(it).i3] += itJ->i2;
            // std::cout << "  UPDATE qualityOfChainEndingIn[" << value(it).i3 << "] == " << qualityOfChainEndingIn[value(it).i3] << std::endl;
            predecessor[value(it).i3] = itJ->i3;
            // std::cout << "         predecessor[" << value(it).i3 << "] == " << itJ->i3 << std::endl;
        } else {  // Is end point.
            // Search for the first triple in intermediateSolutions
            // where the end coordinate in dimension 1 is >= end
            // coordinate in dimension 1 for seed k.  The corresponding
            // seed is seed j.
            //
            // We work with upper_bound here which gives us the first
            // value that is > so we have to work around this to get
            // >= again...
            SEQAN_ASSERT_GT(getEndDim1(seedK), 0u);
            TIntermediateSolution referenceSolution(getEndDim1(seedK), 0, 0);
            TIntermediateSolutionsIterator itSol = intermediateSolutions.upper_bound(referenceSolution);
            if (itSol == intermediateSolutions.end()) {
                // None found.  Insert a new triple for seed k.
                TIntermediateSolution sol(getEndDim1(seedK), qualityOfChainEndingIn[value(it).i3], value(it).i3);
                // std::cout << "  INSERT (" << sol.i1 << ", " << sol.i2 << ", " << sol.i3 << ") " << __LINE__ << std::endl;
                intermediateSolutions.insert(sol);
            } else {
                // Found this intermediate solution.
                SEQAN_ASSERT_GEQ(itSol->i1, getEndDim1(seedK));
                TSeed const & seedJ = value(itSol->i3);
                // Possibly start a new chain at k if the end1 is
                // before the end1 of the chain ending in j or they
                // end at the same coordinate in dim1 but k already
                // has a higher quality than the whole chaing ending
                // at j.
                if (getEndDim1(seedJ) > getEndDim1(seedK) ||
                    (getEndDim1(seedJ) == getEndDim1(seedK) && qualityOfChainEndingIn[value(it).i3] > itSol->i2)) {
                    TIntermediateSolution sol(getEndDim1(seedK), qualityOfChainEndingIn[value(it).i3], value(it).i3);
                    // std::cout << "  INSERT (" << sol.i1 << ", " << sol.i2 << ", " << sol.i3 << ")" << __LINE__  << std::endl;
                    intermediateSolutions.insert(sol);
                }
            }
            
            // Delete all intermediate solutions where end1 >= end1 of k and have a lower quality than k.
            TIntermediateSolutionsIterator itDel = intermediateSolutions.upper_bound(referenceSolution);
            TIntermediateSolutionsIterator itDelEnd = intermediateSolutions.end();
            while (itDel != itDelEnd) {
                TIntermediateSolutionsIterator ptr = itDel;
                ++itDel;
                if (qualityOfChainEndingIn[value(it).i3] > ptr->i2) {
                    // std::cout << "  ERASE (" << ptr->i1 << ", " << ptr->i2 << ", " << ptr->i3 << ")" << std::endl;
                    intermediateSolutions.erase(ptr);
                }
            }
        }
    }

    // std::cout << "Maximal quality: " << intermediateSolutions.rbegin()->i2 << std::endl;

    // -----------------------------------------------------------------------
    // Step 3: Write out the resulting chain.
    // -----------------------------------------------------------------------
    // TODO(holtgrew): We could use two different algorithms for target containers that are strings and those that are lists.
    clear(target);
    TSeed *next = intermediateSolutions.rbegin()->i3;
    while (next != static_cast<TSeed *>(0)) {
        appendValue(target, *next);
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
            SEQAN_ASSERT_LEQ(getEndDim0(*itPrevious), getBeginDim0(*it));
            SEQAN_ASSERT_LEQ(getEndDim1(*itPrevious), getBeginDim1(*it));
            itPrevious = it;
        }
    }
    #endif  // #if SEQAN_ENABLE_DEBUG
}

}  // namespace seqan

#endif  // #ifndef SEQAN_SEEDS_SEEDS_GLOBAL_CHAINING_H_
