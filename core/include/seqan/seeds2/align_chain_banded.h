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
// This header defines the data structure AlignmentChain which supports
// the banded chain alignment algorithm.  It allows to store multiple
// alignments for each seed in the chain and each gap as well as the leading
// and trailing gap.
// 
// The data structure supports the consecutive alignment of seeds and
// the rectangles between them.  The functions encapsulate the
// computation of dimensions for the rectangles and seeds.  They also
// encapsulate the dimension computation for parts to be copied.
// 
// The alignment and especially the copying of the data itself is
// delegated to the alignment algorithms in
// align_dynprog[_banded]_{linear,affine}.h.
// 
// The functions should not have to be adapted to accomodate less
// copying since they are not concerned with the copying themselves.
// 
// bandedChainAlignment() is the public interface to the banded chain
// alignment which calls _bandedChainAlignment() with the appropriate
// alignment algorithm tag (i.e. Needleman-Wunsch or Gotoh).
// ==========================================================================

#ifndef SEQAN_SEEDS_ALIGN_CHAIN_BANDED_LINEAR_H_
#define SEQAN_SEEDS_ALIGN_CHAIN_BANDED_LINEAR_H_

// TODO(holtgrew): Remove includes here, for debug only.
#include <iostream>
#include <fstream>

namespace seqan {

// ===========================================================================
// Enums, Tags, Classes, Specializations
// ===========================================================================

template <typename AlignmentAlgorithm>
struct AlignmentMatrixDimension_;

template <>
struct AlignmentMatrixDimension_<NeedlemanWunsch>
{
    enum { VALUE = 2 };
};

template <>
struct AlignmentMatrixDimension_<Gotoh>
{
    enum { VALUE = 3 };
};

template <typename TSegment, typename TScoringScheme, typename TAlignmentTag>
class AlignmentChain_
{
public:
    // TODO(holtgrew): Underscores should be in front of variable names, right?
    typedef typename Value<TScoringScheme>::Type TScoreValue;
    typedef typename Position<TSegment>::Type TPosition;

    // The bandwidth for the banded alignment.  Actually, this is the
    // delta to the upper and lower diagonals of the seeds.
    // TODO(holtgrew): Rename to reflect this.
    TPosition bandwidth_;

    // The scoring scheme used for the alignment.
    TScoringScheme scoringScheme_;

    // The first sequence / dimension 0 / query sequence, vertical down
    Holder<TSegment> sequence0_;

    // The second sequence / dimension 1 / database sequence, horizontal ltr
    Holder<TSegment> sequence1_;

    // The alignment matrices for the different part.  The first
    // element is the alignment matrix for the rectangle to the upper
    // left right, followed by the one for the first seed in the chain
    // (left-uppermost one), the next rectangle and so on.
    String<Matrix<TScoreValue, AlignmentMatrixDimension_<TAlignmentTag>::VALUE> > alignmentMatrices_;

    // TODO(holtgrew): Default constructor + setSequence{0,1} missing for now.

    AlignmentChain_(TPosition bandwidth, TScoringScheme const & scoringScheme, TSegment /*const*/ & sequence0, TSegment /*const*/ & sequence1)
            : bandwidth_(bandwidth), scoringScheme_(scoringScheme), sequence0_(sequence0), sequence1_(sequence1)
    { SEQAN_CHECKPOINT; }
};

// ===========================================================================
// Metafunctions
// ===========================================================================

template <typename TSegment, typename TScoringScheme, typename TAlignmentTag>
struct ScoringScheme<AlignmentChain_<TSegment, TScoringScheme, TAlignmentTag> >
{
    typedef TScoringScheme Type;
};


template <typename TSegment, typename TScoringScheme, typename TAlignmentTag>
struct ScoringScheme<AlignmentChain_<TSegment, TScoringScheme, TAlignmentTag> const>
        : ScoringScheme<AlignmentChain_<TSegment, TScoringScheme, TAlignmentTag> > {};

// ===========================================================================
// Functions
// ===========================================================================

// Given a seed and an alignment chain, the overlap of the seed with
// the rectangle to the lower right is written to overlap0 and
// overlap1.
template <typename TSize, typename TSequence, typename TScoringScheme, typename TSeedSpec, typename TSeedConfig, typename TAlignmentTag>
inline void
_computeLowerRightOverlap(TSize & overlap0,
                          TSize & overlap1,
                          Seed<TSeedSpec, TSeedConfig> const & seed,
                          AlignmentChain_<TSequence, TScoringScheme, TAlignmentTag> const & alignmentChain)
{
    SEQAN_CHECKPOINT;
    overlap0 = getUpperDiagonal(seed) - getEndDiagonal(seed) + alignmentChain.bandwidth_ + 1;
    overlap1 = getEndDiagonal(seed) - getLowerDiagonal(seed) + alignmentChain.bandwidth_ + 1;
    // Don't overlap more than the seed...
    overlap0 = _min(overlap0, getEndDim0(seed) - getBeginDim0(seed));
    overlap1 = _min(overlap1, getEndDim1(seed) - getBeginDim1(seed));
}


// Given a seed and an alignment chain, the overlap of the seed with
// the rectangle to the upper left is written to overlap0 and
// overlap1.
template <typename TSize, typename TSequence, typename TScoringScheme, typename TSeedSpec, typename TSeedConfig, typename TAlignmentTag>
inline void
_computeUpperLeftOverlap(TSize & overlap0,
                         TSize & overlap1,
                         Seed<TSeedSpec, TSeedConfig> const & seed,
                         AlignmentChain_<TSequence, TScoringScheme, TAlignmentTag> const & alignmentChain)
{
    SEQAN_CHECKPOINT;
    overlap0 = getStartDiagonal(seed) - getLowerDiagonal(seed) + alignmentChain.bandwidth_ + 1;
    overlap1 = getUpperDiagonal(seed) - getStartDiagonal(seed) + alignmentChain.bandwidth_ + 1;
    // Don't overlap more than the seed...
    overlap0 = _min(overlap0, getEndDim0(seed) - getBeginDim0(seed));
    overlap1 = _min(overlap1, getEndDim1(seed) - getBeginDim1(seed));
}


// Performs the alignment matrix filling step for the leading
// rectangle in the alignment chain.
template <typename TSequence, typename TScoringScheme, typename TAlignmentTag, typename TSeedSpec, typename TSeedConfig, bool START1_FREE, bool START0_FREE, bool END1_FREE, bool END0_FREE>
void
_alignLeadingRectangle(
        AlignmentChain_<TSequence, TScoringScheme, TAlignmentTag> & alignmentChain,
        Seed<TSeedSpec, TSeedConfig> const & rightSeed,
        AlignConfig<START1_FREE, START0_FREE, END1_FREE, END0_FREE> const & alignConfig)
{
    SEQAN_CHECKPOINT;
    // std::cout << "_alignLeadingRectangle ---------------------------------------------" << std::endl;

    typedef typename Value<TScoringScheme>::Type TScoreValue;
    typedef Matrix<TScoreValue, AlignmentMatrixDimension_<TAlignmentTag>::VALUE> TMatrix;
    typedef typename Infix<TSequence>::Type TInfix;
    typedef typename Size<TSequence>::Type TSize;

    // Compute overlap of the rectangle with the seed to the lower right.
    TSize rightOverlap0, rightOverlap1;
    _computeUpperLeftOverlap(rightOverlap0, rightOverlap1, rightSeed, alignmentChain);
    // std::cout << "rightSeed == " << rightSeed << std::endl;
    // std::cout << "rightOverlap0 == " << rightOverlap0 << "rightOverlap1 == " << rightOverlap1 << std::endl;

    // Get infixes of the sequences that correspond to this rectangle.
    TInfix prefix0 = prefix(value(alignmentChain.sequence0_), getBeginDim0(rightSeed) + rightOverlap0);
    TInfix prefix1 = prefix(value(alignmentChain.sequence1_), getBeginDim1(rightSeed) + rightOverlap1);

    // TODO(holtgrew): Temporary debug code.
    std::cerr << ",-- _alignLeadingRectangle" << std::endl;
    std::cerr << "| prefix0: '" << prefix0 << "'" << std::endl;
    std::cerr << "| prefix1: '" << prefix1 << "'" << std::endl;
    std::cerr << "`--" << std::endl;

    // Append a new alignment matrix to the chain.
    appendValue(alignmentChain.alignmentMatrices_, TMatrix());
    // Resize the alignment matrix to the appropriate size.
    _alignResizeMatrix(back(alignmentChain.alignmentMatrices_), prefix0, prefix1, TAlignmentTag());
    // Initialize the matrix gutter.
    _alignInitGutter(back(alignmentChain.alignmentMatrices_), alignmentChain.scoringScheme_, alignConfig, TAlignmentTag());
    // Fill the Matrix using standard dynamic programming.
    _alignFillMatrix(back(alignmentChain.alignmentMatrices_), prefix0, prefix1, alignmentChain.scoringScheme_, TAlignmentTag());

    // TODO(holtgrew): Temporary debug code.
    // {
    //     TMatrix & matrix = back(alignmentChain.alignmentMatrices_);
    //     std::cout << ",-- First rectangle, NW Matrix " << length(matrix, 0) << " x " << length(matrix, 1) << std::endl;
    //     for (unsigned i = 0; i < length(matrix, 0); ++i) {
    //         std::cout << "|";
    //         for (unsigned j = 0; j < length(matrix, 1); ++j) {
    //             if (value(matrix, i, j) == MinValue<int>::VALUE / 2)
    //                 std::cout << "\tinf";
    //             else 
    //                 std::cout << "\t" <<value(matrix, i, j);
    //         }
    //         std::cout << std::endl;
    //     }
    //     std::cout << "`--" << std::endl;
    // }
}


template <typename TSequence, typename TScoringScheme, typename TAlignmentTag, typename TSeedSpec, typename TSeedConfig>
void
_alignTrailingRectangle(
        AlignmentChain_<TSequence, TScoringScheme, TAlignmentTag> & alignmentChain,
        Seed<TSeedSpec, TSeedConfig> const & leftSeed)
{
    SEQAN_CHECKPOINT;
    // std::cout << "_alignTrailingRectangle ---------------------------------------------" << std::endl;

    typedef typename Value<TScoringScheme>::Type TScoreValue;
    typedef Matrix<TScoreValue, AlignmentMatrixDimension_<TAlignmentTag>::VALUE> TMatrix;
    typedef typename Suffix<TSequence>::Type TSuffix;
    typedef typename Size<TSequence>::Type TSize;
    typedef Seed<TSeedSpec, TSeedConfig> TSeed;
    typedef typename Diagonal<TSeed>::Type TDiagonal;

    // Compute overlap of the rectangle with the seed to the upper left.
    TSize leftOverlap0, leftOverlap1;
    _computeLowerRightOverlap(leftOverlap0, leftOverlap1, leftSeed, alignmentChain);

    // Get infixes of the sequences that correspond to this rectangle.
    TSuffix suffix0 = suffix(value(alignmentChain.sequence0_), getEndDim0(leftSeed) - leftOverlap0);
    TSuffix suffix1 = suffix(value(alignmentChain.sequence1_), getEndDim1(leftSeed) - leftOverlap1);

    // Compute lower and upper diagonal for the alignment.
    TDiagonal lowerDiagonal = getStartDiagonal(leftSeed) - getLowerDiagonal(leftSeed) - alignmentChain.bandwidth_;
    TDiagonal upperDiagonal = getUpperDiagonal(leftSeed) - getStartDiagonal(leftSeed) + alignmentChain.bandwidth_;

    // // TODO(holtgrew): Temporary debug code.
    std::cerr << ",-- _alignTrailingRectangle" << std::endl;
    std::cerr << "| suffix0: '" << suffix0 << "'" << std::endl;
    std::cerr << "| suffix1: '" << suffix1 << "'" << std::endl;
    std::cerr << "`--" << std::endl;
    std::cerr << "leftOverlap0 = " << leftOverlap0 << std::endl;
    std::cerr << "leftOverlap1 = " << leftOverlap1 << std::endl;

    // Append a new alignment matrix to the chain.
    appendValue(alignmentChain.alignmentMatrices_, TMatrix());
    // Resize the alignment matrix to the appropriate size.
    _alignResizeMatrix(back(alignmentChain.alignmentMatrices_), suffix0, suffix1, TAlignmentTag());
    // Copy over the data from the banded seed alignment matrix into
    // the gutter and initialize the rest of it according to the
    // alignment config object.
    _alignInitGutterFromBanded(back(alignmentChain.alignmentMatrices_), alignmentChain.scoringScheme_, lowerDiagonal, upperDiagonal, value(end(alignmentChain.alignmentMatrices_) - 2), leftOverlap0, leftOverlap1, TAlignmentTag());
    // // TODO(holtgrew): Temporary debug code.
    // {
    //     TMatrix & matrix = back(alignmentChain.alignmentMatrices_);
    //     std::cout << ",-- Trailing rectangle after init from banded " << length(matrix, 0) << " x " << length(matrix, 1) << std::endl;
    //     for (unsigned i = 0; i < length(matrix, 0); ++i) {
    //         std::cout << "|";
    //         for (unsigned j = 0; j < length(matrix, 1); ++j) {
    //             if (value(matrix, i, j) == MinValue<int>::VALUE / 2)
    //                 std::cout << "\tinf";
    //             else 
    //                 std::cout << "\t" <<value(matrix, i, j);
    //         }
    //         std::cout << std::endl;
    //     }
    //     std::cout << "`--" << std::endl;
    // }
    // Fill the Matrix using standard dynamic programming.
    _alignFillMatrix(back(alignmentChain.alignmentMatrices_), suffix0, suffix1, alignmentChain.scoringScheme_, TAlignmentTag());

    // // TODO(holtgrew): Temporary debug code.
    // {
    //     TMatrix & matrix = back(alignmentChain.alignmentMatrices_);
    //     std::cout << ",-- Last rectangle, NW Matrix " << length(matrix, 0) << " x " << length(matrix, 1) << std::endl;
    //     for (unsigned i = 0; i < length(matrix, 0); ++i) {
    //         std::cout << "|";
    //         for (unsigned j = 0; j < length(matrix, 1); ++j) {
    //             if (value(matrix, i, j) == MinValue<int>::VALUE / 2)
    //                 std::cout << "\tinf";
    //             else 
    //                 std::cout << "\t" <<value(matrix, i, j);
    //         }
    //         std::cout << std::endl;
    //     }
    //     std::cout << "`--" << std::endl;
    // }
}


template <typename TSequence, typename TScoringScheme, typename TAlignmentTag, typename TSeedSpec, typename TSeedConfig>
void
_alignRectangle(
        AlignmentChain_<TSequence, TScoringScheme, TAlignmentTag> & alignmentChain,
        Seed<TSeedSpec, TSeedConfig> const & leftSeed,
        Seed<TSeedSpec, TSeedConfig> const & rightSeed)
{
    SEQAN_CHECKPOINT;
    // std::cout << "_alignRectangle ---------------------------------------------" << std::endl;

    typedef typename Value<TScoringScheme>::Type TScoreValue;
    typedef Matrix<TScoreValue, AlignmentMatrixDimension_<TAlignmentTag>::VALUE> TMatrix;
    typedef typename Infix<TSequence>::Type TInfix;
    typedef typename Size<TSequence>::Type TSize;
    typedef Seed<TSeedSpec, TSeedConfig> TSeed;
    typedef typename Diagonal<TSeed>::Type TDiagonal;

    // Compute overlap of the rectangle with the seed to the upper left and to the lower right.
    TSize leftOverlap0, leftOverlap1, rightOverlap0, rightOverlap1;
    _computeLowerRightOverlap(leftOverlap0, leftOverlap1, leftSeed, alignmentChain);
    _computeUpperLeftOverlap(rightOverlap0, rightOverlap1, rightSeed, alignmentChain);

    // Get infixes of the sequences that correspond to this rectangle.
    TInfix infix0 = infix(value(alignmentChain.sequence0_), getEndDim0(leftSeed) - leftOverlap0, getBeginDim0(rightSeed) + rightOverlap0);
    TInfix infix1 = infix(value(alignmentChain.sequence1_), getEndDim1(leftSeed) - leftOverlap1, getBeginDim1(rightSeed) + rightOverlap1);

    // Compute lower and upper diagonal for the alignment.
    TDiagonal lowerDiagonal = getStartDiagonal(leftSeed) - getLowerDiagonal(leftSeed) - alignmentChain.bandwidth_;
    TDiagonal upperDiagonal = getUpperDiagonal(leftSeed) - getStartDiagonal(leftSeed) + alignmentChain.bandwidth_;

    // // TODO(holtgrew): Temporary debug code.
    std::cout << ",-- _alignRectangle" << std::endl;
    std::cout << "| infix0: '" << infix0 << "'" << std::endl;
    std::cout << "| infix1: '" << infix1 << "'" << std::endl;
    std::cout << "`--" << std::endl;
    std::cout << "left seed = " << leftSeed << std::endl;
    std::cout << "right seed = " << rightSeed << std::endl;

    // Append a new alignment matrix to the chain.
    appendValue(alignmentChain.alignmentMatrices_, TMatrix());
    // Resize the alignment matrix to the appropriate size.
    _alignResizeMatrix(back(alignmentChain.alignmentMatrices_), infix0, infix1, TAlignmentTag());
    // Copy over the data from the banded seed alignment matrix into
    // the gutter and initialize the rest of it according to the
    // alignment config object.
    _alignInitGutterFromBanded(back(alignmentChain.alignmentMatrices_), alignmentChain.scoringScheme_, lowerDiagonal, upperDiagonal, value(end(alignmentChain.alignmentMatrices_) - 2), leftOverlap0, leftOverlap1, TAlignmentTag());
    // // TODO(holtgrew): Temporary debug code.
    // {
    //     TMatrix & matrix = back(alignmentChain.alignmentMatrices_);
    //     std::cout << ",-- Rectangle after init from banded " << length(matrix, 0) << " x " << length(matrix, 1) << std::endl;
    //     for (unsigned i = 0; i < length(matrix, 0); ++i) {
    //         std::cout << "|";
    //         for (unsigned j = 0; j < length(matrix, 1); ++j) {
    //             if (value(matrix, i, j) == MinValue<int>::VALUE / 2)
    //                 std::cout << "\tinf";
    //             else 
    //                 std::cout << "\t" <<value(matrix, i, j);
    //         }
    //         std::cout << std::endl;
    //     }
    //     std::cout << "`--" << std::endl;
    // }
    // Fill the Matrix using standard dynamic programming.
    _alignFillMatrix(back(alignmentChain.alignmentMatrices_), infix0, infix1, alignmentChain.scoringScheme_, TAlignmentTag());

    // // TODO(holtgrew): Temporary debug code.
    // {
    //     TMatrix & matrix = back(alignmentChain.alignmentMatrices_);
    //     std::cout << ",-- Middle rectangle, NW Matrix " << length(matrix, 0) << " x " << length(matrix, 1) << std::endl;
    //     for (unsigned i = 0; i < length(matrix, 0); ++i) {
    //         std::cout << "|";
    //         for (unsigned j = 0; j < length(matrix, 1); ++j) {
    //             if (value(matrix, i, j) == MinValue<int>::VALUE / 2)
    //                 std::cout << "\tinf";
    //             else 
    //                 std::cout << "\t" <<value(matrix, i, j);
    //         }
    //         std::cout << std::endl;
    //     }
    //     std::cout << "`--" << std::endl;
    // }
}


template <typename TSequence, typename TScoringScheme, typename TAlignmentTag, typename TSeedSpec, typename TSeedConfig>
void
_alignSeed(
        AlignmentChain_<TSequence, TScoringScheme, TAlignmentTag> & alignmentChain,
        Seed<TSeedSpec, TSeedConfig> const & seed)
{
    SEQAN_CHECKPOINT;
    // std::cout << "_alignSeed ---------------------------------------------" << std::endl;

    typedef typename Value<TScoringScheme>::Type TScoreValue;
    typedef Matrix<TScoreValue, AlignmentMatrixDimension_<TAlignmentTag>::VALUE> TMatrix;
    typedef typename Infix<TSequence>::Type TInfix;
    typedef typename Size<TSequence>::Type TSize;
    typedef Seed<TSeedSpec, TSeedConfig> TSeed;
    typedef typename Diagonal<TSeed>::Type TDiagonal;

    // Compute overlap with the rectangle to the upper left.
    TSize leftOverlap0, leftOverlap1;
    _computeUpperLeftOverlap(leftOverlap0, leftOverlap1, seed, alignmentChain);

    // Get the infixes corresponding to the seed.
    TInfix infix0 = infix(value(alignmentChain.sequence0_), getBeginDim0(seed), getEndDim0(seed));
    TInfix infix1 = infix(value(alignmentChain.sequence1_), getBeginDim1(seed), getEndDim1(seed));

    // Compute lower and upper diagonal for the alignment.
    TDiagonal lowerDiagonal = getLowerDiagonal(seed) - getStartDiagonal(seed) - alignmentChain.bandwidth_;
    TDiagonal upperDiagonal = getUpperDiagonal(seed) - getStartDiagonal(seed) + alignmentChain.bandwidth_;
    // Limit lower and upper diagonal so no more than a whole rectangle is computed.
    TDiagonal maxBandwidth = getEndDim0(seed) - getBeginDim0(seed) - 1;
    // std::cout << "lowerDiagonal = _max(" << lowerDiagonal << ", - " << maxBandwidth << ");" << std::endl;
    // std::cout << "upperDiagonal = _min(" << upperDiagonal << ", " << getEndDiagonal(seed) << " - " << getStartDiagonal(seed) << " + " << maxBandwidth << ");" << std::endl;
    lowerDiagonal = _max(lowerDiagonal, - maxBandwidth);
    upperDiagonal = _min(upperDiagonal, getEndDiagonal(seed) - getStartDiagonal(seed) + maxBandwidth);
    // std::cout << "seed == " << seed << std::endl;
    // std::cout << "leftOverlap0 == " << leftOverlap0 << "leftOverlap1 == " << leftOverlap1 << std::endl;

    // // TODO(holtgrew): Temporary debug code.
    std::cerr << ",-- _alignSeed" << std::endl;
    std::cerr << "| infix0: '" << infix0 << "'" << std::endl;
    std::cerr << "| infix1: '" << infix1 << "'" << std::endl;
    std::cerr << "`--" << std::endl;

    // Append a new alignment matrix to the chain.
    appendValue(alignmentChain.alignmentMatrices_, TMatrix());

    // Resize the alignment matrix to the appropriate size.
    _alignBandedResizeMatrix(back(alignmentChain.alignmentMatrices_), infix0, infix1, lowerDiagonal, upperDiagonal, TAlignmentTag());
    // Initialize the banded DP matrix' gutter.  Some data is copied
    // in from the previous non-banded DP matrix.
    // _alignBandedInitGutter(back(alignmentChain.alignmentMatrices_), alignmentChain.scoringScheme_, lowerDiagonal, upperDiagonal, AlignConfig<false, false, false, false>(), TAlignmentTag());
    _alignBandedInitGutterFromUnbanded(back(alignmentChain.alignmentMatrices_), alignmentChain.scoringScheme_, lowerDiagonal, upperDiagonal, value(end(alignmentChain.alignmentMatrices_) - 2), leftOverlap0, leftOverlap1, TAlignmentTag());
    // // TODO(holtgrew): Debug output, remove when not needed any more.
    // {
    //     TMatrix & matrix = back(alignmentChain.alignmentMatrices_);
    //     std::cout << ",-- matrix after init gutter from unbanded" << std::endl;
    //     for (unsigned i = 0; i < length(matrix, 0); ++i) {
    //         std::cout << "| ";
    //         for (unsigned j = 0; j < i; ++j)
    //             std::cout << "\t";
    //         for (unsigned j = 0; j < length(matrix, 1); ++j) {
    //             if (value(matrix, i, j) == MinValue<int>::VALUE / 2)
    //                 std::cout << "\tinf";
    //             else
    //                 std::cout << "\t" << value(matrix, i, j);
    //         }
    //         std::cout << std::endl;
    //     }
    //     std::cout << "`--" << std::endl;
    // }
    // Fill the Matrix using banded dynamic programming.
    _alignBandedFillMatrix(back(alignmentChain.alignmentMatrices_), infix0, infix1, alignmentChain.scoringScheme_, lowerDiagonal, upperDiagonal, TAlignmentTag());

    // // TODO(holtgrew): Debug output, remove when not needed any more.
    // {
    //     TMatrix & matrix = back(alignmentChain.alignmentMatrices_);
    //     std::cout << ",-- matrix after DP filling" << std::endl;
    //     for (unsigned i = 0; i < length(matrix, 0); ++i) {
    //         std::cout << "| ";
    //         for (unsigned j = 0; j < i; ++j)
    //             std::cout << "\t";
    //         for (unsigned j = 0; j < length(matrix, 1); ++j) {
    //             if (value(matrix, i, j) == MinValue<int>::VALUE / 2)
    //                 std::cout << "\tinf";
    //             else
    //                 std::cout << "\t" << value(matrix, i, j);
    //         }
    //         std::cout << std::endl;
    //     }
    //     std::cout << "`--" << std::endl;
    // }
}


template <typename TAlignment, typename TSequence, typename TScoringScheme, typename TAlignmentTag, typename TSeedChain, bool START1_FREE, bool START0_FREE, bool END1_FREE, bool END0_FREE>
typename Value<typename ScoringScheme<Alignment<TSequence, TScoringScheme, TAlignmentTag> >::Type >::Type
_glueAlignmentChain(
        TAlignment & alignment,
        AlignmentChain_<TSequence, TScoringScheme, TAlignmentTag> const & alignmentChain,
        TSeedChain const & seedChain,
        AlignConfig<START1_FREE, START0_FREE, END1_FREE, END0_FREE> const & alignConfig)
{
    SEQAN_CHECKPOINT;

    typedef typename Value<TScoringScheme>::Type TScoreValue;
    typedef String<Matrix<TScoreValue, AlignmentMatrixDimension_<TAlignmentTag>::VALUE> > const TMatrixString;
    typedef typename Iterator<TMatrixString, Standard>::Type TMatrixStringIterator;
    typedef typename Position<TMatrixString>::Type TPosition;

    typedef typename Value<TSeedChain const>::Type TSeed;
    typedef typename Iterator<TSeedChain const, Standard>::Type TSeedChainIterator;

    // Define type for iterator over alignment rows.
	typedef typename Row<TAlignment>::Type TAlignmentRow;
	typedef typename Iterator<TAlignmentRow, Standard>::Type TTargetIterator;
    typedef typename Infix<TSequence>::Type TSourceInfix;
    typedef typename Iterator<TSourceInfix, Standard>::Type TSourceInfixIterator;

	TSourceInfixIterator sequenceIt0 = end(sourceSegment(row(alignment, 0))) - 1;
	TSourceInfixIterator sequenceIt1 = end(sourceSegment(row(alignment, 1))) - 1;
    TTargetIterator alignmentIt0 = end(row(alignment, 0));
    TTargetIterator alignmentIt1 = end(row(alignment, 1));
    TMatrixStringIterator matricesIt = end(alignmentChain.alignmentMatrices_) - 1;
    TSeedChainIterator seedChainIt = end(seedChain, Standard());
    --seedChainIt;

    // Traceback through trailing rectangle and seed.
    TPosition finalPos0 = 1;
    TPosition finalPos1 = 1;
    std::cout << "trace back through trailing rectangle..." << std::endl;
    TPosition upperLeftOverlap0, upperLeftOverlap1;
    _computeLowerRightOverlap(upperLeftOverlap0, upperLeftOverlap1, value(seedChainIt), alignmentChain);
    std::cerr << "upperLeftOverlap0 = " << upperLeftOverlap0 << " upperLeftOverlap1 = " << upperLeftOverlap1 << std::endl;
    TPosition lowerRightOverlap0 = 1;
    TPosition lowerRightOverlap1 = 1;
    TScoreValue result = _alignTraceback(alignmentIt0, alignmentIt1, sequenceIt0, sequenceIt1, finalPos0, finalPos1, value(matricesIt), alignmentChain.scoringScheme_, lowerRightOverlap0, lowerRightOverlap1, upperLeftOverlap0, upperLeftOverlap1, false, alignConfig, TAlignmentTag());
    // std::cout << "Alignment so far:" << std::endl;
    // std::cout << alignment;
    // std::cout << result << std::endl;
    // std::cout << "finalPos0 = " << finalPos0 << ", finalPos1 = " << finalPos1 << std::endl;
    goPrevious(matricesIt);
    _computeLowerRightOverlap(lowerRightOverlap0, lowerRightOverlap1, value(seedChainIt), alignmentChain);
    TPosition upperTriangleEdgeLength = getStartDiagonal(value(seedChainIt)) - getLowerDiagonal(value(seedChainIt)) + alignmentChain.bandwidth_;
    TPosition lowerTriangleEdgeLength = getUpperDiagonal(value(seedChainIt)) - getEndDiagonal(value(seedChainIt)) + alignmentChain.bandwidth_;
    TPosition maxTriangleEdgeLength = getEndDim0(value(seedChainIt)) - getBeginDim0(value(seedChainIt)) - 1;
    // std::cout << "lowerTriangleEdgeLength == " << lowerTriangleEdgeLength << std::endl;
    // std::cout << "upperTriangleEdgeLength == " << upperTriangleEdgeLength << std::endl;
    // std::cout << "maxTriangleEdgeLength == " << maxTriangleEdgeLength << std::endl;
    lowerTriangleEdgeLength = _min(lowerTriangleEdgeLength, maxTriangleEdgeLength);
    upperTriangleEdgeLength = _min(upperTriangleEdgeLength, maxTriangleEdgeLength);
    std::cout << "trace back through last seed..." << std::endl;
    std::cout << "finalPos0 == " << finalPos0 << ", finalPos1 == " << finalPos1 << std::endl;
    {
        if (length(seedChain) > 1u) {
            TSeedChainIterator seedChainItPrevious = seedChainIt;
            goPrevious(seedChainItPrevious);
            _computeLowerRightOverlap(upperLeftOverlap0, upperLeftOverlap1, value(seedChainItPrevious), alignmentChain);
        } else {
            upperLeftOverlap0 = 0;
            upperLeftOverlap1 = 0;
        }
    }
    _alignBandedTraceback(alignmentIt0, alignmentIt1, sequenceIt0, sequenceIt1, finalPos0, finalPos1, value(matricesIt), alignmentChain.scoringScheme_, lowerRightOverlap0, lowerRightOverlap1, upperTriangleEdgeLength, lowerTriangleEdgeLength, false, AlignConfig<false, false, false, false>(), TAlignmentTag());
    goPrevious(matricesIt);

    // std::cout << "Alignment so far:" << std::endl << alignment;
    
    // Traceback through matrices in reverse order.
    for (TPosition i = 0, iend = length(seedChain) - 1; i < iend; ++i) {
        std::cout << "trace back through rectangle..." << std::endl;
        // std::cout << "finalPos0 == " << finalPos0 << ", finalPos1 == " << finalPos1 << std::endl;
        _computeUpperLeftOverlap(lowerRightOverlap0, lowerRightOverlap1, value(seedChainIt), alignmentChain);
        {
            TSeedChainIterator seedChainItPrevious = seedChainIt;
            goPrevious(seedChainItPrevious);
            _computeLowerRightOverlap(upperLeftOverlap0, upperLeftOverlap1, value(seedChainItPrevious), alignmentChain);
        }
        _alignTraceback(alignmentIt0, alignmentIt1, sequenceIt0, sequenceIt1, finalPos0, finalPos1, value(matricesIt), alignmentChain.scoringScheme_, lowerRightOverlap0, lowerRightOverlap1, upperLeftOverlap0, upperLeftOverlap1, false, AlignConfig<false, false, false, false>(), TAlignmentTag());
        goPrevious(matricesIt);
        goPrevious(seedChainIt);
        // std::cout << "Alignment so far:" << std::endl << alignment;

        std::cout << "trace back through seed..." << std::endl;
        // std::cout << "finalPos0 == " << finalPos0 << ", finalPos1 == " << finalPos1 << std::endl;
        _computeLowerRightOverlap(lowerRightOverlap0, lowerRightOverlap1, value(seedChainIt), alignmentChain);
        upperTriangleEdgeLength = getStartDiagonal(value(seedChainIt)) - getLowerDiagonal(value(seedChainIt)) + alignmentChain.bandwidth_;
        lowerTriangleEdgeLength = getUpperDiagonal(value(seedChainIt)) - getEndDiagonal(value(seedChainIt)) + alignmentChain.bandwidth_;
        maxTriangleEdgeLength = getEndDim0(value(seedChainIt)) - getBeginDim0(value(seedChainIt)) - 1;
        // std::cout << "lowerTriangleEdgeLength == " << lowerTriangleEdgeLength << std::endl;
        // std::cout << "upperTriangleEdgeLength == " << upperTriangleEdgeLength << std::endl;
        // std::cout << "maxTriangleEdgeLength == " << maxTriangleEdgeLength << std::endl;
        lowerTriangleEdgeLength = _min(lowerTriangleEdgeLength, maxTriangleEdgeLength);
        upperTriangleEdgeLength = _min(upperTriangleEdgeLength, maxTriangleEdgeLength);
        _alignBandedTraceback(alignmentIt0, alignmentIt1, sequenceIt0, sequenceIt1, finalPos0, finalPos1, value(matricesIt), alignmentChain.scoringScheme_, lowerRightOverlap0, lowerRightOverlap1, upperTriangleEdgeLength, lowerTriangleEdgeLength, false, AlignConfig<false, false, false, false>(), TAlignmentTag());
        goPrevious(matricesIt);
        // std::cout << "Alignment so far:" << std::endl << alignment;
    }

    std::cout << "trace back through leading..." << std::endl;
    // std::cout << "finalPos0 == " << finalPos0 << ", finalPos1 == " << finalPos1 << std::endl;
    // Traceback through leading rectangle.
    upperLeftOverlap0 = 0;
    upperLeftOverlap1 = 0;
    _computeUpperLeftOverlap(lowerRightOverlap0, lowerRightOverlap1, value(seedChainIt), alignmentChain);
    _alignTraceback(alignmentIt0, alignmentIt1, sequenceIt0, sequenceIt1, finalPos0, finalPos1, value(matricesIt), alignmentChain.scoringScheme_, lowerRightOverlap0, lowerRightOverlap1, upperLeftOverlap0, upperLeftOverlap1, true, alignConfig, TAlignmentTag());
    // std::cout << "Alignment so far:" << std::endl << alignment;
    return result;
}


template <typename TStream, typename TAlignment, typename TAlignmentChain, typename TSeedChain>
void
write(TStream & stream, TAlignment const & alignment, TAlignmentChain const & alignmentChain, TSeedChain const & seedChain, Tikz_ const &)
{
//IOREV _nodoc_ specialization not documented
    typedef typename Row<TAlignment const>::Type TRow;
    typedef typename Position<TAlignment>::Type TPosition;
    typedef typename Iterator<TRow, Standard>::Type TRowIterator;
    typedef typename Iterator<TSeedChain, Standard>::Type TSeedIterator;
  
    stream << "\\documentclass{article}" << std::endl
           << "\\usepackage{tikz,nicefrac,amsmath,pifont}" << std::endl
           << "\\usetikzlibrary{arrows,snakes,backgrounds,patterns,matrix,shapes,fit,calc,shadows,plotmarks}" << std::endl
           << "\\begin{document}" << std::endl
           // << "\\begin{tikzpicture}[scale=.5,core/.style={fill=blue!20, fill opacity=.8},band/.style={fill=blue!10, fill opacity=.8}]" << std::endl;
           << "\\begin{tikzpicture}[scale=.5,core/.style={},band/.style={},rectangle/.style={fill=green!10, fill opacity=.1}]" << std::endl;

    // Grid.
    stream << "\\draw[help lines] (0, -" << length(row(alignment, 0)) << ") grid (" << length(row(alignment, 1)) << ", 0);" << std::endl;
    // Top row of characters.
    stream << "\\draw";
    int i = 0;
    for (TRowIterator it = begin(row(alignment, 0)); it != end(row(alignment, 0)); ++it, ++i) {
        stream << " (-1," << -0.5-i << ") node {" << *it << "}";
    }
    stream << ";" << std::endl;
    // Left row of characters.
    stream << "\\draw";
    i = 0;
    for (TRowIterator it = begin(row(alignment, 1)); it != end(row(alignment, 1)); ++it, ++i) {
        stream << " (" << 0.5+i << ", 1) node {" << *it << "}";
    }
    stream << ";" << std::endl;
    // Draw seeds.
    for (TSeedIterator it = begin(seedChain); it != end(seedChain); ++it) {
        int height = getEndDim0(*it) - getBeginDim0(*it);
        int width = getEndDim1(*it) - getBeginDim1(*it);
        int ext = _min(height, width);
        int s0 = getBeginDim0(*it);
        int s1 = getBeginDim1(*it);
        int e0 = getEndDim0(*it);
        int e1 = getEndDim1(*it);
        int ld = getLowerDiagonal(*it);
        int ud = getUpperDiagonal(*it);
        int sd = getStartDiagonal(*it);
        int ed = getEndDiagonal(*it);
        int b = alignmentChain.bandwidth_;
        // Draw seed "core."
        stream << "\\draw[core] (" << s1 << ", -" << s0 << ") -- (" << s1 + ext << ", -" << s0 + ext << ") -- (" << e1 << ", -" << e0 << ") -- (" << e1 - ext << ", -" << e0 - ext << ") -- cycle;" << std::endl;
        // Draw left band.
        stream << "\\draw[band] (" << s1 - ((ld - sd) + b) << ", -" << s0 << ") -- (" << s1 - ((ld - sd) + b) + ext << ", -" << s0 + ext << ") -- (" << s1 + ext << ", -" << s0 + ext << ") -- (" << s1 << ", -" << s0 << ") -- cycle;" << std::endl;
        // Draw right band.
        stream << "\\draw[band] (" << e1 << ", -" << e0 << ") -- (" << e1 + ((ud - ed) + b) << ", -" << e0 << ") -- (" << e1 + ((ud - ed) + b) - ext << ", -" << e0 - ext << ") -- (" << e1 - ext << ", -" << e0 - ext << ") -- cycle;" << std::endl;
    }

    // Draw rectangles.
    // First rectangle.
    {
        int s0 = getBeginDim0(front(seedChain));
        int s1 = getBeginDim1(front(seedChain));
        int o0, o1;
        _computeLowerRightOverlap(o0, o1, front(seedChain), alignmentChain);
        stream << "\\draw[rectangle] (0, 0) rectangle (" << s1 + o1 << ", -" << s0 + o0 << ");" << std::endl;
    }
    // Middle rectangles.
    if (length(seedChain) > 1u) {
        for (TSeedIterator it = begin(seedChain) + 1; it != end(seedChain); ++it) {
            int e0 = getEndDim0(value(it - 1));
            int e1 = getEndDim1(value(it - 1));
            int s0 = getBeginDim0(value(it));
            int s1 = getBeginDim1(value(it));
            int ou0, ou1, ol0, ol1;
            _computeLowerRightOverlap(ou0, ou1, value(it-1), alignmentChain);
            _computeUpperLeftOverlap(ol0, ol1, value(it), alignmentChain);
            printf("... e0=%d e1=%d s0=%d s1=%d ou0=%d ou1=%d ol0=%d ol1=%d\n", e0, e1, s0, s1, ou0, ou1, ol0, ol1);
            stream << "\\draw[rectangle] (" << e1 - ou1 << ", -" << e0 - ou0 << ") rectangle (" << s1 + ol1 << ", -" << s0 + ol0 << ");" << std::endl;
        }
    }
    // Last rectangle.
    {
        int e0 = getEndDim0(back(seedChain));
        int e1 = getEndDim1(back(seedChain));
        int o0, o1;
        _computeLowerRightOverlap(o0, o1, back(seedChain), alignmentChain);
        stream << "\\draw[rectangle] (" << e1 - o1 << ", -" << e0 - o0 << ") rectangle (" << length(row(alignment, 1)) << ", -" << length(row(alignment, 0)) << ");" << std::endl;
    }

    stream << "\\end{tikzpicture}" << std::endl
           << "\\end{document}" << std::endl;
}

template <typename TContainer, typename TBandwidth, typename TScoreValue, typename TAlign, bool START1_FREE, bool START0_FREE, bool END1_FREE, bool END0_FREE, typename TGlobalAlignmentTag>
TScoreValue
_bandedChainAlignment(
        TAlign & alignment,
        TContainer const & seedChain,
        TBandwidth k,
        Score<TScoreValue, Simple> const & scoringScheme,
        AlignConfig<START1_FREE, START0_FREE, END1_FREE, END0_FREE> const & alignConfig,
        TGlobalAlignmentTag const &)
{
    SEQAN_CHECKPOINT;

    SEQAN_ASSERT_GT(length(seedChain), 0u);
    #if SEQAN_ENABLE_DEBUG
    {
        // Assert that the chain is non-overlapping.
        typedef typename Iterator<TContainer const, Standard>::Type TIterator;
        // std::cout << ".-- Chain (" << __FILE__ << ":" << __LINE__ << "):" << std::endl;
        // for (TIterator it = begin(seedChain, Standard()); it != end(seedChain, Standard()); ++it)
        //     std::cout << "| " << *it << std::endl;
        // std::cout << "`--" << std::endl;
        TIterator itPrevious = begin(seedChain, Standard());
        TIterator it = itPrevious;
        // std::cout << *it << std::endl;
        ++it;
        for (; it != end(seedChain, Standard()); ++it) {
            // std::cout << *it << std::endl;
            SEQAN_ASSERT_LEQ(getEndDim0(*itPrevious), getBeginDim0(*it));
            SEQAN_ASSERT_LEQ(getEndDim1(*itPrevious), getBeginDim1(*it));
            itPrevious = it;
        }
    }
    #endif  // #if SEQAN_ENABLE_DEBUG

    //
    // Function-Wide Typedefs
    //
    typedef typename Source<TAlign>::Type TSequence;
    typedef typename Infix<TSequence>::Type TSegment;
    typedef typename Value<TContainer const>::Type TSeed;
    typedef typename Iterator<TContainer const, Standard>::Type TIterator;
    typedef Score<TScoreValue, Simple> TScoringScheme;
    typedef AlignmentChain_<TSegment, TScoringScheme, TGlobalAlignmentTag> TAlignmentChain;

    //
    // Initialization
    //
	TSegment seq1 = sourceSegment(row(alignment, 0));
	TSegment seq2 = sourceSegment(row(alignment, 1));
    TAlignmentChain alignmentChain(k, scoringScheme, seq1, seq2);

    // std::cerr << "length(seedChain) == " << length(seedChain) << std::endl;

    // // TODO(holtgrew): Temporary debug code.
    // std::cout << ",-- _bandedChainAlignment" << std::endl;
    // std::cout << "| seq1:" << seq1 << std::endl;
    // std::cout << "| seq2:" << seq2 << std::endl;
    // std::cout << "`--" << std::endl;

    //
    // Compute Alignment using Alignment Chain
    //
    // Compute alignment for leading rectangle and the first seed.
    // std::cout << "first seed: " << front(seedChain) << std::endl;
    _alignLeadingRectangle(alignmentChain, front(seedChain), alignConfig);
    _alignSeed(alignmentChain, front(seedChain));

    // For all seeds from the second one from the left to the
    // rightmost one: Align rectangle left of it and then the seed.
    TIterator itPrevious = begin(seedChain, Standard());
    TIterator it = itPrevious;
    ++it;
    for (; it != end(seedChain, Standard()); ++it) {
        _alignRectangle(alignmentChain, value(itPrevious), value(it));
        _alignSeed(alignmentChain, value(it));
        itPrevious = it;
    }

    // Compute alignment for the trailing rectangle.
    _alignTrailingRectangle(alignmentChain, back(seedChain));

    // // Write alignment to file.
    // {
    //     std::ofstream file;
    //     file.open("/tmp/example.tex", std::ios::out);
    //     write(file, alignment, alignmentChain, seedChain, Tikz_());
    // }

    // Glue all alignments together.
    return _glueAlignmentChain(alignment, alignmentChain, seedChain, alignConfig);
}


/**
.Internal.Function.bandedChainAlignment:
..summary:Calculates a banded alignment around a chain of seeds. 
..cat:Seed Handling
..signature:bandedChainAlignment(seedChain, k, alignment, scoreMatrix);
..param.seedChain:A chain of seeds, must be ascendingly sorted in both dimensions.
..param.k:Half of the width of the band.
..param.alignment:The alignment where the result is stored.
...type:Class.Align
..param.scoreMatrix: The score matrix.
...type:Spec.Simple Score
...remarks: Depending on the score matrix the Needleman-Wunsch or the Gotoh algorithm is used. For a description of the algorithm see the masters thesis of C. Kemena, Section 5.3.3 LAGAN Alignment.
..returns: The score of the alignment.
..include:seqan/seeds2.h
*/
// TODO(holtgrew): wholeAlignment is the result and should be the first parameter.
// TODO(holtgrew): Adjust the documentation to the parameter names.
template<typename TContainer, typename TValue, typename TScoreValue, typename TAlign, bool START1_FREE, bool START0_FREE, bool END1_FREE, bool END0_FREE>
TScoreValue
bandedChainAlignment(TAlign & alignment,
                     TContainer const & seedChain, 
					 TValue k,
					 Score<TScoreValue, Simple> const & scoringScheme,
                     AlignConfig<START1_FREE, START0_FREE, END1_FREE, END0_FREE> const & alignConfig)
{
    SEQAN_CHECKPOINT;

	if (scoreGapOpen(scoringScheme) == scoreGapExtend(scoringScheme))
		return _bandedChainAlignment(alignment, seedChain, k, scoringScheme, alignConfig, NeedlemanWunsch());
	else
		return _bandedChainAlignment(alignment, seedChain, k, scoringScheme, alignConfig, Gotoh());
}

}  // namespace seqan

#endif  // #ifndef SEQAN_SEEDS_ALIGN_CHAIN_BANDED_LINEAR_H_
