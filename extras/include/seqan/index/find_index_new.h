// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2013 NVIDIA Corporation
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
//     * Neither the name of NVIDIA Corporation nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL NVIDIA CORPORATION BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_EXTRAS_INDEX_FIND_INDEX_H_
#define SEQAN_EXTRAS_INDEX_FIND_INDEX_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// TODO(esiragusa): Remove this when there will be a base finder class.
template <typename TText, typename TPattern, typename TSpec = void>
struct Finder2;

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction View
// ----------------------------------------------------------------------------
// TODO(esiragusa): move this function in the base finder class.

template <typename TText, typename TPattern, typename TSpec>
struct View<Finder2<TText, TPattern, TSpec> >
{
    typedef Finder2<typename View<TText>::Type, typename View<TPattern>::Type, TSpec>   Type;
};

// ----------------------------------------------------------------------------
// Metafunction RemoveView
// ----------------------------------------------------------------------------
// TODO(esiragusa): move this function in the base finder class.

template <typename TText, typename TPattern, typename TSpec>
struct RemoveView<Finder2<TText, TPattern, TSpec> >
{
    typedef Finder2<typename RemoveView<TText>::Type, typename RemoveView<TPattern>::Type, TSpec>   Type;
};

// ----------------------------------------------------------------------------
// Metafunction IsView
// ----------------------------------------------------------------------------
// TODO(esiragusa): move this function in the base finder class.

template <typename TText, typename TPattern, typename TSpec>
struct IsView<Finder2<TText, TPattern, TSpec> > : IsView<TText> {};

// ----------------------------------------------------------------------------
// Metafunction IsDevice
// ----------------------------------------------------------------------------
// TODO(esiragusa): move this function in the base finder class.

template <typename TText, typename TPattern, typename TSpec>
struct IsDevice<Finder2<TText, TPattern, TSpec> > : IsDevice<TText> {};

// ----------------------------------------------------------------------------
// Metafunction View
// ----------------------------------------------------------------------------
// TODO(esiragusa): move this function in the base finder class.

template <typename TNeedle, typename TSpec>
struct View<Pattern<TNeedle, TSpec> >
{
    typedef Pattern<typename View<TNeedle>::Type, TSpec>    Type;
};

// ----------------------------------------------------------------------------
// Metafunction RemoveView
// ----------------------------------------------------------------------------
// TODO(esiragusa): move this function in the base finder class.

template <typename TNeedle, typename TSpec>
struct RemoveView<Pattern<TNeedle, TSpec> >
{
    typedef Pattern<typename RemoveView<TNeedle>::Type, TSpec>    Type;
};

// ----------------------------------------------------------------------------
// Metafunction IsView
// ----------------------------------------------------------------------------
// TODO(esiragusa): move this function in the base finder class.

template <typename TNeedle, typename TSpec>
struct IsView<Pattern<TNeedle, TSpec> > : IsView<TNeedle> {};

// ----------------------------------------------------------------------------
// Metafunction IsDevice
// ----------------------------------------------------------------------------
// TODO(esiragusa): move this function in the base finder class.

template <typename TNeedle, typename TSpec>
struct IsDevice<Pattern<TNeedle, TSpec> > : IsDevice<TNeedle> {};

// ----------------------------------------------------------------------------
// Metafunction TextIterator_
// ----------------------------------------------------------------------------

// TODO(esiragusa): move this function in the base finder class.
template <typename TText, typename TSpec>
struct TextIterator_
{
    typedef typename Iterator<TText>::Type  Type;
};

// ----------------------------------------------------------------------------
// Metafunction PatternIterator_
// ----------------------------------------------------------------------------

// TODO(esiragusa): move this function in the base finder class.
template <typename TPattern, typename TSpec>
struct PatternIterator_
{
    typedef typename Iterator<TPattern const, Rooted>::Type  Type;
};

// ----------------------------------------------------------------------------
// Metafunction TextIterator_
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TSpec>
struct TextIterator_<Index<TText, TIndexSpec>, TSpec>
{
    typedef typename Iterator<Index<TText, TIndexSpec>, TopDown<> >::Type  Type;
};

template <typename TText, typename TIndexSpec, typename TDistance, typename TSpec>
struct TextIterator_<Index<TText, TIndexSpec>, Backtracking<TDistance, TSpec> >
{
    typedef typename Iterator<Index<TText, TIndexSpec>, TopDown<ParentLinks<> > >::Type  Type;
};

// ----------------------------------------------------------------------------
// Metafunction Score_
// ----------------------------------------------------------------------------

template <typename TSpec>
struct Score_;

template <>
struct Score_<FinderSTree>
{
    typedef Nothing Type;
};

template <typename TSpec>
struct Score_<Backtracking<HammingDistance, TSpec> >
{
    typedef unsigned char   Type;
};

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class Finder
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TPattern, typename TSpec>
struct Finder2<Index<TText, TIndexSpec>, TPattern, TSpec>
{
    typedef Index<TText, TIndexSpec>                                TIndex;
    typedef typename TextIterator_<TIndex, TSpec>::Type             TTextIterator;
    typedef typename PatternIterator_<TPattern, TSpec>::Type        TPatternIterator;
    typedef typename Score_<TSpec>::Type                            TScore;

    TTextIterator       _textIt;
    TPatternIterator    _patternIt;
    TScore              _scoreThreshold;
    TScore              _score;

    SEQAN_HOST_DEVICE
    Finder2() :
        _scoreThreshold(),
        _score()
    {}

    SEQAN_HOST_DEVICE
    Finder2(TIndex /* const */ & index) :
        _textIt(index),
        _scoreThreshold(),
        _score()
    {}

    SEQAN_HOST_DEVICE
    Finder2(TTextIterator const & textIt) :
        _textIt(textIt),
        _scoreThreshold(),
        _score()
    {}
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function textIterator()
// ----------------------------------------------------------------------------
// TODO(esiragusa): move this function in the base finder class.

template <typename TText, typename TPattern, typename TSpec>
SEQAN_HOST_DEVICE inline typename TextIterator_<TText, TSpec>::Type &
textIterator(Finder2<TText, TPattern, TSpec> & finder)
{
    return finder._textIt;
}

template <typename TText, typename TPattern, typename TSpec>
SEQAN_HOST_DEVICE inline typename TextIterator_<TText, TSpec>::Type const &
textIterator(Finder2<TText, TPattern, TSpec> const & finder)
{
    return finder._textIt;
}

// ----------------------------------------------------------------------------
// Function patternIterator()
// ----------------------------------------------------------------------------
// TODO(esiragusa): move this function in the base finder class.

template <typename TText, typename TPattern, typename TSpec>
SEQAN_HOST_DEVICE inline typename PatternIterator_<TPattern, TSpec>::Type &
patternIterator(Finder2<TText, TPattern, TSpec> & finder)
{
    return finder._patternIt;
}

template <typename TText, typename TPattern, typename TSpec>
SEQAN_HOST_DEVICE inline typename PatternIterator_<TPattern, TSpec>::Type const &
patternIterator(Finder2<TText, TPattern, TSpec> const & finder)
{
    return finder._patternIt;
}

// ----------------------------------------------------------------------------
// Function setPatternIterator()
// ----------------------------------------------------------------------------
// TODO(esiragusa): move this function in the base finder class.

template <typename TText, typename TPattern, typename TSpec, typename TPatternIterator>
SEQAN_HOST_DEVICE inline void
setPatternIterator(Finder2<TText, TPattern, TSpec> & finder, TPatternIterator const & patternIt)
{
    finder._patternIt = patternIt;
}

// ----------------------------------------------------------------------------
// Function getScore()
// ----------------------------------------------------------------------------
// TODO(esiragusa): move this function in the base finder class.

template <typename TText, typename TPattern, typename TSpec>
SEQAN_HOST_DEVICE inline typename Score_<TSpec>::Type
getScore(Finder2<TText, TPattern, TSpec> const & finder)
{
    return finder._score;
}

// ----------------------------------------------------------------------------
// Function setScoreThreshold()
// ----------------------------------------------------------------------------
// TODO(esiragusa): move this function in the base finder class.

template <typename TText, typename TPattern, typename TSpec, typename TScore>
SEQAN_HOST_DEVICE inline void
setScoreThreshold(Finder2<TText, TPattern, TSpec> & finder, TScore score)
{
    finder._scoreThreshold = score;
}

// ----------------------------------------------------------------------------
// Function _getVertexScore()
// ----------------------------------------------------------------------------

template <typename TText, typename TPattern, typename TSpec>
SEQAN_HOST_DEVICE inline typename Score_<Backtracking<HammingDistance, TSpec> >::Type
_getVertexScore(Finder2<TText, TPattern, Backtracking<HammingDistance, TSpec> > const & finder)
{
    return parentEdgeLabel(textIterator(finder)) != value(patternIterator(finder));
}

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TPattern, typename TSpec>
SEQAN_HOST_DEVICE inline void
clear(Finder2<Index<TText, TIndexSpec>, TPattern, TSpec> & finder)
{
    // NOTE(esiragusa): should clear() be called on text/patternIt?
    goRoot(textIterator(finder));
    // NOTE(esiragusa): if find wasn't called yet, patternIterator is uninitialized.
//    goBegin(patternIterator(finder));
    finder._score = typename Score_<TSpec>::Type();
}

// ----------------------------------------------------------------------------
// Function printState()
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TPattern, typename TSpec>
SEQAN_HOST_DEVICE inline void
printState(Finder2<Index<TText, TIndexSpec>, TPattern, Backtracking<HammingDistance, TSpec> > & finder)
{
    std::cout << "Text:        " << parentEdgeLabel(textIterator(finder)) << std::endl;
    std::cout << "Pattern:     " << value(patternIterator(finder)) << std::endl;
    std::cout << "Text Len:    " << repLength(textIterator(finder)) << std::endl;
    std::cout << "Pattern Len: " << position(patternIterator(finder)) + 1 << std::endl;
    std::cout << "Errors:      " << static_cast<unsigned>(getScore(finder)) << std::endl;
    std::cout << "Max errors:  " << static_cast<unsigned>(finder._scoreThreshold) << std::endl;
}

// ----------------------------------------------------------------------------
// Function find()
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TPattern, typename TDelegate>
SEQAN_HOST_DEVICE inline void
find(Finder2<Index<TText, TIndexSpec>, TPattern, FinderSTree> & finder,
     TPattern const & pattern,
     TDelegate & delegate)
{
    if (goDown(textIterator(finder), pattern))
    {
        // TODO(esiragusa): update patternIterator.
        delegate(finder);
    }
}

template <typename TText, typename TIndexSpec, typename TPattern, typename TSpec, typename TDelegate>
SEQAN_HOST_DEVICE inline void
find(Finder2<Index<TText, TIndexSpec>, TPattern, Backtracking<HammingDistance, TSpec> > & finder,
     TPattern const & pattern,
     TDelegate & delegate)
{
    typedef Index<TText, TIndexSpec>                                TIndex;
    typedef Backtracking<HammingDistance, TSpec>                    TFinderSpec;
    typedef typename TextIterator_<TIndex, TFinderSpec>::Type       TTextIterator;
    typedef typename PatternIterator_<TPattern, TFinderSpec>::Type  TPatternIterator;

    setPatternIterator(finder, begin(pattern));

    TTextIterator & textIt = textIterator(finder);
    TPatternIterator & patternIt = patternIterator(finder);

    finder._score = 0;
    finder._scoreThreshold = 1;

    do
    {
        // Exact case.
        if (finder._score == finder._scoreThreshold)
        {
            if (goDown(textIt, suffix(pattern, position(patternIt))))
            {
                delegate(finder);
            }

            goUp(textIt);

            // Termination.
            if (isRoot(textIt)) break;
        }

        // Approximate case.
        else if (finder._score < finder._scoreThreshold)
        {
            // Base case.
            if (atEnd(patternIt))
            {
                delegate(finder);
            }

            // Recursive case.
            else if (goDown(textIt))
            {
                finder._score += _getVertexScore(finder);
                goNext(patternIt);
                continue;
            }
        }

        // Backtrack.
        do
        {
            goPrevious(patternIt);
            finder._score -= _getVertexScore(finder);
        }
        while (!goRight(textIt) && goUp(textIt));

        // Termination.
        if (isRoot(textIt)) break;

        finder._score += _getVertexScore(finder);
        goNext(patternIt);
    }
    while (true);
}

}

#endif  // #ifndef SEQAN_EXTRAS_INDEX_FIND_INDEX_H_
