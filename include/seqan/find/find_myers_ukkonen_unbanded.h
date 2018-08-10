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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_HEADER_FIND_MYERS_UKKONEN_UNBANDED_H
#define SEQAN_HEADER_FIND_MYERS_UKKONEN_UNBANDED_H

#include <seqan/find/find_myers_ukkonen_base.h>
#include <seqan/find/find_myers_ukkonen_pattern.h>

namespace seqan
{

//////////////////////////////////////////////////////////////////////////////
//overwrite FindBegin_ to define host member if find begin is switched on

template <typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec, typename TFindBeginPatternSpec2>
struct FindBegin_< Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> >, TFindBeginPatternSpec2>
{
private:
    typedef Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > TPattern;
    typedef typename FindBeginPattern<TPattern>::Type TFindBeginPattern;

public:
    TFindBeginPattern data_findBeginPattern;
//     Holder<TNeedle>    data_host;    //defines the
    typedef False HasHost;
};

template <typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec>
struct FindBegin_< Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> >, void>
{
    typedef False HasHost;
//need no findBegin if FindBeginPatternSpec is void
};


//////////////////////////////////////////////////////////////////////////////
// State Data
//////////////////////////////////////////////////////////////////////////////

// small state
template <typename TNeedle, typename TSpec>
struct MyersSmallState_
{
#ifdef SEQAN_SSE2_INT128
    typedef Sse2Int128 TWord;
#else
    typedef unsigned long TWord;
#endif

    TWord VP0;                    // VP[0] (saves one dereferentiation)
    TWord VN0;                    // VN[0]
    unsigned int errors;        // the current number of errors
    unsigned int maxErrors;        // the maximal number of errors allowed

    MyersSmallState_() : VP0(0), VN0(0), errors(0), maxErrors(0)
    {}
};

// large state
template <typename TNeedle, typename TSpec>
struct MyersLargeState_
{
#ifdef SEQAN_SSE2_INT128
    typedef Sse2Int128 TWord;
#else
    typedef unsigned long TWord;
#endif
    unsigned lastBlock;            // the block containing the last active cell
    String<TWord> VP;
    String<TWord> VN;
    TWord scoreMask;            // the mask with a bit set at the position of the last active cell

    MyersLargeState_() : lastBlock(0), scoreMask(0) {}
};

//////////////////////////////////////////////////////////////////////////////
// Pattern Data
//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle, typename TSpec>
struct MyersSmallPattern_
{
#ifdef SEQAN_SSE2_INT128
    typedef Sse2Int128 TWord;
#else
    typedef unsigned long TWord;
#endif

    Holder<TNeedle> data_host;   // Hold a reference to the needle
    String<TWord> bitMasks;        // encode the needle with bitmasks for each alphabet character
    unsigned needleSize;        // needle size

    MyersSmallPattern_():
        needleSize(0) {}
};

// large basic pattern
template <typename TNeedle, typename TSpec>
struct MyersLargePattern_
{
#ifdef SEQAN_SSE2_INT128
    typedef Sse2Int128 TWord;
#else
    typedef unsigned long TWord;
#endif

    unsigned blockCount;        // the number of blocks
    TWord finalScoreMask;        // a mask with a bit set on the position of the last row
};

//////////////////////////////////////////////////////////////////////////////
// Functions
//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec, typename TNeedle2>
inline void _patternFirstInit(Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > & pattern,
                              TNeedle2 & needle)
{

    typedef typename Pattern<TNeedle, Myers<TSpec, TFindBeginPatternSpec> >::TWord TWord;
    typedef typename Value<TNeedle>::Type TValue;

    pattern.needleSize = length(needle);
    unsigned blockCount = (pattern.needleSize + pattern.MACHINE_WORD_SIZE - 1) / pattern.MACHINE_WORD_SIZE;

    if (blockCount > 1)
    {
        if (pattern.largePattern == NULL)
            pattern.largePattern = new MyersLargePattern_<TNeedle, TSpec>();

        pattern.largePattern->blockCount = blockCount;
        pattern.largePattern->finalScoreMask = (TWord)1 << ((pattern.needleSize + pattern.MACHINE_WORD_SIZE - 1) % pattern.MACHINE_WORD_SIZE);
    }
    else
    {
        delete pattern.largePattern;
        pattern.largePattern = NULL;
    }

    clear(pattern.bitMasks);
    resize(pattern.bitMasks, (ValueSize<TValue>::VALUE + 1) * blockCount, 0, Exact());

    // encoding the letters as bit-vectors
    for (unsigned j = 0; j < pattern.needleSize; j++)
        pattern.bitMasks[
            blockCount * ordValue(convert<typename Value<TNeedle>::Type>(getValue(needle, j)))
            + j / pattern.MACHINE_WORD_SIZE
        ] |= (TWord)1 << (j % pattern.MACHINE_WORD_SIZE);
        //pattern.bitMasks[pattern.blockCount * ordValue((CompareType< Value< TNeedle >::Type, Value< Container< THaystack >::Type >::Type >::Type) needle[j]) + j/pattern.MACHINE_WORD_SIZE] = pattern.bitMasks[pattern.blockCount * ordValue((CompareType< Value< TNeedle >::Type, Value< Container< THaystack >::Type >::Type >::Type) needle[j]) + j/MACHINE_WORD_SIZE] | ((TWord)1 << (j%MACHINE_WORD_SIZE));

    _findBeginInit(pattern, needle);
}

template <typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec>
inline void _patternMatchNOfPatternImpl(Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > & pattern,
                                        bool match)
{
    typedef typename Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> >::TWord TWord;
    unsigned blockCount = (pattern.largePattern == NULL)? 1: pattern.largePattern->blockCount;

    // letters are encoded as bit-vectors
    for (unsigned j = 0; j < pattern.needleSize; j++)
    {
        TWord bit = (TWord)1 << (j % pattern.MACHINE_WORD_SIZE);
        bool allNull = true;
        int idx = j / pattern.MACHINE_WORD_SIZE;

        for (int i = 0; i < 4; ++i, idx += blockCount)
            allNull &= (pattern.bitMasks[idx] & bit) == (TWord)0;

        if (allNull)
        {    // all bits are 0 => this letter must be 'N'
            if (match)
            {
                for (; idx >= 0; idx -= blockCount)
                    pattern.bitMasks[idx] |= bit;
            } else
            {
                for (; idx >= 0; idx -= blockCount)
                    pattern.bitMasks[idx] &= ~bit;
            }
        }
    }
}


template <typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec>
inline void
_patternMatchNOfPattern(Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > & pattern, bool match)
{
    _patternMatchNOfPatternImpl(pattern, match);
    _patternMatchNOfPatternImpl(pattern.data_findBeginPattern, match);
}


template <typename TNeedle, typename TSpec, typename THasState>
inline void
_patternMatchNOfPattern(Pattern<TNeedle, Myers<TSpec, THasState, void> > & pattern, bool match)
{
    _patternMatchNOfPatternImpl(pattern, match);
}


template <typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec>
inline void
_patternMatchNOfFinderImpl(Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > & pattern, bool match)
{

    typedef typename Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> >::TWord TWord;
    unsigned blockCount = (pattern.largePattern == NULL)? 1: pattern.largePattern->blockCount;

    // letters are encoded as bit-vectors
    if (match)
    {
        for (unsigned j = 0; j < pattern.needleSize; j++)
            pattern.bitMasks[blockCount * 4 + j / pattern.MACHINE_WORD_SIZE] |= (TWord)1 << (j % pattern.MACHINE_WORD_SIZE);
    } else {
        for (unsigned j = 0; j < pattern.needleSize; j++)
            pattern.bitMasks[blockCount * 4 + j / pattern.MACHINE_WORD_SIZE] &= ~((TWord)1 << (j % pattern.MACHINE_WORD_SIZE));
    }
}


template <typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec>
inline void
_patternMatchNOfFinder(Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > & pattern, bool match)
{
    _patternMatchNOfFinderImpl(pattern, match);
    _patternMatchNOfFinderImpl(pattern.data_findBeginPattern, match);
}


template <typename TNeedle, typename TSpec, typename THasState>
inline void
_patternMatchNOfFinder(Pattern<TNeedle, Myers<TSpec, THasState, void> > & pattern, bool match)
{
    _patternMatchNOfFinderImpl(pattern, match);
}

template <typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec>
inline void
_reinitPattern(Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > & pattern)
{
    _patternFirstInit(pattern, needle(pattern));
}

template <typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec, typename TFinder>
inline bool
_patternInit(Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > const & pattern,
             PatternState_<TNeedle, Myers<TSpec, True, TFindBeginPatternSpec> > & state,
             TFinder &)
{
    typedef MyersLargeState_<TNeedle, TSpec> TLargeState;
    typedef typename TLargeState::TWord TWord;

    if (pattern.largePattern == NULL)
    {
        state.errors = pattern.needleSize;
        state.VP0 = ~(TWord)0;
        state.VN0 = 0;
        delete state.largeState;
        state.largeState = NULL;
    }
    else
    {
        if (state.largeState == NULL)
            state.largeState = new TLargeState;

        TLargeState &largeState = *state.largeState;
        // localMaxErrors either stores the maximal number of errors (me.maxErrors) or the needle size minus one.
        // It is used for the mask computation and setting the initial score (the minus one is there because of the Ukkonen trick).
        int localMaxErrors = _min(state.maxErrors, pattern.needleSize - 1);
        state.errors = localMaxErrors + 1;
        largeState.scoreMask = (TWord)1 << (localMaxErrors % pattern.MACHINE_WORD_SIZE);
        largeState.lastBlock = localMaxErrors / pattern.MACHINE_WORD_SIZE;
        if (largeState.lastBlock >= pattern.largePattern->blockCount)
            largeState.lastBlock = pattern.largePattern->blockCount - 1;

        clear(largeState.VP);
        resize(largeState.VP, pattern.largePattern->blockCount, ~(TWord)0, Exact());

        clear(largeState.VN);
        resize(largeState.VN, pattern.largePattern->blockCount, 0, Exact());
    }
    return true;
}

template <typename TNeedle, typename TSpec, typename TFindBeginPatternSpec, typename TFinder>
inline bool
_patternInit(Pattern<TNeedle, Myers<TSpec, True, TFindBeginPatternSpec> > & pattern, TFinder & finder)
{
    return _patternInit(pattern, pattern, finder);
}

//////////////////////////////////////////////////////////////////////////////
// Myers-Ukkonen for semi-global edit-distance-alignments
// (version for needles longer than one machineword)
//////////////////////////////////////////////////////////////////////////////

template <typename TFinder, typename TNeedle, typename TSpec, typename THasState, typename THasState2, typename TFindBeginPatternSpec, typename TSize>
inline bool _findMyersLargePatterns (TFinder & finder,
                                     Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > const & pattern,
                                     PatternState_<TNeedle, Myers<TSpec, THasState2, TFindBeginPatternSpec> > & state,
                                     TSize haystack_length)
{
    typedef MyersLargePattern_<TNeedle, TSpec> TLargePattern;
    typedef MyersLargeState_<TNeedle, TSpec> TLargeState;
    typedef typename TLargeState::TWord TWord;

    TWord X, D0, HN, HP, temp;
    TWord carryD0, carryHP, carryHN;
    unsigned shift, limit, currentBlock;

    TLargePattern &largePattern = *pattern.largePattern;
    TLargeState &largeState = *state.largeState;

    while (position(finder) < haystack_length)
    {
        carryD0 = carryHN = 0;
        carryHP = (int)MyersUkkonenHP0_<TSpec>::VALUE; // FIXME: replace Noting with TSpec

        // if the active cell is the last of it's block, one additional block has to be calculated
        limit = largeState.lastBlock + (unsigned)(largeState.scoreMask >> (pattern.MACHINE_WORD_SIZE - 1));

        if (limit == largePattern.blockCount)
            limit--;

        shift = largePattern.blockCount * ordValue((typename Value< TNeedle >::Type) *finder);

        // computing the necessary blocks, carries between blocks following one another are stored
        for (currentBlock = 0; currentBlock <= limit; currentBlock++)
        {
            X = pattern.bitMasks[shift + currentBlock] | largeState.VN[currentBlock];

            temp = largeState.VP[currentBlock] + (X & largeState.VP[currentBlock]) + carryD0;
            if (carryD0 != (TWord)0)
                carryD0 = temp <= largeState.VP[currentBlock];
            else
                carryD0 = temp < largeState.VP[currentBlock];

            D0 = (temp ^ largeState.VP[currentBlock]) | X;
            HN = largeState.VP[currentBlock] & D0;
            HP = largeState.VN[currentBlock] | ~(largeState.VP[currentBlock] | D0);

            X = (HP << 1) | carryHP;
            carryHP = HP >> (pattern.MACHINE_WORD_SIZE - 1);

            largeState.VN[currentBlock] = X & D0;

            temp = (HN << 1) | carryHN;
            carryHN = HN >> (pattern.MACHINE_WORD_SIZE - 1);

             largeState.VP[currentBlock] = temp | ~(X | D0);

            // if the current block is the one containing the last active cell
            // the new number of errors is computed
            if (currentBlock == largeState.lastBlock) {
                if ((HP & largeState.scoreMask) != (TWord)0)
                    state.errors++;
                else if ((HN & largeState.scoreMask) != (TWord)0)
                    state.errors--;
            }
        }

        // updating the last active cell
        while (state.errors > state.maxErrors) {
            if ((largeState.VP[largeState.lastBlock] & largeState.scoreMask) != (TWord)0)
                state.errors--;
            else if ((largeState.VN[largeState.lastBlock] & largeState.scoreMask) != (TWord)0)
                state.errors++;

            largeState.scoreMask >>= 1;
            if (largeState.scoreMask == (TWord)0)
            {
                largeState.lastBlock--;
                if (IsSameType<TSpec, FindPrefix>::VALUE && largeState.lastBlock == (unsigned)-1)
                    break;
                largeState.scoreMask = (TWord)1 << (pattern.MACHINE_WORD_SIZE - 1);
            }
        }

        if ((largeState.scoreMask == largePattern.finalScoreMask) && (largeState.lastBlock == largePattern.blockCount - 1))
        {
            _setFinderEnd(finder);
            SEQAN_IF_CONSTEXPR (IsSameType<TSpec, FindPrefix>::VALUE)
            {
                _setFinderLength(finder, endPosition(finder));
            }
            return true;
        }
        else {
            largeState.scoreMask <<= 1;
            if (!largeState.scoreMask) {
                largeState.scoreMask = 1;
                largeState.lastBlock++;
            }

            if ((largeState.VP[largeState.lastBlock] & largeState.scoreMask) != (TWord)0)
                state.errors++;
            else if ((largeState.VN[largeState.lastBlock] & largeState.scoreMask) != (TWord)0)
                state.errors--;
        }

//        SEQAN_ASSERT (state.errors >= 0);

        goNext(finder);
    }

    return false;
}


template <typename TFinder, typename TNeedle, typename TSpec, typename THasState, typename THasState2, typename TFindBeginPatternSpec, typename TSize>
inline bool
_findMyersSmallPatterns(
    TFinder & finder,
    Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > const & pattern,
    PatternState_<TNeedle, Myers<TSpec, THasState2, TFindBeginPatternSpec> > & state,
    TSize haystack_length)
{

    typedef typename Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> >::TWord TWord;

    TWord X, D0, HN, HP;
    TWord lastBit = (TWord)1 << (pattern.needleSize - 1);

    // computing the blocks
    while (position(finder) < haystack_length)
    {
        X = pattern.bitMasks[ordValue((typename Value<TNeedle>::Type) *finder)] | state.VN0;

        D0 = ((state.VP0 + (X & state.VP0)) ^ state.VP0) | X;
        HN = state.VP0 & D0;
        HP = state.VN0 | ~(state.VP0 | D0);
        X = (HP << 1) | (TWord)(int)MyersUkkonenHP0_<TSpec>::VALUE; // FIXME: replace Nothing by TSpec
        state.VN0 = X & D0;
        state.VP0 = (HN << 1) | ~(X | D0);

        if ((HP & lastBit) != (TWord)0)
            state.errors++;
        else if ((HN & lastBit) != (TWord)0)
            state.errors--;

        if (state.errors <= state.maxErrors)
        {
            _setFinderEnd(finder);
            SEQAN_IF_CONSTEXPR (IsSameType<TSpec, FindPrefix>::VALUE)
            {
                _setFinderLength(finder, endPosition(finder));
            }
            return true;
        }
        //
        // SEQAN_IF_CONSTEXPR (IsSameType<TSpec, FindPrefix>::VALUE)
        // {//limit haystack length during prefix search
        //
        // }

        goNext(finder);
    }

    return false;
}

//////////////////////////////////////////////////////////////////////////////
// find

template <typename TFinder, typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec>
inline bool find (TFinder & finder,
                  Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > const & pattern,
                  PatternState_<TNeedle, Myers<TSpec, True, TFindBeginPatternSpec> > & state)
{
    typedef typename Haystack<TFinder>::Type THaystack;
    typedef typename Size<THaystack>::Type TSize;

    SEQAN_UNUSED TSize prefix_begin_position; //for prefix search: the position where the prefix begins

    if (empty(finder))
    {
        _patternInit(pattern, state, finder);
        _finderSetNonEmpty(finder);

        prefix_begin_position = position(finder);

        //TODO: adapt myers-ukkonnen to dynamically change maxErrors
    }
    else
    {
        if (atEnd(finder)) return false;
        goNext(finder);

        prefix_begin_position = beginPosition(finder);
    }

    TSize haystack_length = length(container(finder));
    // limit search width for prefix search
    SEQAN_IF_CONSTEXPR (IsSameType<TSpec, FindPrefix>::VALUE)
    {
        TSize maxlen = prefix_begin_position + pattern.needleSize - scoreLimit(state) + 1;
        if (haystack_length > maxlen)
            haystack_length = maxlen;
    }

    // distinguish between the version for needles not longer than one machineword and the version for longer needles
    if (pattern.largePattern == NULL)
        return _findMyersSmallPatterns(finder, pattern, state, haystack_length);
    else
        return _findMyersLargePatterns(finder, pattern, state, haystack_length);
}

}// namespace seqan

#endif //#ifndef SEQAN_HEADER_FIND_MYERS_UKKONEN_UNBANDED_H
