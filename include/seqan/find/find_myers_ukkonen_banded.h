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

#ifndef SEQAN_HEADER_FIND_MYERS_UKKONEN_BANDED_H
#define SEQAN_HEADER_FIND_MYERS_UKKONEN_BANDED_H

#include <seqan/find/find_myers_ukkonen_base.h>
#include <seqan/find/find_myers_ukkonen_pattern.h>

namespace seqan
{

//////////////////////////////////////////////////////////////////////////////
// State Data
//////////////////////////////////////////////////////////////////////////////

// small state
template <typename TNeedle, typename TSmallAlphabet>
struct MyersSmallStateBandedShift_ {};
template <typename TNeedle>
struct MyersSmallStateBandedShift_<TNeedle, False> {
    typedef typename Value<TNeedle>::Type TValue;
    unsigned short shift[ValueSize<TValue>::VALUE];
};
template <typename TNeedle, typename TSpec, typename TFinderCSP, typename TPatternCSP>
struct MyersSmallState_<TNeedle, AlignTextBanded<TSpec, TFinderCSP, TPatternCSP> >:
    public MyersSmallStateBandedShift_<TNeedle, typename MyersSmallAlphabet_<typename Value<TNeedle>::Type>::Type>
{
#ifdef SEQAN_SSE2_INT128
    typedef Sse2Int128 TWord;
#else
#ifdef SEQAN_DEBUG_MYERSBITVECTOR
    typedef unsigned char TWord;
#else
    typedef unsigned long TWord;
#endif
#endif
    typedef typename Value<TNeedle>::Type TValue;

    TWord bitMasks[ValueSize<TValue>::VALUE];
    TWord VP0;                    // VP[0] (saves one dereferentiation)
    TWord VN0;                    // VN[0]
    unsigned short errors;      // the current number of errors
    unsigned short maxErrors;   // the maximal number of errors allowed
    unsigned short leftClip;    // clip that many characters from the text begin
//    unsigned short rightClip;   // stop alignment that many characters before the end   <<<< currently unused (autom. determined)

#ifdef SEQAN_DEBUG_MYERSBITVECTOR
    String<int> DPMat;
#endif
    MyersSmallState_() : bitMasks(), VP0(0), VN0(0), errors(0), maxErrors(0), leftClip(0) {}
};

// large state
template <typename TNeedle, typename TSpec, typename TFinderCSP, typename TPatternCSP>
struct MyersLargeState_<TNeedle, AlignTextBanded<TSpec, TFinderCSP, TPatternCSP> >
{
#ifdef SEQAN_SSE2_INT128
    typedef Sse2Int128 TWord;
#else
#ifdef SEQAN_DEBUG_MYERSBITVECTOR
    typedef unsigned char TWord;
#else
    typedef unsigned long TWord;
#endif
#endif
    unsigned lastBlock;            // the block containing the last active cell
    unsigned blockCount;        // the number of blocks
    String<TWord> VP;
    String<TWord> VN;
};


//////////////////////////////////////////////////////////////////////////////
// Pattern Data
//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle, typename TSpec, typename TFinderCSP, typename TPatternCSP>
struct MyersSmallPattern_<TNeedle, AlignTextBanded<TSpec, TFinderCSP, TPatternCSP> >
{
#ifdef SEQAN_SSE2_INT128
    typedef Sse2Int128 TWord;
#else
#ifdef SEQAN_DEBUG_MYERSBITVECTOR
    typedef unsigned char TWord;
#else
    typedef unsigned long TWord;
#endif
#endif

    Holder<TNeedle> data_host;  // needle holder (the banded version needs no preprocessed bitmasks)
};

// large basic pattern
template <typename TNeedle, typename TSpec, typename TFinderCSP, typename TPatternCSP>
struct MyersLargePattern_<TNeedle, AlignTextBanded<TSpec, TFinderCSP, TPatternCSP> > {};

//////////////////////////////////////////////////////////////////////////////
// Functions
//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle, typename TSpec, typename TFinderCSP, typename TPatternCSP, typename THasState, typename TFindBeginPatternSpec, typename TNeedle2>
inline void _patternFirstInit(Pattern<TNeedle, Myers<AlignTextBanded<TSpec, TFinderCSP, TPatternCSP>, THasState, TFindBeginPatternSpec> > & pattern,
                              TNeedle2 & ndl)
{
    _findBeginInit(pattern, ndl);
}

//____________________________________________________________________________
// bitmask operations - small alphabet

template <typename TNeedle, typename TSpec>
finline void
_myersPreInit(PatternState_<TNeedle, TSpec> &state, True)
{
    typedef typename Value<TNeedle>::Type TValue;
    for (unsigned i = 0; i < ValueSize<TValue>::VALUE; ++i)
        state.bitMasks[i] = 0;
}

template <typename TNeedle, typename TSpec>
finline void
_myersPostInit(PatternState_<TNeedle, TSpec> &state, True)
{
    typedef typename Value<TNeedle>::Type TValue;
    for (unsigned i = 0; i < ValueSize<TValue>::VALUE; ++i)
        state.bitMasks[i] >>= 1;
}

template <typename TNeedle, typename TSpec, typename TFinderCSP, typename TPatternCSP, typename TFindBeginPatternSpec, typename TValue, typename TShift>
finline void
_myersAdjustBitmask(PatternState_<TNeedle, Myers<AlignTextBanded<TSpec, TFinderCSP, TPatternCSP>, True, TFindBeginPatternSpec> > &state, TValue const value, TShift, True)
{
    typedef typename PatternState_<TNeedle, Myers<AlignTextBanded<TSpec, TFinderCSP, TPatternCSP>, True, TFindBeginPatternSpec> >::TWord TWord;

    // compiler will optimize that
    if (IsSameType<TPatternCSP, NMatchesAll_>::VALUE && value == unknownValue<TValue>())
    {
        for (unsigned i = 0; i < ValueSize<TValue>::VALUE; ++i)
            state.bitMasks[i] = (state.bitMasks[i] >> 1) | ((TWord)1 << (BitsPerValue<TWord>::VALUE - 1));
    }
    else
    {
        for (unsigned i = 0; i < ValueSize<TValue>::VALUE; ++i)
            state.bitMasks[i] >>= 1;
        if (!(IsSameType<TPatternCSP, NMatchesNone_>::VALUE && value == unknownValue<TValue>()))
            state.bitMasks[ordValue(value)] |= (TWord)1 << (BitsPerValue<TWord>::VALUE - 1);
    }

    SEQAN_IF_CONSTEXPR (IsSameType<TFinderCSP, NMatchesAll_>::VALUE)
        state.bitMasks[ordValue(unknownValue<TValue>())] |= (TWord)1 << (BitsPerValue<TWord>::VALUE - 1);
    SEQAN_IF_CONSTEXPR (IsSameType<TFinderCSP, NMatchesNone_>::VALUE)
        state.bitMasks[ordValue(unknownValue<TValue>())] &= ~((TWord)1 << (BitsPerValue<TWord>::VALUE - 1));
}

template <typename TNeedle, typename TSpec, typename TValue, typename TShift>
finline typename PatternState_<TNeedle, TSpec>::TWord
_myersGetBitmask(PatternState_<TNeedle, TSpec> &state, TValue const value, TShift, True)
{
    return state.bitMasks[ordValue(value)];
}

//____________________________________________________________________________
// bitmask operations - large alphabet

template <typename TNeedle, typename TSpec>
finline void
_myersPreInit(PatternState_<TNeedle, TSpec> &state, False)
{
    typedef typename Value<TNeedle>::Type TValue;
    memset(state.bitMasks, 0, (ValueSize<TValue>::VALUE + 1) * sizeof(state.bitMasks[0]));
    memset(state.shift, 0, ValueSize<TValue>::VALUE * sizeof(state.shift[0]));
}

template <typename TNeedle, typename TSpec>
finline void
_myersPostInit(PatternState_<TNeedle, TSpec> &, False)
{
}

template <typename TNeedle, typename TSpec, typename TFinderCSP, typename TPatternCSP, typename TFindBeginPatternSpec, typename TValue, typename TShift>
finline void
_myersAdjustBitmask(PatternState_<TNeedle, Myers<AlignTextBanded<TSpec, TFinderCSP, TPatternCSP>, True, TFindBeginPatternSpec> > &state, TValue const value, TShift const shift, False)
{
    typedef typename PatternState_<TNeedle, Myers<AlignTextBanded<TSpec, TFinderCSP, TPatternCSP>, True, TFindBeginPatternSpec> >::TWord TWord;

    if (IsSameType<TPatternCSP, NMatchesNone_>::VALUE && value == unknownValue<TValue>())
        return;

    unsigned ord = ordValue(value);
    unsigned short x = shift - state.shift[ord];
    if (x < BitsPerValue<TWord>::VALUE)
        state.bitMasks[ord] = (state.bitMasks[ord] >> x) | ((TWord)1 << (BitsPerValue<TWord>::VALUE - 1));
    else
        state.bitMasks[ord] = (TWord)1 << (BitsPerValue<TWord>::VALUE - 1);
    state.shift[ord] = shift;
}

template <typename TNeedle, typename TSpec, typename TFinderCSP, typename TPatternCSP, typename TFindBeginPatternSpec, typename TValue, typename TShift>
finline typename PatternState_<TNeedle, Myers<AlignTextBanded<TSpec, TFinderCSP, TPatternCSP>, True, TFindBeginPatternSpec> >::TWord
_myersGetBitmask(PatternState_<TNeedle, Myers<AlignTextBanded<TSpec, TFinderCSP, TPatternCSP>, True, TFindBeginPatternSpec> > &state, TValue const value, TShift const shift, False)
{
    typedef typename PatternState_<TNeedle, Myers<AlignTextBanded<TSpec, TFinderCSP, TPatternCSP>, True, TFindBeginPatternSpec> >::TWord TWord;

    if (IsSameType<TFinderCSP, NMatchesNone_>::VALUE && value == unknownValue<TValue>())
        return 0;

    if (IsSameType<TFinderCSP, NMatchesAll_>::VALUE && value == unknownValue<TValue>())
        return static_cast<TWord>((shift < BitsPerValue<TWord>::VALUE) ? (~0u << shift) : ~0u);

    unsigned ord = ordValue(value);
    TWord res;
    TShift x = shift - state.shift[ord];
    if (x < BitsPerValue<TWord>::VALUE)
        res = state.bitMasks[ord] >> x;
    else
        res = 0;

    SEQAN_IF_CONSTEXPR (IsSameType<TPatternCSP, NMatchesAll_>::VALUE)
    {
        ord = ordValue(unknownValue<TValue>());
        x = shift - state.shift[ord];
        if (x < BitsPerValue<TWord>::VALUE)
            res |= state.bitMasks[ord] >> x;
    }
    return res;
}

template <typename TFinder, typename TNeedle, typename TNeedle2, typename TSpec, typename TFinderCSP, typename TPatternCSP, typename TFindBeginPatternSpec>
inline bool
_patternInitSmallStateBanded(
    TFinder &finder,
    TNeedle2 const & needle,
    PatternState_<TNeedle, Myers<AlignTextBanded<TSpec, TFinderCSP, TPatternCSP>, True, TFindBeginPatternSpec> > & state)
{
    typedef Pattern<TNeedle, Myers<AlignTextBanded<TSpec, TFinderCSP, TPatternCSP>, True, TFindBeginPatternSpec> > TPattern;
    typedef typename TPattern::TWord TWord;
    typedef typename Iterator<TNeedle2 const, Standard>::Type TIter;
    typedef typename Value<TNeedle>::Type TValue;

#ifdef SEQAN_DEBUG_MYERSBITVECTOR
    int col = state.leftClip + 1;
#ifdef SEQAN_DEBUG_MYERSBITVECTOR_DUMP
    std::cerr << "     ";
    for (int i = length(needle); i != 0; --i)
        std::cerr << std::setw(5) << needle[i - 1];
    std::cerr << std::endl;
    std::cerr << "     ";
    for (int i = length(needle); i >= 0; --i)
        std::cerr << std::setw(5) << state.DPMat[i];
    std::cerr << std::endl;
#endif
#endif

    _myersPreInit(state, typename MyersSmallAlphabet_<TValue>::Type());

    typename Size<TNeedle>::Type const ndlLength = length(needle);

    // Initialize row 0 either with zeros or increasing numbers
    // This can be realized using the following DP pattern and
    // assuming character mismatches at rows -1, -2,...
    // Thus we initialize the bitmasks and VN with 0.
    // VP depends on global/local alignment
    //
    //  0  1  2  3  4   -2 -2 -2 -2 -2   (row -2)
    //  0  1  2  3  4   -1 -1 -1 -1 -1   (row -1)
    //  0  1  2  3  4    0  0  0  0  0   (row  0)
    //  1                1
    //        global           local
    //
    //  VP = 100...      VP = 111...
    //

    TWord VP = (MyersUkkonenHP0_<TSpec>::VALUE == 1)? (TWord)1 << ((int)BitsPerValue<TWord>::VALUE-1): std::numeric_limits<TWord>::max(); // HP[0]==1 <-> global, HP[0]==0 <-> local
    TWord VN = 0;

    // Errors are counted along the lowest diagonal and the
    // lowest row of the band.
    // The errors in the top-left corner are 0.
    //
    // 0 * * * *
    //   x * * * *
    //     x * * * *
    //       x x x x x
    //
    //       |-------|
    //     diagWidth + 1 = 5
    //
    // diagWidth = length(container(finder)) + state.leftClip + state.rightClip - length(needle)

    unsigned errors = 0;
    TIter ndlIter = begin(needle, Standard());
    TIter ndlEnd;

    // The errors along the diagonal can only increase or stay the same.
    // There is only the last row of length diagWidth where errors can decrease.
    // If errors exceeds cutOff it cannot reach maxErrors again.


    typename Size<TFinder>::Type const columns = length(container(finder)) + state.leftClip;
    unsigned cutOff = state.maxErrors;
    if (columns > ndlLength)
    {
        cutOff += columns - ndlLength;        // clipping case *0
        ndlEnd = end(needle, Standard());
    } else {
        errors += ndlLength - columns;
        ndlEnd = ndlIter + columns;            // clipping case *1
    }

//    std::cerr<<std::hex<<"\t  "<<std::setw(17)<<' '<<"\tVN"<<std::setw(17)<<VN<<"\tVP"<<std::setw(17)<<VP<<std::dec<<std::endl;

    unsigned short shift = 0;

    if (state.leftClip != 0)
    {
        //////////////////////////////////////////////////////////////////
        // PART 0: go down the parallelogram in a empty (clipped) area
        //////////////////////////////////////////////////////////////////

        errors += state.leftClip;
        if (errors > ndlLength) errors = ndlLength;
        if (errors > cutOff) return false;

    // leftClip = 2
    //   |-|
    //
    //   . . * * *
    //     . * * * *
    //       * * * * *
    //         * * * * .
    //           * * * . .
    //             * * . . .
    //
    //                 |---|
    //               rightClip = 3
    //
    // We divide the parallelogam into 3 sections:
    //
    //   A A A A
    //     A A A B
    //       A A B B
    //         A B B C
    //           B B C C
    //             B C C C
    //               C C C C
    //
    // Depending on where the clipping ends we identify 4 different clipping cases:
    //
    //     case 00            case 10            case 01            case 11
    //   . . * *            . . . .            . . * *            . . . .
    //     . * * *            . . . *            . * * *            . . . *
    //       * * * *            . . * *            * * * .            . . * .
    //         * * * *            . * * *            * * . .            . * . .
    //           * * * .            * * * .            * . . .            * . . .
    //             * * . .            * * . .            . . . .            . . . .
    //

        // adjust bitmasks (errors = number of needle chars to preprocess)
        for (; shift < errors; ++ndlIter, ++shift)
            _myersAdjustBitmask(state, *ndlIter, shift, typename MyersSmallAlphabet_<TValue>::Type());

        // initialise left column with
        //
        //  0  1  2  3  4   -2 -2 -2 -2 -2
        //  0  1  2  3  4   -1 -1 -1 -1 -1
        //  0  1  2  3  4    0  0  0  0  0
        //  1                1
        //  2   global       2   local
        //  3                3
        //  4                4
        //
        //  VP = 111100...   VP = 111111...
        if (errors < (unsigned)BitsPerValue<TWord>::VALUE-1)
            VP |= ((TWord) -1) << ((unsigned)BitsPerValue<TWord>::VALUE-1 - errors);
        else
            VP = (TWord)-1;
    }

    for (; ndlIter != ndlEnd; ++ndlIter, goNext(finder), ++shift)
    {
        //////////////////////////////////////////////////////////////////
        // PART 1: go down the parallelogram
        //////////////////////////////////////////////////////////////////

        // adjust bitmask
        _myersAdjustBitmask(state, *ndlIter, shift, typename MyersSmallAlphabet_<TValue>::Type());

        /////////////////////////
        // DIAGONAL MYERS CORE

        // VP/VN --> D0  (original Myers)
        TWord X = _myersGetBitmask(state, ordValue(*finder), shift, typename MyersSmallAlphabet_<TValue>::Type()) | VN;
        TWord D0 = ((VP + (X & VP)) ^ VP) | X;

        // adjust errors corresponding to rightmost bit of D0
        errors += (~D0 >> (BitsPerValue<TWord>::VALUE - 1)) & 1;
        if (errors > cutOff) return false;

        // D0 --> HP/HN  (original Myers)
        TWord HN = VP & D0;
        TWord HP = VN | ~(VP | D0);
    //    const int PADDING = sizeof(TWord)*2 + 1;
    //    std::cerr << std::hex;
    //    std::cerr << "\tD0"<<std::setw(PADDING)<<(uint64_t)D0<<"\tHN"<<std::setw(PADDING)<<(uint64_t)HN<<"\tHP"<<std::setw(PADDING)<<(uint64_t)HP << std::endl;

        // moving register down corresponds to shifting HP/HN up (right shift)
        // HP/HN --> shift --> VP/VN (modified Myers)
        X = D0 >> 1;
        VN = X & HP;
        VP = HN | ~(X | HP);
    //    std::cerr << "\t  "<<std::setw(PADDING)<<' '<<"\tVN"<<std::setw(PADDING)<<(uint64_t)VN<<"\tVP"<<std::setw(PADDING)<<(uint64_t)VP << std::endl;
    //    std::cerr << std::dec;

#ifdef SEQAN_DEBUG_MYERSBITVECTOR
#ifdef SEQAN_DEBUG_MYERSBITVECTOR_DUMP
        std::cerr << "diag ";
#endif
        int val = errors;
        state.DPMat[(col-state.leftClip)*(length(needle)+1)+col] = val;
        for (int i = length(needle); i >=0; --i)
        {
            if (i > col)
            {
#ifdef SEQAN_DEBUG_MYERSBITVECTOR_DUMP
                std::cerr << "     ";
#endif
            } else
            {
                int shft = (int)BitsPerValue<TWord>::VALUE-1 - (col-i);
                if (shft >= 0)
                {
                    if (i < col)
                    {
                        TWord mask = (TWord)1 << (shft);
                        val -= ((VP & mask) != (TWord)0)? 1:0;
                        val += ((VN & mask) != (TWord)0)? 1:0;
                    }
                    state.DPMat[(col-state.leftClip)*(length(needle)+1)+i] = val;
#ifdef SEQAN_DEBUG_MYERSBITVECTOR_DUMP
                    std::cerr << std::setw(5) << val;
                } else
                {
                    std::cerr << "     ";
#endif
                }
            }
        }
#ifdef SEQAN_DEBUG_MYERSBITVECTOR_DUMP
        std::cerr << std::setw(5) << *finder;
        std::cerr << std::setw(5) << errors << std::endl;
#endif
        ++col;
#endif
    }
    state.VP0 = VP;
    state.VN0 = VN;
    state.errors = errors;
    _myersPostInit(state, typename MyersSmallAlphabet_<TValue>::Type());
    return true;
}

template <typename TFinder, typename TNeedle, typename TNeedle2, typename TSpec, typename TFinderCSP, typename TPatternCSP, typename TFindBeginPatternSpec>
bool _stateInit(
    TFinder &finder,
    TNeedle const & needle,
    PatternState_<TNeedle2, Myers<AlignTextBanded<TSpec, TFinderCSP, TPatternCSP>, True, TFindBeginPatternSpec> > & state)
{
    typedef PatternState_<TNeedle2, Myers<AlignTextBanded<TSpec, TFinderCSP, TPatternCSP>, True, TFindBeginPatternSpec> > TState;
    typedef typename TState::TLargeState TLargeState;

//    unsigned diagWidth = length(container(finder)) + state.leftClip - length(needle);
//    unsigned blockCount = diagWidth / state.MACHINE_WORD_SIZE + 1;
    unsigned blockCount = 1;

//    SEQAN_ASSERT_GEQ(length(container(finder)), length(needle));

#ifdef SEQAN_DEBUG_MYERSBITVECTOR
    clear(state.DPMat);
    resize(state.DPMat, (length(container(finder)) + 1) * (length(needle) + 1), -9);
    for (unsigned i = 0; i <= length(needle); ++i)
        state.DPMat[i] = i;
    for (unsigned i = 0; i <= length(container(finder)); ++i)
        state.DPMat[i * (length(needle) + 1)] = 0;
#endif

    if (blockCount <= 1)
    {
        delete state.largeState;
        state.largeState = NULL;
        return _patternInitSmallStateBanded(finder, needle, state);
    }
    else
    {
        // TODO: is that good here?
        if (state.largeState == NULL)
            state.largeState = new TLargeState;

        TLargeState &largeState = *state.largeState;
        largeState.blockCount = blockCount;

        clear(largeState.VP);
        resize(largeState.VP, blockCount, ~0, Exact());

        clear(largeState.VN);
        resize(largeState.VN, blockCount, 0, Exact());

        state.errors = 0;

        return true;
    }
}

//////////////////////////////////////////////////////////////////////////////
// Myers-Ukkonen as a banded alignment
// the band includes the main diagonal and the diagonals above
// the band width is (blockCount * MACHINE_WORD_SIZE)
//////////////////////////////////////////////////////////////////////////////

template <
    typename TFinder,
    typename TNeedle,
    typename TNeedle2,
    typename TSpec,
    typename TFinderCSP,
    typename TPatternCSP,
    typename TFindBeginPatternSpec,
    typename TDoPatternSearch
>
inline bool
_findMyersSmallPatternsBanded(
    TFinder & finder,
    TNeedle const & needle,
    PatternState_<TNeedle2, Myers<AlignTextBanded<TSpec, TFinderCSP, TPatternCSP>, True, TFindBeginPatternSpec> > & state,
    TDoPatternSearch const)
{
    typedef PatternState_<TNeedle, Myers<AlignTextBanded<TSpec, TFinderCSP, TPatternCSP>, True, TFindBeginPatternSpec> > TState;
    typedef typename TState::TWord TWord;
    typedef typename Value<TNeedle>::Type TValue;

    TWord VP = state.VP0;
    TWord VN = state.VN0;
    TWord errors = state.errors;
    TWord const maxErrors = state.maxErrors;
    unsigned short const shift = length(needle);

#ifdef SEQAN_DEBUG_MYERSBITVECTOR
    unsigned col = position(finder) + 1;
#endif

    for (; !atEnd(finder); goNext(finder))
    {
        // PART 2: go right

        // normal Myers
        TWord X = _myersGetBitmask(state, ordValue(*finder), shift, typename MyersSmallAlphabet_<TValue>::Type()) | VN;
        TWord D0 = ((VP + (X & VP)) ^ VP) | X;
        TWord HN = VP & D0;
        TWord HP = VN | ~(VP | D0);
    //    const int PADDING = sizeof(TWord)*2 + 1;
    //    std::cerr << std::hex;
    //    std::cerr << "\tD0"<<std::setw(PADDING)<<(uint64_t)D0<<"\tHN"<<std::setw(PADDING)<<(uint64_t)HN<<"\tHP"<<std::setw(PADDING)<<(uint64_t)HP<<std::endl;
        X = (HP << 1) | 1;
        VN = X & D0;
        VP = (HN << 1) | ~(X | D0);
    //    std::cerr << "\t  "<<std::setw(PADDING)<<' '<<"\tVN"<<std::setw(PADDING)<<(uint64_t)VN<<"\tVP"<<std::setw(PADDING)<<(uint64_t)VP<<std::endl;
    //    std::cerr << std::dec;
        errors += (HP >> (BitsPerValue<TWord>::VALUE - 2)) & 1;
        errors -= (HN >> (BitsPerValue<TWord>::VALUE - 2)) & 1;

        // shift bitmasks and states
#ifdef SEQAN_DEBUG_MYERSBITVECTOR
#ifdef SEQAN_DEBUG_MYERSBITVECTOR_DUMP
        std::cerr << "horiz";
#endif
        int val = errors;
        state.DPMat[col*(length(needle)+1)+length(needle)] = val;
        for (int i = length(needle); i >= 0; --i)
        {
            int shft = (int)BitsPerValue<TWord>::VALUE-1 - (length(needle)-i);
            if (shft >= 0)
            {
                if (i < (int)length(needle))
                {
                    TWord mask = (TWord)1 << (shft);
                    val -= ((VP & mask) != (TWord)0)? 1:0;
                    val += ((VN & mask) != (TWord)0)? 1:0;
                }
                state.DPMat[col*(length(needle)+1)+i] = val;
#ifdef SEQAN_DEBUG_MYERSBITVECTOR_DUMP
                std::cerr << std::setw(5) << val;
            } else {
                std::cerr << "     ";
#endif
            }
        }
#ifdef SEQAN_DEBUG_MYERSBITVECTOR_DUMP
        std::cerr << std::setw(5) << *finder;
        std::cerr << std::setw(5) << errors << std::endl;
#endif
        ++col;
#endif

        if (TDoPatternSearch::VALUE)
        {
            // pattern search
            if (errors <= maxErrors)
            {
                state.VP0 = VP;
                state.VN0 = VN;
                state.errors = errors;
                _setFinderEnd(finder);
                SEQAN_IF_CONSTEXPR (IsSameType<TSpec, FindPrefix>::VALUE)
                {
                    _setFinderLength(finder, endPosition(finder));
                }
                return true;
            }
        }
        else
        {
            // edit distance
        }
    }

    if (!TDoPatternSearch::VALUE)
    {
        // edit distance
        state.errors = errors;
    }
    return false;
}

template <typename TSeq1, typename TSeq2>
inline unsigned
_computeEditDistanceBanded(
    TSeq1 const &seq1,
    TSeq2 const &seq2,
    unsigned maxErrors)
{
    PatternState_<TSeq2, MyersUkkonenGlobalBanded> state;
    state.maxErrors = maxErrors;
    state.leftClip = (length(seq2) - length(seq1) + maxErrors) / 2;
    typename Iterator<TSeq1 const, Rooted>::Type seq1Iter = begin(seq1, Rooted());

    if (!_patternInitSmallStateBanded(seq1Iter, seq2, state))
        return maxErrors + 1;
    _findMyersSmallPatternsBanded(seq1Iter, seq2, state, False());
    return state.errors;
}

//////////////////////////////////////////////////////////////////////////////
// find

template <typename TFinder, typename TNeedle, typename TNeedle2, typename TSpec, typename TFinderCSP, typename TPatternCSP, typename TFindBeginPatternSpec>
inline bool find (TFinder & finder,
                  TNeedle const & needle,
                  PatternState_<TNeedle2, Myers<AlignTextBanded<TSpec, TFinderCSP, TPatternCSP>, True, TFindBeginPatternSpec> > & state)
{
    if (empty(finder))
    {
        _finderSetNonEmpty(finder);
        if (!_stateInit(finder, needle, state))
        {
            goEnd(finder);
            return false;
        }

        if (state.errors <= state.maxErrors)
        {
            goPrevious(finder);
            _setFinderEnd(finder);
            SEQAN_IF_CONSTEXPR (IsSameType<TSpec, FindPrefix>::VALUE)
            {
                _setFinderLength(finder, endPosition(finder));
            }
            return true;
        }
        //TODO: adapt myers-ukkonnen to dynamically change maxErrors
    }
    else
    {
        if (atEnd(finder)) return false;
        goNext(finder);
    }

    // distinguish between the version for needles not longer than one machineword and the version for longer needles
    if (state.largeState == NULL)
        return _findMyersSmallPatternsBanded(finder, needle, state, True());
//    else
//        return _findMyersLargePatterns(finder, needle, state);
    return false;
}

// First two for AlignTextBanded
template <typename TFinder, typename TNeedle, typename TSpec, typename TFinderCSP, typename TPatternCSP, typename THasState, typename TFindBeginPatternSpec>
inline bool find (TFinder & finder,
                  Pattern<TNeedle, Myers<AlignTextBanded<TSpec, TFinderCSP, TPatternCSP>, THasState, TFindBeginPatternSpec> > const & pattern,
                  PatternState_<TNeedle, Myers<AlignTextBanded<TSpec, TFinderCSP, TPatternCSP>, True, TFindBeginPatternSpec> > & state)
{
    return find(finder, host(pattern), state);
}

} // namespace seqan

#endif // #ifndef SEQAN_HEADER_FIND_MYERS_UKKONEN_BANDED_H
