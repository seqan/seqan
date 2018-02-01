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
// Author: Rene Rahn <rene.rahn@fu-berlin.de>
// ==========================================================================

#ifndef INCLUDE_SEQAN_JOURNALED_STRING_TREE_JST_EXTENSION_MYERS_UKKONEN_H_
#define INCLUDE_SEQAN_JOURNALED_STRING_TREE_JST_EXTENSION_MYERS_UKKONEN_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class JstExtension, MyersUkkonen
// ----------------------------------------------------------------------------

template <typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec>
class JstExtension<Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > > :
    public JstExtensionBase<JstExtension<Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > >, ContextEnd>
{
public:

    using TPattern = Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> >;
    using TBase    = JstExtensionBase<JstExtension<TPattern>, ContextEnd>;
    using TWord    = typename TPattern::TWord;

    static constexpr TWord WORD_INDEX_HIGH_BIT = BitsPerValue<TWord>::VALUE - 1;

    TPattern&   _pattern;

    TWord       X, D0, HN, HP;  // Used for standard search.
    TWord       carryD0, carryHP, carryHN, temp;  // Additionally used for long pattern search.
    unsigned    shift, limit;  // Used for long pattern search.
    TWord       lastBit;
    TWord       carry;
    TWord       mask;
    bool        isLongNeedle;

    JstExtension(TPattern & pattern) : TBase(*this), _pattern(pattern)
    {
        // Explicitly force to initialize pattern due to implementation detail.
        state(*this) = _pattern;
        Nothing t;
        _patternInit(_pattern, state(*this), t);
        lastBit = static_cast<TWord>(1) << (_pattern.needleSize - 1);
        isLongNeedle = _pattern.largePattern != nullptr;
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction GetPatternState
// ----------------------------------------------------------------------------

template <typename TNeedle, typename TSpec, typename TFindBeginPatternSpec>
struct GetPatternState<JstExtension<Pattern<TNeedle, Myers<TSpec, True, TFindBeginPatternSpec> > > >
{
    using Type = PatternState_<TNeedle, Myers<TSpec, True, TFindBeginPatternSpec> >;
};

template <typename TNeedle, typename TSpec, typename TFindBeginPatternSpec>
struct GetPatternState<JstExtension<Pattern<TNeedle, Myers<TSpec, True, TFindBeginPatternSpec> > > const>
{
    using Type = PatternState_<TNeedle, Myers<TSpec, True, TFindBeginPatternSpec> > const;
};

// If state is disabled
template <typename TNeedle, typename TSpec, typename TFindBeginPatternSpec>
struct GetPatternState<JstExtension<Pattern<TNeedle, Myers<TSpec, False, TFindBeginPatternSpec> > > >
{
    using Type = Nothing;
};

template <typename TNeedle, typename TSpec, typename TFindBeginPatternSpec>
struct GetPatternState<JstExtension<Pattern<TNeedle, Myers<TSpec, False, TFindBeginPatternSpec> > > const> :
        public GetPatternState<JstExtension<Pattern<TNeedle, Myers<TSpec, False, TFindBeginPatternSpec> > > >{};

// ----------------------------------------------------------------------------
// Metafunction ProxySelectionMethod
// ----------------------------------------------------------------------------

template <typename TNeedle, typename TSpec, typename TState, typename TFindBeginPatternSpec>
struct ProxySelectionMethod<JstExtension<Pattern<TNeedle, Myers<TSpec, TState, TFindBeginPatternSpec> > > >
{
    using Type = SelectFirstProxy;
};

// ============================================================================
// Functions
// ============================================================================

namespace impl
{

// ----------------------------------------------------------------------------
// Function runLongNeedle()
// ----------------------------------------------------------------------------

template <typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec, typename TIterator>
inline std::pair<size_t, bool>
runLongNeedle(JstExtension<Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > > & me,
              TIterator hystkIt)
{
    using TWord = typename MyersLargeState_<TNeedle, TSpec>::TWord;

    auto& largePattern = *me._pattern.largePattern;
    auto& s = state(me);
    auto& largeState = *s.largeState;

    me.carryD0 = me.carryHN = 0;
    me.carryHP = static_cast<int>(MyersUkkonenHP0_<TSpec>::VALUE); // FIXME: replace Noting with TSpec

    // if the active cell is the last of it's block, one additional block has to be calculated
    me.limit = largeState.lastBlock + static_cast<unsigned>(largeState.scoreMask >> (me._pattern.MACHINE_WORD_SIZE - 1));

    if (me.limit == largePattern.blockCount)
        me.limit--;

    me.shift = largePattern.blockCount * ordValue(static_cast<typename Value< TNeedle >::Type>(getValue(hystkIt)));

    // computing the necessary blocks, carries between blocks following one another are stored
    for (TWord currentBlock = 0; currentBlock <= me.limit; currentBlock++)
    {
        me.X = me._pattern.bitMasks[me.shift + currentBlock] | largeState.VN[currentBlock];

        me.temp = largeState.VP[currentBlock] + (me.X & largeState.VP[currentBlock]) + me.carryD0;
        if (me.carryD0 != static_cast<TWord>(0))
            me.carryD0 = me.temp <= largeState.VP[currentBlock];
        else
            me.carryD0 = me.temp < largeState.VP[currentBlock];

        me.D0 = (me.temp ^ largeState.VP[currentBlock]) | me.X;
        me.HN = largeState.VP[currentBlock] & me.D0;
        me.HP = largeState.VN[currentBlock] | ~(largeState.VP[currentBlock] | me.D0);

        me.X = (me.HP << 1) | me.carryHP;
        me.carryHP = me.HP >> (me._pattern.MACHINE_WORD_SIZE - 1);

        largeState.VN[currentBlock] = me.X & me.D0;

        me.temp = (me.HN << 1) | me.carryHN;
        me.carryHN = me.HN >> (me._pattern.MACHINE_WORD_SIZE - 1);

        largeState.VP[currentBlock] = me.temp | ~(me.X | me.D0);

        // if the current block is the one containing the last active cell
        // the new number of errors is computed
        if (currentBlock == largeState.lastBlock)
        {
            if ((me.HP & largeState.scoreMask) != static_cast<TWord>(0))
                s.errors++;
            else if ((me.HN & largeState.scoreMask) != static_cast<TWord>(0))
                s.errors--;
        }
    }

    // updating the last active cell
    while (s.errors > s.maxErrors) {
        if ((largeState.VP[largeState.lastBlock] & largeState.scoreMask) != static_cast<TWord>(0))
            s.errors--;
        else if ((largeState.VN[largeState.lastBlock] & largeState.scoreMask) != static_cast<TWord>(0))
            s.errors++;

        largeState.scoreMask >>= 1;
        if (largeState.scoreMask == static_cast<TWord>(0))
        {
            largeState.lastBlock--;
//                if (IsSameType<TSpec, FindPrefix>::VALUE && largeState.lastBlock == (unsigned)-1)  // TODO(rmaerker): Check influence of PrefixFind -> We will not support this!
//                    break;
            largeState.scoreMask = static_cast<TWord>(1) << (me._pattern.MACHINE_WORD_SIZE - 1);
        }
    }

    if ((largeState.scoreMask == largePattern.finalScoreMask) && (largeState.lastBlock == largePattern.blockCount - 1))
    {
        return std::pair<size_t, bool>(1, true);
    }
    else {
        largeState.scoreMask <<= 1;
        if (!largeState.scoreMask) {
            largeState.scoreMask = 1;
            largeState.lastBlock++;
        }

        if ((largeState.VP[largeState.lastBlock] & largeState.scoreMask) != static_cast<TWord>(0))
            s.errors++;
        else if ((largeState.VN[largeState.lastBlock] & largeState.scoreMask) != static_cast<TWord>(0))
            s.errors--;
    }
    return std::pair<size_t, bool>(1, false);
}

// ----------------------------------------------------------------------------
// Function runShortNeedle()
// ----------------------------------------------------------------------------

template <typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec, typename TIterator>
inline std::pair<size_t, bool>
runShortNeedle(JstExtension<Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > > & me,
               TIterator hystkIt)
{
    using  TWord = typename Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> >::TWord;

    auto& s = state(me);

    // computing the blocks
    me.X = me._pattern.bitMasks[ordValue(static_cast<typename Value<TNeedle>::Type>(getValue(hystkIt)))] | s.VN0;

    me.D0 = ((s.VP0 + (me.X & s.VP0)) ^ s.VP0) | me.X;
    me.HN = s.VP0 & me.D0;
    me.HP = s.VN0 | ~(s.VP0 | me.D0);
    me.X = (me.HP << 1) | static_cast<TWord>(static_cast<int>(MyersUkkonenHP0_<TSpec>::VALUE)); // FIXME: replace Nothing by TSpec
    s.VN0 = me.X & me.D0;
    s.VP0 = (me.HN << 1) | ~(me.X | me.D0);

    if ((me.HP & me.lastBit) != static_cast<TWord>(0))
        s.errors++;
    else if ((me.HN & me.lastBit) != static_cast<TWord>(0))
        s.errors--;

    return std::pair<size_t, bool>(1, s.errors <= s.maxErrors);
}
}  // namespace impl

// ----------------------------------------------------------------------------
// Function run()
// ----------------------------------------------------------------------------

template <typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec, typename TIterator>
inline std::pair<size_t, bool>
run(JstExtension<Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > > & me,
    TIterator hystkIt)
{
    return (me.isLongNeedle) ? impl::runLongNeedle(me, hystkIt) : impl::runShortNeedle(me, hystkIt);
}
}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_JOURNALED_STRING_TREE_JST_EXTENSION_MYERS_UKKONEN_H_
