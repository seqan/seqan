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
// Author: Rene Rahn <rene.rahn@fu-berlin.de>
// ==========================================================================

#ifndef INCLUDE_SEQAN_JOURNALED_STRING_TREE_JST_EXTENSION_SHIFT_AND_H_
#define INCLUDE_SEQAN_JOURNALED_STRING_TREE_JST_EXTENSION_SHIFT_AND_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class PatternStateShiftAnd_
// ----------------------------------------------------------------------------

template <typename TPattern>
struct PatternStateShiftAnd_
{
    using TWord = typename TPattern::TWord;

    String<TWord> prefSufMatch;        // Set of all the prefixes of needle that match a suffix of haystack (called D in "Navarro")
};

// ----------------------------------------------------------------------------
// Class JstExtension<ShiftAnd>
// ----------------------------------------------------------------------------

template <typename TNeedle>
class JstExtension<Pattern<TNeedle, ShiftAnd> > :
    public JstExtensionBase<JstExtension<Pattern<TNeedle, ShiftAnd> >, ContextEnd>
{
public:
    using TPattern = Pattern<TNeedle, ShiftAnd>;
    using TBase    = JstExtensionBase<JstExtension<Pattern<TNeedle, ShiftAnd> >, ContextEnd>;
    using TWord    = typename TPattern::TWord;

    static constexpr TWord WORD_INDEX_HIGH_BIT = BitsPerValue<TWord>::VALUE - 1;

    TPattern&  _pattern;
    TWord      carry;
    TWord      mask;
    bool       isLongNeedle;

    JstExtension(TPattern & pattern) : TBase(*this), _pattern(pattern)
    {
        // Explicitly force to initialize pattern due to implementation detail.
        _patternInit(_pattern);
        state(*this).prefSufMatch = _pattern.prefSufMatch;
        mask = static_cast<TWord>(1) << ((_pattern.needleLength - 1) % BitsPerValue<TWord>::VALUE);
        isLongNeedle = _pattern.blockCount > 1;
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction GetPatternState
// ----------------------------------------------------------------------------

template <typename TNeedle>
struct GetPatternState<JstExtension<Pattern<TNeedle, ShiftAnd> > >
{
    using Type = PatternStateShiftAnd_<Pattern<TNeedle, ShiftAnd> >;
};

template <typename TNeedle>
struct GetPatternState<JstExtension<Pattern<TNeedle, ShiftAnd> > const>
{
    using Type = PatternStateShiftAnd_<Pattern<TNeedle, ShiftAnd> > const;
};

// ----------------------------------------------------------------------------
// Metafunction ProxySelectionMethod
// ----------------------------------------------------------------------------

template <typename TNeedle>
struct ProxySelectionMethod<JstExtension<Pattern<TNeedle, ShiftAnd> > >
{
    using Type = SelectFirstProxy;
};

// ============================================================================
// Functions
// ============================================================================

namespace impl
{

// ----------------------------------------------------------------------------
// Function impl::runLongNeedle()
// ----------------------------------------------------------------------------

template <typename TNeedle, typename TIterator>
inline std::pair<size_t, bool>
runLongNeedle(JstExtension<Pattern<TNeedle, ShiftAnd> > & me,
              TIterator const hystkIt)
{
    using TExt   = JstExtension<Pattern<TNeedle, ShiftAnd> >;
    using TWord  = typename TExt::TWord;
    using TValue = typename Value<TNeedle>::Type;

    me.carry = 1;
    for (TWord block = 0; block < me._pattern.blockCount; ++block)
    {
        bool newCarry = isBitSet(state(me).prefSufMatch[block], me.WORD_INDEX_HIGH_BIT);
        state(me).prefSufMatch[block] = ((state(me).prefSufMatch[block] << 1) | me.carry) &
                                        me._pattern.bitMasks[me._pattern.blockCount * ordValue(convert<TValue>(getValue(hystkIt))) + block];
        me.carry = newCarry;
    }
    return std::pair<size_t, bool>(1, (state(me).prefSufMatch[me._pattern.blockCount-1] & me.mask) != 0);
}

// ----------------------------------------------------------------------------
// Function impl::runShortNeedle()
// ----------------------------------------------------------------------------

template <typename TNeedle, typename TIterator>
inline std::pair<size_t, bool>
runShortNeedle(JstExtension<Pattern<TNeedle, ShiftAnd> > & me,
               TIterator const hystkIt)
{
    using TWord = typename Pattern<TNeedle, ShiftAnd>::TWord;
    using TValue = typename Value<TNeedle>::Type;

    state(me).prefSufMatch[0] = ((state(me).prefSufMatch[0] << 1) | static_cast<TWord>(1)) &
                                me._pattern.bitMasks[ordValue(convert<TValue>(getValue(hystkIt)))];

    return std::pair<size_t, bool>(1, (state(me).prefSufMatch[0] & me.mask) != 0);
}
}  // namespace impl

// ----------------------------------------------------------------------------
// Function run()
// ----------------------------------------------------------------------------

template <typename TNeedle, typename TIterator>
inline std::pair<size_t, bool>
run(JstExtension<Pattern<TNeedle, ShiftAnd> > & me,
    TIterator hystkIt)
{
    return (me.isLongNeedle) ? impl::runLongNeedle(me, hystkIt) : impl::runShortNeedle(me, hystkIt);
}
}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_JOURNALED_STRING_TREE_JST_EXTENSION_SHIFT_AND_H_
