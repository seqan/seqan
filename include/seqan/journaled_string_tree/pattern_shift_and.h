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

#ifndef INCLUDE_SEQAN_JOURNALED_STRING_TREE_PATTERN_SHIFT_AND_H_
#define INCLUDE_SEQAN_JOURNALED_STRING_TREE_PATTERN_SHIFT_AND_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================


template <typename TPattern>
struct PatternStateShiftAnd_
{
    using TWord = unsigned int;
    static constexpr TWord MACHINE_WORD_SIZE = sizeof(TWord) * 8;

    String<TWord> prefSufMatch;        // Set of all the prefixes of needle that match a suffix of haystack (called D in "Navarro")
};

template <typename TNeedle>
class Pattern2<TNeedle, ShiftAnd> :
    public PatternBase<Pattern2<TNeedle, ShiftAnd>, True, ContextBegin>
{
public:
    using TBase  = PatternBase<Pattern2<TNeedle, ShiftAnd>, True, ContextBegin>;
    using TState = typename GetPatternState<Pattern2>::Type;
    using TWord  = typename TState::TWord;

    Holder<TNeedle>     data_host;
    String<TWord>       bitMasks;            // Look up table for each character in the alphabet (called B in "Navarro")
    TWord               needleLength;                // e.g., needleLength=33 --> blockCount=2 (iff w=32 bits)
    TWord               blockCount;                // #unsigned ints required to store needle

    Pattern2() {}

    template <typename TNeedle2>
    Pattern2(TNeedle2 && ndl, SEQAN_CTOR_DISABLE_IF(IsSameType<typename std::remove_reference<TNeedle2>::type const &, Pattern2 const &>)) : TBase(*this)
    {
        setHost(*this, std::forward<TNeedle2>(ndl));
        ignoreUnusedVariableWarning(dummy);
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

template <typename TNeedle>
struct GetPatternState<Pattern2<TNeedle, ShiftAnd> >
{
    using Type = PatternStateShiftAnd_<Pattern2<TNeedle, ShiftAnd> >;
};

template <typename TNeedle>
struct GetPatternState<Pattern2<TNeedle, ShiftAnd> const>
{
    using Type = PatternStateShiftAnd_<Pattern2<TNeedle, ShiftAnd> > const;
};

// ============================================================================
// Functions
// ============================================================================

template <typename TNeedle>
void
_reinitPattern(Pattern2<TNeedle, ShiftAnd> & me)
{
    using TWord = typename Pattern2<TNeedle, ShiftAnd>::TWord;
    using TValue = typename Value<TNeedle>::Type;

    me.needleLength = length(needle(me));
    if (me.needleLength < 1)
        me.blockCount = 1;
    else
        me.blockCount = (me.needleLength - 1) / BitsPerValue<TWord>::VALUE + 1;

    clear(me.bitMasks);
    resize(me.bitMasks, me.blockCount * ValueSize<TValue>::VALUE, 0, Exact());

    clear(state(me).prefSufMatch);
    resize(state(me).prefSufMatch, me.blockCount, 0, Exact());

    for (TWord j = 0; j < me.needleLength; ++j)
        me.bitMasks[me.blockCount * ordValue(convert<TValue>(getValue(needle(me), j))) + j / state(me).MACHINE_WORD_SIZE] |=
                    static_cast<TWord>(1) << (j % state(me).MACHINE_WORD_SIZE);

}

template <typename TNeedle, typename TIterator>
inline std::pair<TIterator, bool>
run(Pattern2<TNeedle, ShiftAnd> & me,
    TIterator hystkIt)
{
    using TWord = typename Pattern2<TNeedle, ShiftAnd>::TWord;
    using TValue = typename Value<TNeedle>::Type;

    state(me).prefSufMatch[0] = ((state(me).prefSufMatch[0] << 1) | (TWord)1) & me.bitMasks[ordValue(convert<TValue>(getValue(hystkIt)))];
    return std::pair<TIterator, bool>(++hystkIt, (state(me).prefSufMatch[0] & (static_cast<TWord>(1) << (me.needleLength - 1))) != 0);
}
}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_JOURNALED_STRING_TREE_PATTERN_SHIFT_AND_H_
