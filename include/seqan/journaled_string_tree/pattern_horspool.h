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

#ifndef INCLUDE_SEQAN_JOURNALED_STRING_TREE_PATTERN_HORSPOOL_H_
#define INCLUDE_SEQAN_JOURNALED_STRING_TREE_PATTERN_HORSPOOL_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TNeedle>
class Pattern2<TNeedle, Horspool> :
    public PatternBase<Pattern2<TNeedle, Horspool>, False, ContextRange>
{
public:
    typedef typename Size<TNeedle>::Type TSize;
    typedef PatternBase<Pattern2<TNeedle, Horspool>, False, ContextRange> TBase;

    Holder<TNeedle>      data_host;
    String<TSize>        data_map;

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

// ============================================================================
// Functions
// ============================================================================

template <typename TNeedle>
void
_reinitPattern(Pattern2<TNeedle, Horspool> & me)
{
    typedef typename Value<TNeedle>::Type TValue;
    typedef typename Size<TNeedle>::Type TSize;
    typedef typename ValueSize<TValue>::Type TValSize;

    TNeedle& ndl = needle(me);
    SEQAN_ASSERT_NOT(empty(ndl));

    TSize value_size = ValueSize<TValue>::VALUE;
    //make room for map
    resize(me.data_map, value_size);

    //fill map
    TValSize jump_width = length(ndl);
    arrayFill(begin(me.data_map, Standard()), begin(me.data_map, Standard()) + value_size, jump_width);

    typename Iterator<TNeedle, Standard>::Type it;
    it = begin(ndl, Standard());
    while (jump_width > 1)
    {
        me.data_map[ordValue(getValue(it))] = --jump_width;
        ++it;
    }
}


template <typename TNeedle, typename TIterator>
inline std::pair<TIterator, bool>
run(Pattern2<TNeedle, Horspool> & me,
    TIterator hystkBegin,
    TIterator hystkEnd)
{
    using TDiff = typename Difference<TIterator>::Type;
    SEQAN_ASSERT(!empty(needle(me)));

    // Sanity check: Range must have same size as needle!
    if (hystkEnd - hystkBegin != static_cast<TDiff>(length(needle(me))))
        return std::pair<TIterator, bool>(hystkBegin + me.data_map[ordValue(getValue(--hystkEnd))], false);

    auto ndlIt = end(needle(me), Standard());
    TIterator hstkIt = hystkEnd;
    while (getValue(--ndlIt) == getValue(--hstkIt) && hstkIt != hystkBegin)
    {}

    return std::pair<TIterator, bool>(hystkBegin + me.data_map[ordValue(getValue(--hystkEnd))],
                                      hstkIt == hystkBegin && getValue(ndlIt) == getValue(hstkIt));
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_JOURNALED_STRING_TREE_PATTERN_HORSPOOL_H_
