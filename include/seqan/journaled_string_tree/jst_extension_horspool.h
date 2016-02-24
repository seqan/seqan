// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2016, Knut Reinert, FU Berlin
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

#ifndef INCLUDE_SEQAN_JOURNALED_STRING_TREE_JST_EXTENSION_HORSPOOL_H_
#define INCLUDE_SEQAN_JOURNALED_STRING_TREE_JST_EXTENSION_HORSPOOL_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class JstExtension; Horspool
// ----------------------------------------------------------------------------

template <typename TNeedle>
class JstExtension<Pattern<TNeedle, Horspool> > :
    public JstExtensionBase<JstExtension<Pattern<TNeedle, Horspool> >, ContextRange>
{
public:
    typedef typename Size<TNeedle>::Type TSize;
    typedef JstExtensionBase<JstExtension<Pattern<TNeedle, Horspool> >, ContextRange> TBase;

    Pattern<TNeedle, Horspool>& _pattern;

    JstExtension(Pattern<TNeedle, Horspool> & pattern) : TBase(*this), _pattern(pattern)
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function run()
// ----------------------------------------------------------------------------

template <typename TNeedle, typename TIterator>
inline std::pair<size_t, bool>
run(JstExtension<Pattern<TNeedle, Horspool> > const & me,
    TIterator hystkBegin,
    TIterator hystkEnd)
{
    using TDiff = typename Difference<TIterator>::Type;
    SEQAN_ASSERT(!empty(needle(me._pattern)));

    // Sanity check: Range must have same size as needle!
    if (hystkEnd - hystkBegin != static_cast<TDiff>(length(needle(me._pattern)) - 1))
        return std::pair<size_t, bool>(me._pattern.data_map[ordValue(getValue(hystkEnd))], false);

    auto ndlIt = end(needle(me._pattern), Standard()) - 1;
    TIterator hstkIt = hystkEnd;
    while (hstkIt != hystkBegin && getValue(ndlIt) == getValue(hstkIt))
    {
        --ndlIt; --hstkIt;
    }

    return std::pair<size_t, bool>(me._pattern.data_map[ordValue(getValue(hystkEnd))],
                                   hstkIt == hystkBegin && getValue(ndlIt) == getValue(hstkIt));
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_JOURNALED_STRING_TREE_JST_EXTENSION_HORSPOOL_H_
