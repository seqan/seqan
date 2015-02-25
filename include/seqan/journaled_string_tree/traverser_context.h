// ==========================================================================
// traverser_context.h
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

#ifndef INCLUDE_SEQAN_JOURNALED_STRING_TREE_TRAVERSER_CONTEXT_H_
#define INCLUDE_SEQAN_JOURNALED_STRING_TREE_TRAVERSER_CONTEXT_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// TODO(rmaerker): Docu!
template <typename TObject, typename TSpec = void>
struct TraverserContext;

template <typename TReference, typename TSpec, typename TConfig, typename TSpec>
struct TraverserContext<JournaledStringTree<TReference, TSpec, TConfig>, TSpec>
{
    typedef JournaleStringTree<TReference, TSpec, TConfig>          TJst;
    typedef typename Member<TJst, JstReference>::Type               TReference;
    typedef Range<typename Iterator<TReference, Standard>::Type >   TRefRange;
    typedef String<TReference, JournaledString<> >                  TContextDataValue;
    typedef StringSet<TContextString, Owner<JournaledSet> >         TContextData;
    typedef typename Position<TContextDataValue>::Type              TPos;
    typedef String<TPos>                                            TOffsetString;

    // Member variables.

    TJst const *    _jstPtr;   // Pointer to the underlying journaled string tree.
    TRefRange       _refRange; // Range over the reference defining the context view.
    TContextData    _data;     // The object holding the data in this context.
    TOffsetString   _offsets;  // Keeps track of offsets in case of dynamic context generation.

    // Member functions.

    TraverserContext() : _jstPtr(nullptr)
    {}

    TraverserContext(TJst const & jst) : _jstPtr(&jst)
    {
        SEQAN_ASSERT_NOT(empty(host(jst)));
        toRange(_refRange, begin(host(jst)), end(host(jst)));

        auto res = impl::fetchContextDataForward(*this);
        SEQAN_ASSERT(res);
    }

    template <typename TSize>
    TraverserContext(TJst const & jst, TSize blockSize) : _jstPtr(&jst)
    {
        SEQAN_ASSERT_NOT(empty(host(jst)));
        toRange(_refRange, begin(host(jst)), begin(host(jst)) + blockSize);

        auto res = impl::fetchContextDataForward(*this);
        SEQAN_ASSERT(res);
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Private Functions
// ============================================================================

namespace impl
{

template <typename THostData, typename TSpec>
inline bool fetchContextDataForward(TraverserContext<THostData, TSpec> & me)
{
    // TODO(rmaerker): Here we load the context data in forward direction.
}

template <typename THostData, typename TSpec>
inline bool fetchContextDataBackward(TraverserContext<THostData, TSpec> & me)  // Only for specific type.
{
    // TODO(rmaerker): Here we load the context data in backward direction.
}

}

// ============================================================================
// Public Functions
// ============================================================================

template <typename TReference, typename TSpec, typename TConfig, typename TSpec, typename TSize>
inline bool fetchNext(TraverserContext<JournaledStringTree<TReference, TSpec, TConfig>, TSpec> & me,
                      TSize blockSize)
{
    // TODO(rmaerker): Implement me!
    return fetchContextDataForward(me);
}

// TODO(rmaerker): Disable for default JST
template <typename TReference, typename TSpec, typename TConfig, typename TSpec, typename TSize>
inline bool fetchPrevious(TraverserContext<JournaledStringTree<TReference, TSpec, TConfig>, TSpec> & me,
                          TSize blockSize)
{
    // TODO(rmaerker): Implement me!
    return fetchContextDataBackward(me);
}

template <typename TReference, typename TSpec, typename TConfig, typename TSpec, typename TOffset, typename TSize>
inline bool fetchUpstream(TraverserContext<JournaledStringTree<TReference, TSpec, TConfig>, TSpec> & me,
                          TOffset offset,
                          TSize blockSize)
{
    // TODO(rmaerker): Implement me!
    // Fetches block upstream to current range end position by offset.
}

// TODO(rmaerker): Disable for default JST
template <typename TReference, typename TSpec, typename TConfig, typename TSpec, typename TOffset, typename TSize>
inline bool fetchDownstream(TraverserContext<JournaledStringTree<TReference, TSpec, TConfig>, TSpec> & me,
                            TOffset offset,
                            TSize blockSize)
{
    // TODO(rmaerker): Implement me!
    // Fetches block downstream to current range begin position by offset.
}

template <typename TReference, typename TSpec, typename TConfig, typename TSpec, typename TPosition>
inline TPosition toGlobalPos(TraverserContext<JournaledStringTree<TReference, TSpec, TConfig>, TSpec> const & me,
                             TPosition localPos)
{
    // TODO(rrahn): Implement me!
}

}

#endif  // #ifndef INCLUDE_SEQAN_JOURNALED_STRING_TREE_TRAVERSER_CONTEXT_H_
