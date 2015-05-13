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
// Basic defintions and forwards used globally for this module.
// ==========================================================================

#ifndef EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TRAVERSER_H_
#define EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TRAVERSER_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// TODO(rrahn): Move to jst_traverser_base.h
// ----------------------------------------------------------------------------
// Metafunction Buffer
// ----------------------------------------------------------------------------

template <typename TObject, typename TConfig, typename TSpec>
struct Buffer<Traverser<TObject, TConfig, TSpec> >;

template <typename TObject, typename TConfig, typename TSpec>
struct Buffer<Traverser<TObject, TConfig, TSpec> const>
{
    typedef typename Buffer<Traverser<TObject, TConfig, TSpec>::Type const Type;
};

// ----------------------------------------------------------------------------
// Metafunction Traits
// ----------------------------------------------------------------------------

// TODO(rrahn): Move to basic meta-function.
// TODO(rrahn): Add Documentation

template <typename TObject>
struct Traits;

// Default const implementation.
template <typename TObject>
struct Traits<TObject const> : Traits<TObject>
{};

// TODO(rrahn): Introduce concepts: BidirectionalTraversableConcept = ForwardTraversableConcept & BackwardTraversableConcept

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Tag JstTraverser
// ----------------------------------------------------------------------------

struct JstTraverser_;
typedef Tag<JstTraverser_> JstTraverser;

// ----------------------------------------------------------------------------
// Config DefaultJstTraverserConfig
// ----------------------------------------------------------------------------

template <typename TObject>
struct DefaultJstTraverserConfig
{
    typedef TraverseForward TDirection;
    typedef StaticJstBuffer TBuffer;
};

// ----------------------------------------------------------------------------
// Class JstTraverser
// ----------------------------------------------------------------------------

/*!
 * @class JstTraverser
 * @headerfile <seqan/journaled_string_tree.h>
 */

template <typename TObject, typename TConfig, typename TSpec>
class Traverser;

template <typename TObject, typename TConfig = DefaultJstTraverserConfig<TObject> >
class Traverser<TObject, TConfig, JstTraverser>
{
public:

    typedef typename Container<Traverser>::Type                                 TDeltaMap;
    typedef typename Buffer<Traverser>::Type                                    TJstBuffer;
    typedef typename Member<TJstBuffer, JstSeqBufferJournaledSetMember>::Type   TJournaledSet;
    typedef typename Value<TJournaledSet>::Type                                 TJournaledSeq;
    typedef typename Size<TDeltaMap>::Type                                      TSize;
    typedef BaseJstTraversalEntry<TDeltaMap, TJournaledSeq>                     TBaseEntry;


    bool                inBranch;
    TBaseEntry *        entryPtr;                   // Pointer to the traverser stack.
    TSize               windowSize;                 // Size of the window.
    TJstBuffer          bufferOwner;                // Owns the sequence buffer.
    TJstBuffer *        bufferPtr;                  // Points to used sequence buffer.

    Traverser() : inBranch(false), entryPtr(nullptr), contextSize(1), bufferOwner(), bufferPtr(&bufferOwner)
    {}

    Traverser(TContainer & jst) :
        inBranch(false),
        entryPtr(nullptr),
        contextSize(1),
        bufferOwner(jst),
        bufferPtr(&bufferOwner)
    {}

    template <typename TPos>
    Traverser(TContainer & jst, TPos baseBegin, TPos baseEnd) :
        inBranch(false),
        entryPtr(nullptr),
        contextSize(1),
        bufferOwner(jst, baseBegin, baseEnd),
        bufferPtr(&bufferOwner)
    {}

    Traverser(TJstBuffer & buffer) :
        inBranch(false),
        entryPtr(nullptr),
        contextSize(1),
        bufferPtr(&buffer)
    {}

};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Container
// ----------------------------------------------------------------------------

template <typename TObject, typename TConfig, typename TSpec>
struct Container<Traverser<TObject, TConfig, TSpec> >
{
    typedef TObject Type;
};

template <typename TObject, typename TConfig, typename TSpec>
struct Container<Traverser<TObject, TConfig, TSpec> const >
{
    typedef TObject const Type;
};

// ----------------------------------------------------------------------------
// Metafunction TraverserContext
// ----------------------------------------------------------------------------

template <typename TObject, typename TConfig>
struct Window<Traverser<TObject, TConfig, JstTraverser> >
{
    typedef typename Buffer<Traverser<TObject, TConfig, JstTraverser> >::Type   TBuffer_;
    typedef typename Member<JstSequenceBuffer, JstSeqBufferJournaledSet>::Type  TJournaledSet_;
    typedef typename Value<TJournaledSet>::Type                                 TJournaledString_;
    typedef typename Iterator<TJournaledString_, Standard>::Type                TJournaledStrIterator_;
    typedef          Range<TJournaledStrIterator_>                              Type;
};

template <typename TObject, typename TConfig>
struct Window<Traverser<TObject, TConfig, JstTraverser> const >
{
    typedef typename TraverserContext<Traverser<TObject, TConfig, JstTraverser> >::Type const Type;
};

// ----------------------------------------------------------------------------
// Metafunction Buffer
// ----------------------------------------------------------------------------

template <typename TObject, typename TConfig>
struct Buffer<Traverser<TObject, TConfig, JstTraverser> >
{
    typedef typename TConfig::TDirection    TDirection_;
    typedef typename TConfig::TBuffer       TBuffer_;
    typedef JstSequenceBuffer<TObject, TDirection_, TBuffer_> Type;
};

// ----------------------------------------------------------------------------
// Metafunction  Traits
// ----------------------------------------------------------------------------

template <typename TObject, typename TConfig, typename TSpec>
struct Traits<Traverser<TObject, TConfig, TSpec> >
{
    typedef TConfig Type;
};

// ============================================================================
// Private Functions
// ============================================================================

namespace impl
{

// ----------------------------------------------------------------------------
// Function impl::minContextSpan()
// ----------------------------------------------------------------------------

template <typename TObject, typename TConfig, typename TSpec>
inline typename Size<typename Container<Traverser<TObject, TConfig, TSpec> >::Type>::Type
minContextSpan(Traverser<TObject, TConfig, TSpec> const & container)
{
    typedef Traverser<TObject, TConfig, TSpec>                          TTraverser;
    typedef typename Size<typename Container<TTraverser>::Type>::Type   TSize;

    // TODO(rrahn): Need to ensure that while dynamic buffering the second block always guarantees the following condition to be not true!
    if (SEQAN_UNLIKELY(traverser.statePtr->currIt - begin(container(traverser.statePtr->currIt), Standard()) <
                       contextSize(traverser)))
        return traverser.statePtr->currIt - begin(container(traverser.statePtr->currIt), Standard())
    return contextSize(traverser)
}

// ----------------------------------------------------------------------------
// Function impl::getContext()
// ----------------------------------------------------------------------------

template <typename TObject, typename TConfig, typename TSpec>
inline typename TraverserContext<Traverser<TObject, TConfig, TSpec> >::Type
getContext(Traverser<TObject, TConfig, TSpec> & traverser)
{
    typedef typename TraverserContext<Traverser<TObject, TConfig, TSpec> >::Type TContext;

    // TODO(rrahn): Add for bidirectional and backward traversal.
    return TContext(traverser.statePtr->currIt - minContextSpan(traverser), traverser.statePtr->currIt);
}

template <typename TObject, typename TConfig, typename TSpec>
inline typename TraverserContext<Traverser<TObject, TConfig, TSpec> const>::Type
getContext(Traverser<TObject, TConfig, TSpec> const & traverser)
{
    typedef typename TraverserContext<Traverser<TObject, TConfig, TSpec> >::Type TContext;

    // TODO(rrahn): Add for bidirectional and backward traversal.
    return TContext(traverser.statePtr->currIt - minContextSpan(traverser), traverser.statePtr->currIt);
}

}  // namespace impl

// ============================================================================
// Public Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function setContainer()
// ----------------------------------------------------------------------------

// JstTraverser option to set the container for the buffer as well.
template <typename TObject, typename TConfig>
inline void
setContainer(Traverser<TObject, TConfig, JstTraverser> & traverser,
             typename Container<Traverser<TObject, TConfig, JstTraverser> >::Type & container)
{
    // TODO(rrahn): Clear the traverser's states.
    setContainer(jstBuffer(traverser), container);
}

// ----------------------------------------------------------------------------
// Function container()
// ----------------------------------------------------------------------------

// RD_ONLY to avoid replacing the container without noticing the internal structures.
template <typename TObject, typename TConfig>
inline typename Container<Traverser<TObject, TConfig, JstTraverser> const>::Type &
container(Traverser<TObject, TConfig, JstTraverser> const & traverser)
{
    return container(jstBuffer(traverser));
}

// ----------------------------------------------------------------------------
// Function setBuffer()
// ----------------------------------------------------------------------------

template <typename TObject, typename TConfig>
inline typename Buffer<Traverser<TObject, TConfig, JstTraverser> >::Type &
setBuffer(Traverser<TObject, TConfig, JstTraverser> & traverser,
          typename Buffer<Traverser<TObject, TConfig, JstTraverser> >::Type & buffer)
{
    // TODO(rrahn): Clear the traverser's states.
    return traverser.bufferPtr = &buffer;
}

// ----------------------------------------------------------------------------
// Function buffer()
// ----------------------------------------------------------------------------

// TODO(rrahn): Document behaviour and possible problem, if multiple traversers depend on the same buffer, e.g. invalid state.
template <typename TObject, typename TConfig>
inline typename Buffer<Traverser<TObject, TConfig, JstTraverser> >::Type &
buffer(Traverser<TObject, TConfig, JstTraverser> & traverser)
{
    return *traverser.bufferPtr;
}

template <typename TObject, typename TConfig>
inline typename Buffer<Traverser<TObject, TConfig, JstTraverser> const>::Type &
buffer(Traverser<TObject, TConfig, JstTraverser> const & traverser)
{
    return *traverser.bufferPtr;
}

// ----------------------------------------------------------------------------
// Function setBeginPosition()
// ----------------------------------------------------------------------------

template <typename TObject, typename TConfig, typename TPos>
inline void
setBeginPosition(Traverser<TObject, TConfig, JstTraverser> & traverser,
                 TPos beginPos)
{
    // TODO(rrahn): Clear the traverser's states, if forward or bidirectional and beginPos right of previous one.
    setBeginPosition(buffer(traverser), beginPos);
}

// ----------------------------------------------------------------------------
// Function setEndPosition()
// ----------------------------------------------------------------------------

//  to avoid unnotified replacement of the jst buffer.
template <typename TObject, typename TConfig, typename TPos>
inline void
setEndPosition(Traverser<TObject, TConfig, JstTraverser> & traverser,
               TPos endPos)
{
    // TODO(rrahn): Clear the traverser's states -> if bidirectional or backward and endPosition left of previous one.
    setEndPosition(buffer(traverser), endPos);
}

// ----------------------------------------------------------------------------
// Function windowSize()
// ----------------------------------------------------------------------------

template <typename TObject, typename TConfig, typename TSpec>
inline typename Size<typename Container<Traverser<TObject, TConfig, TSpec> >::Type>::Type
windowSize(Traverser<TObject, TConfig, TSpec> const & traverser)
{
    return traverser.windowSize;
}

// ----------------------------------------------------------------------------
// Function setWindowSize()
// ----------------------------------------------------------------------------

template <typename TObject, typename TConfig, typename TSpec, typename TSize>
inline void
setWindowSize(Traverser<TObject, TConfig, TSpec> const & traverser, TSize newSize)
{
    traverser.windowSize = newSize;
}

// ----------------------------------------------------------------------------
// Function window()
// ----------------------------------------------------------------------------

template <typename TObject, typename TConfig, typename TSpec>
inline typename Window<Traverser<TObject, TConfig, TSpec> >::Type
window(Traverser<TObject, TConfig, TSpec> & traverser)
{
    return impl::getContext(traverser);
}

template <typename TObject, typename TConfig, typename TSpec>
inline typename Window<Traverser<TObject, TConfig, TSpec> const>::Type
window(Traverser<TObject, TConfig, TSpec> const & traverser)
{
    return impl::getContext(traverser);
}

// ----------------------------------------------------------------------------
// Function windowBegin()
// ----------------------------------------------------------------------------

template <typename TObject, typename TConfig, typename TSpec>
inline typename Iterator<typename Window<Traverser<TObject, TConfig, TSpec> >::Type>::Type
windowBegin(Traverser<TObject, TConfig, TSpec> & traverser)
{
    return impl::getContext(traverser).begin;
}

template <typename TObject, typename TConfig, typename TSpec>
inline typename Iterator<typename Window<Traverser<TObject, TConfig, TSpec> const>::Type>::Type
windowBegin(Traverser<TObject, TConfig, TSpec> const & traverser)
{
    return impl::getContext(traverser).begin;
}

// ----------------------------------------------------------------------------
// Function windowEnd()
// ----------------------------------------------------------------------------

template <typename TObject, typename TConfig, typename TSpec>
inline typename Iterator<typename Window<Traverser<TObject, TConfig, TSpec> >::Type>::Type
windowEnd(Traverser<TObject, TConfig, TSpec> & traverser)
{
    return impl::getContext(traverser).end;
}

template <typename TObject, typename TConfig, typename TSpec>
inline typename Iterator<typename Window<Traverser<TObject, TConfig, TSpec> const>::Type>::Type
windowEnd(Traverser<TObject, TConfig, TSpec> const & traverser)
{
    return impl::getContext(traverser).end;
}

// ----------------------------------------------------------------------------
// Function traverse()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TConfig, typename TExternal, typename TDirection>
inline void
traverse(Traverser<TContainer, TConfig, JstTraverser> & traverser,
         TExternal & external,
         TDirection const & /*dir*/)
{
    typedef Traverser<TContainer, TConfig, JstTraverser>    TTraverser;
    typedef JstTraversalOperator<TTraverser, TExternal>     TOperator;

    TOperator op(traverser, extension);  // Initialize the operator.

    SEQAN_ASSERT(traverser.entryPtr != nullptr);
    // Continue as long as there is work on the stack.
    while (!empty(op.stack))
    {
        // Loop over the current branch.
        while(impl::current(op).cur != impl::current(op).end)
        {
            if (impl::isBase(op))
                if (SEQAN_UNLIKELY(impl::current(op).cur >= chunk(buffer(traverser)).end))
                    advanceChunk(buffer(traverser));  // Buffer new chunk.

            // Check for branch position.
            if (position(impl::current(op).cur) >= position(impl::current(op).bpNextVirtual))
                expand(op);

            advance(op, TDirection());
        }
        pop(op.stack);
    }
}

}  // namespace seqan

#endif  // EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_JOURNALED_STRING_TREE_BASE_H_
