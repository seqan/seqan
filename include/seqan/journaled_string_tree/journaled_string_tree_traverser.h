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

#ifndef INCLUDE_SEQAN_JOURNALED_STRING_TREE_JOURNALED_STRING_TREE_TRAVERSER_H_
#define INCLUDE_SEQAN_JOURNALED_STRING_TREE_JOURNALED_STRING_TREE_TRAVERSER_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TJst, typename TSpec>
class TraverserImpl<TJst, JstTraversalSpec<TSpec> >
{
public:
    typedef typename Member<TraverserImpl, TraverserStackMember>::Type  TStack;
    typedef typename Value<TStack>::Type                                TNode;
    typedef typename Size<TraverserImpl>::Type                          TSize;
    typedef std::shared_ptr<TStack>                                     TStackPtr;
    typedef JstBuffer_<TJst>                                            TBuffer;  // Provides a sequence context.
    typedef std::shared_ptr<TBuffer>                                    TBufferPtr;
    typedef typename Member<TJst, JstDeltaMapMember>::Type              TDeltaMap;
    typedef typename DeltaCoverage<TDeltaMap>::Type                     TCoverage;

    // __ Member Variables ____________________________________________________

    TJst *      _contPtr;  // TODO(rmaerker): Underlying container should be const.
    TSize       _contextSize    = 1;
    TSize       _branchLength   = 1;
    TStackPtr   _stackPtr;
    TBufferPtr  _bufferPtr;
    TCoverage   _baseCov;
    bool        _needInitialization = true;

    // __ Constructors ________________________________________________________

    // Default c'tor.
    TraverserImpl() :
        _contPtr(nullptr),
        _stackPtr(impl::createStack<TStack>()),
        _bufferPtr(impl::createBuffer<TBuffer>())
    {}

    // C'tor with just the jst.
    TraverserImpl(TJst & jst) :
        _contPtr(&jst),
        _stackPtr(impl::createStack<TStack>()),
        _bufferPtr(impl::createBuffer<TBuffer>())
    {}

    // C'tor with the jst and the context or branch size.
    template <typename TSize>
    TraverserImpl(TJst & jst, TSize const contextSize) :
        _contPtr(&jst),
        _contextSize(contextSize),
        _branchLength(contextSize),
        _stackPtr(impl::createStack<TStack>()),
        _bufferPtr(impl::createBuffer<TBuffer>())
    {}

    // C'tor with the jst and the context or branch size.
    template <typename TSize>
    TraverserImpl(TJst & jst, TSize const contextSize, TSize const branchLength) :
        _contPtr(&jst),
        _contextSize(contextSize),
        _branchLength(branchLength),
        _stackPtr(impl::createStack<TStack>()),
        _bufferPtr(impl::createBuffer<TBuffer>())
    {}

    // Copy c'tor.
    template <typename TOtherJst>
    TraverserImpl(TraverserImpl<TOtherJst, JstTraversalSpec<TSpec> > const & other,
                  SEQAN_CTOR_ENABLE_IF(IsConstructible<TJst, TOtherJst>)) :
        _contPtr(other._contPtr),
        _contextSize(other._contextSize),
        _branchLength(other._branchLength),
        _stackPtr(other._stackPtr),
        _bufferPtr(other._bufferPtr),
        _baseCov(other._baseCov),
        _needInitialization(other._needInitialization)

    {
        ignoreUnusedVariableWarning(dummy);
    }

    // Move c'tor.
    template <typename TOtherJst>
    TraverserImpl(TraverserImpl<TOtherJst, JstTraversalSpec<TSpec> > && other,
                  SEQAN_CTOR_ENABLE_IF(IsConstructible<TJst, TOtherJst>)) :
        _contPtr(std::move(other._contPtr)),
        _contextSize(std::move(other._contextSize)),
        _branchLength(std::move(other._branchLength)),
        _stackPtr(std::move(other._stackPtr)),
        _bufferPtr(std::move(other._bufferPtr)),
        _baseCov(std::move(other._baseCov)),
        _needInitialization(std::move(other._needInitialization))
    {
        ignoreUnusedVariableWarning(dummy);
    }

    // __ Member Functions ____________________________________________________

    template <typename TOtherJst>
    inline SEQAN_FUNC_ENABLE_IF(IsConstructible<TJst, TOtherJst>, TraverserImpl &)
    operator=(TraverserImpl<TOtherJst, JstTraversalSpec<TSpec> > other)
    {
        impl::swap(*this, other);
        return *this;
    }

    // Destructor
    ~TraverserImpl() = default;
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Container
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec>
struct Container<TraverserImpl<TJst, JstTraversalSpec<TSpec> > >
{
    typedef TJst Type;
};

template <typename TJst, typename TSpec>
struct Container<TraverserImpl<TJst, JstTraversalSpec<TSpec> > const>
    : Container<TraverserImpl<TJst, JstTraversalSpec<TSpec> > >
{};

// ----------------------------------------------------------------------------
// Metafunction Position
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec>
struct Position<TraverserImpl<TJst, JstTraversalSpec<TSpec> > >
{
    using TPos_ = typename Position<TJst>::Type;
    using TSize_ = typename Size<TJst>::Type;

    using Type = String<Pair<TSize_, TPos_> >;
};

// ----------------------------------------------------------------------------
// Metafunction Size
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec>
struct Size<TraverserImpl<TJst, JstTraversalSpec<TSpec> > >
{
    typedef typename Size<TJst>::Type    Type;
};

// ----------------------------------------------------------------------------
// Metafunction TraverserImpl
// ----------------------------------------------------------------------------

template <typename TSeq, typename TConfig, typename TSpec>
struct Traverser<JournaledStringTree<TSeq, TConfig, TSpec> >
{
    typedef JournaledStringTree<TSeq, TConfig, TSpec> TJst_;
    typedef TraverserImpl<TJst_, JstTraversalSpec<void> > Type;
};

template <typename TSeq, typename TConfig, typename TSpec>
struct Traverser<JournaledStringTree<TSeq, TConfig, TSpec> const>
{
    typedef JournaledStringTree<TSeq, TConfig, TSpec> const TJst_;
    typedef TraverserImpl<TJst_, JstTraversalSpec<void> > Type;
};

// ----------------------------------------------------------------------------
// Metafunction ContextIterator
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec>
struct ContextIterator<TraverserImpl<TJst, JstTraversalSpec<TSpec> > >
{
    typedef TraverserImpl<TJst, JstTraversalSpec<TSpec> >               TThis_;
    typedef typename Container<TThis_>::Type                            TContainer_;
    typedef typename Member<TContainer_, JstSourceMember>::Type         TSource_;

    typedef typename Iterator<TSource_, Standard>::Type                 Type;
};

template <typename TJst, typename TSpec>
struct ContextIterator<TraverserImpl<TJst, JstTraversalSpec<TSpec> > const>
{
    typedef TraverserImpl<TJst, JstTraversalSpec<TSpec> > TTraverser_;

    typedef typename ContextIterator<TTraverser_>::Type  const      Type;
};

// ----------------------------------------------------------------------------
// Metafunction TraverserNode
// ----------------------------------------------------------------------------

template <typename TTraverser>
struct TraverserNode;

template <typename TJst, typename TSpec>
struct TraverserNode<TraverserImpl<TJst, JstTraversalSpec<TSpec> > >
{
    using Type = JstTraversalNode<TJst>;
};

template <typename TJst, typename TSpec>
struct TraverserNode<TraverserImpl<TJst, JstTraversalSpec<TSpec> > const>
{
    using Type = JstTraversalNode<TJst> const;
};

// ----------------------------------------------------------------------------
// Metafunction Member<TraverserStackMember>
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec>
struct Member<TraverserImpl<TJst, JstTraversalSpec<TSpec> >, TraverserStackMember>
{
    using TNode = typename TraverserNode<TraverserImpl<TJst, JstTraversalSpec<TSpec> > >::Type;
    using Type  = String<TNode, Block<> >;
};

// ============================================================================
// Public Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function contextIterator();
// ----------------------------------------------------------------------------

// Returns current sequence iterator.
template <typename TJst, typename TSpec>
inline typename ContextIterator<TraverserImpl<TJst, JstTraversalSpec<TSpec> > >::Type
contextIterator(TraverserImpl<TJst, JstTraversalSpec<TSpec> > const & me)
{
    return impl::getContextIterator(me);
}

// ----------------------------------------------------------------------------
// Function contextBegin();
// ----------------------------------------------------------------------------

// Returns current sequence iterator.
template <typename TJst, typename TSpec>
inline typename ContextIterator<TraverserImpl<TJst, JstTraversalSpec<TSpec> > >::Type
contextBegin(TraverserImpl<TJst, JstTraversalSpec<TSpec> > const & me)
{
    return impl::getContextBegin(me);
}

// Returns current sequence iterator.
template <typename TJst, typename TSpec>
inline typename ContextIterator<TraverserImpl<TJst, JstTraversalSpec<TSpec> > >::Type
contextEnd(TraverserImpl<TJst, JstTraversalSpec<TSpec> > const & me)
{
    return impl::getContextEnd(me);
}

// ----------------------------------------------------------------------------
// Function container();
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec>
inline typename Container<TraverserImpl<TJst, JstTraversalSpec<TSpec> > >::Type &
container(TraverserImpl<TJst, JstTraversalSpec<TSpec> > & me)
{
    return *me._contPtr;
}

template <typename TJst, typename TSpec>
inline typename Container<TraverserImpl<TJst, JstTraversalSpec<TSpec> > const>::Type &
container(TraverserImpl<TJst, JstTraversalSpec<TSpec> > const & me)
{
    return *me._contPtr;
}

// ----------------------------------------------------------------------------
// Function setContainer();
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec>
inline void
setContainer(TraverserImpl<TJst, JstTraversalSpec<TSpec> > & me,
             TJst & jst)
{
    me._needInitialization = true;
    me._contPtr = &jst;
}

// ----------------------------------------------------------------------------
// Function contextSize();
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec>
inline auto
contextSize(TraverserImpl<TJst, JstTraversalSpec<TSpec> > const & me) -> decltype(me._contextSize)
{
    return me._contextSize;
}

// ----------------------------------------------------------------------------
// Function setContextSize();
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec, typename TSize>
inline void
setContextSize(TraverserImpl<TJst, JstTraversalSpec<TSpec> > & me,
               TSize contextSize)
{
    me._needInitialization = true;
    me._contextSize = contextSize;
}

// ----------------------------------------------------------------------------
// Function branchSize();
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec>
inline auto
branchSize(TraverserImpl<TJst, JstTraversalSpec<TSpec> > const & me) -> decltype(me._branchLength)
{
    return me._branchLength;
}

// ----------------------------------------------------------------------------
// Function setBranchSize();
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec, typename TSize>
inline void
setBranchSize(TraverserImpl<TJst, JstTraversalSpec<TSpec> > & me,
              TSize newSize)
{
    me._needInitialization = true;
    me._branchLength = newSize;
}

// ----------------------------------------------------------------------------
// Function position();
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec>
inline typename Position<TraverserImpl<TJst, JstTraversalSpec<TSpec> > >::Type
position(TraverserImpl<TJst, JstTraversalSpec<TSpec> > const & me)
{

    if (length(impl::stack(me)) == 1)
        return impl::positionReference(me);
    else
        return impl::positionBranch(me);
}

// ----------------------------------------------------------------------------
// Function init();
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec,
          typename TObserver,
          typename TProxySelector>
inline void
init(TraverserImpl<TJst, JstTraversalSpec<TSpec> > & me,
     TObserver & observer,
     Tag<TProxySelector> const & /*tag*/)
{
    impl::init(me, observer, Tag<TProxySelector>());
    me._needInitialization = false;
}

template <typename TJst, typename TSpec,
          typename TObserver>
inline void
init(TraverserImpl<TJst, JstTraversalSpec<TSpec> > & me,
     TObserver & observer)
{
    init(me, observer, SelectValidProxy());
}

// ----------------------------------------------------------------------------
// Function goNext();
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec,
          typename TSize,
          typename TObserver,
          typename TProxySelector>
inline void
advance(TraverserImpl<TJst, JstTraversalSpec<TSpec> > & me,
        TSize stepSize,
        TObserver & observer,
        Tag<TProxySelector> const & /*tag*/)
{
    SEQAN_ASSERT(me._stackPtr != nullptr);

    if (SEQAN_UNLIKELY(me._needInitialization))
        init(me, observer);

    auto nodePtr = &back(*me._stackPtr);

#if defined (JST_FIND_DEBUG)
    if (nodePtr->isBase && getDeltaPosition(*nodePtr->nextDelta) == 40)
        std::cout << "     STOP! " << std::endl;

    std::cout << *nodePtr << std::endl;
#endif // defined (JST_FIND_DEBUG)
#if defined (DEBUG_JST_TRAVERSAL)
    std::cout << "\n\nGO NEXT" << std::endl;
    std::cout << "Node: (" << length(*me._stackPtr) << ") " << *nodePtr << std::endl;
#endif  // DEBUG_JST_TRAVERSAL
    stepSize = impl::moveWindow(me, nodePtr, stepSize, observer, Tag<TProxySelector>());

    if ((stepSize > 0 && nodePtr->remainingSize == 0))
    {
        SEQAN_ASSERT_NOT(impl::activeNode(me).isBase);
        impl::popNode(me, observer);
        SEQAN_ASSERT(length(impl::stack(me)) >= 1u);
    }

    // Check if the condition changed.
    nodePtr = &impl::activeNode(me);
    if (SEQAN_UNLIKELY(atEnd(nodePtr->curDelta))) // Current branch is at final end.
    {
        if (nodePtr->isBase)
            return;

        while (atEnd(nodePtr->curDelta) && !(nodePtr->isBase))
        {
            SEQAN_ASSERT(nodePtr->curEdgeIt == nodePtr->endEdgeIt); // edgeIt must be at end of current branch.
            impl::popNode(me, observer);  // Might not be the parent one.
            nodePtr = &impl::activeNode(me);
        }
    }
}

template <typename TJst, typename TSpec,
          typename TSize,
          typename TProxySelector>
inline void
advance(TraverserImpl<TJst, JstTraversalSpec<TSpec> > & me,
        TSize const stepSize,
        Tag<TProxySelector> const & /*tag*/)
{
    auto observer = makeObserverList();
    advance(me, stepSize, observer);
}

template <typename TJst, typename TSpec,
          typename TSize,
          typename TObserver>
inline void
advance(TraverserImpl<TJst, JstTraversalSpec<TSpec> > & me,
        TSize const stepSize,
        TObserver & observer)
{
    advance(me, stepSize, observer, SelectValidProxy());
}

template <typename TJst, typename TSpec,
          typename TSize>
inline void
advance(TraverserImpl<TJst, JstTraversalSpec<TSpec> > & me,
        TSize const stepSize)
{
    auto observer = makeObserverList();
    advance(me, stepSize, observer, SelectValidProxy());
}

// ----------------------------------------------------------------------------
// Function atEnd();
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec>
inline bool
atEnd(TraverserImpl<TJst, JstTraversalSpec<TSpec> > & me)
{
    SEQAN_ASSERT(me._stackPtr != nullptr);
    return length(*me._stackPtr) == 1 && back(*me._stackPtr).curEdgeIt == sourceEnd(impl::buffer(me));
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_JOURNALED_STRING_TREE_JOURNALED_STRING_TREE_TRAVERSER_H_
