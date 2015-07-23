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

template <typename TJst, typename TSpec, typename TObserverList>
class TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserverList> : public Observable<TObserverList>
{
public:
    typedef Observable<TObserverList>                                   TSuper;
    typedef typename Member<TraverserImpl, TraverserStackMember>::Type  TStack;
    typedef typename Value<TStack>::Type                                TNode;
    typedef typename Size<TraverserImpl>::Type                          TSize;
    typedef std::shared_ptr<TStack>                                     TStackPtr;
    typedef JstBuffer_<TJst>                                            TBuffer;  // Provides a sequence context.
    typedef std::shared_ptr<TBuffer>                                    TBufferPtr;



    // __ Member Variables ____________________________________________________

    TJst *      _contPtr;  // TODO(rmaerker): Should always be const.
    TSize       _branchLength   = 1;
    TSize       _contextSize    = 1;
    TStackPtr   _stackPtr;
    TBufferPtr  _bufferPtr;

    // __ Constructors ________________________________________________________

    // Default c'tor.
    TraverserImpl() :
        TSuper(),
        _contPtr(nullptr),
        _stackPtr(impl::createStack<TStack>()),
        _bufferPtr(impl::createBuffer<TBuffer>())
    {}

    // C'tor with just the jst.
    TraverserImpl(TJst & jst) :
        TSuper(),
        _contPtr(nullptr),
        _stackPtr(impl::createStack<TStack>()),
        _bufferPtr(impl::createBuffer<TBuffer>())
    {
        init(*this, jst);
    }

    // C'tor with the jst and the context or branch size.
    template <typename TSize>
    TraverserImpl(TJst & jst, TSize const contextSize) :
    TSuper(),
    _contPtr(nullptr),
    _stackPtr(impl::createStack<TStack>()),
    _bufferPtr(impl::createBuffer<TBuffer>())
    {
        init(*this, jst, contextSize);
    }

    // C'tor with the jst and the context or branch size.
    template <typename TSize>
    TraverserImpl(TJst & jst, TSize const contextSize, TSize const branchLength) :
        TSuper(),
        _contPtr(nullptr),
        _stackPtr(impl::createStack<TStack>()),
        _bufferPtr(impl::createBuffer<TBuffer>())
    {
        init(*this, jst, contextSize, branchLength);
    }

    // C'tor with the jst and list of observers.
    TraverserImpl(TJst & jst, TObserverList && observers) :
        TSuper(std::forward<TObserverList>(observers)),
        _contPtr(nullptr),
        _stackPtr(impl::createStack<TStack>()),
        _bufferPtr(impl::createBuffer<TBuffer>())
    {
        init(*this, jst);
    }

    // C'tor with the jst and list of observers.
    template <typename TSize>
    TraverserImpl(TJst & jst, TSize const contextSize, TObserverList && observers) :
        TSuper(std::forward<TObserverList>(observers)),
        _contPtr(nullptr),
        _stackPtr(impl::createStack<TStack>()),
        _bufferPtr(impl::createBuffer<TBuffer>())
    {
        init(*this, jst, contextSize);
    }

    // C'tor with the jst and list of observers.
    template <typename TSize>
    TraverserImpl(TJst & jst, TSize const contextSize, TSize const branchLength, TObserverList && observers) :
        TSuper(std::forward<TObserverList>(observers)),
        _contPtr(nullptr),
        _stackPtr(impl::createStack<TStack>()),
        _bufferPtr(impl::createBuffer<TBuffer>())
    {
        init(*this, jst, contextSize, branchLength);
    }

    // Copy c'tor.
    template <typename TOtherJst>
    TraverserImpl(TraverserImpl<TOtherJst, JstTraversalSpec<TSpec>, TObserverList> const & other,
                  SEQAN_CTOR_ENABLE_IF(IsConstructible<TJst, TOtherJst>)) :
        _contPtr(other._contPtr),
        _branchLength(other._historySize),
        _contextSize(other._contextSize),
        _stackPtr(other._stackPtr),
        _bufferPtr(other._bufferPtr)
    {
        ignoreUnusedVariableWarning(dummy);
    }

    // __ Member Functions ____________________________________________________

    template <typename TOtherJst>
    inline SEQAN_FUNC_ENABLE_IF(IsConstructible<TJst, TOtherJst>, TraverserImpl &)
    operator=(TraverserImpl<TOtherJst, JstTraversalSpec<TSpec>, TObserverList> const & other)
    {
        if (*this != &other)
        {
            _contPtr = other._contPtr;
            _branchLength =other._historySize;
            _contextSize = other._contextSize;
            _stackPtr = other._stackPtr;
            _bufferPtr = other._bufferPtr;
        }
        return *this;
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Container
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec, typename TObserver>
struct Container<TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> >
{
    typedef TJst Type;
};

template <typename TJst, typename TSpec, typename TObserver>
struct Container<TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> const>
    : Container<TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> >
{};

// ----------------------------------------------------------------------------
// Metafunction Size
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec, typename TObserver>
struct Size<TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> >
{
    typedef typename Size<TJst>::Type    Type;
};

// ----------------------------------------------------------------------------
// Metafunction TraverserImpl
// ----------------------------------------------------------------------------

template <typename TSeq, typename TConfig, typename TSpec,
          typename TObserver>
struct Traverser<JournaledStringTree<TSeq, TConfig, TSpec>,
                 TObserver>
{
    typedef JournaledStringTree<TSeq, TConfig, TSpec> TJst_;
    typedef TraverserImpl<TJst_, JstTraversalSpec<void>, TObserver> Type;
};

template <typename TSeq, typename TConfig, typename TSpec,
          typename TObserver>
struct Traverser<JournaledStringTree<TSeq, TConfig, TSpec> const,
                 TObserver>
{
    typedef JournaledStringTree<TSeq, TConfig, TSpec> const TJst_;
    typedef TraverserImpl<TJst_, JstTraversalSpec<void>, TObserver> Type;
};

// ----------------------------------------------------------------------------
// Metafunction ContextIterator
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec, typename TObserver>
struct ContextIterator<TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> >
{
    typedef TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver>     TThis_;
    typedef typename Container<TThis_>::Type                            TContainer_;
    typedef typename Member<TContainer_, JstSourceMember>::Type         TSource_;

    typedef typename Iterator<TSource_, Standard>::Type                 Type;
};

template <typename TJst, typename TSpec, typename TObserver>
struct ContextIterator<TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> const>
{
    typedef TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> TTraverser_;

    typedef typename ContextIterator<TTraverser_>::Type  const      Type;
};

// ----------------------------------------------------------------------------
// Metafunction Member<TraverserStackMember>
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec, typename TObserver>
struct Member<TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver>, TraverserStackMember>
{
    typedef JstTraversalNode<TJst>              TTraversalNode;
    typedef String<TTraversalNode, Block<> >    Type;
};

// ============================================================================
// Public Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function contextIterator();
// ----------------------------------------------------------------------------

// Returns current sequence iterator.
template <typename TJst, typename TSpec, typename TObserver>
inline typename ContextIterator<TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> >::Type
contextIterator(TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> const & me)
{
    return impl::getContextIterator(me);
}

// ----------------------------------------------------------------------------
// Function contextSize();
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec, typename TObserver>
inline typename Size<TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> >::Type
contextSize(TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> const & me)
{
    return me._contextSize;
}

// ----------------------------------------------------------------------------
// Function container();
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec, typename TObserver>
inline typename Container<TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> >::Type &
container(TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> & me)
{
    return *me._contPtr;
}

template <typename TJst, typename TSpec, typename TObserver>
inline typename Container<TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> const>::Type &
container(TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> const & me)
{
    return *me._contPtr;
}

//// ----------------------------------------------------------------------------
//// Function setContextSize();
//// ----------------------------------------------------------------------------
//
//template <typename TJst, typename TSpec, typename TObserver, typename TSize>
//inline void
//setContextSize(TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> & me,
//               TSize contextSize)
//{
//    me._contextSize = contextSize;
//}
//
//// ----------------------------------------------------------------------------
//// Function setBranchLength();
//// ----------------------------------------------------------------------------
//
//template <typename TJst, typename TSpec, typename TObserver, typename TSize>
//inline void
//setBranchLength(TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> & me,
//                TSize branchLength)
//{
//    me._branchLength = branchLength;
//}

// ----------------------------------------------------------------------------
// Function goNext();
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec, typename TObserver,
          typename TSize>
inline void
goNext(TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> & me,
       TSize stepSize)
{
    SEQAN_ASSERT(me._stackPtr != nullptr);

    auto nodePtr = &back(*me._stackPtr);
#if defined (DEBUG_JST_TRAVERSAL)
    std::cout << "\n\nGO NEXT" << std::endl;
    std::cout << "Node: (" << length(*me._stackPtr) << ") " << *nodePtr << std::endl;
#endif  // DEBUG_JST_TRAVERSAL
    stepSize = impl::moveWindow(me, nodePtr, stepSize);

    if ((stepSize > 0 && nodePtr->remainingSize == 0) || atEnd(nodePtr->curDelta))
    {
        if (!back(*me._stackPtr).fromBase)
        {
            impl::popNode(me);
            SEQAN_ASSERT(length(*me._stackPtr) > 1u);
        }
        else
        {
            eraseBack(*me._stackPtr);  // Remove the old base representing fromBase

            SEQAN_ASSERT(length(*me._stackPtr) == 1u);
            impl::moveWindow(me, &back(*me._stackPtr), stepSize + back(*me._stackPtr).remainingSize);  // Move the base to the next position coming from fromBase node.
        }
    }
    else if (atEnd(nodePtr->curDelta))
    {
        if (!back(*me._stackPtr).fromBase)
        {
            impl::popNode(me);
            SEQAN_ASSERT(length(*me._stackPtr) > 1u);
        }
        else
        {
            eraseBack(*me._stackPtr);  // Remove the old base representing fromBase
            impl::moveWindow(me, &back(*me._stackPtr), stepSize + back(*me._stackPtr).remainingSize);  // Move the base to the next position coming from fromBase node.
        }
    }
}

// Special overload to move by one.
template <typename TJst, typename TSpec, typename TObserver>
inline void
goNext(TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> & it)
{
    goNext(it, 1);
}

// ----------------------------------------------------------------------------
// Function atEnd();
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec, typename TObserver>
inline bool
atEnd(TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> & me)
{
    SEQAN_ASSERT(me._stackPtr != nullptr);
    return length(*me._stackPtr) == 1 && back(*me._stackPtr).endEdgeIt == sourceEnd(impl::buffer(me));
}

// ----------------------------------------------------------------------------
// Function init();
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec, typename TObserver,
          typename TSize>
inline void
init(TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> & me,
     TJst & jst,
     TSize const contextSize,
     TSize const branchLength)
{
    impl::init(me, jst, contextSize, branchLength);
}

template <typename TJst, typename TSpec, typename TObserver, typename TSize>
inline void
init(TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> & me,
     TJst & jst,
     TSize const contextSize)
{
    init(me, jst, contextSize, contextSize);
}

template <typename TJst, typename TSpec, typename TObserver>
inline void
init(TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> & me,
     TJst & jst)
{
    init(me, jst, 1);
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_JOURNALED_STRING_TREE_JOURNALED_STRING_TREE_TRAVERSER_H_
