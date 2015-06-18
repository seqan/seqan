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

template <typename TJst, typename TSpec, typename TObserver>
class TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> : public Observable<TObserver>
{
public:
    typedef Observable<TObserver>                                       TSuper;
    typedef typename Member<TraverserImpl, TraverserStackMember>::Type  TStack;
    typedef typename Value<TStack>::Type                                TNode;
    typedef typename Size<TraverserImpl>::Type                          TSize;
    typedef std::shared_ptr<TStack>                                     TStackPtr;


    // __ Member Variables ____________________________________________________

    TJst *       _contPtr;
    TSize        _historySize;
    TSize        _contextSize;
    TNode        _tmp;      // Better use local vairable.
    TStackPtr    _stackPtr;

    // __ Constructors ________________________________________________________

    TraverserImpl() :
        TSuper(),
        _contPtr(nullptr),
        _historySize(0),
        _contextSize(1),
        _tmp(),
        _stackPtr(impl::createStack<TStack>())
    {}

    TraverserImpl(TJst & jst) :
        TSuper(),
        _contPtr(nullptr),
        _historySize(historySize(jst)),
        _contextSize(1),
        _tmp(),
        _stackPtr(impl::createStack<TStack>())
    {
        impl::init(*this, jst);
    }

    template <typename TObserver_>
    TraverserImpl(TJst & jst, TObserver_ & observer, SEQAN_CTOR_DISABLE_IF(IsSameType<TObserver_, void>)) :
        TSuper(),
        _contPtr(nullptr),
        _historySize(historySize(jst)),
        _contextSize(1),
        _stackPtr(impl::createStack<TStack>())
    {
        addObserver(*this, observer);
        impl::init(*this, jst);
        ignoreUnusedVariableWarning(dummy);
    }

    // Copy Constructor!
    template <typename TOtherJst>
    TraverserImpl(TraverserImpl<TOtherJst, JstTraversalSpec<TSpec>, TObserver> const & other,
                  SEQAN_CTOR_ENABLE_IF(IsConstructible<TJst, TOtherJst>)) :
        _contPtr(other._contPtr),
        _historySize(other._historySize),
        _contextSize(other._contextSize),
        _tmp(other._tmp),
        _stackPtr(other._stackPtr)
    {
        ignoreUnusedVariableWarning(dummy);
    }

    // __ Member Functions ____________________________________________________

    template <typename TOtherJst>
    inline SEQAN_FUNC_ENABLE_IF(IsConstructible<TJst, TOtherJst>, TraverserImpl &)
    operator=(TraverserImpl<TOtherJst, JstTraversalSpec<TSpec>, TObserver> const & other)
    {
        if (*this != &other)
        {
            _contPtr = other._contPtr;
            _historySize =other._historySize;
            _contextSize = other._contextSize;
            _tmp = other._tmp;
            _stackPtr = other._stackPtr;
        }
        return *this;
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

template <typename TJst, typename TSpec, typename TObserver>
struct Container<TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> >
{
    typedef TJst Type;
};

template <typename TJst, typename TSpec, typename TObserver>
struct Size<TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> >
{
    typedef typename Size<TJst>::Type    Type;
};

template <typename TJst, typename TSpec, typename TObserver>
struct StringContext<TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> >
{
    typedef TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver>     TThis_;
    typedef typename Container<TThis_>::Type                            TContainer_;
    typedef typename Source<TContainer_>::Type                          TSource_;
    typedef typename Iterator<TSource_, Standard>::Type                 TSourceIt_;

    typedef Range<TSourceIt_>                                           Type;
};

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
// Function stringContext()
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec, typename TObserver>
inline typename StringContext<TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> >::Type
stringContext(TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> const & traversor)
{
    return impl::getStringContext(traversor);
}

// ----------------------------------------------------------------------------
// Function container
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec, typename TObserver>
inline typename Container<TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> >::Type &
container(TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> & me)
{
    return *me._contPtr;
}

template <typename TJst, typename TSpec, typename TObserver>
inline typename Container<TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> const>::Type &
container(TraverserImpl<TJst, JstTraversalSpec<TSpec> > const & me)
{
    return *me._contPtr;
}

// ----------------------------------------------------------------------------
// Function atBegin()
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec, typename TObserver>
inline bool
atBegin(TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> & me)
{
    SEQAN_ASSERT(resource(me._rm) == nullptr);
    return back(*resource(me._rm)).mappedSrcEndPos == -1;
}


template <typename TJst, typename TSpec, typename TObserver,
          typename TSize>
inline void
goNext(TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> & me,
       TSize stepSize)
{
    SEQAN_ASSERT(me._stackPtr != nullptr);

    auto& node = back(*me._stackPtr);

#if defined (DEBUG_JST_TRAVERSAL)
    std::cout << "goNext: Node ID = " << length(*me._stackPtr) << std::endl;
    std::cout << node << std::endl;
#endif  // DEBUG_JST_TRAVERSAL
    stepSize = impl::moveWindow(me, node, stepSize);

    // Case A) stepSize == 0 && node.remainingSize > 0;  // vaild context
    // Case B) stepSize == 0 && node.remainingSize = 0;  // valid context
    // Case C) stepSize > 0 && node.remaingingSize = 0;  // invalid context -> we have to pop an element from the tree.
    // Case D) stepSize > 0 && node.remainingSize > 0;   // invalid case: ASSERT

    SEQAN_ASSERT(stepSize > static_cast<TSize>(0) && node.remainingSize >static_cast<TSize>(0));

    if (stepSize > 0 && node.remainingSize == 0)
    {
        impl::popNode(me);
    }
}

// ----------------------------------------------------------------------------
// Function goNext()
// ----------------------------------------------------------------------------

// Standard version just iterating over the delta nodes.
template <typename TJst, typename TSpec, typename TObserver>
inline void
goNext(TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> & it)
{
    goNext(it, 1);
}

// ----------------------------------------------------------------------------
// Function atEnd()
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec, typename TObserver>
inline bool
atEnd(TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> & me)
{
    SEQAN_ASSERT(me._stackPtr != nullptr);
    return length(*me._stackPtr) == 1 && back(*me._stackPtr).endEdgeIt == sourceEnd(container(me)._buffer);
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_JOURNALED_STRING_TREE_JOURNALED_STRING_TREE_TRAVERSER_H_
