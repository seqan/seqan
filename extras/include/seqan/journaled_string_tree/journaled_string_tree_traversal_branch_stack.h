// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
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
// Implements the branch stack for journaled string tree traversal.
// ==========================================================================

#ifndef EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_JOURNALED_STRING_TREE_TRAVERSAL_BRANCH_STACK_H_
#define EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_JOURNALED_STRING_TREE_TRAVERSAL_BRANCH_STACK_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class JstBranchStackEntry_
// ----------------------------------------------------------------------------

template <typename TJst, typename TState>
class JstBranchStackEntry_
{
public:

    typedef typename GetStringSet<TJst>::Type TJournaledSet;
    typedef typename Value<TJournaledSet>::Type TJournaledString;
    typedef typename Iterator<TJournaledString, Standard>::Type TJSIterator;

    typedef typename Container<TJst>::Type TDeltaMap;
    typedef typename DeltaCoverage<TDeltaMap>::Type TCoverage;
    typedef typename Iterator<TDeltaMap, Rooted>::Type TBranchNodeIterator;

    typedef typename Position<TJst>::Type TPosition;
    typedef typename Size<TJst>::Type TSize;
    typedef typename MakeSigned<TSize>::Type TDeltaSize;

    // Auxiliary position information.
    TPosition           _branchProxyId;    // The id of the chosen proxy for the current branch.
    TPosition           _mappedHostPos;    // The mapped position to the reference to select split points efficiently.  TODO(rmaerker): Check if we need this?
    TDeltaSize          _proxyEndPosDiff;  // The diff between ref position and virutal position of window end.
    TBranchNodeIterator _firstWindowBranchNode;  // The first branch-node falling into the current window. TODO(rmaerker): Check if we need this?

    // Poxy data.
    TJSIterator _proxyIter;        // The iterator to the journal entry representing the current branch-node.
    TPosition   _proxyEndPos;      // The pruned end of the proxy.
    TDeltaSize  _prefixOffset;     // The offset between branch begin and virtual position of branch-node.
    TCoverage   _branchCoverage;   // The current coverage of the branch.

    // Additional data.
    TState      _externalState;    // The state of the external algorithm.

    JstBranchStackEntry_() : _branchProxyId(-1),
                             _mappedHostPos(-1),
                             _proxyEndPosDiff(0),
                             _firstWindowBranchNode(),
                             _proxyIter(),
                             _proxyEndPos(),
                             _prefixOffset(),
                             _branchCoverage(),
                             _externalState()
    {}

    // Copy constructor.
    JstBranchStackEntry_(JstBranchStackEntry_ const & other)
    {
        _copy(*this, other);
    }

    // Assignment Operator.
    JstBranchStackEntry_ & operator=(JstBranchStackEntry_ const & other)
    {
        if (this != &other)
            _copy(*this, other);
        return *this;
    }
};

// ----------------------------------------------------------------------------
// Class JstBranchStack_
// ----------------------------------------------------------------------------

template <typename TJournaledStringTree, typename TState>
class JstBranchStack_
{
public:

    typedef typename Value<JstBranchStack_>::Type TStackEntry;
    typedef typename MakeSigned<typename Position<TJournaledStringTree>::Type>::Type TPosition;
    typedef String<bool, Packed<> > TStackIndex;

    String<TStackEntry> _stack;
    TStackIndex         _stackIndex;
    TPosition           _activeId;

    JstBranchStack_() : _stack(), _stackIndex(), _activeId(-1)
    {}

    // Copy constructor.
    JstBranchStack_(JstBranchStack_ const & other)
    {
        _copy(*this, other);
    }

    // Assignment Operator.
    JstBranchStack_ & operator=(JstBranchStack_ const & other)
    {
        if (this != &other)
            _copy(*this, other);
        return *this;
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

template <typename TJournaledStringTree, typename TState>
struct Value<JstBranchStack_<TJournaledStringTree, TState> >
{
    typedef JstBranchStackEntry_<TJournaledStringTree, TState> Type;
};

template <typename TJournaledStringTree, typename TState>
struct Value<JstBranchStack_<TJournaledStringTree, TState> const >
{
    typedef JstBranchStackEntry_<TJournaledStringTree, TState> const Type;
};

// ----------------------------------------------------------------------------
// Metafunction Reference
// ----------------------------------------------------------------------------

template <typename TJournaledStringTree, typename TState>
struct Reference<JstBranchStack_<TJournaledStringTree, TState> >
{
    typedef JstBranchStackEntry_<TJournaledStringTree, TState> & Type;
};

template <typename TJournaledStringTree, typename TState>
struct Reference<JstBranchStack_<TJournaledStringTree, TState> const >
{
    typedef JstBranchStackEntry_<TJournaledStringTree, TState> const & Type;
};

// ----------------------------------------------------------------------------
// Metafunction Size
// ----------------------------------------------------------------------------

template <typename TJournaledStringTree, typename TState>
struct Size<JstBranchStack_<TJournaledStringTree, TState> >
{
    typedef typename Size<TJournaledStringTree>::Type Type;
};

template <typename TJournaledStringTree, typename TState>
struct Size<JstBranchStack_<TJournaledStringTree, TState> const > :
    Size<JstBranchStack_<TJournaledStringTree, TState> >{};

// ----------------------------------------------------------------------------
// Metafunction Position
// ----------------------------------------------------------------------------

template <typename TJournaledStringTree, typename TState>
struct Position<JstBranchStack_<TJournaledStringTree, TState> >
{
    typedef typename Position<TJournaledStringTree>::Type Type;
};

template <typename TJournaledStringTree, typename TState>
struct Position<JstBranchStack_<TJournaledStringTree, TState> const > :
    Position<JstBranchStack_<TJournaledStringTree, TState> >{};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _copy()
// ----------------------------------------------------------------------------

template <typename TJournaledStringTree, typename TState>
inline void
_copy(JstBranchStackEntry_<TJournaledStringTree, TState> & entry,
      JstBranchStackEntry_<TJournaledStringTree, TState> const & other)
{
    entry._branchProxyId = other._branchProxyId;
    entry._mappedHostPos = other._mappedHostPos;
    entry._proxyEndPosDiff = other._proxyEndPosDiff;
    entry._firstWindowBranchNode = other._firstWindowBranchNode;

    // Poxy data.
    entry._proxyIter = other._proxyIter;
    entry._proxyEndPos = other._proxyEndPos;
    entry._prefixOffset = other._prefixOffset;
    entry._branchCoverage = other._branchCoverage;

    // Additional data.
    entry._externalState = other._externalState;
}

// ----------------------------------------------------------------------------
// Function _copy()
// ----------------------------------------------------------------------------

template <typename TJournaledStringTree, typename TState>
inline void
_copy(JstBranchStack_<TJournaledStringTree, TState> & stack,
      JstBranchStack_<TJournaledStringTree, TState> const & other)
{
    stack._stack = other._stack;
    stack._stackIndex = other._stackIndex;
    stack._activeId = other._activeId;
}


template <typename TJst, typename TState>
inline typename Reference<JstBranchStack_<TJst, TState> >::Type
createEntry(JstBranchStack_<TJst, TState> & stack)
{
    typedef JstBranchStack_<TJst, TState> TStack;
    typedef typename Position<TStack>::Type TPos;
    typedef typename TStack::TStackIndex TStackIndex;
    typedef typename Position<TStackIndex>::Type TIndPos;

    TPos id = bitScanReverse(stack._stackIndex);
    if (id == MaxValue<TIndPos>::VALUE || ++id >= length(stack._stackIndex))
    {
        id = length(stack._stackIndex); // Set id to begin of new block.
        resize(stack._stackIndex, length(stack._stackIndex) + 64, false, Exact());
    }

    if (id >= length(stack._stack))
        resize(stack._stack, id + 1);

    stack._stackIndex[id] = true;
    return stack._stack[id];
}

template <typename TJst, typename TState>
inline typename Reference<JstBranchStack_<TJst, TState> >::Type
createInitialEntry(JstBranchStack_<TJst, TState> & stack)
{
    clear(stack._stackIndex);
    stack._activeId = 0;
    return createEntry(stack);
}


template <typename TJst, typename TState>
inline bool
pop(JstBranchStack_<TJst, TState> & stack)
{
    SEQAN_ASSERT_GT(length(stack._stack), static_cast<unsigned>(stack._activeId));

    if (stack._activeId > -1 && !empty(stack._stack))
        stack._stackIndex[stack._activeId] = false;

    if (testAllZeros(stack._stackIndex))
    {
        stack._activeId = -1;
        return false;
    }

    stack._activeId = bitScanReverse(stack._stackIndex);
    return true;
}

template <typename TJst, typename TState>
inline typename Reference<JstBranchStack_<TJst, TState> >::Type
top(JstBranchStack_<TJst, TState> & stack)
{
    SEQAN_ASSERT_NOT(stack._activeId == -1);
    SEQAN_ASSERT_GT(length(stack._stack), static_cast<unsigned>(stack._activeId));

    return stack._stack[stack._activeId];
}

}
#endif // EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_JOURNALED_STRING_TREE_TRAVERSAL_BRANCH_STACK_H_
