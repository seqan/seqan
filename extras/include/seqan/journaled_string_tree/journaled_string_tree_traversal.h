// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2014, Knut Reinert, FU Berlin
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
// Implements a traversal interface to efficiently traverse a set of
// journaled sequences simultaneously.
// ==========================================================================

#ifndef EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_JOURNALED_STRING_TREE_TRAVERSAL_H_
#define EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_JOURNALED_STRING_TREE_TRAVERSAL_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================


/*!
 * @defgroup JstTraverserContextPositionTags JST Context Position Tags
 * @brief Tag to specialize the location of the context iterator.
 *
 *
 * @tag JstTraverserContextPositionTags#ContextPositionLeft
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Iterator always points to begin of context.
 *
 * @signature struct ContextPositionLeft_;
 * @signature typedef Tag<ContextPositionLeft_> ContextPositionLeft;
 *
 * @tag JstTraverserContextPositionTags#ContextPositionRight
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Iterator always points to end of context.
 *
 * @signature struct ContextPositionRight_;
 * @signature typedef Tag<ContextPositionRight_> ContextPositionRight;
 */

enum JstTraversalState
{
    JST_TRAVERSAL_STATE_NULL,
    JST_TRAVERSAL_STATE_MASTER,
    JST_TRAVERSAL_STATE_BRANCH
};

/*!
 * @defgroup JstTraversalStates JST Traversal Tags
 * @brief Tags for selecting state of the Journaled String Tree traversal.
 *
 * @tag JstTraversalStates#StateTraverseMaster
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Tag for selecting master strand.
 *
 * @signature struct TraverseStateMaster_;
 * @signature typedef Tag<TraverseStateMaster_> StateTraverseMaster;
 *
 * @tag JstTraversalStates#StateTraverseBranch
 * @headerfile <seqan/journaled_string_tree.h>
 *
 * @brief Tag for selecting branch strand.
 *
 * @signature struct TraverseStateBranch_;
 * @signature typedef Tag<TraverseStateBranch_> StateTraverseBranch;
 */

struct TraverseStateMaster_;
typedef Tag<TraverseStateMaster_> StateTraverseMaster;

struct TraverseStateBranch_;
typedef Tag<TraverseStateBranch_> StateTraverseBranch;

template <typename TContainer, typename TState, typename TSpec>
class JstTraverser;

/*
 * @class JstTraverserCore_
 * @headerfile <seqan/journaled_string_tree.h>
 *
 * @brief State stored on a concurrent queue to allow several tasks run in parallel.
 */

template <typename TMasterIt, typename TCoverage, typename TBranchNodeIt, typename TMergeStack>
class JstTraverserCore_
{
public:
    TMasterIt     _masterIt;
    TCoverage     _activeMasterCoverage;  // Active master coverage.

    TBranchNodeIt _branchNodeIt;
    TBranchNodeIt _branchNodeInContextIt;  // Points to left node within context or behind the context.
    TMergeStack   _mergePointStack;  // Stores merge points, when deletions are connected to the master branch.

    JstTraverserCore_()
    {}

    template <typename TDeltaMap, typename TSpec>
    JstTraverserCore_(JournaledStringTree<TDeltaMap, TSpec> & hystk) : _mergePointStack(container(hystk))
    {}
};

/*!
 * @class JstTraverser
 * @headerfile <seqan/journaled_string_tree.h>
 *
 * @brief Manages traversal over the @link JournaledStringTree @endlink.
 *
 * @signature template <tyepanme TContainer, typename TState, typename TContextPosition, typename TRequireFullContext>
 *            class JstTraverser<TContainer, TState, JstTraverserConfig<TContextPosition, TRequireFullContext>;
 *
 * @tparam TContainer        The type of the container to be traversed. Must be of type @link JournaledStringTree @endlink
 * @tparam TState            The state to be restored for an external algorithm.
 * @tparam TContextPosition  Tag to specialize the location of the context iterator. See @link JstTraverserContextPositionTags @endlink.
 * @tparam TRequireFullContext  Tag to specialize whether the entire context needs to be accessible. One of @link LogicalValuesTags @endlink.
 *
 *  The JstTraverser manages the traversal over a Journaled String Tree. It stores the current state of the traversal
 *  and provides additional interfaces to access the values at the current position. It encapsulates an forward-iterator
 *  to simultaneous scan over the sequences contained in the Journaled String Tree from left to right.
 *
 * Some algorithms that work on a local context are either prefix or suffix based. This means, that they either explore
 * the context in forward direction or in backward direction. An example for prefix based algorithms is the naiive
 * online-search algorithm which scans the context from left to right. On the other hand the horspool algorithm
 * scans the context from right to left. To specialize the iterator position within the current context one can
 * use the tag <tt>TContextPosition</tt> to make the traversal dependent on the preferences of the external algorithm
 * evaluating the sequence context.
 * In addition to the context position it is necessary to distinguish between algorithms that require the entire context
 * to be accessible and those that do not. For example the Myer's bitvector algorithm only needs the last position of
 * the current window, while it uses a state to keep track of the entire sequence context.
 *
 *
 * To traverse a @link JournaledStringTree @endlink (JST) one has to provide an external algorithm that implements the
 * @link JstTraversalConcept @endlink. The function @link JstTraverser#traverse @endlink triggers the iteration
 * over the JST and interrupts the iteration every time a new context is explored.
 *
 */

/*!
 * @fn JstTraverser::JstTraverser
 * @brief constructor
 * @signature JstTraverser::JstTraverser();
 * @signature JstTraverser::JstTraverser(container, w);
 *
 * @param container The container to be traversed. Must be of type @link JournaledStringTree @endlink.
 * @param w         The size of the context.
 */

template <typename TDeltaMap, typename TTreeSpec, typename TState, typename TSpec>
class JstTraverser<JournaledStringTree<TDeltaMap, TTreeSpec>, TState, TSpec>
{
public:
    typedef JstTraverser<JournaledStringTree<TDeltaMap, TTreeSpec>, TState, TSpec> TTraverser;
    typedef JournaledStringTree<TDeltaMap, TTreeSpec> TContainer;
    typedef typename GetStringSet<TContainer>::Type TJournalSet;
    typedef typename Host<TJournalSet>::Type TReference;
    typedef typename Iterator<TReference, Rooted>::Type TMasterBranchIterator;
    typedef typename Container<TContainer>::Type TVariantData;

    typedef typename Value<TJournalSet>::Type TJournalString;
    typedef typename Iterator<TJournalString>::Type TJournalIterator;

    typedef typename Iterator<TDeltaMap, Standard>::Type TBranchNodeIterator;
    typedef typename DeltaCoverage<TDeltaMap>::Type TBitVector;

    typedef MergePointMap_<TVariantData> TMergePointStore;
    typedef JstBranchStack_<TContainer, TState> TBranchStack;

    typedef typename Size<TJournalSet>::Type TSize;
    typedef typename Position<TContainer>::Type TPosition;

    typedef JstTraverserCore_<TMasterBranchIterator, TBitVector, TBranchNodeIterator, TMergePointStore> TJstTraverserCore;

    // Basics.
    JstTraversalState _traversalState;
    TContainer * _haystackPtr;  // Pointer to the underlying data parallel facade.

    // Core Traverser.
    mutable TJstTraverserCore   _traverserCore;

    // Sequence iterators.
    TMasterBranchIterator _masterItEnd;
    TJournalIterator _branchIt;

    // Coverage information.
    mutable TBitVector _activeBranchCoverage;  // Active coverage of the branch.

    // Branch-node information.
    TBranchNodeIterator _branchNodeBlockEnd;
    TBranchNodeIterator _proxyBranchNodeIt;

    // Auxiliary structures.
    mutable TBranchStack     _branchStack;  // Handles the branches of the current tree.
    TSize _contextSize;
    bool _needInit;
    mutable bool _isSynchronized;
    TState _lastMasterState;

    JstTraverser() :
        _traversalState(JST_TRAVERSAL_STATE_NULL),
        _haystackPtr((TContainer*) 0),
        _contextSize(1),
        _needInit(true),
        _isSynchronized(false)
    {}

    JstTraverser(TContainer & haystack, TSize contextSize) :
        _traversalState(JST_TRAVERSAL_STATE_NULL),
        _traverserCore(haystack),
        _contextSize(contextSize),
        _needInit(false),
        _isSynchronized(false)
    {
        init(*this, haystack);
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Container
// ----------------------------------------------------------------------------

/*!
 * @mfn JstTraverser#Container
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Returns the container type.
 *
 * @signature Container<TJstTraverser>::Type;
 * @tparam TJstTraverser The type of the traverser.
 *
 * @return The container type.
 */

template <typename TContainer, typename TState, typename TSpec>
struct Container<JstTraverser<TContainer, TState, TSpec> >
{
    typedef TContainer Type;
};

template <typename TContainer, typename TState, typename TSpec>
struct Container<JstTraverser<TContainer, TState, TSpec> const>
{
    typedef TContainer const Type;
};

// ----------------------------------------------------------------------------
// Metafunction Positions
// ----------------------------------------------------------------------------

/*!
 * @mfn JstTraverser#Position
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Returns the type used to store positions.
 *
 * @signature Position<TJstTraverser>::Type;
 * @tparam TJstTraverser The type of the traverser.
 *
 * @return TPosition The position type. Note this is a string of @link Pair @endlink structures where the first entry refers to
 * the sequence and the second entry to its current position.
 */

template <typename TContainer, typename TState, typename TSpec>
struct Position<JstTraverser<TContainer, TState, TSpec> >
{
    typedef typename GetStringSet<TContainer>::Type TJournalSet_;
    typedef typename Value<TJournalSet_>::Type TJournalString_;
    typedef typename Position<TJournalString_>::Type TPosition_;
    typedef Pair<TPosition_, TPosition_> TValue_;
    typedef String<TValue_> Type;
};

template <typename TContainer, typename TState, typename TSpec>
struct Position<JstTraverser<TContainer, TState, TSpec> const> :
    Position<JstTraverser<TContainer, TState, TSpec> >{};

// ----------------------------------------------------------------------------
// Metafunction Size
// ----------------------------------------------------------------------------

/*!
 * @mfn JstTraverser#Size
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Returns the size type.
 *
 * @signature Size<TJstTraverser>::Type;
 * @tparam TJstTraverser The type of the traverser.
 *
 * @return TSize The position type.
 */

template <typename TContainer, typename TState, typename TSpec>
struct Size<JstTraverser<TContainer, TState, TSpec> >
{
    typedef JstTraverser<TContainer, TState, TSpec> TJstTraverser_;
    typedef typename TJstTraverser_::TReference THost_;
    typedef typename TJstTraverser_::TJournalString TBranch_;

    typedef typename Size<THost_>::Type THostSize_;
    typedef typename Size<TBranch_>::Type TBranchSize_;

    typedef typename IfC<sizeof(THostSize_) >= sizeof(TBranchSize_),
                         THostSize_,
                         TBranchSize_>::Type Type;
};

template <typename TContainer, typename TState, typename TSpec>
struct Size<JstTraverser<TContainer, TState, TSpec> const> :
    Size<JstTraverser<TContainer, TState, TSpec> >{};

// ----------------------------------------------------------------------------
// Metafunction Iterator                                  [StateTraverseMaster]
// ----------------------------------------------------------------------------

/*!
 * @mfn JstTraverser#Iterator
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Returns the type of the iterator of either the master or the branch strand.
 *
 * @signature Iterator<TJstTraverser, TTag>::Type;
 * @tparam TJstTraverser The type of the traverser.
 * @tparam TTag          The tag to specify the correct iterator. Must be one of @link JstTraversalStates @endlink.
 *
 * @return TIterator The type of the iterator.
 */

template <typename TContainer, typename TState, typename TSpec>
struct Iterator<JstTraverser<TContainer, TState, TSpec>,  StateTraverseMaster>
{
    typedef JstTraverser<TContainer, TState, TSpec> TTraverser_;
    typedef typename TTraverser_::TMasterBranchIterator Type;
};

template <typename TContainer, typename TState, typename TSpec>
struct Iterator<JstTraverser<TContainer, TState, TSpec> const,  StateTraverseMaster>
{
    typedef JstTraverser<TContainer, TState, TSpec> const TTraverser_;
    typedef typename TTraverser_::TMasterBranchIterator Type;
};

// ----------------------------------------------------------------------------
// Metafunction Iterator                                  [StateTraverseBranch]
// ----------------------------------------------------------------------------

template <typename TContainer, typename TState, typename TSpec>
struct Iterator<JstTraverser<TContainer, TState, TSpec>,  StateTraverseBranch>
{
    typedef JstTraverser<TContainer, TState, TSpec> TTraverser_;
    typedef typename TTraverser_::TJournalIterator Type;
};

template <typename TContainer, typename TState, typename TSpec>
struct Iterator<JstTraverser<TContainer, TState, TSpec> const,  StateTraverseBranch>
{
    typedef JstTraverser<TContainer, TState, TSpec> const TTraverser_;
    typedef typename TTraverser_::TJournalIterator Type;
};

// ----------------------------------------------------------------------------
// Metafunction BranchNode
// ----------------------------------------------------------------------------

template <typename TJstTraverser>
struct BranchNode;

template <typename TContainer, typename TState, typename TSpec>
struct BranchNode<JstTraverser<TContainer, TState, TSpec> >
{
    typedef typename Container<TContainer>::Type TDeltaMap_;
    typedef typename Iterator<TDeltaMap_, Rooted>::Type Type;
};

template <typename TContainer, typename TState, typename TSpec>
struct BranchNode<JstTraverser<TContainer, TState, TSpec> const>
{
    typedef typename Container<TContainer>::Type TDeltaMap_;
    typedef typename Iterator<TDeltaMap_ const, Rooted>::Type Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _contextBeginPosition()                        [StateTraverseMaster]
// ----------------------------------------------------------------------------

template <typename TContainer, typename TState, typename TContextPos, typename TRequireFullContext>
inline typename MakeSigned<
       typename Size<JstTraverser<TContainer, TState, JstTraverserConfig<TContextPos, TRequireFullContext> > >::Type
       >::Type
_contextBeginPosition(JstTraverser<TContainer, TState, JstTraverserConfig<TContextPos, TRequireFullContext> > const & traverser,
                      StateTraverseMaster const & /*tag*/)
{
    typedef typename Position<TContainer>::Type TPosition;
    typedef typename MakeSigned<TPosition>::Type TPos_;

    if (IsSameType<TContextPos, ContextPositionLeft>::VALUE)
        return position(traverser._traverserCore._masterIt);
    else
        return static_cast<TPos_>(position(traverser._traverserCore._masterIt))
            - static_cast<TPos_>(traverser._contextSize - 1);
}

// ----------------------------------------------------------------------------
// Function _contextBeginPosition()                        [StateTraverseBranch]
// ----------------------------------------------------------------------------

template <typename TContainer, typename TState, typename TContextPos, typename TRequireFullContext>
inline typename MakeSigned<
       typename Size<JstTraverser<TContainer, TState, JstTraverserConfig<TContextPos, TRequireFullContext> > >::Type
       >::Type
_contextBeginPosition(JstTraverser<TContainer, TState, JstTraverserConfig<TContextPos, TRequireFullContext> > const & traverser,
                      StateTraverseBranch const & /*tag*/)
{
    typedef typename Position<TContainer>::Type TPosition;
    typedef typename MakeSigned<TPosition>::Type TPos_;

    if (IsSameType<TContextPos, ContextPositionLeft>::VALUE)
        return position(traverser._branchIt);
    else
        return static_cast<TPos_>(position(traverser._branchIt)) - static_cast<TPos_>(traverser._contextSize - 1);
}

// ----------------------------------------------------------------------------
// Function _clippedContextBeginPosition()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TState, typename TSpec, typename TTraversalState>
inline typename MakeSigned<typename Size<JstTraverser<TContainer, TState, TSpec> >::Type>::Type
_clippedContextBeginPosition(JstTraverser<TContainer, TState, TSpec> const & traverser,
                             TTraversalState const & /*tag*/)
{
    return _max(0, _contextBeginPosition(traverser, TTraversalState()));
}

// ----------------------------------------------------------------------------
// Function _contextEndPosition()                          [StateTraverseMaster]
// ----------------------------------------------------------------------------

template <typename TContainer, typename TState, typename TContextPos, typename TRequireFullContext>
inline typename MakeSigned<
       typename Size<JstTraverser<TContainer, TState, JstTraverserConfig<TContextPos, TRequireFullContext> > >::Type
       >::Type
_contextEndPosition(JstTraverser<TContainer, TState, JstTraverserConfig<TContextPos, TRequireFullContext> > const & traverser,
                    StateTraverseMaster const & /*tag*/)
{
    if (IsSameType<TContextPos, ContextPositionLeft>::VALUE)
        return position(traverser._traverserCore._masterIt) + (traverser._contextSize - 1);
    else
        return position(traverser._traverserCore._masterIt);
}

// ----------------------------------------------------------------------------
// Function _contextEndPosition()                          [StateTraverseBranch]
// ----------------------------------------------------------------------------

template <typename TContainer, typename TState, typename TContextPos, typename TRequireFullContext>
inline typename MakeSigned<
       typename Size<JstTraverser<TContainer, TState, JstTraverserConfig<TContextPos, TRequireFullContext> > >::Type
       >::Type
_contextEndPosition(JstTraverser<TContainer, TState, JstTraverserConfig<TContextPos, TRequireFullContext> > const & traverser,
                    StateTraverseBranch const & /*tag*/)
{
    if (IsSameType<TContextPos, ContextPositionLeft>::VALUE)
        return position(traverser._branchIt) + (traverser._contextSize - 1);
    else
        return position(traverser._branchIt);
}

// ----------------------------------------------------------------------------
// Function _clippedContextEndPosition()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TState, typename TSpec, typename TTraversalState>
inline typename MakeSigned<typename Size<JstTraverser<TContainer, TState, TSpec> >::Type>::Type
_clippedContextEndPosition(JstTraverser<TContainer, TState, TSpec> const & traverser,
                           TTraversalState const & /*tag*/)
{
    if (IsSameType<TTraversalState, StateTraverseMaster>::VALUE)
        return _min(length(host(container(traverser))), _contextEndPosition(traverser, TTraversalState()));
    else
        return _min(length(*traverser._branchIt._journalStringPtr), _contextEndPosition(traverser, TTraversalState()));
}

// ----------------------------------------------------------------------------
// Function contextBegin()                                 [StateTraverseMaster]
// ----------------------------------------------------------------------------

/*!
 * @fn JstTraverser#contextBegin
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Returns an iterator pointing to the begin of the current context.
 *
 * @signature TIterator contextBegin(traverser, tag);
 * @param[in] traverser The traverser to query the iterator pointing to the begin of the context for.
 * @param[in] tag       Tag indicating the source of the iterator. Must be one of @link JstTraversalStates @endlink.
 *
 * @return TIterator An iterator of type @link JstTraverser#Iterator @endlink pointing to the begin of the current
 * context in either the master or the current branch segment.
 *
 * @note The result is undefined if the iterator for the branch segment is requested and the current traversal is not
 * in the branch mode. For the master segment the iterator is always valid. To check for the current state of the
 * traverser use the functions @link JstTraverser#isMasterState @endlink and @link JstTraverser#isBranchState @endlink.
 *
 * @see JstTraverser#contextEnd
 * @see JstTraverser#contextIterator
 */

// ContextPositionLeft.
template <typename TContainer, typename TState, typename TFullContext>
inline typename Iterator<JstTraverser<TContainer, TState, JstTraverserConfig<ContextPositionLeft, TFullContext> >,  StateTraverseMaster>::Type
contextBegin(JstTraverser<TContainer, TState, JstTraverserConfig<ContextPositionLeft, TFullContext> > const & traverser,
            StateTraverseMaster const & /*tag*/)
{
    return traverser._traverserCore._masterIt;
}

// ContextPositionRight.
template <typename TContainer, typename TState, typename TFullContext>
inline typename Iterator<JstTraverser<TContainer, TState, JstTraverserConfig<ContextPositionRight, TFullContext> >,  StateTraverseMaster>::Type
contextBegin(JstTraverser<TContainer, TState, JstTraverserConfig<ContextPositionRight, TFullContext> > const & traverser,
            StateTraverseMaster const & /*tag*/)
{
    return traverser._traverserCore._masterIt - (traverser._contextSize - 1);
}

// ----------------------------------------------------------------------------
// Function contextBegin()                                 [StateTraverseBranch]
// ----------------------------------------------------------------------------

// ContextPositionLeft.
template <typename TContainer, typename TState, typename TFullContext>
inline typename Iterator<JstTraverser<TContainer, TState, JstTraverserConfig<ContextPositionLeft, TFullContext> >,  StateTraverseBranch>::Type
contextBegin(JstTraverser<TContainer, TState, JstTraverserConfig<ContextPositionLeft, TFullContext> > const & traverser,
            StateTraverseBranch const & /*tag*/)
{
    return traverser._branchIt;
}

// ContextPositionRight.
template <typename TContainer, typename TState, typename TFullContext>
inline typename Iterator<JstTraverser<TContainer, TState, JstTraverserConfig<ContextPositionRight, TFullContext> >,  StateTraverseBranch>::Type
contextBegin(JstTraverser<TContainer, TState, JstTraverserConfig<ContextPositionRight, TFullContext> > const & traverser,
            StateTraverseBranch const & /*tag*/)
{
    return traverser._branchIt - (traverser._contextSize - 1);
}

// ----------------------------------------------------------------------------
// Function contextEnd()                                   [StateTraverseMaster]
// ----------------------------------------------------------------------------

/*!
 * @fn JstTraverser#contextEnd
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Returns an iterator pointing to the end of the current context.
 *
 * @signature TIterator contextEnd(traverser, tag);
 * @param[in] traverser The traverser to query the iterator pointing to the end of the context for.
 * @param[in] tag       Tag indicating the source of the iterator. Must be one of @link JstTraversalStates @endlink.
 *
 * @return TIterator An iterator of type @link JstTraverser#Iterator @endlink pointing to the end of the current
 * context in either the master or the current branch segment.
 *
 * @note The result is undefined if the iterator for the branch segment is requested and the current traversal is not
 * in the branch mode. For the master segment the iterator is always valid. To check for the current state of the
 * traverser use the functions @link JstTraverser#isMasterState @endlink and @link JstTraverser#isBranchState @endlink.
 *
 * @see JstTraverser#contextBegin
 * @see JstTraverser#contextIterator
 */

// ContextPositionLeft.
template <typename TContainer, typename TState, typename TFullContext>
inline typename Iterator<JstTraverser<TContainer, TState, JstTraverserConfig<ContextPositionLeft, TFullContext> >,  StateTraverseMaster>::Type
contextEnd(JstTraverser<TContainer, TState, JstTraverserConfig<ContextPositionLeft, TFullContext> > const & traverser,
          StateTraverseMaster const & /*tag*/)
{
    return traverser._traverserCore._masterIt + (traverser._contextSize - 1);
}

// ContextPositionRight.
template <typename TContainer, typename TState, typename TFullContext>
inline typename Iterator<JstTraverser<TContainer, TState, JstTraverserConfig<ContextPositionRight, TFullContext> >,  StateTraverseMaster>::Type
contextEnd(JstTraverser<TContainer, TState, JstTraverserConfig<ContextPositionRight, TFullContext> > const & traverser,
          StateTraverseMaster const & /*tag*/)
{
    return traverser._traverserCore._masterIt;
}

// ----------------------------------------------------------------------------
// Function contextEnd()                                   [StateTraverseBranch]
// ----------------------------------------------------------------------------

template <typename TContainer, typename TState, typename TFullContext>
inline typename Iterator<JstTraverser<TContainer, TState, JstTraverserConfig<ContextPositionLeft, TFullContext> >,  StateTraverseBranch>::Type
contextEnd(JstTraverser<TContainer, TState, JstTraverserConfig<ContextPositionLeft, TFullContext> > const & traverser,
          StateTraverseBranch const & /*tag*/)
{
    return traverser._branchIt + (traverser._contextSize - 1);
}

template <typename TContainer, typename TState, typename TFullContext>
inline typename Iterator<JstTraverser<TContainer, TState, JstTraverserConfig<ContextPositionRight, TFullContext> >,  StateTraverseBranch>::Type
contextEnd(JstTraverser<TContainer, TState, JstTraverserConfig<ContextPositionRight, TFullContext> > const & traverser,
          StateTraverseBranch const & /*tag*/)
{
    return traverser._branchIt;
}

// ----------------------------------------------------------------------------
// Function contextIterator()
// ----------------------------------------------------------------------------

/*!
 * @fn JstTraverser#contextIterator
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Returns the iterator to the current context.
 *
 * @signature TIterator contextIterator(traverser, tag);
 * @param[in] traverser The traverser to query current context iterator for.
 * @param[in] tag       Tag indicating the source of the iterator. Must be one of @link JstTraversalStates @endlink.
 *
  * @return TIterator An iterator of type @link JstTraverser#Iterator @endlink pointing to the current of the current
 * context in either the master or the current branch segment.
 *
 * @note The result is undefined if the iterator for the branch segment is requested and the current traversal is not
 * in the branch mode. For the master segment the iterator is always valid. To check for the current state of the
 * traverser use the functions @link JstTraverser#isMasterState @endlink and @link JstTraverser#isBranchState @endlink.
 *
 * @see JstTraverser#contextBegin
 * @see JstTraverser#contextEnd
 */

// StateTraverseMaster.
template <typename TContainer, typename TState, typename TSpec>
inline typename Iterator<JstTraverser<TContainer, TState, TSpec>, StateTraverseMaster>::Type &
contextIterator(JstTraverser<TContainer, TState, TSpec> & traverser,
                StateTraverseMaster const & /*tag*/)
{
    return traverser._traverserCore._masterIt;
}

template <typename TContainer, typename TState, typename TSpec>
inline typename Iterator<JstTraverser<TContainer, TState, TSpec> const, StateTraverseMaster>::Type &
contextIterator(JstTraverser<TContainer, TState, TSpec> const & traverser,
                StateTraverseMaster const & /*tag*/)
{
    return traverser._traverserCore._masterIt;
}

// StateTraverseBranch.
template <typename TContainer, typename TState, typename TSpec>
inline typename Iterator<JstTraverser<TContainer, TState, TSpec>, StateTraverseBranch>::Type &
contextIterator(JstTraverser<TContainer, TState, TSpec> & traverser,
                StateTraverseBranch const & /*tag*/)
{
    return traverser._branchIt;
}

template <typename TContainer, typename TState, typename TSpec>
inline typename Iterator<JstTraverser<TContainer, TState, TSpec> const, StateTraverseBranch>::Type &
contextIterator(JstTraverser<TContainer, TState, TSpec> const & traverser,
                StateTraverseBranch const & /*tag*/)
{
    return traverser._branchIt;
}

// ----------------------------------------------------------------------------
// Function _globalInit()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TState, typename TSpec>
inline void
_globalInit(JstTraverser<TContainer, TState, TSpec> & traverser)
{
    resize(traverser._traverserCore._activeMasterCoverage, getCoverageSize(container(container(traverser))), true, Exact());  // Set all seq active.
    push(traverser._traverserCore._mergePointStack, length(host(container(traverser))) + 1, end(container(container(traverser)), Standard()));
}

// ----------------------------------------------------------------------------
// Function position()
// ----------------------------------------------------------------------------

/*!
 * @fn JstTraverser#position
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Returns a string containing all the sequences and their virtual positions that are active for the current
 * context.
 *
 * @signature TPosition position(traverser);
 * @param[in,out] traverser The traverser to query the position vector for.
 *
 * @return TPosition The vector containing the positions of all active sequences for the current context  of type
 * @link JstTraverser#Position @endlink. This string might also be empty if no sequences are valid.
 *
 * @note This function internally updates some states of the traverser despite of it's constness qualification.
 * To ensure that no internal states are updated call the function @link JstTraverser#sync @endlink explicitly before
 * calling this function.
 */

template <typename TContainer, typename TState, typename TContextPos, typename TContextBegin>
inline typename Position<JstTraverser<TContainer, TState, JstTraverserConfig<TContextPos, TContextBegin> > >::Type
position(JstTraverser<TContainer, TState,  JstTraverserConfig<TContextPos, TContextBegin> > const & traverser)
{
    typedef JstTraverser<TContainer, TState, JstTraverserConfig<TContextPos, TContextBegin> > TJstTraverser;
    typedef typename Position<TJstTraverser>::Type TPositions;
    typedef typename Value<TPositions>::Type TPositionValue;

    typedef typename GetStringSet<TContainer>::Type TJournalSet;
    typedef typename Value<TJournalSet>::Type TJournalString;
    typedef typename Iterator<TJournalString, Standard>::Type TJournalStringIt;

    typedef typename Container<TContainer>::Type TDeltaMap;
    typedef typename DeltaCoverage<TDeltaMap>::Type TCoverage;
    typedef typename Iterator<TCoverage, Standard>::Type TIterator;

    typedef typename Size<TJstTraverser>::Type TSize;

    TPositions posVec;

    if (traverser._traversalState == JST_TRAVERSAL_STATE_MASTER)
    {  // Easy case, since all sequences must be in a original node.
        if (!traverser._isSynchronized)
            _syncAndUpdateCoverage(traverser, StateTraverseMaster());

        TSize hostPos = position(traverser._traverserCore._masterIt); //position(contextBegin(traverser, StateTraverseMaster()));  // Current position within the host.
        TIterator itBegin = begin(traverser._traverserCore._activeMasterCoverage, Standard());
        TIterator itEnd = end(traverser._traverserCore._activeMasterCoverage, Standard());

        for (TIterator it = itBegin; it != itEnd; ++it)
        {
            if (getValue(it))
            {
                TSize seqId = it - itBegin;
                appendValue(posVec,
                            TPositionValue(seqId,
                                           localToGlobalPos(hostToVirtualPosition(value(stringSet(container(traverser)),
                                                            seqId), hostPos), seqId, container(traverser))));
            }
        }
    }
    else
    {  // Harder case, the virtual positions cannot easily be re-based to a common host position.
        SEQAN_ASSERT_EQ(traverser._traversalState, JST_TRAVERSAL_STATE_BRANCH);

        if (!traverser._isSynchronized)
            _syncAndUpdateCoverage(traverser, StateTraverseBranch());

        TIterator itBegin = begin(traverser._activeBranchCoverage, Standard());
        TIterator itEnd = end(traverser._activeBranchCoverage, Standard());

        for (TIterator it = itBegin; it != itEnd; ++it)
        {
            if (getValue(it))
            {
                TSize seqId = it - itBegin;
                TJournalStringIt journalIt;
                // This is ok, since we do not change the underlying string.
                journalIt._journalStringPtr = &traverser._haystackPtr->_journalSet[seqId];
                journalIt = begin(*journalIt._journalStringPtr, Standard());
                _mapVirtualToVirtual(journalIt, traverser._branchIt, top(traverser._branchStack)._proxyIter,
                                     traverser._traverserCore._branchNodeIt, container(container(traverser)), seqId);
//                if (IsSameType<TContextPos, ContextPositionRight>::VALUE)
//
//                else
//                    _mapVirtualToVirtual(journalIt, contextBegin(traverser, StateTraverseBranch()),
//                                         traverser._traverserCore._branchNodeIt, container(container(traverser)), seqId);
                appendValue(posVec, TPositionValue(seqId, localToGlobalPos(position(journalIt), seqId,
                                                                           container(traverser))));
            }
        }
    }
    return posVec;
}

/*!
 * @fn JstTraverser#sync
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Synchronizes the tarversal state to the current context position.
 *
 * @signature sync(traverser);
 * @param traverser The object managing the traversal, whose state is to be synchronized.
 *
 * This function updates the internal state of the traversal and synchronizes to the current window.
 * This includes synchronizing the active coverage as well as to set some pointers.
 */

template <typename TContainer, typename TState, typename TSpec>
inline void
sync(JstTraverser<TContainer, TState, TSpec> & traverser)
{
    if (isMasterState(traverser))
        _syncAndUpdateCoverage(traverser, StateTraverseMaster());
    else
        _syncAndUpdateCoverage(traverser, StateTraverseBranch());
}

// ----------------------------------------------------------------------------
// Function state()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TState, typename TSpec>
inline JstTraversalState
state(JstTraverser<TContainer, TState, TSpec> const & traverser)
{
    return traverser._traversalState;
}

// ----------------------------------------------------------------------------
// Function isMasterState()
// ----------------------------------------------------------------------------

/*!
 * @fn JstTraverser#isMasterState
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Returns a <tt>bool</tt> indicating if the traverser is operating on the master strand.
 *
 * @signature bool isMasterState(traverser);
 * @param[in] traverser The traverser to check the current state for.
 *
 * @return <tt>true</tt> if the current context is in the master strand, otherwise <tt>false</tt>.
 *
 * @see JstTraverser#isBranchState
 */

template <typename TContainer, typename TState, typename TSpec>
inline bool
isMasterState(JstTraverser<TContainer, TState, TSpec> const & traverser)
{
    return state(traverser) == JST_TRAVERSAL_STATE_MASTER;
}

// ----------------------------------------------------------------------------
// Function isBranchState()
// ----------------------------------------------------------------------------

/*!
 * @fn JstTraverser#isBranchState
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Returns a <tt>bool</tt> indicating if the traverser is operating on the branch strand.
 *
 * @signature bool isBranchState(traverser);
 * @param[in] traverser The traverser to check the current state for.
 *
 * @return <tt>true</tt> if the current context is in the branch strand, otherwise <tt>false</tt>.
 *
 * @see JstTraverser#isBranchState
 */

template <typename TContainer, typename TState, typename TSpec>
inline bool
isBranchState(JstTraverser<TContainer, TState, TSpec> const & traverser)
{
    return state(traverser) == JST_TRAVERSAL_STATE_BRANCH;
}

// ----------------------------------------------------------------------------
// Function coverage()
// ----------------------------------------------------------------------------

/*!
 * @fn JstTraverser#coverage
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Returns the coverage for the current window.
 *
 * @signature TCoverage coverage(traverser);
 * @param[in]   traverser The traverser to get the active coverage for.
 *
 * @return TCoverage The coverage for the current context of type @link DeltaMap#DeltaCoverage @endlink.
 */

template <typename TContainer, typename TState, typename TSpec>
inline typename DeltaCoverage<typename Container<TContainer const>::Type> ::Type &
coverage(JstTraverser<TContainer, TState, TSpec> const & traverser)
{
    if (isMasterState(traverser))
        return traverser._traverserCore._activeMasterCoverage;
    return traverser._activeBranchCoverage;
}

template <typename TContainer, typename TState, typename TSpec>
inline typename DeltaCoverage<typename Container<TContainer>::Type> ::Type &
coverage(JstTraverser<TContainer, TState, TSpec> & traverser)
{
    if (isMasterState(traverser))
        return traverser._traverserCore._activeMasterCoverage;
    return traverser._activeBranchCoverage;
}

// ----------------------------------------------------------------------------
// Function branchNode()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TState, typename TSpec>
inline typename BranchNode<JstTraverser<TContainer, TState, TSpec> >::Type &
branchNode(JstTraverser<TContainer, TState, TSpec> & traverser)
{
    return traverser._traverserCore._branchNodeIt;
}

template <typename TContainer, typename TState, typename TSpec>
inline typename BranchNode<JstTraverser<TContainer, TState, TSpec> const>::Type &
branchNode(JstTraverser<TContainer, TState, TSpec> const & traverser)
{
    return traverser._traverserCore._branchNodeIt;
}

// ----------------------------------------------------------------------------
// Function _selectNextSplitPoint()
// ----------------------------------------------------------------------------

template <typename TBranchStackNode, typename TNodeIt>
inline typename Size<TNodeIt>::Type
_selectNextSplitPoint(TBranchStackNode const & proxyWindow,
                      TNodeIt const & nodeIt,
                      TNodeIt const & branchNode)
{
#ifdef DEBUG_DATA_PARALLEL
//    int virtualMapping = (deltaPosition(nodeIt) - deltaPosition(branchNode));
    unsigned splitPointPos = position(proxyWindow._proxyIter) + (deltaPosition(nodeIt) - deltaPosition(branchNode)) -
                             proxyWindow._proxyEndPosDiff;
    std::cerr << "Virtual Positions: " << position(proxyWindow._proxyIter) << " to " << splitPointPos <<  std::endl;
    std::cerr << "Physical Positions: " << *branchNode << " to " << *nodeIt << std::endl;
    return splitPointPos;
#else
    return position(proxyWindow._proxyIter) + (deltaPosition(nodeIt) - deltaPosition(branchNode)) -
           proxyWindow._proxyEndPosDiff;
#endif
}

// ----------------------------------------------------------------------------
// Function _updateAuxiliaryBranchStructures()
// ----------------------------------------------------------------------------

template <typename TBranchStackEntry, typename TMapIter>
inline void
_updateAuxiliaryBranchStructures(TBranchStackEntry & branchEntry,
                                 TMapIter & mapIter)
{
    if (deltaType(mapIter) == DELTA_TYPE_DEL)
    {
        branchEntry._proxyEndPosDiff += deltaValue(mapIter, DeltaTypeDel());
        branchEntry._mappedHostPos += deltaValue(mapIter, DeltaTypeDel()) - 1;
    }
    else if (deltaType(mapIter) == DELTA_TYPE_INS)
    {
        branchEntry._proxyEndPosDiff -= static_cast<int>(length(deltaValue(mapIter, DeltaTypeIns())));
    }
    else if (deltaType(mapIter) == DELTA_TYPE_SV)
    {
        branchEntry._proxyEndPosDiff += deltaValue(mapIter, DeltaTypeSV()).i1;
        branchEntry._mappedHostPos += deltaValue(mapIter, DeltaTypeSV()).i1;
        branchEntry._proxyEndPosDiff -= static_cast<int>(length(deltaValue(mapIter, DeltaTypeSV()).i2));
    }
}

// ----------------------------------------------------------------------------
// Function _traverseBranch()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TState, typename TContextPosition, typename TRequireFullContext,
          typename TExternal, typename TDelegate>
inline void
_traverseBranch(JstTraverser<TContainer, TState, JstTraverserConfig<TContextPosition, TRequireFullContext> > & traverser,
                TExternal & externalAlg,
                TDelegate & delegate)
{
    typedef JstTraverserConfig<TContextPosition, TRequireFullContext> TJstTraverserConfig_;
    typedef JstTraverser<TContainer, TState, TJstTraverserConfig_ > TJstTraverser;
    typedef typename TJstTraverser::TBranchNodeIterator TBranchNodeIter;

    typedef typename Size<TJstTraverser>::Type TSize;
    typedef typename MakeSigned<TSize>::Type TSignedSize;

    typedef typename TJstTraverser::TBranchStack TBranchStack;
    typedef typename Value<TBranchStack>::Type TBranchStackEntry;

    typedef typename Container<TContainer>::Type TVariantMap;
    typedef typename DeltaCoverage<TVariantMap>::Type TBitVector;

#ifdef DEBUG_DATA_PARALLEL
    bool check = false;
    std::cerr << "o Search Branch: " <<  top(traverser._branchStack)._branchProxyId << std::endl;
    std::cerr << "Virtual Space: " << position(top(traverser._branchStack)._proxyIter) - top(traverser._branchStack)._prefixOffset << " - " << top(traverser._branchStack)._proxyEndPos << std::endl;
    std::cerr << "Break Point: " << position(top(traverser._branchStack)._proxyIter) << std::endl;
//        std::cerr << "Coverage: " << traverser._activeBranchCoverage << std::endl;
#endif

    // Select the next node.
    traverser._proxyBranchNodeIt = traverser._traverserCore._branchNodeIt;
    TBranchNodeIter nodeItEnd = end(container(container(traverser)), Standard()) - 1;

    TSize splitPointPos = 0;
    if (traverser._proxyBranchNodeIt != nodeItEnd)
    {
        // Move right until first node whose host pos is equal or greater than the current mappedHostPos.
        while (traverser._proxyBranchNodeIt != nodeItEnd && deltaPosition(traverser._proxyBranchNodeIt) < top(traverser._branchStack)._mappedHostPos)
            ++traverser._proxyBranchNodeIt;

        if (deltaPosition(traverser._proxyBranchNodeIt) < top(traverser._branchStack)._mappedHostPos)
        {
            splitPointPos = top(traverser._branchStack)._proxyEndPos;
            ++traverser._proxyBranchNodeIt; // Set to the end.
        }
        else
        {
            splitPointPos = _selectNextSplitPoint(top(traverser._branchStack), traverser._proxyBranchNodeIt, traverser._traverserCore._branchNodeIt);
        }
    }
    else
    {
        splitPointPos = top(traverser._branchStack)._proxyEndPos;
        ++traverser._proxyBranchNodeIt; // Set to the end.
    }

#ifdef DEBUG_DATA_PARALLEL
    std::cerr << "split point: " << splitPointPos << " ("<< deltaPosition(traverser._proxyBranchNodeIt) << ")" << std::endl;
#endif

    // Set the correct iterator here.
    if (IsSameType<TContextPosition, ContextPositionLeft>::VALUE)
    {
        if (top(traverser._branchStack)._prefixOffset < 0)
            traverser._branchIt = top(traverser._branchStack)._proxyIter + -(top(traverser._branchStack)._prefixOffset);  // Needed because of the iterator implementation.
        else
            traverser._branchIt = top(traverser._branchStack)._proxyIter - top(traverser._branchStack)._prefixOffset;
    }
    else
    {
        traverser._branchIt = top(traverser._branchStack)._proxyIter + ((traverser._contextSize -1) - top(traverser._branchStack)._prefixOffset);
    }

    // Note, we need the position of the iterator, because the journal string iterator never
    // exceeds the end of the journal string. And it might be faster than evaluating the iterator.
    TSize branchPos = _contextEndPosition(traverser, StateTraverseBranch());
#ifdef DEBUG_DATA_PARALLEL
    if (check == true)
    {
        std::cerr << "The character: " << static_cast<int>(*traverser._branchIt) << " with ordValue: " << ordValue(*traverser._branchIt) << " at position: " << *traverser._branchIt._journalEntriesIterator << std::endl;
        std::cerr << "The infix: " << top(traverser._branchStack)._proxyIter._journalStringPtr->_insertionBuffer << std::endl;
    }
#endif

    TBitVector splitVec;

    while(branchPos < top(traverser._branchStack)._proxyEndPos && branchPos < length(*(top(traverser._branchStack)._proxyIter._journalStringPtr)))  // Check end condition.
    {
        if (branchPos >= splitPointPos)
        {
            // First check if we need to update the positions first.
            // Get the variant coverage for the split point.
            TBitVector & varCoverage = deltaCoverage(traverser._proxyBranchNodeIt);
            transform(splitVec, varCoverage, top(traverser._branchStack)._branchCoverage, FunctorBitwiseAnd());

            // Check if found split point effects the current branch.
            if (!testAllZeros(splitVec) && splitVec != top(traverser._branchStack)._branchCoverage)
            {
                // Appending to the branch stack might invalidate the current pointer.
                TBranchStackEntry & splitProxyWindow = createEntry(traverser._branchStack);
                splitProxyWindow._mappedHostPos = deltaPosition(traverser._proxyBranchNodeIt);  // The host position of the split point, where we spawn the search tree.
                splitProxyWindow._proxyEndPosDiff = top(traverser._branchStack)._proxyEndPosDiff;
                splitProxyWindow._externalState = getState(externalAlg); // Current state before the split.
                splitProxyWindow._firstWindowBranchNode = top(traverser._branchStack)._firstWindowBranchNode;  // TODO(rmaerker): Do we need this?

                splitProxyWindow._prefixOffset = static_cast<TSignedSize>(position(top(traverser._branchStack)._proxyIter)) -
                                                 static_cast<TSignedSize>(_contextBeginPosition(traverser, StateTraverseBranch()));
                splitProxyWindow._proxyIter = top(traverser._branchStack)._proxyIter;

                // Check if current branch covers the new variant.
                if (splitVec[top(traverser._branchStack)._branchProxyId])
                {  // Split branch for all sequences that do not share the new delta.
                    // Set the coverage of the split branch.
                    transform(splitProxyWindow._branchCoverage, top(traverser._branchStack)._branchCoverage, splitVec,
                              FunctorNested<FunctorBitwiseAnd, FunctorIdentity, FunctorBitwiseNot>());
                    // Update the coverage of the current branch.
                    top(traverser._branchStack)._branchCoverage = splitVec;  // TODO(rmaerker): Use swap/move.

                    // Update auxiliary branch structures.
                    _updateAuxiliaryBranchStructures(top(traverser._branchStack), traverser._proxyBranchNodeIt);
                }
                else
                {
                    // Udpate the split branch coverage.
                    splitProxyWindow._branchCoverage = splitVec; // TODO(rmaerker): swap/move here.
                    // Update the current branch coverage.
                    transform(top(traverser._branchStack)._branchCoverage, top(traverser._branchStack)._branchCoverage, splitVec,
                              FunctorNested<FunctorBitwiseAnd, FunctorIdentity, FunctorBitwiseNot>());
                    // Get the type of the next variant.
                    ++splitProxyWindow._mappedHostPos;
                    _updateAuxiliaryBranchStructures(splitProxyWindow, traverser._proxyBranchNodeIt);
                }

#ifdef DEBUG_DATA_PARALLEL
                std::cerr << "-> split branch proxy: " << splitProxyWindow._branchProxyId << std::endl;
                std::cerr << "-> split branch vp: " << position(splitProxyWindow._proxyIter) - splitProxyWindow._prefixOffset << " - " << splitProxyWindow._proxyEndPos<< std::endl;
                std::cerr << "-> Original branch point: " << position(splitProxyWindow._proxyIter) << std::endl;
#endif
            }
            else
            {
                if (deltaCoverage(traverser._proxyBranchNodeIt)[top(traverser._branchStack)._branchProxyId])
                {
                    _updateAuxiliaryBranchStructures(top(traverser._branchStack), traverser._proxyBranchNodeIt);
                }
            }
            if (deltaType(traverser._proxyBranchNodeIt) == DELTA_TYPE_DEL ||
                deltaType(traverser._proxyBranchNodeIt) == DELTA_TYPE_SV)
                while(nodeItEnd != traverser._proxyBranchNodeIt && deltaPosition(traverser._proxyBranchNodeIt + 1) < top(traverser._branchStack)._mappedHostPos)
                    ++traverser._proxyBranchNodeIt;

            if (traverser._proxyBranchNodeIt != nodeItEnd)
            {
                splitPointPos = _selectNextSplitPoint(top(traverser._branchStack), ++traverser._proxyBranchNodeIt, traverser._traverserCore._branchNodeIt);
            }
            else
            {
                splitPointPos = top(traverser._branchStack)._proxyEndPos;
                ++traverser._proxyBranchNodeIt;
            }

#ifdef DEBUG_DATA_PARALLEL
                std::cerr << "-> split branch split point: " << splitPointPos << " ("<< deltaPosition(traverser._proxyBranchNodeIt)  << ")" << std::endl;
#endif
            continue;
        }
        traverser._isSynchronized = false;
        TSize shiftSize = notify(externalAlg, delegate, traverser, StateTraverseBranch());
        branchPos += shiftSize;
        traverser._branchIt += shiftSize;

#ifdef DEBUG_DATA_PARALLEL
        std::cerr << "--- position: " << position(contextBegin(traverser, StateTraverseBranch())) << std::endl;
#endif
    }
}

// ----------------------------------------------------------------------------
// Function _updateBranchBegin()
// ----------------------------------------------------------------------------

template <typename TProxyId, typename TJstTraverser, typename TPosition, typename TCoverage>
inline typename Size<TJstTraverser>::Type
_selectValidBeginAndProxy(TProxyId & proxyId,
                          TJstTraverser & traverser,
                          TPosition const & contextBeginPosHost,
                          TCoverage const & branchCoverage)
{
    typedef typename TJstTraverser::TBranchNodeIterator TBranchNodeIterator;

    typedef typename TJstTraverser::TMergePointStore TMergePointStack;
    typedef typename TMergePointStack::TMergePoints TMergePoints;
    typedef typename Iterator<TMergePoints>::Type TMergePointIterator;

    typedef typename Size<TJstTraverser>::Type TSize;

    TCoverage tmp;
    transform(tmp, branchCoverage, traverser._traverserCore._activeMasterCoverage, FunctorBitwiseAnd());
    if (!testAllZeros(tmp))  // There exist a valid begin at the current mapped master position.
    {
        proxyId = bitScanForward(tmp);
        return deltaPosition(traverser._traverserCore._branchNodeIt) - contextBeginPosHost;
    }

    // Harder case: We need to parse all possible positions to find the first valid proxy, which is furthest away
    // from the current branch point.

    // Initialization
    TCoverage seenVariants;
    resize(seenVariants, length(stringSet(container(traverser))), false, Exact());

    // The merge point iter is always valid.
    SEQAN_ASSERT_NOT(empty(traverser._traverserCore._mergePointStack._mergePoints));
    TMergePointIterator beginMp = begin(traverser._traverserCore._mergePointStack._mergePoints, Standard());
    TMergePointIterator endMp = end(traverser._traverserCore._mergePointStack._mergePoints, Standard());
    TMergePointIterator itMpLeft = endMp -1;
    // TODO(rmaerker): Check if binary search could be faster?

    // NOTE(rmaerker): The merge points are sorted in decreasing order from left to right.
    for (; itMpLeft != beginMp && itMpLeft->i1 < static_cast<TSize>(contextBeginPosHost); --itMpLeft);
    TMergePointIterator itMp = itMpLeft;
    for (; itMp != beginMp && itMp->i1 < deltaPosition(traverser._traverserCore._branchNodeIt); --itMp);
    ++itMp;  // Either it points to the first no

    TBranchNodeIterator itBp = traverser._traverserCore._branchNodeIt - 1;  // The first node within or before the current context.
    TBranchNodeIterator itBpBegin = traverser._traverserCore._branchNodeIt;

    // As long as the value is greater or equal to the contecxtBeginPosHost we decrement the branch-node pointer.
    for (; static_cast<TPosition>(deltaPosition(itBpBegin)) >= contextBeginPosHost; --itBpBegin)
    {
        if (itBpBegin == begin(container(container(traverser)), Standard()))
        {
            --itBpBegin;
            break;
        }
    }
    ++itBpBegin;  // Points to the first node within context.

    TSize newOffset = 0;
    proxyId = bitScanForward(branchCoverage);

    // Linear scan over merge and branch points to find the branch begin point which is furthest to the left.
    while(true)
    {
        // One of the iterators is at the end.
        if (itBp < itBpBegin || itMp > itMpLeft)
            break;

        // What happens if one is empty and the other not!
        if (itMp->i1 > deltaPosition(itBp))
        {
            transform(tmp, getMergeCoverage(traverser._traverserCore._mergePointStack, itMp - beginMp), seenVariants,
                      FunctorNested<FunctorBitwiseAnd, FunctorIdentity, FunctorBitwiseNot>());
            if (!testAllZeros(tmp))  // Some of this variants are new.
            {
                transform(seenVariants, seenVariants, tmp, FunctorBitwiseOr());
                transform(tmp, tmp, branchCoverage, FunctorBitwiseAnd());
                if (!testAllZeros(tmp))  // There is a coverage containing this position.
                {
                    proxyId = bitScanForward(tmp);
                    newOffset = deltaPosition(traverser._traverserCore._branchNodeIt) - itMp->i1;
                }
            }
            ++itMp;
        }
        else if (deltaPosition(itBp) >= itMp->i1)
        {
            transform(tmp, deltaCoverage(itBp), seenVariants, FunctorNested<FunctorBitwiseAnd, FunctorIdentity,
                      FunctorBitwiseNot>());
            if (!testAllZeros(tmp))  // Some of this variants are new.
            {
                transform(seenVariants, seenVariants, tmp, FunctorBitwiseOr());
                transform(tmp, tmp, branchCoverage, FunctorBitwiseAnd());
                if (!testAllZeros(tmp))  // There is a coverage containing this position.
                {
                    newOffset = (deltaPosition(traverser._traverserCore._branchNodeIt) - deltaPosition(itBp)) - 1;
                    proxyId = bitScanForward(tmp);
                }
            }
            --itBp;
        }
    }

    if (itMp > itMpLeft)
    {
        while (itBp >= itBpBegin)
        {
            transform(tmp, deltaCoverage(itBp), seenVariants, FunctorNested<FunctorBitwiseAnd, FunctorIdentity,
                      FunctorBitwiseNot>());
            if (!testAllZeros(tmp))  // Some of this variants are new.
            {
                transform(seenVariants, seenVariants, tmp, FunctorBitwiseOr());
                transform(tmp, tmp, branchCoverage, FunctorBitwiseAnd());
                if (!testAllZeros(tmp))  // There is a coverage containing this position.
                {
                    newOffset = (deltaPosition(traverser._traverserCore._branchNodeIt) - deltaPosition(itBp)) - 1;
                    proxyId = bitScanForward(tmp);
                }
            }
            --itBp;
        }
    }

    if (itBp < itBpBegin)
    {
        while (itMp <= itMpLeft)
        {
            transform(tmp, getMergeCoverage(traverser._traverserCore._mergePointStack, itMp - beginMp), seenVariants,
                      FunctorNested<FunctorBitwiseAnd, FunctorIdentity, FunctorBitwiseNot>());
            if (!testAllZeros(tmp))  // Some of this variants are new.
            {
                transform(seenVariants, seenVariants, tmp, FunctorBitwiseOr());
                transform(tmp, tmp, branchCoverage, FunctorBitwiseAnd());
                if (!testAllZeros(tmp))  // There is a coverage containing this position.
                {
                    proxyId = bitScanForward(tmp);
                    newOffset = deltaPosition(traverser._traverserCore._branchNodeIt) - itMp->i1;
                }
            }
            ++itMp;
        }
    }
    return newOffset;
}

// ----------------------------------------------------------------------------
// Function _traverseBranchWithAlt()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TState, typename TContextPosition, typename TRequireFullContext,
          typename TExternal, typename TDelegate>
void
_traverseBranchWithAlt(JstTraverser<TContainer, TState, JstTraverserConfig<TContextPosition, TRequireFullContext> > & traverser,
                       TExternal & externalAlg,
                       TDelegate & delegate)
{
    typedef JstTraverser<TContainer, TState, JstTraverserConfig<TContextPosition, TRequireFullContext> > TJstTraverser;
    typedef typename GetStringSet<TContainer>::Type TJournalSet;
    typedef typename Value<TJournalSet>::Type TJournalString;
    typedef typename Iterator<TJournalString, Standard>::Type TJournalStringIterator;
    typedef typename Container<TContainer>::Type TDeltaMap;
    typedef typename DeltaValue<TDeltaMap, DeltaTypeIns>::Type TIns;
    typedef typename DeltaValue<TDeltaMap, DeltaTypeSV>::Type TIndel;

    typedef typename Size<TJstTraverser>::Type TSize;
    typedef typename MakeSigned<TSize>::Type TSignedSize;

    createInitialEntry(traverser._branchStack);

    // We start by initializing the branch coverage and the active branch coverage.
    top(traverser._branchStack)._branchCoverage = deltaCoverage(traverser._traverserCore._branchNodeIt);

    if (TRequireFullContext::VALUE)
        top(traverser._branchStack)._prefixOffset =
            _selectValidBeginAndProxy(top(traverser._branchStack)._branchProxyId, traverser,
                                      _contextBeginPosition(traverser, StateTraverseMaster()),
                                      top(traverser._branchStack)._branchCoverage);
#ifdef DEBUG_DATA_PARALLEL
        std::cerr << "Selected Proxy: " << top(traverser._branchStack)._branchProxyId << std::endl;
#endif

    // The coverage cannot be empty.
    SEQAN_ASSERT_GEQ(top(traverser._branchStack)._branchProxyId, 0u);
    SEQAN_ASSERT_LT(top(traverser._branchStack)._branchProxyId, length(stringSet(container(traverser))));

//    TJournalString & proxySeq = value(stringSet(container(traverser)), top(traverser._branchStack)._branchProxyId);
    top(traverser._branchStack)._proxyIter = begin(value(stringSet(container(traverser)), top(traverser._branchStack)._branchProxyId), Standard());

    _mapHostToVirtual(top(traverser._branchStack)._proxyIter, container(container(traverser)), top(traverser._branchStack)._branchProxyId,
                      deltaPosition(traverser._traverserCore._branchNodeIt));  // Maybe this can be simplified.

    top(traverser._branchStack)._mappedHostPos = deltaPosition(traverser._traverserCore._branchNodeIt) + 1;
    TSize contextSizeRight = contextSize(traverser);
    switch(deltaType(traverser._traverserCore._branchNodeIt))
    {
        case DELTA_TYPE_DEL:
        {
            top(traverser._branchStack)._proxyEndPosDiff = deltaValue(traverser._traverserCore._branchNodeIt, DeltaTypeDel());
            if (top(traverser._branchStack)._proxyEndPosDiff > 1)
                top(traverser._branchStack)._mappedHostPos += top(traverser._branchStack)._proxyEndPosDiff - 1;  // Moves right by the size of the deletion.

            if (top(traverser._branchStack)._prefixOffset == 0)
            {
#ifdef DEBUG_DATA_PARALLEL
                std::cerr << "Points directly into deleted area." << std::endl;
#endif
                return;  // Points directly into deleted area!
            }
            --contextSizeRight;
            break;
        }
        case DELTA_TYPE_INS:
        {
            TIns & insVar = deltaValue(traverser._traverserCore._branchNodeIt, DeltaTypeIns());
            top(traverser._branchStack)._proxyEndPosDiff = -static_cast<int>(length(insVar));
            contextSizeRight += length(insVar);
            break;
        }
        case DELTA_TYPE_SV:
        {
            TIndel & indel = deltaValue(traverser._traverserCore._branchNodeIt, DeltaTypeSV());
            top(traverser._branchStack)._proxyEndPosDiff = indel.i1;
            top(traverser._branchStack)._proxyEndPosDiff -= static_cast<int>(length(indel.i2));
            contextSizeRight += length(indel.i2) - 1;
            if (indel.i1 > 1)
                top(traverser._branchStack)._mappedHostPos += indel.i1 - 1;  // Moves right by the size of the deletion.
            break;
        }
        default:
        {
            break;
        }
    }

    top(traverser._branchStack)._externalState = traverser._lastMasterState;
    top(traverser._branchStack)._proxyEndPos = position(top(traverser._branchStack)._proxyIter) + contextSizeRight;
    
    _traverseBranch(traverser, externalAlg, delegate);
 	pop(traverser._branchStack);
    while (!empty(traverser._branchStack))
    {
        // We first check if we need to update the current branch split.
        if (top(traverser._branchStack)._prefixOffset <= 0)  // The branch starts directly within the variant.
        {
            top(traverser._branchStack)._branchProxyId = bitScanForward(top(traverser._branchStack)._branchCoverage);
        }
        else
        {
            if (TRequireFullContext::VALUE)
            {
                TSignedSize newOffset = _selectValidBeginAndProxy(top(traverser._branchStack)._branchProxyId, traverser,
                                                          _contextBeginPosition(traverser, StateTraverseMaster()),
                                                          top(traverser._branchStack)._branchCoverage);

                if (newOffset < top(traverser._branchStack)._prefixOffset)
                    top(traverser._branchStack)._prefixOffset = newOffset;
            }

            if (top(traverser._branchStack)._prefixOffset == 0 && deltaType(traverser._traverserCore._branchNodeIt) == DELTA_TYPE_DEL)
                continue;
        }

        SEQAN_ASSERT_GEQ(top(traverser._branchStack)._branchProxyId, 0u);
        SEQAN_ASSERT_LT(top(traverser._branchStack)._branchProxyId, length(stringSet(container(traverser))));

        TJournalStringIterator targetIt;
        targetIt._journalStringPtr = &value(stringSet(container(traverser)), top(traverser._branchStack)._branchProxyId);

        // It might be necessary to select from the host position instead of the virtual mapping.
        if (deltaType(traverser._traverserCore._branchNodeIt) == DELTA_TYPE_DEL)
            _mapHostToVirtual(targetIt,
                              container(container(traverser)), top(traverser._branchStack)._branchProxyId,
                              deltaPosition(traverser._traverserCore._branchNodeIt));
        else
            _mapVirtualToVirtual(targetIt, top(traverser._branchStack)._proxyIter, top(traverser._branchStack)._proxyIter,
                                 traverser._traverserCore._branchNodeIt,
                                 container(container(traverser)), top(traverser._branchStack)._branchProxyId);

        top(traverser._branchStack)._proxyIter = targetIt;  // NOTE(rmaerker): We could do a swap here, because we don't need the target iter no more.
        top(traverser._branchStack)._proxyEndPos = position(top(traverser._branchStack)._proxyIter) + contextSizeRight;

        setState(externalAlg, top(traverser._branchStack)._externalState);
        _traverseBranch(traverser, externalAlg, delegate);
        pop(traverser._branchStack);
    }
}

template <typename TContainer, typename TState, typename TExternal, typename TDelegate>
void _traverseBranchWithAlt(JstTraverser<TContainer, TState, JstTraverserConfig<ContextPositionRight, False> > & traverser,
                            TExternal & externalAlg,
                            TDelegate & delegate)
{
    typedef typename GetStringSet<TContainer>::Type TJournalSet;
    typedef typename Value<TJournalSet>::Type TJournalString;
    typedef typename Iterator<TJournalString, Standard>::Type TJournalStringIterator;
    typedef typename Container<TContainer>::Type TDeltaMap;
    typedef typename DeltaValue<TDeltaMap, DeltaTypeIns>::Type TIns;
    typedef typename DeltaValue<TDeltaMap, DeltaTypeSV>::Type TIndel;

    typedef typename Size<TContainer>::Type TSize;

    createInitialEntry(traverser._branchStack);

    top(traverser._branchStack)._branchCoverage = deltaCoverage(traverser._traverserCore._branchNodeIt);
    top(traverser._branchStack)._branchProxyId = bitScanForward(top(traverser._branchStack)._branchCoverage);
    top(traverser._branchStack)._proxyEndPosDiff = 0;

#ifdef DEBUG_DATA_PARALLEL
        std::cerr << "Selected Proxy: " << top(traverser._branchStack)._branchProxyId << std::endl;
#endif

    // Stop looking into the branch because it is invalid.
    if (top(traverser._branchStack)._branchProxyId >= length(stringSet(container(traverser))))
        return;

    // The coverage cannot be empty.
    SEQAN_ASSERT_GEQ(top(traverser._branchStack)._branchProxyId, 0u);
    SEQAN_ASSERT_LT(top(traverser._branchStack)._branchProxyId, length(stringSet(container(traverser))));

//    TJournalString & proxySeq = value(stringSet(container(traverser)), top(traverser._branchStack)._branchProxyId);
    top(traverser._branchStack)._proxyIter = begin(value(stringSet(container(traverser)), top(traverser._branchStack)._branchProxyId), Standard());

    _mapHostToVirtual(top(traverser._branchStack)._proxyIter, container(container(traverser)), top(traverser._branchStack)._branchProxyId,
                      deltaPosition(traverser._traverserCore._branchNodeIt));

    top(traverser._branchStack)._mappedHostPos = deltaPosition(traverser._traverserCore._branchNodeIt) + 1;
    TSize contextSizeRight =  contextSize(traverser);
    switch(deltaType(traverser._traverserCore._branchNodeIt))
    {
        case DELTA_TYPE_DEL:
        {
            top(traverser._branchStack)._proxyEndPosDiff = deltaValue(traverser._traverserCore._branchNodeIt, DeltaTypeDel());
            // We only consider deletions bigger than 1.
            if (top(traverser._branchStack)._proxyEndPosDiff > 1)
                top(traverser._branchStack)._mappedHostPos += top(traverser._branchStack)._proxyEndPosDiff - 1;  // Moves right by the size of the deletion.
            --contextSizeRight;  // We update the right size.
            break;
        }
        case DELTA_TYPE_INS:
        {
            TIns & insVar = deltaValue(traverser._traverserCore._branchNodeIt, DeltaTypeIns());
            top(traverser._branchStack)._proxyEndPosDiff = -static_cast<int>(length(insVar));
            contextSizeRight += length(insVar);
            break;
        }
        case DELTA_TYPE_SV:
        {
            TIndel & indel = deltaValue(traverser._traverserCore._branchNodeIt, DeltaTypeSV());
            top(traverser._branchStack)._proxyEndPosDiff = indel.i1;
            top(traverser._branchStack)._proxyEndPosDiff -= static_cast<int>(length(indel.i2));
            contextSizeRight += length(indel.i2);
            if (indel.i1 > 1)
                top(traverser._branchStack)._mappedHostPos += indel.i1 - 1;  // Moves right by the size of the deletion.
            break;
        }
        default:
            break;
    }

    // Fill the state with the selected proxy until the end.
    top(traverser._branchStack)._externalState = traverser._lastMasterState;
    top(traverser._branchStack)._proxyEndPos  = position(top(traverser._branchStack)._proxyIter) + contextSizeRight;
    top(traverser._branchStack)._prefixOffset = contextSize(traverser) - 1; // We assume single base movement. -> This is a strong assumption maybe not always valid!

    _traverseBranch(traverser, externalAlg, delegate);
    pop(traverser._branchStack);
    while (!empty(traverser._branchStack))
    {
        // Now we come back with the reduced branch coverage.
        SEQAN_ASSERT(!testAllZeros(top(traverser._branchStack)._branchCoverage));  // The coverage cannot be empty.

        top(traverser._branchStack)._branchProxyId = bitScanForward(top(traverser._branchStack)._branchCoverage);
        if (top(traverser._branchStack)._branchProxyId > length(stringSet(container(traverser))))  // No valid proxy can be found.
            continue;

        TJournalStringIterator targetIt;
        targetIt._journalStringPtr = &value(stringSet(container(traverser)), top(traverser._branchStack)._branchProxyId);

        // It might be necessary to select from the host position instead of the virtual mapping.
        if (deltaType(traverser._traverserCore._branchNodeIt) == DELTA_TYPE_DEL)
            _mapHostToVirtual(targetIt,
                              container(container(traverser)), top(traverser._branchStack)._branchProxyId,
                              deltaPosition(traverser._traverserCore._branchNodeIt));
        else
            _mapVirtualToVirtual(targetIt, top(traverser._branchStack)._proxyIter, top(traverser._branchStack)._proxyIter,
                                 traverser._traverserCore._branchNodeIt, container(container(traverser)),
                                 top(traverser._branchStack)._branchProxyId);

        // Now targetIt points to the position of
        top(traverser._branchStack)._proxyIter = targetIt;  // NOTE(rmaerker): We could do a swap here, because we don't need the target iter any more.
        top(traverser._branchStack)._proxyEndPos = position(top(traverser._branchStack)._proxyIter) + contextSizeRight;

        // We need to be careful, because, when the current variant is a deletion this might point directly into the deletion.
        setState(externalAlg, top(traverser._branchStack)._externalState);
        _traverseBranch(traverser, externalAlg, delegate);
        pop(traverser._branchStack);
    }
}

// ----------------------------------------------------------------------------
// Function _syncAndUpdateCoverage()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TState, typename TSpec>
inline void
_syncAndUpdateCoverage(JstTraverser<TContainer, TState, TSpec> const & traverser,
                       StateTraverseMaster const & /*tag*/)
{
    typedef typename Container<TContainer const>::Type TDeltaMap;
    typedef typename Iterator<TDeltaMap, Standard>::Type TBranchNodeIterator;

    typedef typename Position<TContainer const>::Type TPosition;
    typedef typename MakeSigned<TPosition const>::Type TPos_;

    traverser._isSynchronized = true;
    // First update the merge points.
    bool isUpdated = _updateMergePoints(traverser._traverserCore._mergePointStack,
                                        _clippedContextBeginPosition(traverser, StateTraverseMaster()));

    // First check if new context begin position is after the left most node within the context.
    if (!atEnd(traverser._traverserCore._branchNodeInContextIt, container(container(traverser))) &&
        _contextBeginPosition(traverser, StateTraverseMaster()) > static_cast<TPos_>(deltaPosition(traverser._traverserCore._branchNodeInContextIt)))
    {
        // Set correct context branch node iterator.
         while (!atEnd(++traverser._traverserCore._branchNodeInContextIt, container(container(traverser))) &&
                _contextBeginPosition(traverser, StateTraverseMaster()) > static_cast<TPos_>(deltaPosition(traverser._traverserCore._branchNodeInContextIt)))
         {}
         isUpdated = true;
    }
     // The context node must be smaller or equal the current branch node.
     SEQAN_ASSERT(traverser._traverserCore._branchNodeInContextIt <= traverser._traverserCore._branchNodeIt);

    if (isUpdated)
    {
        // Now we update the current master coverage, while considering all deltas within the current context.
        traverser._traverserCore._activeMasterCoverage = ~traverser._traverserCore._mergePointStack._mergeCoverage;
        for(TBranchNodeIterator it = traverser._traverserCore._branchNodeInContextIt; it != traverser._traverserCore._branchNodeIt; ++it)
            transform(traverser._traverserCore._activeMasterCoverage, traverser._traverserCore._activeMasterCoverage, deltaCoverage(it),
                      FunctorNested<FunctorBitwiseAnd, FunctorIdentity, FunctorBitwiseNot>());
    }
}

template <typename TContainer, typename TState, typename TSpec>
inline void
_syncAndUpdateCoverage(JstTraverser<TContainer, TState, TSpec> const & traverser,
                       StateTraverseBranch const & /*tag*/)
{
    typedef JstTraverser<TContainer, TState, TSpec> TJstTraverser;
    typedef typename Container<TContainer const>::Type TDeltaMap;
    typedef typename Iterator<TDeltaMap, Standard>::Type TBranchNodeIterator;

    typedef typename Position<TContainer const>::Type TPos;
    typedef typename MakeSigned<TPos>::Type TPos_;

    typedef typename Size<TJstTraverser>::Type TSize;

    typedef typename TJstTraverser::TMergePointStore TMergePointStack;
    typedef typename TMergePointStack::TMergePoints TMergePoints;
    typedef typename Iterator<TMergePoints const, Standard>::Type TMergePointIterator;

    typedef typename TJstTraverser::TBranchStack const TBranchStack;
    typedef typename Value<TBranchStack>::Type TBranchStackEntry;

    traverser._isSynchronized = true;
    // At the current branch position there might be a deletion and we merge the values.
    // register unsigned oldLength = length(traverser._traverserCore._mergePointStack._mergePoints);

    // We need the correct begin position within the master branch.
    // Now we assume the context can begin in any prefix that has a path that includes the current branch node.
    // Since, we already searched all previous branch nodes that include the current node on their path,
    // we need to find only those who come directly from the master branch.

    // To find only the valid coverage we move left through the previous nodes and merge points that can affect the
    // current prefix.

    // Simple case, the current context begin maps directly onto the current delta.
    TBranchStackEntry & branchStackEntry = top(traverser._branchStack);
    traverser._activeBranchCoverage = branchStackEntry._branchCoverage;
    if (_contextBeginPosition(traverser, StateTraverseBranch()) >= static_cast<TPos_>(position(branchStackEntry._proxyIter)))
        return;  // Nothing to do, we can keep the current coverage.

    // Harder case, we need to determine all sequences that had not been destroyed by the previous operations.
    // But we don't know, from which path the current proxy comes. -> Can be the correct one or not.
    // But we know:
    TSize sDelta = _contextEndPosition(traverser, StateTraverseBranch()) - position(branchStackEntry._proxyIter);
    SEQAN_ASSERT_LEQ(sDelta, contextSize(traverser) -1);
    TSize pDelta = contextSize(traverser) - 1 - sDelta;

    // In some rare cases, where there is a delta direct at beginning, it might be possible that we
    // reach over the begin position of the delta.
    TSize beginHostPos = 0;
    if (pDelta <= deltaPosition(traverser._traverserCore._branchNodeIt))
        beginHostPos = deltaPosition(traverser._traverserCore._branchNodeIt) - pDelta;

    SEQAN_ASSERT_LEQ(beginHostPos, deltaPosition(traverser._traverserCore._branchNodeIt));

    // Considering merge points.
    // Only the intersection of *mp and cov is valid for all: beginHostPos <= *mp <= *branchNodeIt;
    TMergePointIterator it = end(traverser._traverserCore._mergePointStack._mergePoints) - 1;
    TMergePointIterator itBegin = begin(traverser._traverserCore._mergePointStack._mergePoints);

    for(;it != itBegin; --it)  // Moves over all merge points! Only merge points within range are interesting.
    {
        if (it->i1 > beginHostPos && it->i1 <= deltaPosition(traverser._traverserCore._branchNodeIt))
            transform(traverser._activeBranchCoverage, traverser._activeBranchCoverage,
                      getMergeCoverage(traverser._traverserCore._mergePointStack, it - itBegin),
                      FunctorNested<FunctorBitwiseAnd, FunctorIdentity, FunctorBitwiseNot>());
    }

    // Considering previous deltas.
    // For all bp: beginHostPos <= ~*bp < *branchNodeIt
    TBranchNodeIterator tmpIt = traverser._traverserCore._branchNodeIt;
    TBranchNodeIterator tmpItBegin = begin(container(container(traverser)), Standard());
    for (; tmpIt != tmpItBegin && deltaPosition(tmpIt) == deltaPosition(traverser._traverserCore._branchNodeIt); --tmpIt)
    {}  // Only nodes that are not equal to the current.

    for (; tmpIt != tmpItBegin && deltaPosition(tmpIt) >= beginHostPos; --tmpIt)
        transform(traverser._activeBranchCoverage, traverser._activeBranchCoverage, deltaCoverage(tmpIt),
                  FunctorNested<FunctorBitwiseAnd, FunctorIdentity, FunctorBitwiseNot>());

    if (tmpIt == tmpItBegin && deltaPosition(tmpIt) >= beginHostPos  &&
        deltaPosition(tmpIt) < deltaPosition(traverser._traverserCore._branchNodeIt))
        transform(traverser._activeBranchCoverage, traverser._activeBranchCoverage, deltaCoverage(tmpIt),
                  FunctorNested<FunctorBitwiseAnd, FunctorIdentity, FunctorBitwiseNot>());
}

// ----------------------------------------------------------------------------
// Function _execProducerThread()
// ----------------------------------------------------------------------------

template <typename TConcurrentQueue, typename TContainer, typename TState, typename TContextPosition,
          typename TRequireFullContext, typename TExternal, typename TDelegate>
inline void
_execProducerThread(TConcurrentQueue & queue,
                    JstTraverser<TContainer, TState, JstTraverserConfig<TContextPosition, TRequireFullContext> > & traverser,
                    TExternal & externalAlg,
                    TDelegate & delegate,
                    Parallel const & /*tag*/)
{
    typedef JstTraverser<TContainer, TState, JstTraverserConfig<TContextPosition, TRequireFullContext> > TTraverser;
    typedef typename TTraverser::TBranchNodeIterator TBranchNodeIt;

    traverser._lastMasterState = getState(externalAlg);

    // Loop over the branch nodes.
    while (traverser._traverserCore._branchNodeIt != traverser._branchNodeBlockEnd)
    {
        // We add the first node including the iterator initialized to the begin of the current block.
        appendValue(queue, traverser._traverserCore);
        // Before we go to the next branch node.
        if (IsSameType<TContextPosition, ContextPositionLeft>::VALUE)
            setPosition(traverser._traverserCore._masterIt, _max(0, static_cast<int>(deltaPosition(traverser._traverserCore._branchNodeIt)) - static_cast<int>(contextSize(traverser) - 1)));
        else
            setPosition(traverser._traverserCore._masterIt, deltaPosition(traverser._traverserCore._branchNodeIt));

        // Update the current coverage if necessary.
        _syncAndUpdateCoverage(traverser, StateTraverseMaster());
        // Go to next node.
        TBranchNodeIt currNode = traverser._traverserCore._branchNodeIt;
        while(currNode != traverser._branchNodeBlockEnd && deltaPosition(traverser._traverserCore._branchNodeIt) == deltaPosition(currNode))
        {
            _recordMergePointEnds(traverser, currNode);
            transform(traverser._traverserCore._activeMasterCoverage, traverser._traverserCore._activeMasterCoverage, deltaCoverage(currNode),
                      FunctorNested<FunctorBitwiseAnd, FunctorIdentity, FunctorBitwiseNot>());
            ++currNode;
        }
        // TODO(rmaerker): Check if the consumer can do this steps here.
        traverser._traverserCore._branchNodeIt = currNode;
    }

    traverser._traversalState = JST_TRAVERSAL_STATE_MASTER;
    initState(externalAlg);  // Reactivate the last state.

    // Continue with last part until end of current block or sequence.
    while (contextEnd(traverser, StateTraverseMaster()) < traverser._masterItEnd)
    {
        traverser._isSynchronized = false;
        traverser._traverserCore._masterIt += notify(externalAlg, delegate, traverser, StateTraverseMaster());
    }
    // Synchronize master coverage in the end.
    _updateMergePoints(traverser._traverserCore._mergePointStack, position(contextBegin(traverser, StateTraverseMaster())));
    transform(traverser._traverserCore._activeMasterCoverage, traverser._traverserCore._activeMasterCoverage,
              traverser._traverserCore._mergePointStack._mergeCoverage,
              FunctorNested<FunctorBitwiseAnd, FunctorIdentity, FunctorBitwiseNot>());
}

// ----------------------------------------------------------------------------
// Function _internallyExecuteConsumerThread()
// ----------------------------------------------------------------------------

template <typename TJstTraverserState, typename TExternal, typename TDelegate>
inline void
_internallyExecuteConsumerThread(TJstTraverserState & traverser,
                                 TExternal & externalAlg,
                                 TDelegate & delegate)
{
    typedef typename TJstTraverserState::TBitVector TCoverage;
    typedef typename Size<TJstTraverserState>::Type TSize;

    // We need to initialize the state here.
    traverser._traversalState = JST_TRAVERSAL_STATE_MASTER;
    initState(externalAlg);
    // Search along the master strand.
    while (_contextEndPosition(traverser, StateTraverseMaster()) < deltaPosition(traverser._traverserCore._branchNodeIt))
    {
        traverser._isSynchronized = false;
        traverser._traverserCore._masterIt += notify(externalAlg, delegate, traverser, StateTraverseMaster());
    }

    traverser._traversalState = JST_TRAVERSAL_STATE_BRANCH;

    // Processing the current node.
    TSize branchPosition = deltaPosition(traverser._traverserCore._branchNodeIt);

    _syncAndUpdateCoverage(traverser, StateTraverseMaster());
    traverser._lastMasterState = getState(externalAlg);  // Keep the last active caller state.

    // Search all haplotypes with the alternative allel at this position.
    while (traverser._traverserCore._branchNodeIt != traverser._branchNodeBlockEnd &&
           deltaPosition(traverser._traverserCore._branchNodeIt) == branchPosition)
    {
        TCoverage& mappedCov = deltaCoverage(traverser._traverserCore._branchNodeIt); //mappedCoverage(container(container(traverser)), position(traverser._traverserCore._branchNodeIt));
        if (!testAllZeros(mappedCov))
        {
            _traverseBranchWithAlt(traverser, externalAlg, delegate);
            // Remove the coverage from the current delta from the active master coverage.
            transform(traverser._traverserCore._activeMasterCoverage, traverser._traverserCore._activeMasterCoverage, mappedCov,
                      FunctorNested<FunctorBitwiseAnd, FunctorIdentity, FunctorBitwiseNot>());
        }
        // We increase the delta
        ++traverser._traverserCore._branchNodeIt;
    }
}

// ----------------------------------------------------------------------------
// Function _execConsumerThread()
// ----------------------------------------------------------------------------

template <typename TValue, typename TJstTraverserState, typename TExternal, typename TDelegate>
inline void
_execConsumerThread(ConcurrentQueue<TValue> & queue,
                    TJstTraverserState & traverser,
                    TExternal & externalAlg,
                    TDelegate & delegate,
                    Parallel /*tag*/)
{
    while (popFront(traverser._traverserCore, queue))
        _internallyExecuteConsumerThread(traverser, externalAlg, delegate);

    // TODO(rmaerker): Check if this can be removed.
    while (!empty(queue))
    {
        if (tryPopFront(traverser._traverserCore, queue, Parallel()))
            _internallyExecuteConsumerThread(traverser, externalAlg, delegate);
    }
}

// ----------------------------------------------------------------------------
// Function _recordMergePointEnds()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TState, typename TSpec, typename TBranchNode>
inline void
_recordMergePointEnds(JstTraverser<TContainer, TState, TSpec> & traverser, TBranchNode const & branchNodeIt)
{
    DeltaType dType = deltaType(branchNodeIt);
    if (dType == DELTA_TYPE_DEL)
        if (deltaValue(branchNodeIt, DeltaTypeDel()) > 1)
            push(traverser._traverserCore._mergePointStack, deltaPosition(branchNodeIt) +
                 deltaValue(branchNodeIt, DeltaTypeDel()), branchNodeIt);
    if (dType == DELTA_TYPE_SV)
        if (deltaValue(branchNodeIt, DeltaTypeSV()).i1 > 1)
            push(traverser._traverserCore._mergePointStack, deltaPosition(branchNodeIt) +
                 deltaValue(branchNodeIt, DeltaTypeSV()).i1, branchNodeIt);
}

// ----------------------------------------------------------------------------
// Function _execProducerThread()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TState, typename TContextPosition, typename TRequireFullContext,
          typename TExternal, typename TDelegate>
inline void
_execTraversal(JstTraverser<TContainer, TState, JstTraverserConfig<TContextPosition, TRequireFullContext> > & traverser,
               TExternal & externalAlg,
               TDelegate & delegate,
               Serial const & /*tag*/)
{
    typedef typename Container<TContainer>::Type TDeltaMap;
    typedef typename Position<TDeltaMap>::Type TPosition;
    typedef typename DeltaCoverage<TDeltaMap>::Type TBitVector;

#ifdef PROFILE_JST_TRAVERSAL
    String<double> timeTable;
    resize(timeTable, 3, 0.0, Exact());
    double timeAll = sysTime();
    unsigned counter = 0;
    unsigned currentPercentage = 0;
    unsigned fivePercentInterval = ((traverser._branchNodeBlockEnd - traverser._traverserCore._branchNodeIt) * 5) / 100;
    std::cerr << currentPercentage << "% " << std::flush;
#endif //PROFILE_JST_TRAVERSAL

    traverser._lastMasterState = getState(externalAlg);

    // Loop over the branch nodes.
    while (traverser._traverserCore._branchNodeIt != traverser._branchNodeBlockEnd)
    {
#ifdef PROFILE_JST_TRAVERSAL
        double timeMaster = sysTime();
#endif

#ifdef DEBUG_DATA_PARALLEL
        std::cerr << "\n" << "#####################" << std::endl;
        std::cerr << "Search Master Segment: " << position(contextBegin(traverser, StateTraverseMaster())) << " - " << deltaPosition(traverser._traverserCore._branchNodeIt) << std::endl;
        std::cerr << "Breakpoint: " << deltaPosition(traverser._traverserCore._branchNodeIt) << std::endl;
        std::cerr << "Coverage: " << traverser._traverserCore._activeMasterCoverage << std::endl;
#endif

        traverser._traversalState = JST_TRAVERSAL_STATE_MASTER;
        setState(externalAlg, traverser._lastMasterState);  // Reactivate the last state.

        while (_contextEndPosition(traverser, StateTraverseMaster()) <
               deltaPosition(traverser._traverserCore._branchNodeIt))
        {
            traverser._isSynchronized = false;
            traverser._traverserCore._masterIt += notify(externalAlg, delegate, traverser, StateTraverseMaster());
        }

#ifdef PROFILE_JST_TRAVERSAL
        timeTable[0] += sysTime() - timeMaster;
#endif

        traverser._traversalState = JST_TRAVERSAL_STATE_BRANCH;

        // Processing the current node.
        TPosition branchPosition = deltaPosition(traverser._traverserCore._branchNodeIt);

#ifdef DEBUG_DATA_PARALLEL
        std::cerr << "#####################" << std::endl;
        std::cerr << "Search Branch Segment: " << std::endl;
        std::cerr << "Master Branch Coverage: " << traverser._traverserCore._activeMasterCoverage << std::endl;
#endif

#ifdef PROFILE_JST_TRAVERSAL
        double timeBranchAll = sysTime();
#endif
        _syncAndUpdateCoverage(traverser, StateTraverseMaster());
        traverser._lastMasterState = getState(externalAlg);  // Keep the last active caller state.

        // Search all haplotypes with the alternative allel at this position.
        while (traverser._traverserCore._branchNodeIt != traverser._branchNodeBlockEnd &&
               deltaPosition(traverser._traverserCore._branchNodeIt) == branchPosition)
        {
#ifdef DEBUG_DATA_PARALLEL
            std::cerr << "Coverage: " << traverser._activeBranchCoverage << std::endl;
#endif
            TBitVector& mappedCov = deltaCoverage(traverser._traverserCore._branchNodeIt); //mappedCoverage(container(container(traverser)), position(traverser._traverserCore._branchNodeIt));
            if (!testAllZeros(mappedCov))  // TODO(rmaerker): Should be removed.
            {
#ifdef PROFILE_JST_TRAVERSAL
                double timeBranch1 = sysTime();
#endif
                std::cerr << "[LOG] coverage: " << mappedCov << std::endl;
                switch(deltaType(traverser._traverserCore._branchNodeIt))
                {
                case DELTA_TYPE_SNP: std::cerr << "[LOG] deltaValue: " << deltaValue(traverser._traverserCore._branchNodeIt, DeltaTypeSnp()) << std::endl; break;
                case DELTA_TYPE_DEL: std::cerr << "[LOG] deltaValue: " << deltaValue(traverser._traverserCore._branchNodeIt, DeltaTypeDel()) << std::endl; break;
                case DELTA_TYPE_INS: std::cerr << "[LOG] deltaValue: " << deltaValue(traverser._traverserCore._branchNodeIt, DeltaTypeIns()) << std::endl; break;
                case DELTA_TYPE_SV: std::cerr << "[LOG] deltaValue: " << deltaValue(traverser._traverserCore._branchNodeIt, DeltaTypeSV()) << std::endl; break;
                }
                std::cerr << "[LOG] deltaPosition: " << deltaPosition(traverser._traverserCore._branchNodeIt) << std::endl;
                
                _recordMergePointEnds(traverser, traverser._traverserCore._branchNodeIt);
                _traverseBranchWithAlt(traverser, externalAlg, delegate);
#ifdef PROFILE_JST_TRAVERSAL
                timeTable[1] += sysTime() - timeBranch1;
#endif
                // Remove the coverage from the current delta from the active master coverage.
                transform(traverser._traverserCore._activeMasterCoverage, traverser._traverserCore._activeMasterCoverage, mappedCov,
                          FunctorNested<FunctorBitwiseAnd, FunctorIdentity, FunctorBitwiseNot>());
            }
            // We increase the delta
            ++traverser._traverserCore._branchNodeIt;

#ifdef PROFILE_JST_TRAVERSAL
            if (++counter == fivePercentInterval)
            {
                currentPercentage += 5;
                std::cerr << currentPercentage << "% " << std::flush;
                counter = 0;
            }
#endif //PROFILE_DATA_PARALLEL
        }
#ifdef PROFILE_JST_TRAVERSAL
        timeTable[2] += sysTime() - timeBranchAll;
#endif
    }
#ifdef PROFILE_JST_TRAVERSAL
    double timeMaster = sysTime();
#endif
    traverser._traversalState = JST_TRAVERSAL_STATE_MASTER;
    setState(externalAlg, traverser._lastMasterState);  // Reactivate the last state.

    // Set end of segment to next breakpoint.
#ifdef DEBUG_DATA_PARALLEL
    std::cerr << "#####################" << std::endl;
    std::cerr << "Search Master Segment: " << position(contextBegin(traverser, StateTraverseMaster())) << " - " << position(traverser._masterItEnd) << std::endl;
#endif

    while (contextEnd(traverser, StateTraverseMaster()) < traverser._masterItEnd)
    {
        traverser._isSynchronized = false;
        traverser._traverserCore._masterIt += notify(externalAlg, delegate, traverser, StateTraverseMaster());
#ifdef DEBUG_DATA_PARALLEL
    std::cerr << "--- position: " << position(contextBegin(traverser, StateTraverseMaster())) << std::endl;
#endif
    }
    // Synchronize master coverage in the end.
    _updateMergePoints(traverser._traverserCore._mergePointStack, position(contextBegin(traverser, StateTraverseMaster())));
    transform(traverser._traverserCore._activeMasterCoverage, traverser._traverserCore._activeMasterCoverage,
              traverser._traverserCore._mergePointStack._mergeCoverage,
              FunctorNested<FunctorBitwiseAnd, FunctorIdentity, FunctorBitwiseNot>());

#ifdef PROFILE_JST_TRAVERSAL
    timeTable[0] += sysTime() - timeMaster;
#endif

#ifdef PROFILE_JST_TRAVERSAL
    std::cerr <<  std::endl;
    std::cerr << "Time Master: " << timeTable[0] << " s." << std::endl;
    std::cerr << "Time Branch iterate: " << timeTable[1] << " s." << std::endl;
    std::cerr << "Time Branch all: " << timeTable[2] << " s." << std::endl;
    std::cerr << "Time total: " << sysTime() - timeAll << " s." <<std::endl;
#endif // PROFILE_DATA_PARALLEL
}

// ----------------------------------------------------------------------------
// Function _execTraversal()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TState, typename TContextPosition, typename TRequireFullContext,
          typename TExternal, typename TDelegate>
inline void
_execTraversal(JstTraverser<TContainer, TState, JstTraverserConfig<TContextPosition, TRequireFullContext> > & traverser,
               TExternal externalAlg,
               TDelegate & delegate,
               Parallel const & parallelTag)
{
    typedef JstTraverser<TContainer, TState, JstTraverserConfig<TContextPosition, TRequireFullContext> > TTraverserState;
    typedef typename TTraverserState::TJstTraverserCore TJstTraverserCore;
    typedef ConcurrentQueue<TJstTraverserCore> TQueue;

    // Use for each thread an own traverser state.
    String<TTraverserState> jobs;
    resize(jobs, omp_get_max_threads(), traverser, Exact());

#ifdef PROFILE_JST_TRAVERSAL
    String<double> timeTable;
    resize(timeTable, 3, 0.0, Exact());
    double timeAll = sysTime();
    unsigned counter = 0;
    unsigned currentPercentage = 0;
    unsigned fivePercentInterval = ((traverser._branchNodeBlockEnd - traverser._traverserCore._branchNodeIt) * 5) / 100;
    std::cerr << currentPercentage << "% " << std::flush;
#endif //PROFILE_DATA_PARALLEL

    TQueue queue(0u);  // Concurrently scheduling the jobs.

    // Parallelize with SPMC-model.
    // Everyone works on its own external Algorithm.
    SEQAN_OMP_PRAGMA(parallel firstprivate(externalAlg))
    {
        // Call the function for the producer
        SEQAN_OMP_PRAGMA(master)
        {
            ScopedWriteLock<TQueue> writeLock(queue);
            waitForWriters(queue, 1);  // Barrier for writers until all are registered to the queue.

            _execProducerThread(queue, jobs[omp_get_thread_num()], externalAlg, delegate, parallelTag);
            traverser = jobs[0];
        }

//        SEQAN_OMP_PRAGMA(critical(cout))
//        {
//            printf("Thread: %i registered for popping.\n", omp_get_thread_num());
//        }
        ScopedReadLock<TQueue> readLock(queue);
        waitForFirstValue(queue); // Barrier to wait for all writers to set up.

        _execConsumerThread(queue, jobs[omp_get_thread_num()], externalAlg, delegate, parallelTag);
    }

    SEQAN_ASSERT(empty(queue));
}

// ----------------------------------------------------------------------------
// Function _initSegment()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TState, typename TContextPosition, typename TRequireFullContext,
          typename TDeltaNodeIterator, typename THostSegmentBegin, typename THostSegmentEnd>
inline void
_initSegment(JstTraverser<TContainer, TState, JstTraverserConfig<TContextPosition, TRequireFullContext> > & traverser,
             TDeltaNodeIterator const & nodeItBegin,
             TDeltaNodeIterator const & nodeItEnd,
             THostSegmentBegin const & hostSegmentBeginPosition,
             THostSegmentEnd const & hostSegmentEndPosition)
{
    traverser._traverserCore._masterIt = begin(host(container(traverser)), Rooted()) + hostSegmentBeginPosition;
    if (IsSameType<TContextPosition, ContextPositionRight>::VALUE && TRequireFullContext::VALUE)
        traverser._traverserCore._masterIt +=  (contextSize(traverser) - 1);
    traverser._masterItEnd = begin(host(container(traverser)), Rooted()) + hostSegmentEndPosition;
    traverser._traverserCore._branchNodeInContextIt = traverser._traverserCore._branchNodeIt = nodeItBegin;
    traverser._branchNodeBlockEnd = nodeItEnd;
    _globalInit(traverser);
}

// ----------------------------------------------------------------------------
// Function _reinitBlockEnd()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TState, typename TContextPosition, typename TRequireFullContext>
inline void
_reinitBlockEnd(JstTraverser<TContainer, TState, JstTraverserConfig<TContextPosition, TRequireFullContext> > & traverser)
{
    // We do not need to update the end if the full tree is journaled.
    if (fullJournalRequired(container(traverser)))
        return;

    traverser._branchNodeBlockEnd = container(traverser)._mapBlockEnd;
    if (traverser._branchNodeBlockEnd == end(container(container(traverser)), Standard()))  // Last block.
        traverser._masterItEnd = end(host(container(traverser)), Rooted());
    else
        traverser._masterItEnd = begin(host(container(traverser)), Rooted()) + deltaPosition(traverser._branchNodeBlockEnd);
}

// ----------------------------------------------------------------------------
// Function init()
// ----------------------------------------------------------------------------

/*!
 * @fn JstTraverser#init
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Initializes the @link JstTraverser @endlink.
 *
 * @signature init(traverser);
 * @signature init(traverser, cont);
 * @signature init(traverser, cont, w);
 * @param[in,out]  traverser   The traverser to be initialized.
 * @param[in]      cont        The container to be set.
 * @param[in]      w           The context size to be set.
 *
 * Before using the init function one has to make sure, that the context size and the container of the
 * @link JstTraverser @endlink instance are set properly. Otherwise the behavior is undefined.
 */

template <typename TContainer, typename TState, typename TContextPosition, typename TRequireFullContext>
inline void
init(JstTraverser<TContainer, TState, JstTraverserConfig<TContextPosition, TRequireFullContext> > & traverser)
{
    _initSegment(traverser, begin(container(container(traverser)), Rooted()),
                 end(container(container(traverser)), Rooted()), 0, length(host(container(traverser))));
}

template <typename TContainer, typename TState, typename TContextPosition, typename TRequireFullContext>
inline void
init(JstTraverser<TContainer, TState, JstTraverserConfig<TContextPosition, TRequireFullContext> > & traverser,
     TContainer & obj)
{
    setContainer(traverser, obj);
    init(traverser);
}

template <typename TContainer, typename TState, typename TContextPosition, typename TRequireFullContext, typename TSize>
inline void
init(JstTraverser<TContainer, TState, JstTraverserConfig<TContextPosition, TRequireFullContext> > & traverser,
     TContainer & obj,
     TSize const & contextSize)
{
    setContainer(traverser, obj);
    setContextSize(traverser, contextSize);
    init(traverser);
}

// ----------------------------------------------------------------------------
// Function traverse()
// ----------------------------------------------------------------------------

/*!
 * @fn JstTraverser#traverse
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Triggers the traversal.
 *
 * @signature traverse(ext, delegate, traverser[, tag]);
 *
 * @param[in]      ext       An external algorithm. Has to implement the @link JstTraversalConcept @endlink.
 * @param[in,out]  delegate  A functor which is called by the external algorithm.
 * @param[in,out]  traverser The traverser that manages the traverser. Has to be of type @link JstTraverser @endlink.
 * @param[in]      tag       Tag to enable parallel traversal. One of @link ParallelismTags @endlink.
 */

template <typename TOperator, typename TDelegate, typename TContainer, typename TState, typename TSpec,
          typename TParallelSpec>
inline
SEQAN_FUNC_ENABLE_IF(Is<JstTraversalConcept<TOperator> >, void)
traverse(TOperator & traversalCaller,
         TDelegate & delegate,
         JstTraverser<TContainer, TState, TSpec> & traverser,
         Tag<TParallelSpec> const & tag)
{
    create(container(traverser), contextSize(traverser), tag);

    _reinitBlockEnd(traverser);
    _execTraversal(traverser, traversalCaller, delegate, tag);

    while(createNext(container(traverser), contextSize(traverser), tag))
    {
        _reinitBlockEnd(traverser);
        _execTraversal(traverser, traversalCaller, delegate, tag);
    }
}

template <typename TOperator, typename TDelegate, typename TContainer, typename TState, typename TSpec>
inline
SEQAN_FUNC_ENABLE_IF(Is<JstTraversalConcept<TOperator> >, void)
traverse(TOperator & traversalCaller,
         TDelegate & delegate,
         JstTraverser<TContainer, TState, TSpec> & traverser)
{
    traverse(traversalCaller, delegate, traverser, Serial());
}

// ----------------------------------------------------------------------------
// Function setContextSize()
// ----------------------------------------------------------------------------

/*!
 * @fn JstTraverser#setContextSize
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Sets the context size of the context.
 *
 * @signature TContainer setContextSize(traverser, w);
 * @param[in,out]  traverser The traverser to set the context size to.
 * @param[in]       w         The context size.
 *
 * @see JstTraverser#contextSize
 */

template <typename TContainer, typename TState, typename TSpec, typename TSize>
inline void
setContextSize(JstTraverser<TContainer, TState, TSpec> & traverser,
               TSize const & newWindowSize)
{
    traverser._contextSize = newWindowSize;
    
}

// ----------------------------------------------------------------------------
// Function contextSize()
// ----------------------------------------------------------------------------

/*!
 * @fn JstTraverser#contextSize
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Returns the context size of the context.
 *
 * @signature TContainer contextSize(traverser);
 * @param[in] traverser The taverser to query the contextSize for.
 *
 * @return TSize The context size of type @link JstTraverser#Size @endlink.
 *
 * @see JstTraverser#setContextSize
 */

template <typename TContainer, typename TState, typename TSpec>
inline typename Size<JstTraverser<TContainer, TState, TSpec> >::Type
contextSize(JstTraverser<TContainer, TState, TSpec> const & traverser)
{
    return traverser._contextSize;
}

// ----------------------------------------------------------------------------
// Function setContainer()
// ----------------------------------------------------------------------------

/*!
 * @fn JstTraverser#setContainer
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Sets the container.
 *
 * @signature TContainer setContainer(traverser, container);
 * @param[in,out]  traverser The traverser to set the container to.
 * @param[in]       container The container to be set.
 *
 * @see JstTraverser#container
 */

template <typename TContainer, typename TState, typename TSpec>
inline void
setContainer(JstTraverser<TContainer, TState, TSpec> & traverser,
             TContainer & container)
{
    traverser._haystackPtr = &container;
}

// ----------------------------------------------------------------------------
// Function container()
// ----------------------------------------------------------------------------

/*!
 * @fn JstTraverser#container
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Returns a reference to the underlying container.
 *
 * @signature TContainer container(traverser);
 * @param[in] traverser The traverser to query the container for.
 *
 * @return TContainer A reference to the container of type @link JstTraverser#Container @endlink the traverser is operating on.
 *
 * @see JstTraverser#setContainer
 */

template <typename TContainer, typename TState, typename TSpec>
inline typename Container<JstTraverser<TContainer, TState, TSpec> >::Type &
container(JstTraverser<TContainer, TState, TSpec> & traverser)
{
    return *traverser._haystackPtr;
}

template <typename TContainer, typename TState, typename TSpec>
inline typename Container<JstTraverser<TContainer, TState, TSpec> const>::Type &
container(JstTraverser<TContainer, TState, TSpec> const & traverser)
{
    return *traverser._haystackPtr;
}

// ----------------------------------------------------------------------------
// Function getObjId()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TState, typename TSpec>
inline void const *
getObjectId(JstTraverser<TContainer, TState, TSpec> const & obj)
{
    return static_cast<void const *>(&obj);
}

#ifdef DEBUG_DATA_PARALLEL
template <typename TContainer, typename TState, typename TContextPos, typename TFullContext>
void
_printContext(JstTraverser<TContainer, TState, JstTraverserConfig<TContextPos, TFullContext> > & traverser)
{
    typedef JstTraverser<TContainer, TState, JstTraverserConfig<TContextPos, TFullContext> > TTraverser;
    typedef typename TTraverser::TMasterBranchIterator THostIterator;
    typedef typename TTraverser::TJournalIterator TJournalIterator;

    THostIterator it;
    THostIterator itEnd;

    TJournalIterator itJ;
    TJournalIterator itJEnd;

    if (IsSameType<TContextPos, ContextPositionLeft>::VALUE)
    {
        if (state(traverser) == JST_TRAVERSAL_STATE_MASTER)
        {
            it = contextIterator(traverser, StateTraverseMaster());
            itEnd = it + contextSize(traverser);
        }
        else
        {
            itJ = contextIterator(traverser, StateTraverseBranch());
            itJEnd = itJ + contextSize(traverser);
        }
    }
    else
    {
        if (state(traverser) == JST_TRAVERSAL_STATE_MASTER)
        {
            itEnd = contextIterator(traverser, StateTraverseMaster()) + 1;
            it = itEnd - contextSize(traverser);
        }
        else
        {
            itJEnd = contextIterator(traverser, StateTraverseBranch()) + 1;
            itJ = itJEnd - contextSize(traverser);
        }
    }
    if (state(traverser) == JST_TRAVERSAL_STATE_MASTER)
    {
        std::cerr << "Context-M: ";
        for (; it != itEnd; ++it)
            std::cerr << *it;
        std::cerr << std::endl;
    }
    else
    {
        std::cerr << "Context-B: ";
        for (; itJ != itJEnd; ++itJ)
            std::cerr << *itJ;
        std::cerr << std::endl;
    }
}

#endif

}  // namespace seqan

#endif  // EXTRAS_INCLUDE_SEQAN_DATA_PARALLEL_DATA_PARALLEL_TRAVERSAL_H_
