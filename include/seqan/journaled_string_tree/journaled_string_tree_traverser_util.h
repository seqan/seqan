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

#ifndef INCLUDE_SEQAN_JOURNALED_STRING_TREE_JOURNALED_STRING_TREE_TRAVERSER_UTIL_H_
#define INCLUDE_SEQAN_JOURNALED_STRING_TREE_JOURNALED_STRING_TREE_TRAVERSER_UTIL_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct TraverserStackMember_;
typedef Tag<TraverserStackMember_> TraverserStackMember;

template <typename TSpec = void>
struct JstTraversalSpec
{};

template <typename TContainer, typename TSpec, typename TObserver = void>
class TraverserImpl;

// ============================================================================
// Metafunctions
// ============================================================================

template <typename TObject>
struct StringContext;

// ============================================================================
// Private Functions
// ============================================================================

namespace impl
{

// ----------------------------------------------------------------------------
// Function impl::mapSourceToVirtual()
// ----------------------------------------------------------------------------

template <typename TIterator, typename TDeltaMap, typename TProxyId, typename THostPos>
inline void
mapBranchPointToVirtual(TIterator & resultIt,
                        TDeltaMap const & variantStore,
                        TProxyId const & proxyId,
                        THostPos const & hostPos)
{
    typedef typename Container<TIterator>::Type TJournalString;
    typedef typename JournalType<TJournalString>::Type TJournalEntries;
    typedef typename Value<TJournalEntries>::Type TCargo;
    typedef typename Iterator<TJournalEntries>::Type TEntriesIterator;
    typedef typename Position<TCargo>::Type TCargoPos;
    typedef typename Size<TCargo>::Type TCargoSize;
    typedef typename Value<TDeltaMap>::Type TMapEntry;
    typedef JournalEntryLtByPhysicalOriginPos<TCargoPos, TCargoSize> TComp;

    typedef typename Iterator<TDeltaMap const, Standard>::Type TVarIterator;
    typedef typename Position<TDeltaMap const>::Type TDeltaMapPos;

    // We need to set the iterator to the correct position within the proxy sequence given the host pos.
    TJournalEntries & journalEntries = _journalEntries(*resultIt._journalStringPtr);

    TCargo refCargo;
    refCargo.physicalOriginPosition = hostPos;
    TEntriesIterator it = std::lower_bound(begin(journalEntries._journalNodes, Standard()),
                                           end(journalEntries._journalNodes, Standard()), refCargo, TComp());

    // This is now the first position whose var is equal or greater to the host pos.
    // Since this is either a position that is deleted
    // or a position after the insertion made -> Even for a SNP
    // We have to go backwards.
    if (it != begin(journalEntries, Standard()))
        --it;

    while (it != begin(journalEntries, Standard()) && it->segmentSource == SOURCE_PATCH)
        --it;

    if (it->segmentSource == SOURCE_PATCH)  // The iterator has to be at the beginning.
    {
        TVarIterator itVar = begin(variantStore, Standard());  // TODO(rrahn): Optimize!
        SEQAN_ASSERT_LEQ(getDeltaPosition(*itVar), static_cast<TDeltaMapPos const>(hostPos));

        TDeltaMapPos virtualOffset = 0;
        // Now we move to the right until we find the node that we are looking for and reconstruct the offset of the virtual positions.
        while(getDeltaPosition(*itVar) != static_cast<TDeltaMapPos const>(hostPos) && !atEnd(itVar, variantStore))
        {
            if (getDeltaCoverage(*itVar)[proxyId] != true)  // irrelevant variant.
            {
                ++itVar;
                continue;
            }

            if (getDeltaType(*itVar) == DELTA_TYPE_INS)
                virtualOffset += length(deltaValue(itVar, DeltaTypeIns()));
            else if (getDeltaType(*itVar) == DELTA_TYPE_SNP)
                ++virtualOffset;
            else if (getDeltaType(*itVar) == DELTA_TYPE_SV)
                virtualOffset += length(deltaValue(itVar, DeltaTypeSV()).i2);
            ++itVar;
        }
        resultIt += virtualOffset;  // Set the iterator to the beginning of the variant.
        return;
    }

    SEQAN_ASSERT_EQ(it->segmentSource, SOURCE_ORIGINAL);

    // We assume that the operation begins here!
    resultIt._journalEntriesIterator = it;
    if (it->physicalOriginPosition + it->length > static_cast<TDeltaMapPos const>(hostPos))
    {
        _updateSegmentIterators(resultIt);
        if (it->physicalOriginPosition < hostPos)
            resultIt += hostPos - it->physicalOriginPosition;
        return;
    }

    _updateSegmentIteratorsLeft(resultIt);  // Set the iterator to the end of the current original node.
    if (_physicalPosition(resultIt) + 1 == static_cast<TDeltaMapPos const>(hostPos))
    {
        ++resultIt;
        return;
    }

    // TODO(rmaerker): Can remove the binary Search here!
    // Find the first node that is left or equal to the current physical position!
    TMapEntry child;
    child.deltaPosition = _physicalPosition(resultIt);
    TVarIterator itVar = std::upper_bound(begin(variantStore, Standard()), end(variantStore, Standard()), child,
                                          DeltaMapEntryPosLessThanComparator_());

    SEQAN_ASSERT_LEQ(getDeltaPosition(*itVar), static_cast<TDeltaMapPos const>(hostPos));

    TDeltaMapPos virtualOffset = 0;
    // Now we move to the right until we find the node that we are looking for and reconstruct the offset of the virtual positions.
    while (getDeltaPosition(*itVar) != static_cast<TDeltaMapPos const>(hostPos) && !atEnd(itVar, variantStore))
    {
        if (getDeltaCoverage(*itVar)[proxyId] != true)  // irrelevant variant.
        {
            ++itVar;
            continue;
        }

        if (getDeltaType(*itVar) == DELTA_TYPE_INS)
            virtualOffset += length(deltaValue(itVar, DeltaTypeIns()));
        else if (getDeltaType(*itVar) == DELTA_TYPE_SNP)
            ++virtualOffset;
        else if (getDeltaType(*itVar) == DELTA_TYPE_SV)
            virtualOffset += length(deltaValue(itVar, DeltaTypeSV()).i2);
        ++itVar;
    }
    resultIt += virtualOffset + 1;  // Set the iterator to the beginning of the variant.
}


// ----------------------------------------------------------------------------
// Function impl::inBranch()
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec, typename TObserver>
inline bool
inBranch(TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> const & me)
{
    return length(*resource(me._rm)) > 1;
}

// ----------------------------------------------------------------------------
// Function impl::isMergePoint()
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec, typename TObserver>
inline bool
isMergePoint(TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> const & me)
{
    return !impl::inBranch(me) && back(*resource(me._rm)).mappedSrcEndPos != deltaPosition(back(*resource(me._rm)).curDelta);
}

// ----------------------------------------------------------------------------
// Function impl::nextChild()
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec, typename TObserver,
          typename TTraversalNode>
inline bool
nextChild(TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> & it,
          TTraversalNode & parent,
          TTraversalNode & child)
{
    // We set the coverage of the left child to be the one of the parent & coverage(curDelta);
    transform(child.coverage, parent.coverage, getDeltaCoverage(*impl::hostIter(parent.nextDelta)), FunctorBitwiseAnd());
    transform(parent.coverage, parent.coverage, getDeltaCoverage(*impl::hostIter(parent.nextDelta)),
              FunctorNested<FunctorBitwiseAnd, FunctorIdentity, FunctorBitwiseNot>());

    // Set the current delta of the child to the next delta of the parent.
    child.curDelta = parent.nextDelta;

    // A) First get the proxyId of the updated child coverage.
    auto proxyId = bitScanForward(child.coverage);
    if (proxyId >= dimension(container(it)))
        return false;  // We can skip this child, since it has an empty coverage.

    // B) Remap the sequence iterators according to the new positions.
    if (&container(it)._buffer._journaledSet[proxyId] == &container(parent.begEdgeIt))
    {
        // remap parent node to new journal sequence.
        proxyId = bitScanForward(parent.coverage);
        if (proxyId < dimension(container(it)) && !parent.fromBase)
        {
            auto endToBegDist = parent.endEdgeIt - parent.begEdgeIt;
            impl::mapBranchPointToVirtual(parent.endEdgeIt, host(container(it)), proxyId,
                                          getDeltaPosition(*impl::hostIter(parent.nextDelta)));
            parent.begEdgeIt = parent.endEdgeIt - endToBegDist;
            parent.curEdgeIt = parent.endEdgeIt;
            // We have the sequence iterators updated and the coverage.
            // So far the remainingLength and the curNode remain the same.
        }
    }
    else
    {
        // map child node to new journal sequence.
        // remap child node to new journal sequence.
        auto endToBegDist = child.endEdgeIt - child.begEdgeIt;
        impl::mapBranchPointToVirtual(child.endEdgeIt, host(container(it)), proxyId,
                                      getDeltaPosition(*impl::hostIter(child.curDelta)));
        child.begEdgeIt = child.endEdgeIt - endToBegDist;
        child.curEdgeIt = child.endEdgeIt;
    }

    // Move to next valid delta position, note that we skip all end points and all deltas that occur at the same position.
    while (getDeltaPosition(*(impl::hostIter(++child.nextDelta))) == getDeltaPosition(*impl::hostIter(child.curDelta))
           || impl::isEndPoint(child.nextDelta))
    {}

    // C) Update remaining length if we come directly from the base sequence.
    if (parent.isBase)  // TODO(rmaerker): Make sure that the base node always has the complete window size as remaining Size.
    {
        // We should keep the mappedSrcEndPos of the base node already in sync.
        child.mappedSrcEndPos += historySize(container(it));
        switch (getDeltaType(*impl::hostIter(child.curDelta)))
        {
            case DELTA_TYPE_DEL:
            {
                child.mappedSrcEndPos += deletionSize(host(container(it))._deltaStore,
                                                      getStorePosition(*impl::hostIter(child.curDelta)), DeltaTypeDel());
                --child.remainingSize;
                break;
            }
            case DELTA_TYPE_INS:
            {
                child.remainingSize += insertionSize(host(container(it))._deltaStore,
                                                     getStorePosition(*impl::hostIter(child.curDelta)), DeltaTypeIns());
                break;
            }
            case DELTA_TYPE_SV:
            {
                child.mappedSrcEndPos += deletionSize(host(container(it))._deltaStore,
                                                      getStorePosition(*impl::hostIter(child.curDelta)), DeltaTypeSV());
                child.remainingSize += insertionSize(host(container(it))._deltaStore,
                                                     getStorePosition(*impl::hostIter(child.curDelta)), DeltaTypeSV()) - 1;
                break;
            }
            default: break;
        }
    }
    else
    {
        // Now we are down somewhere in an internal subtree.
        switch (getDeltaType(*impl::hostIter(child.curDelta)))
        {
            case DELTA_TYPE_DEL:
            {
                child.mappedSrcEndPos += deletionSize(host(container(it))._deltaStore,
                                                      getStorePosition(*impl::hostIter(child.curDelta)), DeltaTypeDel());
                break;
            }
            case DELTA_TYPE_SV:
            {
                child.mappedSrcEndPos += deletionSize(host(container(it))._deltaStore,
                                                      getStorePosition(*impl::hostIter(child.curDelta)), DeltaTypeSV());
                break;
            }
            default: break;
        }
    }

    child.begEdgeIt = child.endEdgeIt;
    child.curEdgeIt = child.begEdgeIt;
    child.endEdgeIt = child.begEdgeIt + (getDeltaPosition(*impl::hostIter(child.nextDelta)) -
                                         getDeltaPosition(*impl::hostIter(child.curDelta)));
    child.isBase = false;
    child.fromBase = false;
    return true;
}


// ----------------------------------------------------------------------------
// Function impl::advanceBaseParent()
// ----------------------------------------------------------------------------

template <typename TTraversalNode>
inline void
advanceBaseParent(TTraversalNode & base)
{
    SEQAN_ASSERT(base.curDelta == base.nextDelta);

#if defined(DEBUG_JST_TRAVERSAL)
    std::cout << "       BASE: " << base << std::endl;
#endif // defined(DEBUG_JST_TRAVERSAL)
    // TODO(rrahn): Check if we need this?
    arrayFill(begin(base.coverage, Standard()), end(base.coverage, Standard()), true);  // Make sure the coverage is set to 1.
    while (getDeltaPosition(*impl::hostIter(base.curDelta)) == getDeltaPosition(*impl::hostIter(base.nextDelta)))
    {
        if (SEQAN_UNLIKELY(getDeltaType(*impl::hostIter(base.curDelta)) == DELTA_TYPE_DEL ||
                           getDeltaType(*impl::hostIter(base.curDelta)) == DELTA_TYPE_SV))
            transform(base.coverage, base.coverage, getDeltaCoverage(*impl::hostIter(base.curDelta)),
                      FunctorNested<FunctorBitwiseAnd, FunctorIdentity, FunctorBitwiseNot>());
        ++base.nextDelta;
    }
    // base.nextDelta might point to an merge point.
    base.begEdgeIt += getDeltaPosition(*impl::hostIter(base.curDelta)) - position(base.begEdgeIt);
    base.curEdgeIt = base.begEdgeIt + 1;
    base.endEdgeIt += getDeltaPosition(*impl::hostIter(base.nextDelta)) - base.mappedSrcEndPos;
    base.mappedSrcEndPos = getDeltaPosition(*impl::hostIter(base.nextDelta));
}

// ----------------------------------------------------------------------------
// Function impl::updateParent()
// ----------------------------------------------------------------------------

template <typename TTraversalNode>
inline bool
updateParent(TTraversalNode & parent)
{
    SEQAN_ASSERT_NOT(parent.isBase);

    if (testAllZeros(parent.coverage))
        return false;
    // Update:
    parent.begEdgeIt = parent.endEdgeIt;
    parent.curEdgeIt = parent.begEdgeIt;

    // Make sure that for a right subtree we don't consider merge points, because we cannot be within a branch
    // coming from a delta while at the same time there is a deletion to be covered.
    while(impl::isEndPoint(++parent.nextDelta))  // Move the nextDelta iterator to the next branch point.
    {}
    parent.mappedSrcEndPos = getDeltaPosition(*impl::hostIter(parent.nextDelta));
    parent.endEdgeIt += parent.mappedSrcEndPos - getDeltaPosition(*impl::hostIter(parent.curDelta));
    return true;
}


// ----------------------------------------------------------------------------
// Function impl::moveWindow()
// ----------------------------------------------------------------------------

template <typename TNode, typename TSize>
inline TSize
moveWindow(TNode & node, TSize stepSize)
{
    SEQAN_ASSERT_GEQ(stepSize, 0);
    
    if (node.isBase)
    {
        return stepSize - _min(node.endEdgeIt - node.curEdgeIt, stepSize);
    }
    
    TSize remainingEdgeSize = node.endEdgeIt - node.curEdgeIt;
    if (stepSize < remainingEdgeSize)
    {
        if (stepSize > static_cast<TSize>(node.remainingSize))
        {
            node.remainingSize = 0;
            node.curEdgeIt += node.remainingSize;
            return stepSize - node.remainingSize;
        }
        node.curEdgeIt += stepSize;
        node.remainingSize -= stepSize;
        return 0;
    }
    
    if (stepSize > static_cast<TSize>(node.remainingSize))
    {
        node.remainingSize = 0;
        node.curEdgeIt += node.remainingSize;
        return stepSize - node.remainingSize;
    }
    node.curEdgeIt = node.endEdgeIt;
    node.remainingSize -= remainingEdgeSize;
    return stepSize - remainingEdgeSize;
}

template <typename TJst, typename TSpec, typename TObserver,
          typename TTraversalNode>
inline void
pushNode(TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> & traversor,
         TTraversalNode SEQAN_FORWARD_CARG child)
{
    appendValue(*traversor._stackPtr, SEQAN_FORWARD(TTraversalNode, child));
    notify(traversor, PushEvent());
}

template <typename TJst, typename TSpec, typename TObserver>
inline void
popNode(TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> & traversor)
{
    eraseBack(*traversor._stackPtr);
    notify(traversor, PopEvent());
}

template <typename TStack>
inline TStack*
createStack()
{
    return new TStack;
}


template <typename TJst, typename TSpec, typename TObserver>
inline typename Member<TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver>, TraverserStackMember>::Type &
stack(TraverserImpl<TJst, JstTraversalSpec<TSpec> > & me)
{
    return *me._stackPtr;
}

template <typename TJst, typename TSpec, typename TObserver>
inline typename Member<TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> const, TraverserStackMember>::Type &
stack(TraverserImpl<TJst, JstTraversalSpec<TSpec> > const & me)
{
    return *me._stackPtr;
}

// ----------------------------------------------------------------------------
// Function impl::getStringContext()
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec, typename TObserver>
inline typename StringContext<TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> >::Type
getStringContext(TraverserImpl<TJst, JstTraversalSpec<TSpec> > const & traverser)
{
    // TODO(rrahn): Write me!
}


template <typename TJst, typename TSpec, typename TObserver,
          typename TTraversalNode,
          typename TSize>
inline TSize
moveWindow(TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> &, TTraversalNode &, TSize);

// ----------------------------------------------------------------------------
// Function impl::expandSubtree()
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec, typename TObserver,
          typename TTraversalNode,
          typename TSize>
inline void
expandSubtree(TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> & it,
              TTraversalNode & parent,
              TSize stepSize)
{
    parent.curDelta = parent.nextDelta;
    if (parent.isBase)
        it._tmp = parent;  // Store copy of the original base node.

    while (getDeltaPosition(*(*parent.nextDelta).hostIter) == getDeltaPosition(*(*parent.curDelta).hostIter))
    {
        auto child = parent;
        if (impl::nextChild(it, parent, child) && impl::moveWindow(it, child, stepSize) == 0 &&
            child.remainingSize >= 0)
        {
#if defined(DEBUG_JST_TRAVERSAL)
            std::cout << "     PUSH: " << child << std::endl;
#endif //defined(DEBUG_JST_TRAVERSAL)
            impl::pushNode(it, SEQAN_MOVE(child));
        }

        ++parent.nextDelta;  // Move to the next position.
    }

    if (parent.isBase)
    {
        // We need to move the base parent to the next segment redarding the branch information.
        impl::advanceBaseParent(it._tmp);
        swap(parent, it._tmp);  // Set advanced base parent back to stack at the correct position.
        it._tmp.isBase = false;
        it._tmp.fromBase = true;
        impl::updateParent(it._tmp);
        if (impl::moveWindow(it, it._tmp, stepSize) == 0 && it._tmp.remainingSize >= 0)
            impl::pushNode(it, it._tmp);
    }
    else
    {
        if (impl::updateParent(parent))
            return impl::moveWindow(it, parent, stepSize);
    }
}

// ----------------------------------------------------------------------------
// Function impl::moveWindow()
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec, typename TObserver,
          typename TTraversalNode,
          typename TSize>
inline TSize
moveWindow(TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> & it,
           TTraversalNode & parent,
           TSize stepSize)
{
#if defined(DEBUG_JST_TRAVERSAL)
    std::cout << "MOVE: " << parent << std::endl;
    std::cout << "STEP: " << stepSize << std::endl;
#endif //defined(DEBUG_JST_TRAVERSAL)
    stepSize = impl::moveWindow(parent, stepSize);
    if (stepSize > 0)
    {
        if (parent.remainingSize > 0)
        {
            impl::expandSubtree(it, parent, stepSize);
        }
    }
    return stepSize;
}

// ----------------------------------------------------------------------------
// Function init()
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec, typename TObserver>
inline void
init(TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> & me,
     TJst & jst)
{
    typedef typename TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver>::TNode TNode;

    me._contPtr = &jst;
    TNode node;
    resize(node.coverage, dimension(jst), true, Exact());

    node.curDelta = jst._buffer._deltaRangeBegin - 1;
    node.nextDelta = jst._buffer._deltaRangeBegin;

    node.begEdgeIt = begin(source(jst), Standard());  // This points to some value already -> what could this position be?
    node.curEdgeIt = node.begEdgeIt;
    node.endEdgeIt = node.begEdgeIt + (getDeltaPosition(*(*node.nextDelta).hostIter) - position(node.begEdgeIt));
    node.mappedSrcEndPos = position(node.endEdgeIt);
    node.remainingSize = me._historySize;
    node.isBase = true;
    node.fromBase = true;
    appendValue(*me._stackPtr, SEQAN_MOVE(node));  // Push onto stack.

    // After we realized this.
    TNode & base = back(*me._stackPtr);
    if ((*base.nextDelta).deltaPos == position(base.begEdgeIt))
    {
        expandSubtree(me, base, 0);
    }
}



}  // namespace impl

}  // namespace seqan.

#endif  // #ifndef INCLUDE_SEQAN_JOURNALED_STRING_TREE_JOURNALED_STRING_TREE_TRAVERSER_UTIL_H_
