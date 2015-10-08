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

template <typename TObject>
struct StringContext;

template <typename TContainer, typename TSpec, typename TObserverList = ObserverList<> >
class TraverserImpl;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct TraverserStackMember_;
typedef Tag<TraverserStackMember_> TraverserStackMember;

template <typename TSpec = void>
struct JstTraversalSpec
{};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Member<TJst, JstBufferMember>
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec, typename TObserver>
struct Member<TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver>, JstBufferMember>
{
    typedef JstBuffer_<TJst> Type;
};

// ============================================================================
// Private Functions
// ============================================================================

namespace impl
{

// ----------------------------------------------------------------------------
// Function impl::createBuffer();
// ----------------------------------------------------------------------------

template <typename TBuffer>
inline TBuffer*
createBuffer()
{
    return new TBuffer;
}

// ----------------------------------------------------------------------------
// Function impl::buffer();
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec, typename TObserver>
inline typename Member<TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver>, JstBufferMember>::Type &
buffer(TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> & me)
{
    return *me._bufferPtr;
}

template <typename TJst, typename TSpec, typename TObserver>
inline typename Member<TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> const, JstBufferMember>::Type &
buffer(TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> const & me)
{
    return *me._bufferPtr;
}

// ----------------------------------------------------------------------------
// Function impl::createStack();
// ----------------------------------------------------------------------------

template <typename TStack>
inline TStack*
createStack()
{
    return new TStack;
}

// ----------------------------------------------------------------------------
// Function impl::stack();
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec, typename TObserver>
inline typename Member<TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver>, TraverserStackMember>::Type &
stack(TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> & me)
{
    return *me._stackPtr;
}

template <typename TJst, typename TSpec, typename TObserver>
inline typename Member<TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> const, TraverserStackMember>::Type &
stack(TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> const & me)
{
    return *me._stackPtr;
}

// ----------------------------------------------------------------------------
// Function impl::pushNode();
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec, typename TObserver,
          typename TTraversalNode>
inline void
pushNode(TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> & me,
         TTraversalNode SEQAN_FORWARD_CARG node)
{
#if defined(DEBUG_JST_TRAVERSAL)
    //    std::cout << "-----> Journal " << container(node.endEdgeIt) << std::endl;
    std::cout << "        PUSH: (" << node << ")" << std::endl;
#endif //defined(DEBUG_JST_TRAVERSAL)
    appendValue(impl::stack(me), SEQAN_FORWARD(TTraversalNode, node));
    notify(me, PushEvent());
}

// ----------------------------------------------------------------------------
// Function impl::popNode();
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec, typename TObserver>
inline void
popNode(TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> & me)
{

#if defined(DEBUG_JST_TRAVERSAL)
    std::cout << "        POP: (" << impl::activeNode(me) << ")" << std::endl;
#endif //defined(DEBUG_JST_TRAVERSAL)
    eraseBack(impl::stack(me));
    notify(me, PopEvent());
}

// ----------------------------------------------------------------------------
// Function impl::activeNode()
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec, typename TObserver>
inline JstTraversalNode<TJst> const &
activeNode(TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> const & me)
{
    return back(impl::stack(me));
}

template <typename TJst, typename TSpec, typename TObserver>
inline JstTraversalNode<TJst> &
activeNode(TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> & me)
{
    return back(impl::stack(me));
}

// ----------------------------------------------------------------------------
// Function impl::getContextIterator()
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec, typename TObserver>
inline typename ContextIterator<TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> const>::Type
getContextIterator(TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> const & me)
{
    return impl::activeNode(me).curEdgeIt;
}


// ----------------------------------------------------------------------------
// Function impl::mapSourceToVirtual()
// ----------------------------------------------------------------------------

template <typename TDeltaMapIter>
struct MapSourceToVirtualHelper_
{
    typedef typename Size<TDeltaMapIter>::Type TSize;


    TSize           virtOffset;
    TDeltaMapIter   iter;

    MapSourceToVirtualHelper_() : virtOffset(0)
    {}

    template <typename TTag>
    inline void
    operator()(TTag const & /*deltaType*/)
    {
        virtOffset += insertionSize(container(iter)._deltaStore, getStorePosition(*iter), TTag());
    }

};

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
    MapSourceToVirtualHelper_<TVarIterator>  f;
    TMapEntry child;
    child.deltaPosition = _physicalPosition(resultIt);
    f.iter = std::upper_bound(begin(variantStore, Standard()), end(variantStore, Standard()), child,
                              DeltaMapEntryPosLessThanComparator_());

    SEQAN_ASSERT_LEQ(getDeltaPosition(*f.iter), static_cast<TDeltaMapPos const>(hostPos));

//    TDeltaMapPos virtualOffset = 0;
    // Now we move to the right until we find the node that we are looking for and reconstruct the offset of the virtual positions.
    while (getDeltaPosition(*f.iter) != static_cast<TDeltaMapPos const>(hostPos) && !atEnd(f.iter))
    {
        if (getDeltaCoverage(*f.iter)[proxyId] != true || isRightEnd(*f.iter))  // irrelevant variant.
        {
            ++f.iter;
            continue;
        }
        DeltaTypeSelector selector;

        applyOnDelta(f, getDeltaType(*f.iter), selector);
        ++f.iter;
    }
    resultIt += f.virtOffset + 1;  // Set the iterator to the beginning of the variant.
}

// ----------------------------------------------------------------------------
// Function impl::getPos()
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec, typename TObserver,
          typename THostIter>
inline typename Position<THostIter>::Type
getPos(TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> const & me,
       THostIter const & it)
{
    if (SEQAN_UNLIKELY(atEnd(it)))
        return length(host(container(me)));
    return getDeltaPosition(*it);
}

// ----------------------------------------------------------------------------
// Function impl::toNextDeltaBehindDeletion();
// ----------------------------------------------------------------------------

template <typename TJst, typename TSize>
inline void
toNextDeltaBehindDeletion(JstTraversalNode<TJst> & node,
                          TSize const delSize)
{
    // Move to next valid delta position. Skip all nodes that lie within the range of the deletion.
    while (!atEnd(node.nextDelta) && (isRightEnd(*node.nextDelta) ||
           (getDeltaPosition(*node.nextDelta)) < getDeltaPosition(*node.curDelta) + delSize))
    {
        ++node.nextDelta;
    }
}

// ----------------------------------------------------------------------------
// Function impl::branchOut()
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec, typename TObserver,
          typename TTraversalNode>
inline bool
createBranch(TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> & me,
             TTraversalNode & parent,
             TTraversalNode & child)
{
    typedef typename TTraversalNode::TDeltaIterator TDeltaIt   SEQAN_TYPEDEF_FOR_DEBUG;
    typedef typename Position<TDeltaIt>::Type       TPos       SEQAN_TYPEDEF_FOR_DEBUG;
    // We set the coverage of the left child to be the one of the parent & coverage(curDelta);
    transform(child.coverage, parent.coverage, getDeltaCoverage(*parent.nextDelta), FunctorBitwiseAnd());
    transform(parent.coverage, parent.coverage, getDeltaCoverage(*parent.nextDelta),
              FunctorNested<FunctorBitwiseAnd, FunctorIdentity, FunctorBitwiseNot>());

//#if defined(DEBUG_JST_TRAVERSAL)
//    std::cout << "Coverage Child: " << _printCoverage(child.coverage) << std::endl;
//    std::cout << "Coverage Paren: " << _printCoverage(parent.coverage) << std::endl;
//#endif // defined(DEBUG_JST_TRAVERSAL)
    // Set the current delta of the child to the next delta of the parent.
    child.curDelta = parent.nextDelta;

    // A) First get the proxyId of the updated child coverage.
    auto proxyId = bitScanForward(child.coverage);
    if (proxyId >= length(container(me)))
        return false;  // We can skip this child, since it has an empty coverage.

    // B) Remap the sequence iterators according to the new positions.
    if (&buffer(me)._journaledSet[proxyId] == &container(parent.endEdgeIt))
    {
        // remap parent node to new journal sequence.
        proxyId = bitScanForward(parent.coverage);
        if (proxyId < length(container(me)) && !parent.fromBase)
        {
            parent.endEdgeIt = begin(buffer(me)._journaledSet[proxyId], Standard());
            impl::mapBranchPointToVirtual(parent.endEdgeIt, impl::member(container(me), JstDeltaMapMember()), proxyId,
                                          getDeltaPosition(*parent.nextDelta));
            // We have the sequence iterators updated and the coverage.
            // So far the remainingLength and the curNode remain the same.
        }
        else  // We move to the base sequence.
        {
            parent.endEdgeIt = iter(impl::member((container(me)), JstSourceMember()), getDeltaPosition(*parent.nextDelta));
        }
    }
    else
    {
        // map child node to new journal sequence.
        // remap child node to new journal sequence.
        child.endEdgeIt = begin(buffer(me)._journaledSet[proxyId], Standard());
        impl::mapBranchPointToVirtual(child.endEdgeIt, impl::member(container(me), JstDeltaMapMember()), proxyId,
                                      getDeltaPosition(*child.curDelta));
    }

    // Move to next valid delta position, note that we skip all end points and all deltas that occur at the same position.
    while (!atEnd(child.nextDelta) && (isRightEnd(*child.nextDelta) || getPos(me, child.nextDelta) == getPos(me, child.curDelta)))
    {
        ++child.nextDelta;
    }

    // C) Update remaining length if we come directly from the base sequence.
    auto iSize = insSize(child.curDelta);
    auto dSize = delSize(child.curDelta);
    if (parent.isBase)
    {
        // We should keep the mappedSrcEndPos of the base node already in sync.
        child.mappedSrcEndPos += me._branchLength;
        switch (getDeltaType(*child.curDelta))
        {
            case DELTA_TYPE_DEL:
            {
                if (child.remainingSize == 0)
                    return false;
                child.mappedSrcEndPos += dSize;
                // We need to skip all deleted
                impl::toNextDeltaBehindDeletion(child, dSize);
                --child.remainingSize;
                break;
            }
            case DELTA_TYPE_INS:
            {
                child.remainingSize += iSize;
                break;
            }
            case DELTA_TYPE_SV:
            {
                child.mappedSrcEndPos += dSize;
                child.remainingSize += iSize - 1;
                impl::toNextDeltaBehindDeletion(child, dSize);
                break;
            }
            default: break;
        }
    }
    else
    {
        // Now we are down somewhere in an internal subtree.
        switch (getDeltaType(*child.curDelta))
        {
            case DELTA_TYPE_DEL:
            {
                child.mappedSrcEndPos += dSize;
                impl::toNextDeltaBehindDeletion(child, dSize);
                break;
            }
            case DELTA_TYPE_SV:
            {
                child.mappedSrcEndPos += dSize;
                impl::toNextDeltaBehindDeletion(child, dSize);
                break;
            }
            default: break;
        }
    }

    SEQAN_ASSERT(static_cast<TPos>(dSize) <= impl::getPos(me, child.nextDelta) - impl::getPos(me, child.curDelta));

    child.begEdgeIt = child.endEdgeIt;
    child.curEdgeIt = child.begEdgeIt;
    child.endEdgeIt += (iSize + (getPos(me, child.nextDelta) - getPos(me, child.curDelta)) - dSize);
    child.isBase = false;
    child.fromBase = false;
    return true;
}

// ----------------------------------------------------------------------------
// Function impl::updateOnDeletion()
// ----------------------------------------------------------------------------

template <typename TTraversalNode,
          typename THostIter>
inline void
updateOnDeletion(TTraversalNode & base, THostIter const & hostIt)
{
    switch (getDeltaType(*hostIt))
    {
        case DELTA_TYPE_DEL:
        {
            if (deletionSize(container(hostIt)._deltaStore, getStorePosition(*hostIt), DeltaTypeDel()) > 1)
                transform(base.coverage, base.coverage, getDeltaCoverage(*hostIt),
                          FunctorNested<FunctorBitwiseAnd, FunctorIdentity, FunctorBitwiseNot>());
            break;
        }
        case DELTA_TYPE_SV:
        {
            if (deletionSize(container(hostIt)._deltaStore, getStorePosition(*hostIt), DeltaTypeSV()) > 1)
                transform(base.coverage, base.coverage, getDeltaCoverage(*hostIt),
                          FunctorNested<FunctorBitwiseAnd, FunctorIdentity, FunctorBitwiseNot>());
            break;
        }
        default: break;
    }
}

// ----------------------------------------------------------------------------
// Function impl::advanceBase()
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec, typename TObserver,
          typename TTraversalNode>
inline void
advanceBase(TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> const & me,
            TTraversalNode & base)
{
    SEQAN_ASSERT(base.curDelta == base.nextDelta);

#if defined(DEBUG_JST_TRAVERSAL)
    std::cout << "        BASE before --> " << base << std::endl;
#endif // defined(DEBUG_JST_TRAVERSAL)

    while (!atEnd(base.nextDelta) && impl::getPos(me, base.curDelta) == impl::getPos(me, base.nextDelta))
    {
        if (isRightEnd(*base.nextDelta))
            transform(base.coverage, base.coverage, getDeltaCoverage(*base.nextDelta), FunctorBitwiseOr());
        else
            impl::updateOnDeletion(base, base.nextDelta);  // Update the coverage of the base node if the current delta represents a deletion.
        ++base.nextDelta;
    }
    // base.nextDelta might point to an merge point.
    base.begEdgeIt = base.endEdgeIt;
    //getDeltaPosition(*impl::hostIter(base.curDelta)) - position(base.begEdgeIt);
    base.curEdgeIt = base.begEdgeIt;
    auto posN = impl::getPos(me, base.nextDelta);
    auto posC = impl::getPos(me, base.curDelta);
    base.endEdgeIt += posN - posC;
    base.mappedSrcEndPos = posN;  // TODO(rrahn): Remove.
#if defined(DEBUG_JST_TRAVERSAL)
    std::cout << "        BASE after ---> " << base << std::endl;
#endif // defined(DEBUG_JST_TRAVERSAL)
}

// ----------------------------------------------------------------------------
// Function impl::updateParent()
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec, typename TObserver,
          typename TTraversalNode>
inline bool
updateParent(TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> const & me,
             TTraversalNode & parent)
{
    SEQAN_ASSERT_NOT(parent.isBase);        // Should never be the base node.

#if defined(DEBUG_JST_TRAVERSAL)
    std::cout << "        EXPAND Parent before --> " << parent << std::endl;
#endif // defined(DEBUG_JST_TRAVERSAL)

//    if (testAllZeros(parent.coverage))
//        return false;
    // Update:
    parent.begEdgeIt = parent.endEdgeIt;
    parent.curEdgeIt = parent.begEdgeIt;

    // Make sure that for a right subtree we don't consider merge points, because we cannot be within a branch
    // coming from a delta while at the same time there is a deletion to be covered.
    while(!atEnd(parent.nextDelta) && isRightEnd(*parent.nextDelta))  // Move the nextDelta iterator to the next branch point, in case the current one is an endPoint.
    {
        ++parent.nextDelta;
    }
    parent.mappedSrcEndPos += impl::getPos(me, parent.nextDelta) - impl::getPos(me, parent.curDelta);
    auto posC = impl::getPos(me, parent.curDelta);
    parent.endEdgeIt += impl::getPos(me, parent.nextDelta) - posC;

#if defined(DEBUG_JST_TRAVERSAL)
    std::cout << "        EXPAND Parent after ---> " << parent << std::endl;
#endif // defined(DEBUG_JST_TRAVERSAL)
    return true;
}

// ----------------------------------------------------------------------------
// Function impl::shiftWindowBy()
// ----------------------------------------------------------------------------

template <typename TNode, typename TSize>
inline TSize
shiftWindowBy(TNode & node, TSize stepSize)
{
    typedef typename TNode::TSeqIterator TSeqIter;
    typedef typename Difference<TSeqIter>::Type TDiff;

    SEQAN_ASSERT_GEQ(stepSize, static_cast<TSize>(0));

    if (stepSize < static_cast<TSize>(node.remainingSize) || node.isBase)
    {
        auto minSteps = _min(static_cast<TDiff>(stepSize), node.endEdgeIt - node.curEdgeIt);
        node.curEdgeIt += minSteps;
        if (!node.isBase)
            node.remainingSize -= minSteps;
        return stepSize - minSteps;
    }

    auto minSteps = _min(static_cast<TDiff>(node.remainingSize), node.endEdgeIt - node.curEdgeIt);
    node.curEdgeIt += minSteps;
    if (!node.isBase)
        node.remainingSize -= minSteps;
    return stepSize - minSteps;
}

// ----------------------------------------------------------------------------
// Function impl::moveWindow();
// ----------------------------------------------------------------------------

// Forward declaration for recursive call.
template <typename TJst, typename TSpec, typename TObserver,
          typename TTraversalNode,
          typename TSize>
inline TSize
moveWindow(TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> &, TTraversalNode*, TSize);

// ----------------------------------------------------------------------------
// Function impl::expandNode();
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec, typename TObserver,
          typename TTraversalNode,
          typename TSize>
inline TSize
expandNode(TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> & it,
           TTraversalNode * parentPtr,
           TSize stepSize)
{

#if defined(DEBUG_JST_TRAVERSAL)
    std::cout << "    Expand (" << &parent << ") " << std::endl;
#endif //defined(DEBUG_JST_TRAVERSAL)
    parentPtr->curDelta = parentPtr->nextDelta;

    if (SEQAN_UNLIKELY(atEnd(parentPtr->curDelta)))  // Reached end of tree.
    {
        arrayFill(begin(parentPtr->coverage), end(parentPtr->coverage), false);  // Set coverage to 0.
        return stepSize;
    }

    if (parentPtr->isBase)  // At a branching point coming from base.
    {
        SEQAN_ASSERT(length(*it._stackPtr) == 1u);
        appendValue(*it._stackPtr, *parentPtr);  // Copy the base node and mark it from base. This represents all contexts without the delta at the current position.
        impl::advanceBase(it, *parentPtr);
        parentPtr = &back(*it._stackPtr);
        parentPtr->fromBase = true;
    }

    while (!atEnd(parentPtr->nextDelta) && impl::getPos(it, parentPtr->nextDelta) == impl::getPos(it, parentPtr->curDelta))
    {
        if (SEQAN_LIKELY(!isRightEnd(*parentPtr->nextDelta)))  // Skip points, where we merge a deletion.
        {
            auto child = *parentPtr;
            if (impl::createBranch(it, *parentPtr, child) && impl::moveWindow(it, &child, stepSize) == 0 &&
                child.remainingSize >= 0)
            {
                if (SEQAN_LIKELY(!atEnd(child.curDelta)))  // Skip the node in case we reached the end already.
                    impl::pushNode(it, SEQAN_MOVE(child));
            }
        }
        ++parentPtr->nextDelta;  // Move to the next branch point.
    }
    parentPtr->isBase = false;
    impl::updateParent(it, *parentPtr);
    return impl::moveWindow(it, parentPtr, stepSize);  // Recursive call to move as long as stepSize is greater 0.
}

// ----------------------------------------------------------------------------
// Function impl::moveWindow()
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec, typename TObserver,
          typename TTraversalNode,
          typename TSize>
inline TSize
moveWindow(TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> & it,
           TTraversalNode* parentPtr,
           TSize stepSize)
{
#if defined(DEBUG_JST_TRAVERSAL)
    std::cout << "############################## BEGIN ########################################" << std::endl;
    std::cout << "MOVE by (" << stepSize << ") -> " << parentPtr << std::endl;
#endif //defined(DEBUG_JST_TRAVERSAL)
    stepSize = impl::shiftWindowBy(*parentPtr, stepSize);
    if (parentPtr->curEdgeIt == parentPtr->endEdgeIt)  // Reached branching point => expand node.
        stepSize = impl::expandNode(it, parentPtr, stepSize);

#if defined(DEBUG_JST_TRAVERSAL)
    std::cout << "Remaining (" << stepSize << ") -> " << parentPtr << std::endl;
    std::cout << "================================ END =========================================" << std::endl;
#endif //defined(DEBUG_JST_TRAVERSAL)
    return stepSize;
}

// ----------------------------------------------------------------------------
// Function init()
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec, typename TObserver,
          typename TSize>
inline void
init(TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> & me,
     TJst & jst,
     TSize const contextSize,
     TSize const branchLength)
{
    typedef typename TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver>::TNode TNode;

    SEQAN_ASSERT(contextSize > 0);
    SEQAN_ASSERT(branchLength > 0);

    me._contextSize = contextSize;
    me._branchLength = _max(contextSize, branchLength);
    me._contPtr = &jst;

    SEQAN_ASSERT(me._bufferPtr);
    SEQAN_ASSERT(me._stackPtr);

    clear(buffer(me));
    clear(stack(me));

    init(buffer(me), jst);

    TNode node;
    resize(node.coverage, length(jst), true, Exact());

    node.curDelta = buffer(me)._deltaRangeBegin - 1;
    node.nextDelta = buffer(me)._deltaRangeBegin;

    node.begEdgeIt = begin(impl::member(jst, JstSourceMember()), Standard());  // This points to some value already -> what could this position be?
    node.curEdgeIt = node.begEdgeIt;
    node.endEdgeIt = node.begEdgeIt + (getDeltaPosition(*node.nextDelta) - position(node.begEdgeIt));
    node.mappedSrcEndPos = position(node.endEdgeIt);
    node.remainingSize = me._branchLength - 1;
    node.isBase = true;
    node.fromBase = false;
    appendValue(*me._stackPtr, SEQAN_MOVE(node));  // Push onto stack.

    // After we realized this.
    TNode* basePtr = &impl::activeNode(me);

    SEQAN_ASSERT_GEQ(me._contextSize, 1u);
    impl::moveWindow(me, basePtr, me._contextSize - 1);   // We move the traverser to the beginning of the
}

}  // namespace impl

}  // namespace seqan.

#endif  // #ifndef INCLUDE_SEQAN_JOURNALED_STRING_TREE_JOURNALED_STRING_TREE_TRAVERSER_UTIL_H_
