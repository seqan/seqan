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
// Function impl::getPos()
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec, typename TObserver,
          typename THostIter>
inline typename Position<THostIter>::Type
getPos(TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> const & me,
       THostIter const & it)
{
    if (SEQAN_UNLIKELY(atEnd(it)))
        return length(source(container(me)));
    return (*it).deltaPos;
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
// Function impl::toNextDeltaBehindDeletion();
// ----------------------------------------------------------------------------

template <typename TJst, typename TSize>
inline void
toNextDeltaBehindDeletion(JstTraversalNode<TJst> & node,
                          TSize const delSize)
{
    // Move to next valid delta position. Skip all nodes that lie within the range of the deletion.
    while (!atEnd(node.nextDelta) && (impl::isEndPoint(node.nextDelta) ||
           (getDeltaPosition(*(impl::hostIter(node.nextDelta))) <
           getDeltaPosition(*impl::hostIter(node.curDelta)) + delSize)))
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
branchOut(TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> & it,
          TTraversalNode & parent,
          TTraversalNode & child)
{
    // We set the coverage of the left child to be the one of the parent & coverage(curDelta);
    transform(child.coverage, parent.coverage, getDeltaCoverage(*impl::hostIter(parent.nextDelta)), FunctorBitwiseAnd());
    transform(parent.coverage, parent.coverage, getDeltaCoverage(*impl::hostIter(parent.nextDelta)),
              FunctorNested<FunctorBitwiseAnd, FunctorIdentity, FunctorBitwiseNot>());

//#if defined(DEBUG_JST_TRAVERSAL)
    std::cout << "Coverage Child: " << _printCoverage(child.coverage) << std::endl;
    std::cout << "Coverage Paren: " << _printCoverage(parent.coverage) << std::endl;
//#endif // defined(DEBUG_JST_TRAVERSAL)
    // Set the current delta of the child to the next delta of the parent.
    child.curDelta = parent.nextDelta;

    // A) First get the proxyId of the updated child coverage.
    auto proxyId = bitScanForward(child.coverage);
    if (proxyId >= dimension(container(it)))
        return false;  // We can skip this child, since it has an empty coverage.

    // B) Remap the sequence iterators according to the new positions.
    bool inBegOfDelta = (parent.endEdgeIt - parent.begEdgeIt) <= 1;
    if (&container(it)._buffer._journaledSet[proxyId] == &container(parent.begEdgeIt))
    {
        // remap parent node to new journal sequence.
        proxyId = bitScanForward(parent.coverage);
        if (proxyId < dimension(container(it)) && !parent.fromBase)
        {
            parent.endEdgeIt = begin(container(it)._buffer._journaledSet[proxyId], Standard());
            impl::mapBranchPointToVirtual(parent.endEdgeIt, host(container(it)), proxyId,
                                          getDeltaPosition(*impl::hostIter(parent.nextDelta)));
            // We have the sequence iterators updated and the coverage.
            // So far the remainingLength and the curNode remain the same.
        }
        else  // We move to the base sequence.
        {
            parent.endEdgeIt = iter(container(it)._source, getDeltaPosition(*impl::hostIter(parent.nextDelta)));
        }
//        parent.begEdgeIt = parent.endEdgeIt - endToBegDist;
//        parent.curEdgeIt = parent.endEdgeIt;
    }
    else
    {
        // map child node to new journal sequence.
        // remap child node to new journal sequence.
        child.endEdgeIt = begin(container(it)._buffer._journaledSet[proxyId], Standard());
        impl::mapBranchPointToVirtual(child.endEdgeIt, host(container(it)), proxyId,
                                      getDeltaPosition(*impl::hostIter(child.curDelta)));
    }

    // Move to next valid delta position, note that we skip all end points and all deltas that occur at the same position.
    while ((!atEnd(child.nextDelta)) && (impl::isEndPoint(child.nextDelta) ||
           getPos(it, child.nextDelta) == getPos(it, child.curDelta)))
//    while (getDeltaPosition(*(impl::hostIter(++child.nextDelta))) == getDeltaPosition(*impl::hostIter(child.curDelta))
//           || impl::isEndPoint(child.nextDelta))
    {
        ++child.nextDelta;
    }

    // C) Update remaining length if we come directly from the base sequence.
    auto insSize = 0;
    auto delSize = 0;
    if (parent.isBase)
    {
        // We should keep the mappedSrcEndPos of the base node already in sync.
        child.mappedSrcEndPos += historySize(container(it));
        switch (getDeltaType(*impl::hostIter(child.curDelta)))
        {
            case DELTA_TYPE_DEL:
            {
                if (inBegOfDelta)
                    return false;
                delSize = deletionSize(host(container(it))._deltaStore,
                                       getStorePosition(*impl::hostIter(child.curDelta)), DeltaTypeDel());
                child.mappedSrcEndPos += delSize;
                // We need to skip all deleted
                impl::toNextDeltaBehindDeletion(child, delSize);
                --child.remainingSize;
                break;
            }
            case DELTA_TYPE_INS:
            {
                insSize = insertionSize(host(container(it))._deltaStore,
                                        getStorePosition(*impl::hostIter(child.curDelta)), DeltaTypeIns());
                child.remainingSize += insSize;
                break;
            }
            case DELTA_TYPE_SV:
            {
                delSize = deletionSize(host(container(it))._deltaStore,
                                       getStorePosition(*impl::hostIter(child.curDelta)), DeltaTypeSV());
                child.mappedSrcEndPos += delSize;
                insSize = insertionSize(host(container(it))._deltaStore,
                                        getStorePosition(*impl::hostIter(child.curDelta)), DeltaTypeSV());
                child.remainingSize += insSize - 1;
                impl::toNextDeltaBehindDeletion(child, delSize);
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
                delSize = deletionSize(host(container(it))._deltaStore,
                                       getStorePosition(*impl::hostIter(child.curDelta)), DeltaTypeDel());
                child.mappedSrcEndPos += delSize;
                impl::toNextDeltaBehindDeletion(child, delSize);
                break;
            }
            case DELTA_TYPE_INS:
            {
                insSize = insertionSize(host(container(it))._deltaStore,
                                        getStorePosition(*impl::hostIter(child.curDelta)), DeltaTypeIns());
                break;
            }
            case DELTA_TYPE_SV:
            {
                delSize = deletionSize(host(container(it))._deltaStore,
                                       getStorePosition(*impl::hostIter(child.curDelta)), DeltaTypeSV());

                child.mappedSrcEndPos += delSize;
                insSize = insertionSize(host(container(it))._deltaStore,
                                        getStorePosition(*impl::hostIter(child.curDelta)), DeltaTypeSV());
                impl::toNextDeltaBehindDeletion(child, delSize);
                break;
            }
            default: break;
        }
    }

    SEQAN_ASSERT(delSize <= impl::getPos(it, child.nextDelta) - impl::getPos(it, child.curDelta));

    child.begEdgeIt = child.endEdgeIt;
    child.curEdgeIt = child.begEdgeIt;
    child.endEdgeIt += (insSize + (getPos(it, child.nextDelta) - getPos(it, child.curDelta)) - delSize);
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
    // TODO(rrahn): Check if we need this?
//    arrayFill(begin(base.coverage, Standard()), end(base.coverage, Standard()), true);  // Make sure the coverage is set to 1.
    while (!atEnd(base.nextDelta) && impl::getPos(me, base.curDelta) == impl::getPos(me, base.nextDelta))
    {
        if (impl::isEndPoint(base.nextDelta))
            transform(base.coverage, base.coverage, getDeltaCoverage(*impl::hostIter(base.nextDelta)),
                      FunctorBitwiseOr());
        else
            impl::updateOnDeletion(base, impl::hostIter(base.nextDelta));  // Update the coverage of the base node if the current delta represents a deletion.
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
    while(!atEnd(parent.nextDelta) && impl::isEndPoint(parent.nextDelta))  // Move the nextDelta iterator to the next branch point, in case the current one is an endPoint.
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
    SEQAN_ASSERT_GEQ(stepSize, 0);

    if (stepSize < node.remainingSize || node.isBase)
    {
        auto minSteps = _min(stepSize, node.endEdgeIt - node.curEdgeIt);
        node.curEdgeIt += minSteps;
        if (!node.isBase)
            node.remainingSize -= minSteps;
        return stepSize - minSteps;
    }

    auto minSteps = _min(node.remainingSize, node.endEdgeIt - node.curEdgeIt);
    node.curEdgeIt += minSteps;
    if (!node.isBase)
        node.remainingSize -= minSteps;
    return stepSize - minSteps;
}

template <typename TJst, typename TSpec, typename TObserver,
          typename TTraversalNode>
inline void
pushNode(TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> & traverser,
         TTraversalNode SEQAN_FORWARD_CARG node)
{
#if defined(DEBUG_JST_TRAVERSAL)
//    std::cout << "-----> Journal " << container(node.endEdgeIt) << std::endl;
    std::cout << "        PUSH: (" << node << ")" << std::endl;
#endif //defined(DEBUG_JST_TRAVERSAL)
    appendValue(*traverser._stackPtr, SEQAN_FORWARD(TTraversalNode, node));
    notify(traverser, PushEvent());
}

template <typename TJst, typename TSpec, typename TObserver>
inline void
popNode(TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> & traverser)
{

#if defined(DEBUG_JST_TRAVERSAL)
    std::cout << "        POP: (" << back(*traverser._stackPtr) << ")" << std::endl;
#endif //defined(DEBUG_JST_TRAVERSAL)
    eraseBack(*traverser._stackPtr);
    notify(traverser, PopEvent());
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

template <typename TJst>
inline bool skipNode(JstTraversalNode<TJst> const & node)
{
    if (impl::isEndPoint(node.nextDelta))
        return true;
    if (node.isBase && getDeltaType(*impl::hostIter(node.nextDelta)) == DELTA_TYPE_DEL)
        return node.endEdgeIt - node.begEdgeIt <= 1;
    return false;
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
    if (parentPtr->isBase)
    {
        if ((*parentPtr->curDelta).deltaPos == 20)
        {
            std::cout << "###############\n"
                      << "Last Range: " << (*(parentPtr->curDelta - 1)).deltaPos << " to " << (*parentPtr->nextDelta).deltaPos << std::endl;
        }
        SEQAN_ASSERT(length(*it._stackPtr) == 1u);
        appendValue(*it._stackPtr, *parentPtr);  // Copy the base node and mark it from base. This represents all contexts without the delta at the current position.
        impl::advanceBase(it, *parentPtr);
//        back(*it._stackPtr).coverage = parentPtr->coverage;
        parentPtr = &back(*it._stackPtr);
        parentPtr->fromBase = true;
    }

    while (!atEnd(parentPtr->nextDelta) && impl::getPos(it, parentPtr->nextDelta) == impl::getPos(it, parentPtr->curDelta))
    {
        if (SEQAN_LIKELY(!impl::isEndPoint(parentPtr->nextDelta)))
        {
            auto child = *parentPtr;
            if (impl::branchOut(it, *parentPtr, child) && impl::moveWindow(it, &child, stepSize) == 0 &&
                child.remainingSize >= 0)
            {
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
    if (parentPtr->curEdgeIt == parentPtr->endEdgeIt && !atEnd(parentPtr->nextDelta))  // Reached branching point => expand node.
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

template <typename TJst, typename TSpec, typename TObserver>
inline void
init(TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver> & me,
     TJst & jst)
{
    typedef typename TraverserImpl<TJst, JstTraversalSpec<TSpec>, TObserver>::TNode TNode;

    SEQAN_ASSERT(me._historySize > 0);

    me._contPtr = &jst;
    TNode node;
    resize(node.coverage, dimension(jst), true, Exact());

    node.curDelta = jst._buffer._deltaRangeBegin - 1;
    node.nextDelta = jst._buffer._deltaRangeBegin;

    node.begEdgeIt = begin(source(jst), Standard());  // This points to some value already -> what could this position be?
    node.curEdgeIt = node.begEdgeIt;
    node.endEdgeIt = node.begEdgeIt + (getDeltaPosition(*(*node.nextDelta).hostIter) - position(node.begEdgeIt));
    node.mappedSrcEndPos = position(node.endEdgeIt);
    node.remainingSize = me._historySize - 1;
    node.isBase = true;
    node.fromBase = false;
    appendValue(*me._stackPtr, SEQAN_MOVE(node));  // Push onto stack.

    // After we realized this.
    TNode* base = &back(*me._stackPtr);

    // Record deletions at beginning.
    auto tmpDelta = base->nextDelta;
    while ((*tmpDelta).deltaPos == position(base->begEdgeIt))
    {
        impl::updateOnDeletion(*base, impl::hostIter(tmpDelta));
        ++tmpDelta;
    }

    // Initialize the traverser if there are events at the beginning of the tree.
    if ((*(base->nextDelta)).deltaPos == position(base->begEdgeIt))
    {
        impl::expandNode(me, base, 0);
    }
}


}  // namespace impl

}  // namespace seqan.

#endif  // #ifndef INCLUDE_SEQAN_JOURNALED_STRING_TREE_JOURNALED_STRING_TREE_TRAVERSER_UTIL_H_
