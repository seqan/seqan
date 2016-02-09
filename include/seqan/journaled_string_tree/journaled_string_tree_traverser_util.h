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

#if defined(JST_FIND_DEBUG)
    StringSet<DnaString> __testSet;
#endif // JST_FIND_DEBUG

// ============================================================================
// Forwards
// ============================================================================

template <typename TObject>
struct StringContext;

template <typename TContainer, typename TSpecList = ObserverList<> >
class TraverserImpl;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct TraverserStackMember_;
typedef Tag<TraverserStackMember_> TraverserStackMember;

template <typename TSpec = void>
struct JstTraversalSpec
{};

// ----------------------------------------------------------------------------
// Tag SlectFirstProxy
// ----------------------------------------------------------------------------

struct SelectFirstProxy_;
typedef Tag<SelectFirstProxy_> SelectFirstProxy;

// ----------------------------------------------------------------------------
// Tag SlectValidProxy
// ----------------------------------------------------------------------------

struct SelectValidProxy_;
typedef Tag<SelectValidProxy_> SelectValidProxy;

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Member<TJst, JstBufferMember>
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec>
struct Member<TraverserImpl<TJst, JstTraversalSpec<TSpec> >, JstBufferMember>
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

template <typename TJst, typename TSpec>
inline typename Member<TraverserImpl<TJst, JstTraversalSpec<TSpec> >, JstBufferMember>::Type &
buffer(TraverserImpl<TJst, JstTraversalSpec<TSpec> > & me)
{
    return *me._bufferPtr;
}

template <typename TJst, typename TSpec>
inline typename Member<TraverserImpl<TJst, JstTraversalSpec<TSpec> > const, JstBufferMember>::Type &
buffer(TraverserImpl<TJst, JstTraversalSpec<TSpec> > const & me)
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

template <typename TJst, typename TSpec>
inline typename Member<TraverserImpl<TJst, JstTraversalSpec<TSpec> >, TraverserStackMember>::Type &
stack(TraverserImpl<TJst, JstTraversalSpec<TSpec> > & me)
{
    return *me._stackPtr;
}

template <typename TJst, typename TSpec>
inline typename Member<TraverserImpl<TJst, JstTraversalSpec<TSpec> > const, TraverserStackMember>::Type &
stack(TraverserImpl<TJst, JstTraversalSpec<TSpec> > const & me)
{
    return *me._stackPtr;
}

// ----------------------------------------------------------------------------
// Function impl::pushNode();
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec,
          typename TTraversalNode,
          typename TObserver>
inline void
pushNode(TraverserImpl<TJst, JstTraversalSpec<TSpec> > & me,
         TTraversalNode && node,
         TObserver & observer)
{
#if defined(DEBUG_JST_TRAVERSAL)
    //    std::cout << "-----> Journal " << container(node.endEdgeIt) << std::endl;
    std::cout << "        PUSH: (" << node << ")" << std::endl;
#endif //defined(DEBUG_JST_TRAVERSAL)
    appendValue(impl::stack(me), std::forward<TTraversalNode>(node));
    notify(observer, PushEvent());
}

// ----------------------------------------------------------------------------
// Function impl::popNode();
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec,
          typename TObserver>
inline void
popNode(TraverserImpl<TJst, JstTraversalSpec<TSpec> > & me,
        TObserver & observer)
{

#if defined(DEBUG_JST_TRAVERSAL)
    std::cout << "        POP: (" << impl::activeNode(me) << ")" << std::endl;
#endif //defined(DEBUG_JST_TRAVERSAL)
    eraseBack(impl::stack(me));
    notify(observer, PopEvent());
}

// ----------------------------------------------------------------------------
// Function impl::activeNode()
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec>
inline JstTraversalNode<TJst> const &
activeNode(TraverserImpl<TJst, JstTraversalSpec<TSpec> > const & me)
{
    return back(impl::stack(me));
}

template <typename TJst, typename TSpec>
inline JstTraversalNode<TJst> &
activeNode(TraverserImpl<TJst, JstTraversalSpec<TSpec> > & me)
{
    return back(impl::stack(me));
}

// ----------------------------------------------------------------------------
// Function impl::baseNode()
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec>
inline JstTraversalNode<TJst> const &
baseNode(TraverserImpl<TJst, JstTraversalSpec<TSpec> > const & me)
{
    return front(impl::stack(me));
}

template <typename TJst, typename TSpec>
inline JstTraversalNode<TJst> &
baseNode(TraverserImpl<TJst, JstTraversalSpec<TSpec> > & me)
{
    return front(impl::stack(me));
}

// ----------------------------------------------------------------------------
// Function impl::getContextIterator()
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec>
inline typename ContextIterator<TraverserImpl<TJst, JstTraversalSpec<TSpec> > const>::Type
getContextIterator(TraverserImpl<TJst, JstTraversalSpec<TSpec> > const & me)
{
    return impl::activeNode(me).curEdgeIt;
}

// ----------------------------------------------------------------------------
// Function impl::getContextBegin()
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec>
inline typename ContextIterator<TraverserImpl<TJst, JstTraversalSpec<TSpec> > const>::Type
getContextBegin(TraverserImpl<TJst, JstTraversalSpec<TSpec> > const & me)
{
    if (length(*me._stackPtr) == 1)
        return  (sourceBegin(buffer(me)) + contextSize(me) > activeNode(me).curEdgeIt) ? sourceBegin(buffer(me)) :
            activeNode(me).curEdgeIt - (contextSize(me) - 1);
    else
        return activeNode(me).curEdgeIt  - (contextSize(me) - 1);
}

// ----------------------------------------------------------------------------
// Function impl::getContextEnd()
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec>
inline typename ContextIterator<TraverserImpl<TJst, JstTraversalSpec<TSpec> > const>::Type
getContextEnd(TraverserImpl<TJst, JstTraversalSpec<TSpec> > const & me)
{
    return getContextIterator(me);
}

// ----------------------------------------------------------------------------
// Function impl::refineCoverage()
// ----------------------------------------------------------------------------

template <typename TCoverage,
          typename TNode>
inline void
refineCoverage(TCoverage & cov,
               TNode const & node)
{
    auto tmp = node.headDelta;
    while (!atEnd(tmp) && tmp != node.branchRoot)
    {
#if defined(JST_FIND_DEBUG)
        std::cout << "(" << cov << ") & ~(" << getDeltaCoverage(*tmp) << ") = " << std::flush;
#endif  // defi ned(JST_FIND_DEBUG)
        if (!isRightEnd(*tmp))
            transform(cov, cov, getDeltaCoverage(*tmp), FunctorNested<FunctorBitwiseAnd, FunctorIdentity, FunctorBitwiseNot>());
#if defined(JST_FIND_DEBUG)
        std::cout << "(" << cov << ")" << std::endl;
#endif  // defined(JST_FIND_DEBUG)
        ++tmp;
    }
}

// ----------------------------------------------------------------------------
// Function impl::refineCoverageBranch()
// ----------------------------------------------------------------------------

template <typename TCoverage, typename TNode>
inline void
refineCoverageBranch(TCoverage & cov, TNode const & activeNode)
{
    // We go from the head of the active node to the branch root.
    // Now we can simply switch of all sequences that still affect the region between the head of the
    // active node and the root of the branch. Since those candidates have been searched already previously.
    auto tmp = activeNode.headDelta;
    while (!atEnd(tmp) && tmp != activeNode.branchRoot)
    {
        transform(cov, cov, getDeltaCoverage(*tmp), FunctorNested<FunctorBitwiseAnd, FunctorIdentity, FunctorBitwiseNot>());
        ++tmp;
    }
}

// ----------------------------------------------------------------------------
// Function impl::mapBranchPointToVirtual()
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

template <typename TJst, typename TSpec,
          typename THostIter>
inline typename Position<THostIter>::Type
getPos(TraverserImpl<TJst, JstTraversalSpec<TSpec> > const & me,
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
// Function impl::selectProxy()
// ----------------------------------------------------------------------------

template <typename TTraverserNode>
inline auto
selectProxy(TTraverserNode const & node, SelectFirstProxy const & /*tag*/) -> decltype(bitScanForward(node.coverage))
{
    return bitScanForward(node.coverage);
}

template <typename TTraverserNode>
inline auto
selectProxy(TTraverserNode const & node, SelectValidProxy const & /*tag*/) -> decltype(bitScanForward(node.coverage))
{
    // Simply take the first proxy, if there are no deltas reaching into
    // the region between head and branch root.
    if (node.headDelta == node.branchRoot)
        return (node.coverage[node.proxyId]) ? node.proxyId : bitScanForward(node.coverage);

    auto tmp1 = node.coverage;
    decltype(tmp1) tmp2;
    auto deltaIt = node.branchRoot;
    while(deltaIt != node.headDelta)
    {
        --deltaIt;
        transform(tmp2, tmp1, getDeltaCoverage(*deltaIt),
                  FunctorNested<FunctorBitwiseAnd, FunctorIdentity, FunctorBitwiseNot>());
        if (testAllZeros(tmp2))
            return (tmp1[node.proxyId]) ? node.proxyId : bitScanForward(tmp1);
        tmp1 = tmp2;
    }
    return (tmp2[node.proxyId]) ? node.proxyId : bitScanForward(tmp2);
}

// ----------------------------------------------------------------------------
// Function impl::createBranch()
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec,
          typename TTraversalNode,
          typename TProxySelection>
inline bool
createBranch(TraverserImpl<TJst, JstTraversalSpec<TSpec> > & me,
             TTraversalNode & parent,
             TTraversalNode & child,
             Tag<TProxySelection> const & /*tag*/)
{
    typedef typename TTraversalNode::TDeltaIterator TDeltaIt   SEQAN_TYPEDEF_FOR_DEBUG;
    typedef typename Position<TDeltaIt>::Type       TPos       SEQAN_TYPEDEF_FOR_DEBUG;

    // We set the coverage of the left child to be the one of the parent & coverage(curDelta);
    // If coming from the base we use the base coverage, i.e. all bits true, to avoid
    // cicumstantial updating of the coverage when the context head passes previously
    // visited branch nodes, who still affect the current base node.
    if (parent.isBase)
    {
        child.branchRoot = parent.nextDelta;
        updateOnDeletion(parent, parent.nextDelta);  // Only update base coverage if current delta is Del|SV.
        transform(child.coverage, me._baseCov, getDeltaCoverage(*parent.nextDelta), FunctorBitwiseAnd());
    }
    else  // Update the parent coverage.
    {
        transform(child.coverage, parent.coverage, getDeltaCoverage(*parent.nextDelta), FunctorBitwiseAnd());
        transform(parent.coverage, parent.coverage, getDeltaCoverage(*parent.nextDelta),
                  FunctorNested<FunctorBitwiseAnd, FunctorIdentity, FunctorBitwiseNot>());
    }

    // Set the current delta of the child to the next delta of the parent.
    child.curDelta = parent.nextDelta;

    // A) First get the proxyId of the updated child coverage.
    child.proxyId = impl::selectProxy(child, Tag<TProxySelection>());
    if (child.proxyId >= length(container(me)))
        return false;  // We can skip this child, since it has an empty coverage.

    // B) Remap the sequence iterators according to the new positions.
    if (&buffer(me)._journaledSet[child.proxyId] == &container(parent.endEdgeIt))
    {
        // remap parent node to new journal sequence.
        parent.proxyId = impl::selectProxy(parent, Tag<TProxySelection>());
        if (parent.proxyId < length(container(me)) && !parent.isBase)
        {
            parent.endEdgeIt = begin(buffer(me)._journaledSet[parent.proxyId], Standard());
            impl::mapBranchPointToVirtual(parent.endEdgeIt, impl::member(container(me), JstDeltaMapMember()), parent.proxyId,
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
        child.endEdgeIt = begin(buffer(me)._journaledSet[child.proxyId], Standard());
        impl::mapBranchPointToVirtual(child.endEdgeIt, impl::member(container(me), JstDeltaMapMember()), child.proxyId,
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
        switch (getDeltaType(*child.curDelta))
        {
            case DELTA_TYPE_DEL:
            {
                if (child.remainingSize == 0)
                    return false;
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
                impl::toNextDeltaBehindDeletion(child, dSize);
                break;
            }
            case DELTA_TYPE_SV:
            {
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
    return true;
}

// ----------------------------------------------------------------------------
// Function updateContextHead();
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec,
          typename TTraversalNode>
inline void
updateContextHead(TraverserImpl<TJst, JstTraversalSpec<TSpec> > const & me,
                  TTraversalNode & parent)
{
    using TPos = decltype(parent.headSrcPos);
    // headSrcPos is ahead of headDelta
    while (!atEnd(parent.headDelta) && parent.headDelta < parent.branchRoot &&
           static_cast<TPos>(impl::getPos(me, parent.headDelta)) < parent.headSrcPos)
    {
        if (parent.isBase)
        {
            if (isRightEnd(*parent.headDelta))  // If we passed the head of the context beyond a merge point we update the coverage accordingly.
                transform(parent.coverage, parent.coverage, getDeltaCoverage(*parent.headDelta), FunctorBitwiseOr());
            else  // If we reach a deletion again, we need to record this.
                updateOnDeletion(parent, parent.headDelta);
        }
        ++parent.headDelta;
    }

    // Now we have to parse all nodes between current head and root for the base node.
    // If we find a deletion or SV we make sure the corresponding sequences are
    // deleted from the base coverage.
    if (parent.isBase)
    {
        auto tmp = parent.headDelta;
        while (!atEnd(tmp) && tmp < parent.branchRoot)
        {
            updateOnDeletion(parent, tmp);
            ++tmp;
        }

    }
}

// ----------------------------------------------------------------------------
// Function impl::advanceParent()
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec,
          typename TTraversalNode>
inline bool
advanceParent(TraverserImpl<TJst, JstTraversalSpec<TSpec> > const & me,
             TTraversalNode & parent)
{
//    SEQAN_ASSERT_NOT(parent.isBase);        // Should never be the base node.

#if defined(DEBUG_JST_TRAVERSAL)
    std::cout << "        EXPAND Parent before --> " << parent << std::endl;
#endif // defined(DEBUG_JST_TRAVERSAL)

    // Update:
    parent.begEdgeIt = parent.endEdgeIt;
    parent.curEdgeIt = parent.begEdgeIt;

    // Make sure that for an inner subtree we don't consider merge points, because we cannot be within a branch
    // coming from a delta while at the same time there is a deletion to be covered.
    // However when coming from the base the next node might as well be a mergepoint to be considered.
    // Move the nextDelta iterator to the next branch point, in case the current one is an endPoint.
    while(!atEnd(parent.nextDelta) && isRightEnd(*parent.nextDelta))
        ++parent.nextDelta;

    if (parent.isBase)
        parent.branchRoot = parent.nextDelta;

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
        node.headSrcPos += minSteps;
        if (!node.isBase)
            node.remainingSize -= minSteps;
        return stepSize - minSteps;
    }

    auto minSteps = _min(static_cast<TDiff>(node.remainingSize), node.endEdgeIt - node.curEdgeIt);
    node.curEdgeIt += minSteps;
    node.headSrcPos += minSteps;
    if (!node.isBase)
        node.remainingSize -= minSteps;
    return stepSize - minSteps;
}

// ----------------------------------------------------------------------------
// Function impl::moveWindow();
// ----------------------------------------------------------------------------

// Forward declaration for recursive call.
template <typename TJst, typename TSpec,
          typename TTraversalNode,
          typename TSize,
          typename TObserver,
          typename TProxySelector>
inline TSize
moveWindow(TraverserImpl<TJst, JstTraversalSpec<TSpec> > &, TTraversalNode*, TSize, TObserver&, TProxySelector const &);

// ----------------------------------------------------------------------------
// Function impl::expandNode();
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec,
          typename TTraversalNode,
          typename TSize,
          typename TObserver,
          typename TProxySelector>
inline TSize
expandNode(TraverserImpl<TJst, JstTraversalSpec<TSpec> > & it,
           TTraversalNode * parentPtr,
           TSize stepSize,
           TObserver & observer,
           TProxySelector const & /*tag*/)
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

    while (!atEnd(parentPtr->nextDelta) && impl::getPos(it, parentPtr->nextDelta) == impl::getPos(it, parentPtr->curDelta))
    {
        if (SEQAN_LIKELY(!isRightEnd(*parentPtr->nextDelta)))  // Skip points, where we merge a deletion.
        {
            auto child = *parentPtr;
            if (impl::createBranch(it, *parentPtr, child, TProxySelector()) &&
                impl::moveWindow(it, &child, stepSize, observer, TProxySelector()) == 0 && child.remainingSize >= 0)
            {
                if (SEQAN_LIKELY(!atEnd(child.curDelta)))  // Skip the node in case we reached the end already.
                    impl::pushNode(it, SEQAN_MOVE(child), observer);
            }
        }
        ++parentPtr->nextDelta;  // Move to the next branch point.
    }
    impl::advanceParent(it, *parentPtr);
    return impl::moveWindow(it, parentPtr, stepSize, observer, TProxySelector());  // Recursive call to move as long as stepSize is greater 0.
}

// ----------------------------------------------------------------------------
// Function impl::moveWindow()
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec,
          typename TTraversalNode,
          typename TSize,
          typename TObserver,
          typename TProxySelector>
inline TSize
moveWindow(TraverserImpl<TJst, JstTraversalSpec<TSpec> > & me,
           TTraversalNode* parentPtr,
           TSize stepSize,
           TObserver & observer,
           TProxySelector const & /*tag*/)
{
#if defined(DEBUG_JST_TRAVERSAL)
    std::cout << "############################## BEGIN ########################################" << std::endl;
    std::cout << "MOVE by (" << stepSize << ") -> " << parentPtr << std::endl;
#endif //defined(DEBUG_JST_TRAVERSAL)
    stepSize = impl::shiftWindowBy(*parentPtr, stepSize);
    // we moved the context head -> might need to update the coverage and the
    impl::updateContextHead(me, *parentPtr);

    if (parentPtr->curEdgeIt == parentPtr->endEdgeIt)  // Reached branching point => expand node.
        stepSize = impl::expandNode(me, parentPtr, stepSize, observer, TProxySelector());

#if defined(DEBUG_JST_TRAVERSAL)
    std::cout << "Remaining (" << stepSize << ") -> " << parentPtr << std::endl;
    std::cout << "================================ END =========================================" << std::endl;
#endif //defined(DEBUG_JST_TRAVERSAL)
    return stepSize;
}

// ----------------------------------------------------------------------------
// Function init()
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec,
          typename TObserver,
          typename TProxySelector>
inline void
init(TraverserImpl<TJst, JstTraversalSpec<TSpec> > & me,
     TObserver & observer,
     Tag<TProxySelector> const & /*tag*/)
{
    typedef typename TraverserImpl<TJst, JstTraversalSpec<TSpec> >::TNode TNode;

    SEQAN_ASSERT(me._contPtr != nullptr);
    SEQAN_ASSERT(me._contextSize > 0);
    SEQAN_ASSERT(me._branchLength >= me._contextSize);

    SEQAN_ASSERT(me._bufferPtr);
    SEQAN_ASSERT(me._stackPtr);

    clear(buffer(me));
    clear(stack(me));

    init(buffer(me), container(me));

    resize(me._baseCov, length(container(me)), true, Exact());

    TNode node;
    node.coverage = me._baseCov;

    node.headSrcPos = -static_cast<decltype(node.headSrcPos)>(me._contextSize) + 1;
    node.headDelta = buffer(me)._deltaRangeBegin - 1;
    node.branchRoot = node.headDelta;

    node.curDelta = buffer(me)._deltaRangeBegin - 1;
    node.nextDelta = buffer(me)._deltaRangeBegin;
    node.headDelta = node.nextDelta;
    node.branchRoot = node.headDelta;

    node.proxyId = 0;

    node.begEdgeIt = begin(impl::member(container(me), JstSourceMember()), Standard());  // This points to some value already -> what could this position be?
    node.curEdgeIt = node.begEdgeIt;
    node.endEdgeIt = node.begEdgeIt + (getDeltaPosition(*node.nextDelta) - position(node.begEdgeIt));
    node.remainingSize = me._branchLength - 1;
    node.isBase = true;

    impl::pushNode(me, std::move(node), observer);  // Push onto stack.
    // After we realized this.
    TNode* basePtr = &impl::activeNode(me);

    SEQAN_ASSERT_GEQ(me._contextSize, 1u);
    if (IsSameType<Tag<TProxySelector>, SelectFirstProxy>::VALUE)
        impl::moveWindow(me, basePtr, 0, observer, Tag<TProxySelector>());   // We move the traverser to the first position.
    else
        impl::moveWindow(me, basePtr, me._contextSize - 1, observer, Tag<TProxySelector>());
}

// ----------------------------------------------------------------------------
// Function impl::swap()
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec,
          typename TOtherJst>
inline void
swap(TraverserImpl<TJst, JstTraversalSpec<TSpec> > & me,
     TraverserImpl<TOtherJst, JstTraversalSpec<TSpec> > & other)
{
    std::swap(me._contPtr, other._contPtr);
    std::swap(me._branchLength, other._branchLength);
    std::swap(me._contextSize, other._contextSize);
    std::swap(me._stackPtr, other._stackPtr);
    std::swap(me._bufferPtr, other._bufferPtr);
    swap(me._baseCov, other._baseCov);
    std::swap(me._needInitialization, other._needInitialization);
}

// ----------------------------------------------------------------------------
// Function impl::positionReference();
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec>
inline typename Position<TraverserImpl<TJst, JstTraversalSpec<TSpec> > >::Type
positionReference(TraverserImpl<TJst, JstTraversalSpec<TSpec> > const & me)
{
    // Easy case, since all sequences must be in a original node.
    using TPosVec = typename Position<TraverserImpl<TJst, JstTraversalSpec<TSpec> > >::Type;
    using TPosVal = typename Value<TPosVec>::Type;

    auto tmp = impl::baseNode(me).coverage;
    // Refine the coverage.
    impl::refineCoverage(tmp, impl::baseNode(me));

    TPosVec posVec;
    auto hostPos = position(impl::baseNode(me).curEdgeIt, host(container(me)));
    auto covBegin = begin(tmp, Standard());
    auto covEnd = end(tmp, Standard());

    for (auto it = covBegin; it != covEnd; ++it)
    {
        if (getValue(it))
        {
            auto seqId = it - covBegin;
            appendValue(posVec, TPosVal(seqId, hostToVirtualPosition(impl::buffer(me)._journaledSet[seqId], hostPos)));
        }
    }
    return posVec;
}

// ----------------------------------------------------------------------------
// Function impl::positionBranch();
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec>
inline typename Position<TraverserImpl<TJst, JstTraversalSpec<TSpec> > >::Type
positionBranch(TraverserImpl<TJst, JstTraversalSpec<TSpec> > const & me)
{
    using TPosVec = typename Position<TraverserImpl<TJst, JstTraversalSpec<TSpec> > >::Type;
    using TPosVal = typename Value<TPosVec>::Type;

    auto tmp = impl::activeNode(me).coverage;
    impl::refineCoverageBranch(tmp, impl::activeNode(me));

    TPosVec posVec;
    auto dist = impl::activeNode(me).curEdgeIt - impl::activeNode(me).begEdgeIt;
    auto covBegin = begin(tmp, Standard());
    auto covEnd = end(tmp, Standard());

    for (auto it = covBegin; it != covEnd; ++it)
    {
        if (getValue(it))
        {
            auto seqId = it - covBegin;
            auto tmpJournalIt = begin(buffer(me)._journaledSet[seqId], Standard());

            impl::mapBranchPointToVirtual(tmpJournalIt, impl::member(container(me), JstDeltaMapMember()), seqId,
                                          getDeltaPosition(*impl::activeNode(me).curDelta));
            appendValue(posVec, TPosVal(seqId, position(tmpJournalIt) + dist));
        }
    }
    return posVec;
}

}  // namespace impl

// Helper functions:

#if defined(JST_FIND_DEBUG)
    template <typename TTraverser>
    inline void _fillTestSet(TTraverser const & trav)
    {
        auto pos = position(trav);
        for (auto p : pos)
            appendValue(__testSet[p.i1], impl::buffer(trav)._journaledSet[p.i1][p.i2]);
    }

    inline void _printTestSet()
    {
        auto count = 0;
        for (auto seq : __testSet)
            std::cout << "Seq" << count++ << ": " << seq << std::endl;
    }
#endif // JST_FIND_DEBUG

}  // namespace seqan.

#endif  // #ifndef INCLUDE_SEQAN_JOURNALED_STRING_TREE_JOURNALED_STRING_TREE_TRAVERSER_UTIL_H_
