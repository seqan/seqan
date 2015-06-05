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

#ifndef INCLUDE_SEQAN_JOURNALED_STRING_TREE_JOURNALED_STRING_TREE_TRAVERSOR_H_
#define INCLUDE_SEQAN_JOURNALED_STRING_TREE_JOURNALED_STRING_TREE_TRAVERSOR_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TSpec = void>
struct JstTraversalSpec
{};

template <typename TJst, typename TSpec>
class TraversorImpl<TJst, JstTraversalSpec<TSpec> >
{
public:
    typedef JstTraversalNode<TJst>              TTraversalNode;
    typedef String<TTraversalNode, Block<> >    TTraversalStack;
    typedef typename Size<TraversorImpl>::Type  TSize;

    // __ Member Variables ____________________________________________________

    TJst *          _contPtr;
    TSize           _contextSize;
    TTraversalNode  _tmp;
    TTraversalStack _stack;           // Parent stack used for iterative traversal.

    // __ Constructors ________________________________________________________

    TraversorImpl() : _contPtr(nullptr), _contextSize(1)
    {}

    TraversorImpl(TJst & jst) : _contPtr(nullptr), _contextSize(1), _stack()
    {
        _init(*this, jst);
    }

    // Copy Constructor!
    template <typename TOtherJst>
    TraversorImpl(TraversorImpl<TOtherJst, JstTraversalSpec<TSpec> > const & other,
         SEQAN_CTOR_ENABLE_IF(IsConstructible<TJst, TOtherJst>)) :
            _contPtr(other._contPtr),
            _contextSize(other._contextSize),
            _tmp(other._tmp),
            _stack(other._stack)
    {
        ignoreUnusedVariableWarning(dummy);
    }

    // __ Member Functions ____________________________________________________

    template <typename TOtherJst>
    inline SEQAN_FUNC_ENABLE_IF(IsConstructible<TJst, TOtherJst>, TraversorImpl &)
    operator=(TraversorImpl<TOtherJst, JstTraversalSpec<TSpec> > const & other)
    {
        _contPtr = other._contPtr;
        _contextSize = other._contextSize;
        _tmp = other._tmp;
        _stack = other._stack;
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

template <typename TJst, typename TSpec>
struct Container<TraversorImpl<TJst, JstTraversalSpec<TSpec> > >
{
    typedef TJst Type;
};

template <typename TJst, typename TSpec>
struct Size<TraversorImpl<TJst, JstTraversalSpec<TSpec> > >
{
    typedef typename Size<TJst>::Type    Type;
};


template <typename TObject>
struct StringContext;

template <typename TJst, typename TSpec>
struct StringContext<TraversorImpl<TJst, JstTraversalSpec<TSpec> > >
{
    typedef TraversorImpl<TJst, JstTraversalSpec<TSpec> >   TThis_;
    typedef typename Container<TThis_>::Type                TContainer_;
    typedef typename Source<TContainer_>::Type              TSource_;
    typedef typename Iterator<TSource_, Standard>::Type     TSourceIt_;

    typedef Range<TSourceIt_>                               Type;
};

// ============================================================================
// Private Functions
// ============================================================================

namespace impl
{

// ----------------------------------------------------------------------------
// Function init()
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec>
inline void
init(Iter<TJst, JstTraversalSpec<TSpec> > & me,
     TJst & jst)
{
    typedef typename Iter<TJst, JstTraversalSpec<TSpec> >::TTraversalNode TNode;

    me._contPtr = &jst;
    TNode node;
    node.seqBegin = sourceBegin(jst);
    node.seqEnd = node.seqBegin;
    node.mappedSrcPos = -1;
    resize(node.coverage, dimension(jst), true, Exact());
    node.curDelta = jst._buffer._deltaRangeBegin - 1;  // Point before beginning.
    appendValue(me._history, SEQAN_MOVE(node));  // Push onto stack.
}

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
        SEQAN_ASSERT_LEQ(deltaPosition(itVar), static_cast<TDeltaMapPos const>(hostPos));

        TDeltaMapPos virtualOffset = 0;
        // Now we move to the right until we find the node that we are looking for and reconstruct the offset of the virtual positions.
        while(deltaPosition(itVar) != static_cast<TDeltaMapPos const>(hostPos) && !atEnd(itVar, variantStore))
        {
            if (deltaCoverage(itVar)[proxyId] != true)  // irrelevant variant.
            {
                ++itVar;
                continue;
            }

            if (deltaType(itVar) == DELTA_TYPE_INS)
                virtualOffset += length(deltaValue(itVar, DeltaTypeIns()));
            else if (deltaType(itVar) == DELTA_TYPE_SNP)
                ++virtualOffset;
            else if (deltaType(itVar) == DELTA_TYPE_SV)
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

    SEQAN_ASSERT_LEQ(deltaPosition(itVar), static_cast<TDeltaMapPos const>(hostPos));

    TDeltaMapPos virtualOffset = 0;
    // Now we move to the right until we find the node that we are looking for and reconstruct the offset of the virtual positions.
    while (deltaPosition(itVar) != static_cast<TDeltaMapPos const>(hostPos) && !atEnd(itVar, variantStore))
    {
        if (deltaCoverage(itVar)[proxyId] != true)  // irrelevant variant.
        {
            ++itVar;
            continue;
        }
        
        if (deltaType(itVar) == DELTA_TYPE_INS)
            virtualOffset += length(deltaValue(itVar, DeltaTypeIns()));
        else if (deltaType(itVar) == DELTA_TYPE_SNP)
            ++virtualOffset;
        else if (deltaType(itVar) == DELTA_TYPE_SV)
            virtualOffset += length(deltaValue(itVar, DeltaTypeSV()).i2);
        ++itVar;
    }
    resultIt += virtualOffset + 1;  // Set the iterator to the beginning of the variant.
}


// ----------------------------------------------------------------------------
// Function impl::inBranch()
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec>
inline bool
inBranch(Iter<TJst, JstTraversalSpec<TSpec> > const & me)
{
    return length(me._history) > 1;
}

// ----------------------------------------------------------------------------
// Function impl::isMergePoint()
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec>
inline bool
isMergePoint(Iter<TJst, JstTraversalSpec<TSpec> > const & me)
{
    return !impl::inBranch(me) && back(me._history).mappedSrcEndPos != deltaPosition(back(me._history).curDelta);
}

// ----------------------------------------------------------------------------
// Function impl::nextChild()
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec,
          typename TTraversalNode>
inline bool
nextChild(Iter<TJst, JstTraversalSpec<TSpec> > & it,
         TTraversalNode & parent,
         TTraversalNode & child)
{
    // We set the coverage of the left child to be the one of the parent & coverage(curDelta);
    transform(child.coverage, parent.coverage, getDeltaCoverage(*parent.nextDelta), FunctorBitwiseAnd());
    transform(parent.coverage, parent.coverage, getDeltaCoverage(*parent.nextDelta),
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
            impl::mapBranchPointToVirtual(parent.endEdgeIt, host(container(it)), proxyId, getDeltaPosition(*parent.nextDelta));
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
        impl::mapBranchPointToVirtual(child.endEdgeIt, host(container(it)), proxyId, getDeltaPosition(*child.curDelta));
        child.begEdgeIt = child.endEdgeIt - endToBegDist;
        child.curEdgeIt = child.endEdgeIt;
    }

    // Move to next valid delta position, note that we skip all end points and all deltas that occur at the same position.
    while (getDeltaPosition(*(++child.nextDelta)) == getDeltaPosition(*child.curDelta) || child.nextDelta.isEndPoint())
    {}

    // C) Update remaining length if we come directly from the base sequence.
    if (parent.isBase)  // TODO(rmaerker): Make sure that the base node always has the complete window size as remaining Size.
    {
        // We should keep the mappedSrcEndPos of the base node already in sync.
        child.mappedSrcEndPos += historySize(container(it));
        switch (getDeltaType(*child.curDelta))
        {
            case DELTA_TYPE_DEL:
            {
                child.mappedSrcEndPos += deletionSize(host(container(it))._deltaStore, getStorePosition(*child.curDelta), DeltaTypeDel());
                --child.remainingSize;
                break;
            }
            case DELTA_TYPE_INS:
            {
                child.remainingSize += insertionSize(host(container(it))._deltaStore, getStorePosition(*child.curDelta), DeltaTypeIns());
                break;
            }
            case DELTA_TYPE_SV:
            {
                child.mappedSrcEndPos += deletionSize(host(container(it))._deltaStore, getStorePosition(*child.curDelta), DeltaTypeSV());
                child.remainingSize += insertionSize(host(container(it))._deltaStore, getStorePosition(*child.curDelta), DeltaTypeSV()) - 1;
                break;
            }
        }
    }
    else
    {
        // Now we are down somewhere in an internal subtree.
        switch (getDeltaType(*child.curDelta))
        {
            case DELTA_TYPE_DEL:
            {
                child.mappedSrcEndPos += deletionSize(host(container(it))._deltaStore, getStorePosition(*child.curDelta), DeltaTypeDel());
                break;
            }
            case DELTA_TYPE_SV:
            {
                child.mappedSrcEndPos += deletionSize(host(container(it))._deltaStore, getStorePosition(*child.curDelta), DeltaTypeSV());
                break;
            }
        }
    }

    child.begEdgeIt = child.endEdgeIt;
    child.curEdgeIt = child.begEdgeIt;
    child.endEdgeIt = getDeltaPosition(*child.nextDelta) - getDeltaPosition(*child.curDelta);
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

    // TODO(rrahn): Check if we need this?
    arrayFill(begin(base.cover, Standard()), end(base.coverage, Standard()), true);  // Make sure the coverage is set to 1.
    while (base.curDelta.getDeltaPosition() == base.nextDelta.getDeltaPosition())
    {
        if (SEQAN_UNLIKELY(getDeltaType(*base.curDelta) == DELTA_TYPE_DEL || getDeltaType(*base.curDelta) == DELTA_TYPE_SV))
            transform(base.coverage, base.coverage, getDeltaCoverage(*base.curDelta),
                      FunctorNested<FunctorBitwiseAnd, FunctorIdentity, FunctorBitwiseAnd>());
        ++base.nextDelta;
    }
    // base.nextDelta might point to an merge point.
    base.begEdgeIt += base.curDelta.getDeltaPosition() - position(base.begEdgeIt);
    base.curEdgeIt = base.begEdgeIt + 1;
    base.endEdgeIt += base.nextDelta.getDeltaPosition() - base.mappedSrcEndPos;
    base.mappedSrcEndPos = base.nextDelta.getDeltaPosition();
}

// ----------------------------------------------------------------------------
// Function impl::updateParent()
// ----------------------------------------------------------------------------

template <typename TTraversalNode>
inline bool
updateParent(TTraversalNode & parent)
{
    SEQAN_ASSERT_NOT(parent.isBase);

    if (setAllZeros(parent.coverage))
        return false;
    // Update:
    parent.begEdgeIt = parent.endEdgeIt;
    parent.curEdgeIt = parent.begEdgeIt;

    // Make sure that for a right subtree we don't consider merge points, because we cannot be within a branch
    // coming from a delta while at the same time there is a deletion to be covered.
    while(isMergePoint(++parent.nextDelta))  // Move the nextDelta iterator to the next branch point.
    {}
    parent.mappedSrcEndPos = getDeltaPosition(parent.nextDelta);
    parent.endEdgeIt += parent.mappedSrcEndPos - getDeltaPosition(*parent.curNode);
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
        if (stepSize > node.remainingSize)
        {
            node.remainingSize = 0;
            node.curEdgeIt += node.remainingSize;
            return stepSize - node.remainingSize;
        }
        node.curEdgeIt += stepSize;
        node.remainingSize -= stepSize;
        return 0;
    }

    if (stepSize > node.remainingSize)
    {
        node.remainingSize = 0;
        node.curEdgeIt += node.remainingSize;
        return stepSize - node.remainingSize;
    }
    node.curEdgeIt = node.endEdgeIt;
    node.remainingSize -= remainingEdgeSize;
    return stepSize - remainingEdgeSize;
}

template <typename TJst, typename TSpec,
          typename TTraversalNode>
inline void
pushNode(Iter<TJst, JstTraversalSpec<TSpec> > & traversor,
         TTraversalNode SEQAN_FORWARD_CARG child)
{
    appendValue(traversor._stack, SEQAN_FORWARD(TTraversalNode, child));
}

template <typename TJst, typename TSpec>
inline void
popNode(Iter<TJst, JstTraversalSpec<TSpec> > & traversor)
{
    eraseBack(traversor._stack);
}

// ----------------------------------------------------------------------------
// Function impl::getStringContext()
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec>
inline typename StringContext<TraversorImpl<TJst, JstTraversalSpec<TSpec> > >::Type
getStringContext(TraversorImpl<TJst, JstTraversalSpec<TSpec> > const & traversor)
{
    // TODO(rrahn): Write me!
}

}  // namespace impl

template <typename TJst, typename TSpec>
inline void
_init(Iter<TJst, JstTraversalSpec<TSpec> > & me,
      TJst & jst)
{
    impl::init(me, jst);
}

// ============================================================================
// Public Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function stringContext()
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec>
inline typename StringContext<TraversorImpl<TJst, JstTraversalSpec<TSpec> > >::Type
stringContext(TraversorImpl<TJst, JstTraversalSpec<TSpec> > const & traversor)
{
    return impl::getStringContext(traversor);
}

// ----------------------------------------------------------------------------
// Function container
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec>
inline typename Container<TraversorImpl<TJst, JstTraversalSpec<TSpec> > >::Type &
container(TraversorImpl<TJst, JstTraversalSpec<TSpec> > & me)
{
    return *me._contPtr;
}

template <typename TJst, typename TSpec>
inline typename Container<TraversorImpl<TJst, JstTraversalSpec<TSpec> > const>::Type &
container(TraversorImpl<TJst, JstTraversalSpec<TSpec> > const & me)
{
    return *me._contPtr;
}

// ----------------------------------------------------------------------------
// Function atBegin()
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec>
inline bool
atBegin(Iter<TJst, JstTraversalSpec<TSpec> > & me)
{
    return back(me._stack).mappedSrcEndPos == -1;
}

template <typename TJst, typename TSpec,
          typename TTraversalNode,
          typename TSize>
inline TSize
moveWindow(Iter<TJst, JstTraversalSpec<TSpec> > & it,
           TTraversalNode & parent,
           TSize stepSize)
{
    stepSize = impl::moveWindow(parent, stepSize);
    if (stepSize > 0)
    {
        if (parent.remainingSize > 0)
        {
            parent.curDelta = parent.nextDelta;
            if (parent.isBase)
                it._tmp = parent;  // Store copy of the original base node.

            while (getDeltaPosition(*parent.nextDelta) == getDeltaPosition(*parent.curDelta))
            {
                auto child = parent;
                if (impl::nextChild(it, parent, child) && moveWindow(it, child, stepSize) == 0 &&
                    child.remainingSize >= 0)
                    impl::pushNode(it, SEQAN_MOVE(child));
                ++parent.nextDelta;  // Move to the next position.
            }

            if (parent.isBase)
            {
                // We need to move the base parent to the next segment redarding the branch information.
                impl::advanceBaseParent(it, it._tmp);
                swap(parent, it._tmp);  // Set advanced base parent back to stack at the correct position.
                it._tmp.isBase = false;
                it._tmp.fromBase = true;
                impl::updateParent(it, it._tmp);
                if (moveWindow(it, it._tmp, stepSize) == 0 && it._tmp.remainingSize >= 0)
                    impl::pushNode(it, it._tmp);
            }
            else
            {
                if (impl::updateParent(it, parent))
                    return moveWindow(it, parent, stepSize);
            }
        }
    }
    return stepSize;
}

template <typename TJst, typename TSpec, typename TSize>
inline void
goNext(Iter<TJst, JstTraversalSpec<TSpec> > & it,
       TSize stepSize)
{
    auto& node = back(it._stack);
    stepSize = moveWindow(it, node, stepSize);

    // Case A) stepSize == 0 && node.remainingSize > 0;  // vaild context
    // Case B) stepSize == 0 && node.remainingSize = 0;  // valid context
    // Case C) stepSize > 0 && node.remaingingSize = 0;  // invalid context -> we have to pop an element from the tree.
    // Case D) stepSize > 0 && node.remainingSize > 0;   // invalid case: ASSERT

    SEQAN_ASSERT(stepSize > static_cast<TSize>(0) && node.remainingSize >static_cast<TSize>(0));

    if (stepSize > 0 && node.remainingSize == 0)
    {
        impl::popNode(it);
    }
}

// ----------------------------------------------------------------------------
// Function goNext()
// ----------------------------------------------------------------------------

// Standard version just iterating over the delta nodes.
template <typename TJst, typename TSpec>
inline void
goNext(Iter<TJst, JstTraversalSpec<TSpec> > & it)
{
    goNext(it, 1);
}

// ----------------------------------------------------------------------------
// Function atEnd()
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec>
inline bool
atEnd(Iter<TJst, JstTraversalSpec<TSpec> > & me)
{
    return length(me._stack) == 1 && back(me._stack).srcEnd == sourceEnd(container(me));
}



}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_JOURNALED_STRING_TREE_JOURNALED_STRING_TREE_TRAVERSOR_H_
