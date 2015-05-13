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

#ifndef INCLUDE_SEQAN_JOURNALED_STRING_TREE_TRAVERSAL_OPERATOR_H_
#define INCLUDE_SEQAN_JOURNALED_STRING_TREE_TRAVERSAL_OPERATOR_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================


namespace impl
{
// ----------------------------------------------------------------------------
// Struct DeltaMapIterLessByEndPosition
// ----------------------------------------------------------------------------

struct DeltaMapIterLessByEndPosition
{
    template <typename TIterator>
    bool operator()(TIterator itL, TIterator itR)
    {
        return deltaEndPosition(itL) < deltaEndPosition(itR);
    }
};
}

// ----------------------------------------------------------------------------
// Class JstTraversalOperator
// ----------------------------------------------------------------------------

template <typename TTraverser, typename TExtenal>
class JstTraversalOperator
{
public:
    typedef typename Container<TTraverser>::Type                            TDeltaMap;
    typedef typename Iterator<TDeltaMap, Standard>::Type                    TDeltaIter;
    typedef typename DeltaCoverage<TDeltaMap>::Type                         TCoverage;
    typedef typename Buffer<TTraverser>::Type                               TBuffer;
    typedef typename Member<TBuffer, JstSeqBufferJournaledSet>::Type        TJournaledSet;
    typedef typename Value<TJournaledSet>::Type                             TJournaledSeq;
    typedef JstTraversalEntry<TDeltaMap, TJournaledSeq, TState, TCoverage>  TEntry;

    typedef std::stack<TEntry>                                              TStack;
    typedef std::deque<TDeltaIter>                                          TAuxTable;

    TTraverser& traverser;          // Registered traverser.
    TExtenal&   extFunctor;         // Registered extension functor, called for every new position.

    TCoverage   backUpCoverge;
    TDeque      stack;              // State deque.
    TAuxTable   auxDeltaPoints;     // Current state over tracked deletions and SVs.

    // Custom c'tor.
    JstTraversalOperator(TTraverser & obj, TExternal & alg) :
        traverser(obj),
        externalFunctor(alg)
    {  // Create initial base entry.
        TEntry tmp;
        tmp.bpNext        = buffer(traverser).deltaRange.begin;
        tmp.bp            = tmp.bpNext;
        tmp.cur           = buffer(traverser).sourceRange.begin;  // need to get buffered range.
        tmp.end           = buffer(traverser).sourceRange.end;
        tmp.bpNextVirtual = deltaPosition(tmp.bpNext);
        arrayFill(begin(tmp.coverage), end(tmp.coverage), true);
        push(stack, SEQAN_MOVE(tmp));
        traverser.entryPtr = &top(stack);
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Private Functions
// ============================================================================

namespace impl
{

// ----------------------------------------------------------------------------
// Function impl::current()
// ----------------------------------------------------------------------------

template <typename TTraverser, typename TExternal>
inline typename Reference<Member<JstTraversalOperator<TTraverser, TExternal>, JstTraversalStackMember>::Type>::Type
current(JstTraversalOperator<TTraverser, TExternal> & op)
{
    return top(op.stack);
}

// ----------------------------------------------------------------------------
// Function impl::isBase()
// ----------------------------------------------------------------------------

template <typename TTraverser, typename TExternal>
inline bool impl::isBase(JstTraversalOperator<TTraverser, TExternal> const & op)
{
    return length(op.stack) == 1;
}

// ----------------------------------------------------------------------------
// Function impl::recordMergePoint()
// ----------------------------------------------------------------------------

template <typename TTraverser, typename TExternal, typename TDeltaIter, typename TDeltaType>
inline void
recordMergePoint(JstTraversalOperator<TTraverser, TExternal> & op,
                 TDeltaIter const & branchPoint)
{
    typedef JstTraversalOperator<TTraverser, TExternal>     TOperator;
    typedef typename TOperator::TAuxTable                   TAuxTable;
    typedef typename Iterator<TAuxTable, Standard>::Type    TIter;

    if (SEQAN_UNLIKELY(empty(op.auxDeltaPoints)))
        op.auxDeltaPoints.push_back(branchPoint);

    TIter it = std::lower_bound(begin(op.auxDeltaPoints, Standard()), end(op.auxDeltaPoints, Standard()),
                                branchPoint, DeltaMapIterLessByEndPosition());
    op.auxDeltaPoints.insert(it, branchPoint);
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
        TVarIterator itVar = begin(variantStore, Standard());
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
                                          DeltaMapEntryCompareLessByDeltaPosition_());

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
// Function impl::selectBPByDeltaType
// ----------------------------------------------------------------------------

template <typename TEntry, typename TIter, typename TTraverser, typename TExternal, typename TDeltaType>
inline void
selectBPByDeltaType(TEntry & /*entry*/,
                    TIter const & /*oldBp*/,
                    JstTraversalOperator<TTraverser, TExternal> const & /*op*/,
                    TDeltaType const & /*tag*/)
{
    // no-op.
}

template <typename TEntry, typename TIter, typename TTraverser, typename TExternal>
inline void
selectBPByDeltaType(TEntry & entry,
                    TIter const & oldBp,
                    JstTraversalOperator<TTraverser, TExternal> const & op,
                    DeltaTypeIns const & /*tag*/)
{
    entry.bpNextVirtual += deltaPosition(entry.bpNext) - deltaPosition(oldBp) + length(deltaValue(oldBp, DeltaTypeIns()));
}

template <typename TEntry, typename TIter, typename TTraverser, typename TExternal>
inline void
selectBPByDeltaType(TEntry & entry,
                    TIter const & oldBp,
                    JstTraversalOperator<TTraverser, TExternal> const & op,
                    DeltaTypeDel const & /*tag*/)
{
    typedef typename Container<TTraverser>::Type                TContainer;
    typedef typename Iterator<TContainer>::Type                 TIter;
    typedef typename DeltaValue<TContainer, DeltaTypeDel>::Type TDel;

    TDel del = deltaValue(oldBp, DeltaTypeDel());
    // Find first branch point which comes after the deletion of the current branch point.
    while (entry.bpNext != impl::branchPointEnd(op) && deltaPosition(entry.bpNext) < (deltaPosition(oldBp) + del))
        ++entry.bpNext;

    if (SEQAN_UNLIKELY(entry.bpNext == impl::branchPointEnd(op)))
    {
        entry.bpNextVirtual = length(container(entry.cur));
        return;
    }

    SEQAN_ASSERT_GEQ(deltaPosition(entry.bpNext) - deltaPosition(oldBp), del);
    entry.bpNextVirtual += deltaPosition(entry.bpNext) - deltaPosition(oldBp) - del;
}

template <typename TEntry, typename TIter, typename TTraverser, typename TExternal>
inline void
selectBPByDeltaType(TEntry & entry,
                    TIter const & oldBp,
                    JstTraversalOperator<TTraverser, TExternal> const & op,
                    DeltaTypeSV const & /*tag*/)
{
    typedef typename Container<TTraverser const>::Type          TContainer;
    typedef typename Iterator<TContainer>::Type                 TIter;
    typedef typename DeltaValue<TContainer, DeltaTypeSV>::Type  TSV;

    TSV & sv = deltaValue(oldBp, DeltaTypeSV());
    // Find first branch point which comes after the deletion of the current branch point.
    while (entry.bpNext != impl::branchPointEnd(op) && deltaPosition(entry.bpNext) < (deltaPosition(oldBp) + sv.i1))
        ++entry.bpNext;

    if (SEQAN_UNLIKELY(entry.bpNext == impl::branchPointEnd(op)))
    {
        entry.bpNextVirtual = length(container(entry.cur));
        return;
    }

    SEQAN_ASSERT_GEQ(deltaPosition(entry.bpNext) - deltaPosition(oldBp), sv.i1);
    entry.bpNextVirtual += deltaPosition(entry.bpNext) - deltaPosition(oldBp) - sv.i1 + length(sv.i2);
}

// ----------------------------------------------------------------------------
// Function impl::selectNextBranchPoint
// ----------------------------------------------------------------------------

template <typename TEntry, typename TTraverser, typename TExternal, typename TDeltaType>
inline void
selectNextBranchPoint(TEntry & entry, JstTraversalOperator<TTraverser, TExternal> & op, TDeltaType const & /*tag*/)
{
    typedef typename Container<TTraverser>::Type    TContainer;
    typedef typename Iterator<TContainer>::Type     TIter;

    TIter oldBp = entry.bpNext;
    // Find first branch point which is not at the same position.
    while (++entry.bpNext != impl::branchPointEnd(op) && deltaPosition(entry.bpNext) == deltaPosition(oldBp))
    {}

    if (SEQAN_UNLIKELY(entry.bpNext == impl::branchPointEnd(op)))
    {
        entry.bpNextVirtual = length(container(entry.cur));
        return;
    }

    // Assert if the next branch point is not larger than the previous one.
    SEQAN_ASSERT_GT(deltaPosition(entry.bpNext), deltaPosition(oldBp));
    selectBpByDeltaType(entry, oldBp, op, TDeltaType());

    // TODO(rrahn): Probably don't need it.
    //    if (SEQAN_UNLIKELY(entry.bpNext != impl::branchPointEnd(op) && deltaType(entry.bpNext) == DELTA_TYPE_DEL))
    //        ++entry.bpNextVirtual;  // Need to increase by one in case the next bp is a deletion.
}

// ----------------------------------------------------------------------------
// Function imp::selectProxy
// ----------------------------------------------------------------------------

template <typename TEntry, typename TTraverser, typename TExternal>
inline bool
selectProxy(TEntry & target,
            TEntry & source,
            JstTraversalOperator<TTraverser, TExternal> & op,
            TBool const & /*flag*/)
{
    if (testAllZero(target.coverage))
        return false;

    SEQAN_ASSERT_LT(bitScanForward(target.coverage), length(buffer(op.traverser).journaledSet));

    impl::mapBranchPointToVirtual(target.begBp, container(op.traverser), bitScanForward(target.coverage),
                                  deltaPosition(target.bp));  // Map the iterator to the beginning of this branch.
    target.cur = target.begBp + (source.cur - source.begBp);

    if (IsSameType<TBool, True>::VALUE)
    {
        target.bpNextVirtual = position(target.begBp);
        target.end = target.cur + windowSize(op.traverser);
    }
    else
    {
        target.bpNextVirtual = position(target.begBp) + (source.bpNextVirtual - position(source.begBp));
        target.end = target.cur + (source.end - source.cur);
    }
    return true;
}

// ----------------------------------------------------------------------------
// Function impl::updateBranch()
// ----------------------------------------------------------------------------

template <typename TEntry, typename TTraverser, typename TExternal, typename TBool>
inline void
impl::updateBranch(TEntry & target,
                   JstTraversalOperator<TTraverser, TExternal> & op,
                   TBool const & /*flag*/)
{
    if (SEQAN_LIKELEY(deltaType(child.bp) == DELTA_TYPE_SNP))  // Expect mostly SNPs.
    {
        impl::selectNextBranchPoint(child, op, DeltaTypeSnp());
    }
    else  // In case the branch point is another variant.
    {
        switch (deltaType(child.bp))
        {
            case DELTA_TYPE_DEL:
            {
                if (IsSameType<TBool, True>::VALUE)
                {
                    --child.end; // The iterator points behind the deletion, so we reduce the end by one.
                    impl::recordMergePoint(op, child.bp));
                }
                impl::selectNextBranchPoint(child, op, DeltaTypeDel());
                break;
            }
            case DELTA_TYPE_INS:
            {
                if (IsSameType<TBool, True>::VALUE)
                    child.end += length(deltaValue(child.bp, DeltaTypeIns()));  // The window searches the entire insertion including the source at which the insertion is recorded.
                impl::selectNextBranchPoint(child, op, DeltaTypeIns());
                break;
            }
            case DELTA_TYPE_SV:
            {
                if (IsSameType<TBool, True>::VALUE)
                {
                    child.end += (length(deltaValue(child.bp), DeltaTypeSV()).i2) - 1);  // The window searches the entire insertion part not including the position after the insertion which is apparently deleted.
                    impl::recordMergePoint(op, child.bp);
                }
                impl::selectNextBranchPoint(child, op, DeltaTypeSV());
                break;
            }
        }
    }
}

// ----------------------------------------------------------------------------
// Function impl::mappedSourcePosition()
// ----------------------------------------------------------------------------

typename TTraverser, typename TExternal>
inline typename Position<TTraverser>::Type
mappedSourcePosition(JstTraversalOperator<TTraverser, TExternal> const & op)
{


}

// ----------------------------------------------------------------------------
// Function impl::syncCoverageWithMergePoints()
// ----------------------------------------------------------------------------

template <typename TTraverser, typename TExternal, typename TBool>
inline void
syncCoverageWithMergePoints(JstTraversalOperator<TTraverser, TExternal> & op, TBool const & /*flag*/)
{
    typedef typename Container<TTraverser>::Type            TDeltaMap;
    typedef typename Value<TDeltaMap>::Type                 TEntry;
    typedef typename Iterator<TDeltaMap, Standard>::Type    TMapIterator;
    typedef typename Position<TTraverser>::Type TPos;
    // Implement logic.
    // We wan't to find the
    if (empty(op.auxDeltaPoints))
        return;  // Nothing to do.

    TPos pos;
    if (impl::isBase(op))
        pos = position(windowBegin(op.traverser));
    else
    {
        TPos pos = position(windowBegin(op.traverser));
        if (pos >= position(impl::current(op).begBp))
            // In this case the coverage is at least the coverage of the branch point.
    }


    // Otherwise check all stored merge points that overlap with the current window and update the coverage.
    // If flag is true, remove all previous ones.

    TEntry entry;
    entry.deltaPosition = impl::mappedSourcePosition(op);
    TMapIterator it = std::lower_bound(buffer(op.traverser).deltaRange.begin, buffer(op.traverser).deltaRange.end,
                                       entry, f);
    // A) find first branchNode within window.
    // B) find first mergePoint that overlaps.
    // C) create a coverage that contains all recorded ones
    // D) Update current coverage with this one.
    // E) If flag true -> remove all previous ones.
}
    
}  // namespace impl

// ============================================================================
// Public Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function advance()
// ----------------------------------------------------------------------------

template <typename TTraverser, typename TExternal>
inline void
advance(JstTraversalOperator<TTraverser, TExternal> & op, Forward const & /*dir*/)
{
    typedef JstTraversalOperator<TTraverser, TExternal> TOperator;
    typedef typename Member<TOperator, JstTraversalStackMember>::Type   TStack;
    typedef typename Value<TStack>::Type                                TEntry;
    typedef typename Size<TTraverser>::Type                             TSize;

    TEntry & entry = impl::current(op);

    op.backUpCoverge = entry.coverage;
    if (impl::isBase(op))
        impl::syncCoverageWithMergePoints(op, True());
    else
        impl::syncCoverageWithMergePoints(op, False());
    entry.currentIt += apply(op.extension, op.traverser);  // Move window by returned step size from the external algorithm.
    // sync2 if we left previous delta nodes.
    swap(entry.coverage, op.backUpCoverge);
}

// ----------------------------------------------------------------------------
// Function expand()
// ----------------------------------------------------------------------------

template <typename TTraverser, typename TExternal>
inline void
expand(JstTraversalOperator<TTraverser, TExternal> & op)
{
    // TODO(rrahn): Need to cover empty coverages due to deletions and SVs, as well as consecutive delta operations.
    // In this case one cannot find a proxy target so one has to continue the traversal until the coverage is not empty.
    // TODO(rrahn): Check end of strings.

    // Add new branch.
    typedef JstTraversalOperator<TTraverser, TExternal> TOperator;
    typedef typename Member<TOperator, JstTraversalStackMember>::Type   TStack;
    typedef typename Value<TStack>::Type                                TEntry;

    TEntry& parent = impl::current(op);
    parent.state = getState(op.extFunctor);  // Save last state to continue traversal when coming back.

    do  // Parse all branch points that are at the same location within the parent branch.
    {
        TEntry child;
        child.bpNext = parent.bpNext;
        // The new branch results from the current cov & cov'branchPoint.
        transform(child.coverage, parent.coverage, getDeltaCoverage(parent.bpNext), FunctorBitwiseAnd());
        // Update parent coverage to represent all but the sequences that support the current branch point.
        transform(parent.coverage, parent.coverage, getDeltaCoverage(parent.bpNext),
                  FunctorNested<FunctorBitwiseAnd, FunctorIdentity, FunctorBitwiseNot>());

        if (impl::isBase(op))  // Spawn off an intial branch from the source.
        {
            child.bp = child.bpNext;  // Set the initial branchPoint for this branch.
            SEQAN_ASSERT_GEQ(position(parent.cur), deltaPosition(parent.bpNext));
            parent.begBp = parent.cur - (position(parent.cur) - deltaPosition(parent.bpNext));
            impl::selectProxy(child, parent, op, True());
            impl::updateProxy(child, op, True());
        }
        else  // Spawn off another branch within the current subtree.
        {
            child.bp = parent.bp;

            unsigned proxyId = bitScanForward(child.coverage);
            if (&buffer(op.traverser).journaledSet[proxyId] == &container(parent.cur))
            {  // The current proxy of the parent is the new proxy of the child. In that case, we need to find a new proxy for the parent and only update the child iterators.

                // Switch child and parent.
                child.bpNextVirtual = parent.bpNextVirtual;
                child.cur = parent.cur;
                child.end = parent.end;
                child.begBp = parent.begBp;

                impl::selectProxy(parent, child, op, False());
            }
            else
            {  // parent proxy does not support new delta.
                impl::selectProxy(child, parent, op, False());
            }
        }
        impl::updateProxy(child, op, False());
        child.state = parent.state;  // Transfer the current state to the child branch.
        push(op.stack, SEQAN_MOVE(child));
    } while (deltaPosition(parent.bpNext++) == deltaPosition(parent.bpNext))
        // Now they are different and we know that the parent does not include a variant.
        parent.bpNextVirtual += deltaPosition(parent.bpNext) - deltaPosition(parent.bpNext - 1);
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_JOURNALED_STRING_TREE_TRAVERSAL_OPERATOR_H_
