// ==========================================================================
// traverser_context.h
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

#ifndef INCLUDE_SEQAN_JOURNALED_STRING_TREE_TRAVERSER_CONTEXT_H_
#define INCLUDE_SEQAN_JOURNALED_STRING_TREE_TRAVERSER_CONTEXT_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// TODO(rmaerker): Docu!
template <typename TObject, typename TSpec = void>
struct TraverserContext;

template <typename TReference, typename TSpec, typename TConfig, typename TSpec>
struct TraverserContext<JournaledStringTree<TReference, TSpec, TConfig>, TSpec>
{
    typedef JournaleStringTree<TReference, TSpec, TConfig>          TJst;
    typedef typename Member<TJst, JstReference>::Type               TReference;
    typedef Range<typename Iterator<TReference, Standard>::Type >   TRefRange;
    typedef String<TReference, JournaledString<> >                  TContextDataValue;
    typedef StringSet<TContextString, Owner<JournaledSet> >         TContextData;
    typedef typename Position<TContextDataValue>::Type              TPos;
    typedef String<TPos>                                            TOffsetString;

    // Member variables.

    TJst const *    _jstPtr;   // Pointer to the underlying journaled string tree.
    TRefRange       _refRange; // Range over the reference defining the context view.
    TContextData    _data;     // The object holding the data in this context.
    TOffsetString   _offsets;  // Keeps track of offsets in case of dynamic context generation.

    // Member functions.

    TraverserContext() : _jstPtr(nullptr)
    {}

    TraverserContext(TJst const & jst) : _jstPtr(&jst)
    {
        SEQAN_ASSERT_NOT(empty(host(jst)));
        toRange(_refRange, begin(host(jst)), end(host(jst)));

        auto res = impl::fetchContextDataForward(*this);
        SEQAN_ASSERT(res);
    }

    template <typename TSize>
    TraverserContext(TJst const & jst, TSize blockSize) : _jstPtr(&jst)
    {
        SEQAN_ASSERT_NOT(empty(host(jst)));
        toRange(_refRange, begin(host(jst)), begin(host(jst)) + blockSize);

        auto res = impl::fetchContextDataForward(*this);
        SEQAN_ASSERT(res);
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

template <typename THostData, typename TSpec>
inline bool fetchContextDataForward(TraverserContext<THostData, TSpec> & me)
{
    // TODO(rmaerker): Here we load the context data in forward direction.
}

template <typename THostData, typename TSpec>
inline bool fetchContextDataBackward(TraverserContext<THostData, TSpec> & me)  // Only for specific type.
{
    // TODO(rmaerker): Here we load the context data in backward direction.
}

}

// ============================================================================
// Public Functions
// ============================================================================

template <typename TReference, typename TSpec, typename TConfig, typename TSpec, typename TSize>
inline bool fetchNext(TraverserContext<JournaledStringTree<TReference, TSpec, TConfig>, TSpec> & me,
                      TSize blockSize)
{
    // TODO(rmaerker): Implement me!
    return fetchContextDataForward(me);
}

// TODO(rmaerker): Disable for default JST
template <typename TReference, typename TSpec, typename TConfig, typename TSpec, typename TSize>
inline bool fetchPrevious(TraverserContext<JournaledStringTree<TReference, TSpec, TConfig>, TSpec> & me,
                          TSize blockSize)
{
    // TODO(rmaerker): Implement me!
    return fetchContextDataBackward(me);
}

template <typename TReference, typename TSpec, typename TConfig, typename TSpec, typename TOffset, typename TSize>
inline bool fetchUpstream(TraverserContext<JournaledStringTree<TReference, TSpec, TConfig>, TSpec> & me,
                          TOffset offset,
                          TSize blockSize)
{
    // TODO(rmaerker): Implement me!
    // Fetches block upstream to current range end position by offset.
}

// TODO(rmaerker): Disable for default JST
template <typename TReference, typename TSpec, typename TConfig, typename TSpec, typename TOffset, typename TSize>
inline bool fetchDownstream(TraverserContext<JournaledStringTree<TReference, TSpec, TConfig>, TSpec> & me,
                            TOffset offset,
                            TSize blockSize)
{
    // TODO(rmaerker): Implement me!
    // Fetches block downstream to current range begin position by offset.
}

template <typename TReference, typename TSpec, typename TConfig, typename TSpec, typename TPosition>
inline TPosition toGlobalPos(TraverserContext<JournaledStringTree<TReference, TSpec, TConfig>, TSpec> const & me,
                             TPosition localPos)
{
    // TODO(rrahn): Implement me!
}

    // TODO(rrahn): This is going to be part of the traversal buffer interface.
    // ----------------------------------------------------------------------------
    // Function _bufferJournaledStringEnds()
    // ----------------------------------------------------------------------------

    template <typename TDeltaMap, typename TSpec, typename TSize, typename TPos>
    inline void
    _bufferJournaledStringEnds(JournaledStringTree<TDeltaMap, TSpec> const & jst,
                               TSize contextSize,
                               TPos jobBegin,
                               TPos jobEnd)
    {
        typedef JournaledStringTree<TDeltaMap, TSpec> const TJst;
        typedef BufferJournaledStringsEndHelper_<TJst, TSize> THelper;

        THelper helper(jst, contextSize);
        // Buffer the next branch nodes depending on the context size.
        for (helper.jsPos = jobBegin; helper.jsPos < jobEnd; ++helper.jsPos)
        {
            // Step 1: Store the current offset of all virtual positions.
            jst._activeBlockVPOffset[helper.jsPos] = length(jst._journalSet[helper.jsPos]) - length(host(jst._journalSet));

            // Step 2: Check if deltas need to be buffered at the end of the current block for some j'strings.
            helper.localContextEndPos =
            getValue(end(value(stringSet(jst), helper.jsPos)._journalEntries, Standard()) - 1).physicalOriginPosition +
            contextSize - 1;
            if (jst._mapBlockEnd != end(container(jst), Standard()) &&
                helper.localContextEndPos > deltaPosition(jst._mapBlockEnd))
            {
                helper.mapIter = jst._mapBlockEnd;
                while (helper.mapIter != end(container(jst), Standard()) &&
                       deltaPosition(helper.mapIter) < helper.localContextEndPos)
                {
                    tagApply(helper, DeltaTypes());
                    ++helper.mapIter;
                }
            }
        }
    }

    // ----------------------------------------------------------------------------
    // Function _doJournalBlock()
    // ----------------------------------------------------------------------------

    template <typename TDeltaMap, typename TSpec, typename TContextSize, typename TParallelTag>
    inline bool
    _doJournalBlock(JournaledStringTree<TDeltaMap, TSpec> const & jst,
                    TContextSize contextSize,
                    Tag<TParallelTag> parallelTag = Serial())
    {
        typedef JournaledStringTree<TDeltaMap, TSpec > TJst;
        typedef typename Iterator<TDeltaMap, Standard>::Type TMapIterator;
        typedef typename GetStringSet<TJst>::Type TJournalSet;
        typedef typename Iterator<TJournalSet, Standard>::Type TJournalSetIter;
        typedef typename Size<TJournalSet>::Type TSize;
        typedef typename MakeSigned<TSize>::Type TSignedSize;

        typedef JournalDeltaContext_<TJst const, TMapIterator> TJournalDeltaContext;

        // Define the block limits.
        if ((jst._activeBlock * jst._blockSize) >= length(container(jst)))
            return false;

        unsigned lastBlock = jst._activeBlock * jst._blockSize;
        TSize blockJump = _min(length(container(jst)), lastBlock + jst._blockSize) - lastBlock;

        // Auxiliary variables.
        TDeltaMap const & variantMap = container(jst);

        String<int> _lastVisitedNodes;
        if (!fullJournalRequired(jst))
            blockJump -= _max(0, (static_cast<TSignedSize>(jst._mapBlockEnd - jst._mapBlockBegin) -
                                  static_cast<TSignedSize>(jst._blockSize)));

        jst._mapBlockBegin = jst._mapBlockEnd;
        jst._mapBlockEnd += blockJump;

        for (; jst._mapBlockEnd != end(variantMap, Standard()) &&
             deltaPosition(jst._mapBlockEnd) == deltaPosition(jst._mapBlockEnd - 1); ++jst._mapBlockEnd)
        {}

        if (jst._mapBlockBegin == jst._mapBlockEnd)
            return false;

        // Use parallel processing.
        Splitter<TJournalSetIter> jSetSplitter(begin(jst._journalSet, Standard()), end(jst._journalSet, Standard()),
                                               parallelTag);

        SEQAN_OMP_PRAGMA(parallel for)
        for (int jobId = 0; jobId < static_cast<int>(length(jSetSplitter)); ++jobId)
        {
            //        printf("Thread: %i of %i\n", jobId, omp_get_num_threads());
            unsigned jobBegin = jSetSplitter[jobId] - begin(jst._journalSet, Standard());
            unsigned jobEnd = jSetSplitter[jobId + 1] - begin(jst._journalSet, Standard());

            // Pre-processing: Update VPs for last block.
            if (!fullJournalRequired(jst))
            {
                for (unsigned i = jobBegin; i < jobEnd; ++i)
                {
                    clear(jst._journalSet[i]);  // Reinitialize the journal strings.
                    jst._blockVPOffset[i] += jst._activeBlockVPOffset[i];
                }
            }
            
            TJournalDeltaContext context(jst, jst._mapBlockBegin, jobBegin, jobEnd);
            for (; context.deltaIterator != jst._mapBlockEnd; ++context.deltaIterator)
                tagApply(context, DeltaTypes());
            
            // Post-processing: Store VPs for current block.
            if (!fullJournalRequired(jst))
                _bufferJournaledStringEnds(jst, contextSize, jobBegin, jobEnd);
        }
        return true;
    }


    // ----------------------------------------------------------------------------
    // Function localToGlobalPos()
    // ----------------------------------------------------------------------------

    /*!
     * @fn JournaledStringTree#localToGlobalPos
     * @headerfile <seqan/journaled_string_tree.h>
     * @brief Returns the global position of the given local position for the corresponding
     *        @link JournaledStringTree @endlink.
     *
     * @signature TSize localToGlobalPos(pos, id, jst);
     *
     * @param[in] pos    The local position.
     * @param[in] id     The id of the sequence to get the offset for.
     * @param[in] jst    The Journal String Tree.
     *
     * If the JST is constructed block wise, than the journaled strings only represent the block that is currently
     * active. To get the corresponding global position of the @link JournaledString @endlink respectively, one has
     * to call this function with the local position the sequence id and the @link JournaledStringTree @endlink (JST).
     * If the JST is constructed entirely, than the returned position corresponds to the local position.
     *
     * @return TSize The global position for the given @link JournaledString @endlink of
     type @link JournaledStringTree#Size @endlink.
     */

    template <typename TPos, typename TSeqId, typename TDeltaMap, typename TSpec>
    inline typename Size<JournaledStringTree<TDeltaMap, TSpec> const>::Type
    localToGlobalPos(TPos pos,
                     TSeqId seqId,
                     JournaledStringTree<TDeltaMap, TSpec> const & stringTree)
    {
        SEQAN_ASSERT_LT(seqId, static_cast<TSeqId>(length(stringSet(stringTree))));

        return pos + stringTree._blockVPOffset[seqId];
    }

    // ----------------------------------------------------------------------------
    // Function fullJournalRequired()
    // ----------------------------------------------------------------------------

    /*!
     * @fn JournaledStringTree#fullJournalRequired
     * @headerfile <seqan/journaled_string_tree.h>
     * @brief Checks whether the sequences are constructed block-wise or not.
     *
     * @signature bool fullJournalRequired(jst);
     *
     * @param[in] jst    The Journal String Tree.
     *
     * @return bool <tt>false</tt> if constructed block-wise, <tt>true</tt> otherwise.
     */

    template <typename TDeltaMap, typename TSpec>
    inline bool
    fullJournalRequired(JournaledStringTree<TDeltaMap, TSpec> const & stringTree)
    {
        typedef JournaledStringTree<TDeltaMap, TSpec> TJst;

        return stringTree._blockSize == TJst::REQUIRE_FULL_JOURNAL;
    }

    // ----------------------------------------------------------------------------
    // Function create()
    // ----------------------------------------------------------------------------

    /*!
     * @fn JournaledStringTree#create
     * @headerfile <seqan/journaled_string_tree.h>
     * @brief Creates the journaled strings or the first block if block-wise construction is enabled.
     *
     * @signature bool create(jst, cs[, parallel);
     *
     * @param[in,out] jst       The Journal String Tree to create the journaled strings for.
     * @param[in] cs            The context size.
     * @param[in] parallel      Tag to specialize parallel construction. One of @link ParallelismTags @endlink.
     *                          Defaults to @link ParallelismTags#Serial @endlink.
     *
     * This function creates either the complete journaled strings if @link JournaledStringTree#fullJournalRequired @endlink
     * returns <tt>true<\tt> or the first block if block-wise construction is enabled. To generate the following blocks
     * use the function @link JournaledStringTree#createNext @endlink.
     *
     * @see JournaledStringTree#createNext
     * @see JournaledStringTree#fullJournalRequired
     * @see JournaledStringTree#setBlockSize
     */

    template <typename TDeltaMap, typename TSize, typename TParallelTag>
    void create(JournaledStringTree<TDeltaMap, StringTreeDefault> & stringTree,
                TSize contextSize,
                Tag<TParallelTag> parallelTag)
    {
        typedef JournaledStringTree<TDeltaMap, StringTreeDefault> TJst;

        if (stringTree._blockSize != TJst::REQUIRE_FULL_JOURNAL)  // Require full creation.
        {
            createNext(stringTree, contextSize, parallelTag);
            return;
        }

        if (!stringTree._emptyJournal)  // String tree already created.
            return;
        reinit(stringTree);
        _doJournalBlock(stringTree, contextSize, parallelTag);
        ++stringTree._activeBlock;
        stringTree._emptyJournal = false;
    }

    template <typename TDeltaMap, typename TSize>
    void create(JournaledStringTree<TDeltaMap, StringTreeDefault> const & stringTree,
                TSize contextSize)
    {
        create(stringTree, contextSize, Serial());
    }

    // ----------------------------------------------------------------------------
    // Function createNext()
    // ----------------------------------------------------------------------------

    /*!
     * @fn JournaledStringTree#createNext
     * @headerfile <seqan/journaled_string_tree.h>
     * @brief Constructs the sequence content for the next block if available.
     *
     * @signature bool createNext(jst, cs[, tag]);
     *
     * @param[in,out] jst   The Journal String Tree to create the journaled strings for.
     * @param[in] cs        The context size.
     * @param[in] tag       Tag to enable parallel processing. One of @link ParallelismTags @endlink.
     *                      Defaults to @link ParallelismTags#Serial @endlink.
     *
     *
     * @return bool <tt>true</tt> if a new block was generated, <tt>false</tt> otherwise.
     *
     * Before calling this function the function @link JournaledStringTree#create @endlink must be called.
     *
     * @see JournaledStringTree#create
     */

    template <typename TDeltaMap, typename TSize, typename TParallelTag>
    bool createNext(JournaledStringTree<TDeltaMap, StringTreeDefault> const & stringTree,
                    TSize contextSize,
                    Tag<TParallelTag> tag)
    {
        typedef JournaledStringTree<TDeltaMap, StringTreeDefault> TJst;

        if (stringTree._blockSize == TJst::REQUIRE_FULL_JOURNAL)
            return false;

        bool res = _doJournalBlock(stringTree, contextSize, tag);
        stringTree._emptyJournal = false;
        ++stringTree._activeBlock;
        return res;
    }

    template <typename TDeltaMap, typename TSpec, typename TSize>
    bool createNext(JournaledStringTree<TDeltaMap, TSpec> & stringTree,
                    TSize contextSize)
    {
        return createNext(stringTree, contextSize, Serial());
    }

    // ----------------------------------------------------------------------------
    // Function reinit()
    // ----------------------------------------------------------------------------

    /*!
     * @fn JournaledStringTree#reinit
     * @headerfile <seqan/journaled_string_tree.h>
     * @brief Reinitializes the journaled string tree.
     *
     * @signature void reinit(jst);
     *
     * @param[in,out] jst    The Journal String Tree.
     *
     * This function resets the tree to the first block. Note, that the journaled strings are not cleard such that in case
     * of no block-wise generation the journaled strings are not reconstructed.
     */

    // Keeps the journal strings alive, if the full block has loaded.
    template <typename TDeltaMap, typename TSpec>
    inline void
    reinit(JournaledStringTree<TDeltaMap, TSpec> & jst)
    {
        // Reset all positions to 0.
        jst._activeBlock = 0;
        if (!fullJournalRequired(jst))
        {
            jst._emptyJournal = true;
            arrayFill(begin(jst._activeBlockVPOffset, Standard()), end(jst._activeBlockVPOffset, Standard()), 0);
            arrayFill(begin(jst._blockVPOffset, Standard()), end(jst._blockVPOffset, Standard()), 0);
        }
        jst._mapBlockBegin = jst._mapBlockEnd = begin(container(jst), Standard());
    }

    // ----------------------------------------------------------------------------
    // Function init()
    // ----------------------------------------------------------------------------

    /*!
     * @fn JournaledStringTree#init
     * @headerfile <seqan/journaled_string_tree.h>
     * @brief Initializes the journaled string tree.
     *
     * @signature void init(jst, host, delta);
     *
     * @param[in,out] jst    The Journal String Tree.
     * @param[in] host       The reference sequence.
     * @param[in] delta      The object containing the delta information.
     *
     * This function does not construct the sequences. Use the function @link JournaledStringTree#create @endlink and
     * @link JournaledStringTree#createNext @endlink to dynamically generate the sequences on-demand.
     */
    
    template <typename TDeltaMap, typename TSpec, typename THost>
    inline void
    init(JournaledStringTree<TDeltaMap, TSpec> & jst,
         THost & referenceSeq,
         TDeltaMap & varData)
    {
        typedef JournaledStringTree<TDeltaMap, TSpec> TJst;
        typedef typename GetStringSet<TJst>::Type TStringSet;
        typedef typename Value<TStringSet>::Type TString;
        
        setValue(jst._container, varData);
        setHost(stringSet(jst), referenceSeq);
        
        TString tmp;
        setHost(tmp, host(stringSet(jst)));
        resize(stringSet(jst), getCoverageSize(varData), tmp, Exact());
        resize(jst._blockVPOffset, length(stringSet(jst)), 0, Exact());
        resize(jst._activeBlockVPOffset, length(stringSet(jst)), 0, Exact());
        
        jst._mapBlockBegin = jst._mapBlockEnd = begin(varData, Standard());
    }

    // ----------------------------------------------------------------------------
    // Function setBlockSize()
    // ----------------------------------------------------------------------------

    /*!
     * @fn JournaledStringTree#setBlockSize
     * @headerfile <seqan/journaled_string_tree.h>
     * @brief Sets the number of variants processed at a time.
     *
     * @signature void setBlockSize(jst, block);
     *
     * @param[in]   jst     The Journal String Tree.
     * @param[in]   block   The number of deltas contained in one block.
     *
     * Per default all deltas are processed in a single block.
     *
     * @see JournaledStringTree#getBlockSize
     */

    template <typename TDeltaMap, typename TSpec, typename TPosition>
    inline void
    setBlockSize(JournaledStringTree<TDeltaMap, TSpec> & stringTree,
                 TPosition const & newBlockSize)
    {
        SEQAN_ASSERT_NOT(empty(stringTree._container));  // The delta map needs to be set before.

        stringTree._blockSize = newBlockSize;
        stringTree._activeBlock = 0;
        stringTree._numBlocks = (length(container(stringTree)) + newBlockSize - 1) / newBlockSize;
        resize(stringTree._blockVPOffset, length(stringSet(stringTree)), 0, Exact());
    }

    // ----------------------------------------------------------------------------
    // Function getBlockSize()
    // ----------------------------------------------------------------------------

    /*!
     * @fn JournaledStringTree#getBlockSize
     * @headerfile <seqan/journaled_string_tree.h>
     * @brief Returns the size of the blocks.
     *
     * @signature TSize getBlockSize(jst);
     *
     * @param[in]   jst     The Journal String Tree.
     *
     * @return TSize The size of the blocks of type @link JournaledStringTree#Size @endlink.
     *
     * @see JournaledStringTree#setBlockSize
     */

    template <typename TDeltaMap, typename TSpec>
    inline typename Size<JournaledStringTree<TDeltaMap, TSpec> >::Type
    getBlockSize(JournaledStringTree<TDeltaMap, TSpec> const & stringTree)
    {
        return stringTree._blockSize;
    }

    // ----------------------------------------------------------------------------
    // Function container()
    // ----------------------------------------------------------------------------

    /*!
     * @fn JournaledStringTree#container
     * @headerfile <seqan/journaled_string_tree.h>
     * @brief Returns a reference to the object holding the delta information for a set of sequences.
     *
     * @signature TContainer container(jst);
     *
     * @param[in]   jst     The Journal String Tree.
     *
     * @return TContainer A reference to the object holding the delta information of type @link JournaledStringTree#Container @endlink.
     *
     * @see JournaledStringTree#stringSet
     */

    template <typename TDeltaMap, typename TSpec>
    inline typename Container<JournaledStringTree<TDeltaMap, TSpec> >::Type &
    container(JournaledStringTree<TDeltaMap, TSpec> & stringTree)
    {
        return value(stringTree._container);
    }

    template <typename TDeltaMap, typename TSpec>
    inline typename Container<JournaledStringTree<TDeltaMap, TSpec> const>::Type &
    container(JournaledStringTree<TDeltaMap, TSpec> const & stringTree)
    {
        return value(stringTree._container);
    }

    // ----------------------------------------------------------------------------
    // Function stringSet()
    // ----------------------------------------------------------------------------

    /*!
     * @fn JournaledStringTree#stringSet
     * @headerfile <seqan/journaled_string_tree.h>
     * @brief Returns a reference to the constructed sequences.
     *
     * @signature TStringSet stringSet(jst);
     *
     * @param[in]   jst     The Journal String Tree.
     *
     * @return TStringSet A reference to the object containing the compressed strings @link JournaledStringTree#GetStringSet @endlink.
     *
     * @see JournaledStringTree#container
     */

    template <typename TDeltaMap, typename TSpec>
    inline typename GetStringSet<JournaledStringTree<TDeltaMap, TSpec> >::Type &
    stringSet(JournaledStringTree<TDeltaMap, TSpec> & stringTree)
    {
        return stringTree._journalSet;
    }

    template <typename TDeltaMap, typename TSpec>
    inline typename GetStringSet<JournaledStringTree<TDeltaMap, TSpec> const>::Type &
    stringSet(JournaledStringTree<TDeltaMap, TSpec> const & stringTree)
    {
        return stringTree._journalSet;
    }

}

#endif  // #ifndef INCLUDE_SEQAN_JOURNALED_STRING_TREE_TRAVERSER_CONTEXT_H_
