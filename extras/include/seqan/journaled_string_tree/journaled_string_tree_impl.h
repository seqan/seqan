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
// Implements the default virtual string tree to traverse multiple sequences
// in parallel. This is a facade combining the variant store with the
// journal set.
// ==========================================================================

#ifndef EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_JOURNALED_STRING_TREE_IMPL_H_
#define EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_JOURNALED_STRING_TREE_IMPL_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename T>
struct GetStringSet{};

// ----------------------------------------------------------------------------
// Class JournaledStringTree                                [StringTreeDefault]
// ----------------------------------------------------------------------------

/*!
 * @class JournaledStringTree Journaled String Tree
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief A dynamic data structure representing a set of similar sequences in succinct form.
 *
 * @signature template <typename TDelta[, typename TSpec]>
 *            class JournaledStringTree<TDeltaStore, TSpec>;
 *
 * @tparam TDelta Type of the object storing delta information for a set of sequences. Of type @link DeltaMap @endlink.
 * @tparam TSpec The journaled string tree type. Defaults to @link Default @endlink.
 *
 * This data structure stores delta information between a set of sequences and a common reference sequence. The
 * @link DeltaMap @endlink is used to store these information efficiently.
 *
 * In addition, the journaled string tree manages a set of referentially compressed strings @link JournaledSet @endlink
 * to represent the sequences encoded by the the delta map. The sequence information is built on-demand and can
 * be constructed dynamically in blocks.
 *
 * Use the function @link JournaledStringTree#create @endlink and @link JournaledStringTree#createNext @endlink to force
 * the creation of the sequence information for the next block.
 * Per default the entire sequence set is build. To build the sequences in blocks
 * over the mapped variants use the function @link JournaledStringTree#setBlockSize @endlink. The block size
 * is the number of variants that should be integrated in the current sequence content.
 * Note that when generating the sequences block-wise, the actual position of each sequence depends only on the
 * currently incoporated variants. Use the function @link JournaledStringTree#localToGlobalPos @endlink to map
 * the local position to its global position for a selected sequence.
 *
 * The sequences represented by the delta map and the journaled strings can be traversed in forward direction.
 * To do so one need to instantiate a @link JstTraverser @endlink object, which manages the traversal over the
 * Journaled String Tree. To communicate with the Journaled String Tree during the traversal the caller of the
 * @link JstTraverser#traverse @endlink has to implement the @link JstTraversalConcept @endlink.
 */

template <typename TDeltaMap>
class JournaledStringTree<TDeltaMap, StringTreeDefault>
{
public:
    typedef typename GetStringSet<JournaledStringTree>::Type    TJournalData;
    typedef typename Value<TJournalData>::Type                  TJournaledString;
    typedef typename Size<TJournalData>::Type                   TSize;
    typedef typename MakeSigned<TSize>::Type                    TSignedSize;
    typedef typename Container<JournaledStringTree>::Type       TContainer;
    typedef typename Iterator<TContainer, Standard>::Type       TMapIterator;

    // TODO(rmaerker): Maybe no holder.
    Holder<TDeltaMap>               _container;

    mutable TJournalData            _journalSet;
    mutable String<TSignedSize>     _blockVPOffset;  // Virtual position offset for each sequence visited so far.
    mutable String<TSignedSize>     _activeBlockVPOffset;  // Virtual position offset for each sequence visited so far.
    mutable TSize                   _activeBlock;
    mutable TMapIterator            _mapBlockBegin;
    mutable TMapIterator            _mapBlockEnd;
    mutable bool                    _emptyJournal;

    TSize                           _blockSize;
    TSize                           _numBlocks;
    static const TSize REQUIRE_FULL_JOURNAL = MaxValue<TSize>::VALUE;

    JournaledStringTree() :
        _activeBlock(0),
        _emptyJournal(true),
        _blockSize(REQUIRE_FULL_JOURNAL),
        _numBlocks(1)
    {
        create(_container);
    }

    template <typename THost>
    JournaledStringTree(THost & reference, TDeltaMap & varData) :
        _activeBlock(0),
        _mapBlockBegin(),
        _mapBlockEnd(),
        _emptyJournal(true),
        _blockSize(REQUIRE_FULL_JOURNAL),
        _numBlocks(1)
    {
        init(*this, reference, varData);
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Spec
// ----------------------------------------------------------------------------

/*!
 * @mfn JournaledStringTree#Spec
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Returns the specialization type.
 *
 * @signature Spec<TJst>::Type;
 *
 * @tparam TJst The type of the journaled string tree to get the specialization type for.
 *
 * @return TSpec The specialization type.
 */

template <typename TDeltaMap, typename TSpec>
struct Spec<JournaledStringTree<TDeltaMap, TSpec> >
{
    typedef TSpec Type;
};

template <typename TDeltaMap, typename TSpec>
struct Spec<JournaledStringTree<TDeltaMap, TSpec> const> :
    Spec<JournaledStringTree<TDeltaMap, TSpec> >{};

// ----------------------------------------------------------------------------
// Metafunction Position
// ----------------------------------------------------------------------------

/*!
 * @mfn JournaledStringTree#Position
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief The position type.
 *
 * @signature Position<TJst>::Type;
 *
 * @tparam TJst The type of the journal string tree to get the position type for.
 *
 * @return TPosition The position type.
 */

template <typename TDeltaMap, typename TSpec>
struct Position<JournaledStringTree<TDeltaMap, TSpec> > :
    Position<TDeltaMap>{};

template <typename TDeltaMap, typename TSpec>
struct Position<JournaledStringTree<TDeltaMap, TSpec> const> :
    Position<TDeltaMap const>{};

// ----------------------------------------------------------------------------
// Metafunction Size
// ----------------------------------------------------------------------------

/*!
 * @mfn JournaledStringTree#Size
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief The size type.
 *
 * @signature Size<TJst>::Type;
 *
 * @tparam TJst The type of the journal string tree to get the size type for.
 *
 * @return TSize The size type.
 */

template <typename TDeltaMap, typename TSpec>
struct Size<JournaledStringTree<TDeltaMap, TSpec> > :
    Size<TDeltaMap>{};

template <typename TDeltaMap, typename TSpec>
struct Size<JournaledStringTree<TDeltaMap, TSpec> const> :
    Size<TDeltaMap const>{};

// ----------------------------------------------------------------------------
// Metafunction Host
// ----------------------------------------------------------------------------

/*!
 * @mfn JournaledStringTree#Host
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Returns the type of the global reference.
 *
 * @signature Host<TJst>::Type;
 *
 * @tparam TJst The type of the journal string tree to get the host type for.
 *
 * @return THost The type of the global reference sequence.
 */

template <typename TDeltaMap, typename TSpec>
struct Host<JournaledStringTree<TDeltaMap, TSpec> >
{
    typedef JournaledStringTree<TDeltaMap, TSpec> TVStringTree_;
    typedef typename GetStringSet<TVStringTree_>::Type TJournalData_;
    typedef typename Host<TJournalData_>::Type Type;
};

template <typename TDeltaMap, typename TSpec>
struct Host<JournaledStringTree<TDeltaMap, TSpec> const>
{
    typedef JournaledStringTree<TDeltaMap, TSpec> TVStringTree_;
    typedef typename Host<TVStringTree_>::Type const Type;
};

// ----------------------------------------------------------------------------
// Metafunction Container
// ----------------------------------------------------------------------------

/*!
 * @mfn JournaledStringTree#Container
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Returns the type of the branch node map holding the delta information.
 *
 * @signature Container<TJst>::Type;
 *
 * @tparam TJst The type of the journal string tree to get the branch node map type for.
 *
 * @return TContainer The branch node map containing the delta information.
 */

template <typename TDeltaMap, typename TSpec>
struct Container<JournaledStringTree<TDeltaMap, TSpec> >
{
    typedef TDeltaMap Type;
};

template <typename TDeltaMap, typename TSpec>
struct Container<JournaledStringTree<TDeltaMap, TSpec> const>
{
    typedef TDeltaMap const Type;
};

// ----------------------------------------------------------------------------
// Metafunction GetStringSet
// ----------------------------------------------------------------------------

/*!
 * @mfn JournaledStringTree#GetStringSet
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Returns the type of the string set holding the compressed sequences.
 *
 * @signature GetStringSet<TJst>::Type;
 *
 * @tparam TJst The type of the journaled string tree to get the string set type for.
 *
 * @return TStringSet The string set containing the constructed sequences.
 */

template <typename TDeltaMap, typename TSpec>
struct GetStringSet<JournaledStringTree<TDeltaMap, TSpec> >
{
    typedef typename DeltaValue<TDeltaMap, DeltaTypeSnp>::Type TAlphabet_;

    typedef String<TAlphabet_, Journaled<Alloc<>, SortedArray> > TString_;
    typedef StringSet<TString_, Owner<JournaledSet> > Type;
};

template <typename TDeltaMap, typename TSpec>
struct GetStringSet<JournaledStringTree<TDeltaMap, TSpec> const>
{
    typedef typename DeltaValue<TDeltaMap, DeltaTypeSnp>::Type TAlphabet_;

    typedef String<TAlphabet_, Journaled<Alloc<>, SortedArray> > TString_;
    typedef StringSet<TString_, Owner<JournaledSet> > const Type;
};


// ============================================================================
// Functions
// ============================================================================

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
    typedef typename Iterator<TDeltaMap const, Standard>::Type TMapIterator;
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
// Function host()
// ----------------------------------------------------------------------------

/*!
 * @fn JournaledStringTree#host
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Returns a reference to the global reference sequence.
 *
 * @signature THost host(jst);
 *
 * @param[in] jst    The Journal String Tree.
 *
 * @return THost A reference to the @link JournaledStringTree#Host @endlink object.
 */

template <typename TDeltaMap, typename TSpec>
inline typename Host<JournaledStringTree<TDeltaMap, TSpec> >::Type &
host(JournaledStringTree<TDeltaMap, TSpec> & stringTree)
{
    return host(stringSet(stringTree));
}

template <typename TDeltaMap, typename TSpec>
inline typename Host<JournaledStringTree<TDeltaMap, TSpec> const>::Type &
host(JournaledStringTree<TDeltaMap, TSpec> const & stringTree)
{
    return host(stringSet(stringTree));
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
// Function clear()
// ----------------------------------------------------------------------------

/*!
 * @fn JournaledStringTree#clear
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Clears the container.
 *
 * @signature void clear(jst);
 *
 * @param[in,out] jst    The Journal String Tree to be cleared.
 *
 * @see JournaledStringTree#reinit
 * @see JournaledStringTree#init

 */

template <typename TDeltaMap, typename TSpec>
inline void
clear(JournaledStringTree<TDeltaMap, TSpec> & jst)
{

    clear(jst._container);
    clear(jst._journalSet);
    clear(jst._activeBlockVPOffset);
    clear(jst._blockVPOffset);
    jst._numBlocks = 1;
    jst._activeBlock = 0;
    jst._blockSize = JournaledStringTree<TDeltaMap, TSpec>::REQUIRE_FULL_JOURNAL;
    jst._emptyJournal = true;
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

// ----------------------------------------------------------------------------
// Function _isValiJstdReference()
// ----------------------------------------------------------------------------

template <typename THost, typename TRefId, typename TConfig>
inline bool _isValiJstdReference(THost const & host,
                                 TRefId const & refId,
                                 GdfHeader const & gdfHeader,
                                 GdfFileConfiguration<TConfig> const & config)
{
    unsigned int crc = computeCrc(host);
    if (crc == 0)
        return refId == gdfHeader.referenceId;
    return crc == config.refHash;
}

// ----------------------------------------------------------------------------
// Function open()
// ----------------------------------------------------------------------------

/*!
 * @fn JournaledStringTree#open
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Opens the journaled string tree from a gdf file.
 *
 * @signature int open(jst, gdfHeader, filename);
 *
 * @param[out] jst          The journaled string tree to be loaded from disk.
 * @param[out] gdfHeader    The header information of the gdf file of type @link GdfHeader @endlink.
 * @param[in]  filename     The name of the file used to load the data from.
 *
 * The journaled string tree is stored and loaded as GDF file. For more information see @link GdfIO @endlink.
 *
 * @throw GdfIOException The type of the exception thrown if unexpected error occur. See @link GdfIOException @endlink.
 *
 * @see JournaledStringTree#save
 */

template <typename TDeltaMap, typename TSpec, typename TFilename>
inline void open(JournaledStringTree<TDeltaMap, TSpec> & jst,
                 GdfHeader & gdfHeader,
                 TFilename const & filename)
{
    typedef JournaledStringTree<TDeltaMap, TSpec> TJst;
    typedef typename DeltaValue<TDeltaMap, DeltaTypeSnp>::Type TSnp;
    typedef typename Host<TJst>::Type THost;
    typedef typename GetStringSet<TJst>::Type TStringSet;
    typedef typename Value<TStringSet>::Type TJString;

    GdfFileConfiguration<TSnp> config(0u);

    // Step 1) Read the delta map.
    std::ifstream inputFile;
    inputFile.open(toCString(filename), std::ios_base::in | std::ios_base::binary);
    if (!inputFile.good())
    {
        inputFile.close();
        std::stringstream errMessage;
        errMessage << "Unknown file <" << filename << ">!";
        SEQAN_THROW(GdfIOException(errMessage.str()));
    }
    read(container(jst), gdfHeader, config, inputFile, Gdf());
    inputFile.close();

    // Step 2) Read the reference file.
    std::ifstream refFile;
    refFile.open(toCString(gdfHeader.referenceFilename), std::ios_base::in);
    if (!refFile.good())
    {
        refFile.close();
        std::stringstream errMessage;
        errMessage << "Unknown file <" << gdfHeader.referenceFilename << ">!";
        SEQAN_THROW(GdfIOException(errMessage.str()));
    }
    THost tmpHost = "";
    CharString refId;
    createHost(stringSet(jst), tmpHost);
    RecordReader<std::ifstream, SinglePass<> > reader(refFile);
    if (readRecord(refId, host(jst), reader, Fasta()) != 0)
    {
        refFile.close();
        std::stringstream errMessage;
        errMessage << "Error while parsing <" << gdfHeader.referenceFilename << ">!";
        SEQAN_THROW(GdfIOException(errMessage.str()));
    }
    refFile.close();

    if (!_isValiJstdReference(host(jst), refId, gdfHeader, config))
        SEQAN_THROW(GdfIOWrongReferenceException(refId, gdfHeader.referenceId, config.refHash, computeCrc(host(jst))));
    
    // Initialize the JST.
    resize(stringSet(jst), getCoverageSize(container(jst)), TJString(host(stringSet(jst))), Exact());
    resize(jst._blockVPOffset, length(stringSet(jst)), 0, Exact());
    resize(jst._activeBlockVPOffset, length(stringSet(jst)), 0, Exact());
    
    jst._mapBlockBegin = jst._mapBlockEnd = begin(container(jst), Standard());
}

// ----------------------------------------------------------------------------
// Function save()
// ----------------------------------------------------------------------------

/*!
 * @fn JournaledStringTree#save
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Saves the journaled string tree at the given location in GDF format.
 *
 * @signature int save(jst, gdfHeader, filename);
 *
 * @param[in] jst The   The journaled string tree to be saved.
 * @param[in] gdfHeader The header file of type @link GdfHeader @endlink.
 * @param[in] filename  The name of the file used to save the journaled string tree.
 *
 * Stores the journaled string tree in gdf format at the specified file location.
 *
 * @throw GdfIOException The exception type thrown if an unexpected error occurs. See @link GdfIOException @endlink.
 *
 * @see JournaledStringTree#open
 */

template <typename TDeltaMap, typename TSpec, typename TFilename>
inline void save(JournaledStringTree<TDeltaMap, TSpec> const & jst,
                 GdfHeader const & gdfHeader,
                 TFilename const & filename)
{
    typedef typename DeltaValue<TDeltaMap, DeltaTypeSnp>::Type TSnpType;

    // Check if sequence names are available.
    if (length(gdfHeader.nameStore) < getCoverageSize(container(jst)))
    {
        std::stringstream errMessage;
        errMessage << "Too few sequence names - Needed " << getCoverageSize(container(jst)) << " but " << length(gdfHeader.nameStore) << " were provided!";
        SEQAN_THROW(GdfIOException(errMessage.str()));
    }

    // Write reference.
    if (gdfHeader.referenceMode == GdfIOMode::SAVE_REFERENCE_MODE_ENABLED)
    {
        std::ofstream refStream;
        refStream.open(toCString(gdfHeader.referenceFilename), std::ios_base::out);
        if (!refStream.good())
        {
            std::stringstream errMessage;
            errMessage << "Cannot open file: "<< gdfHeader.referenceFilename << "!";
            SEQAN_THROW(GdfIOException(errMessage.str()));
        }

        if (writeRecord(refStream, gdfHeader.referenceId, host(jst), Fasta()) != 0)
            SEQAN_THROW(GdfIOException("While writing reference to disk!"));
        refStream.close();
    }

    // Compute hash for reference sequence.
    GdfFileConfiguration<TSnpType> gdfConfig(getCoverageSize(container(jst)));
    gdfConfig.refHash = computeCrc(host(jst));

    std::ofstream fileStream;
    fileStream.open(toCString(filename), std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);

    if (!fileStream.good())
    {
        std::stringstream errMessage;
        errMessage << "Cannot open file: " << filename << "!";
        SEQAN_THROW(GdfIOException(errMessage.str()));
    }
    write(fileStream, container(jst), gdfHeader, gdfConfig, Gdf());
    fileStream.close();
}

}

#endif // EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_JOURNALED_STRING_TREE_IMPL_H_
