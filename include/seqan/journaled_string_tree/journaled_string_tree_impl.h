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

struct JstBaseSequenceMember_;
typedef Tag<JstBaseSequenceMember_> JstBaseSequenceMember;

struct JstDeltaMapMember_;
typedef Tag<JstDeltaMapMember_> JstDeltaMapMember;

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

template <typename TSequence, typename TConfig = DefaultJstConfig<TSequence>, typename TSpec = Default>
class JournaledStringTree<TSequence, TConfig, TSpec>
{
public:

    typedef typename Member<JournaledStringTree, JstDeltaMapMember>::Type       TDeltaMap;
    typedef typename Member<JournaledStringTree, JstBaseSequenceMember>::Type   TBaseSeq;

    TBaseSeq   _baseSeq;
    TDeltaMap  _deltaMap;

    JournaledStringTree(TSize newDepth) : _deltaMap(newDepth)
    {}

    template <typename TSeq>
    JournaledStringTree(TSize newDepth, TSeq & seq) : _deltaMap(newDepth), _baseSeq(seq)
    {}
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

template <typename TSequence, typename TConfig, typename TSpec>
struct Spec<JournaledStringTree<TSequence, TConfig, TSpec> const>
{
    typedef TSpec Type;
};

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

template <typename TSequence, typename TConfig, typename TSpec>
struct Size<JournaledStringTree<TSequence, TConfig, TSpec> >
{
    typedef size_t Type;
};

template <typename TSequence, typename TConfig, typename TSpec>
struct Size<JournaledStringTree<TSequence, TConfig, TSpec> const>
{
    typedef size_t const Type;
};

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

template <typename TSequence, typename TConfig, typename TSpec>
struct Position<JournaledStringTree<TSequence, TConfig, TSpec> >
{
    typedef typename TConfig::TDeltaPos Type;
};

template <typename TSequence, typename TConfig, typename TSpec>
struct Position<JournaledStringTree<TSequence, TConfig, TSpec> const>
{
    typedef typename TConfig::TDeltaPos const Type;
};

// ----------------------------------------------------------------------------
// Metafunction Member                                           [BaseSequence]
// ----------------------------------------------------------------------------

/*!
 * @mfn JournaledStringTree#JstBaseSequence
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Returns the type of the base sequence.
 *
 * @signature Member<TJst, JstBaseSequenceMember>::Type;
 *
 * @tparam TJst The type of the journal string tree to get the host type for.
 *
 * @return TMember The type of the base sequence.
 */

template <typename TSequence, typename TConfig, typename TSpec>
struct Member<JournaledStringTree<TSequence, TConfig, TSpec>, JstDeltaMapMember const>
{
    typedef typename DeltaMap<typename TConfig::TDeltaPos, typename TConfig::TSnpValue> Type;
};

template <typename TSequence, typename TConfig, typename TSpec>
struct Member<JournaledStringTree<TSequence, TConfig, TSpec> const, JstDeltaMapMember const>
{
    typedef typename DeltaMap<typename TConfig::TDeltaPos, typename TConfig::TSnpValue> const Type;
};

// ----------------------------------------------------------------------------
// Metafunction Member                                           [BaseSequence]
// ----------------------------------------------------------------------------

/*!
 * @mfn JournaledStringTree#JstBaseSequence
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Returns the type of the base sequence.
 *
 * @signature Member<TJst, JstBaseSequenceMember>::Type;
 *
 * @tparam TJst The type of the journal string tree to get the host type for.
 *
 * @return TMember The type of the base sequence.
 */

template <typename TSequence, typename TConfig, typename TSpec>
struct Member<JournaledStringTree<TSequence, TConfig, TSpec>, JstBaseSequenceMember const>
{
    typedef typename TSequence Type;
};

template <typename TSequence, typename TConfig, typename TSpec>
struct Member<JournaledStringTree<TSequence, TConfig, TSpec> const, JstBaseSequenceMember const>
{
    typedef typename TSequence const Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function baseSeq()
// ----------------------------------------------------------------------------

/*!
 * @fn JournaledStringTree#baseSeq
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Returns a reference to the base sequence used for the journaled string tree.
 *
 * @signature TSequence basSeq(jst);
 *
 * @param[in] jst    The Journal String Tree.
 *
 * @return TSequence A reference to the base sequence of type @link JournaledStringTree#JstBaseSequence @endlink.
 *
 * @see JournaledStringTree#setBaseSeq
 */

template <typename TSequence, typename TConfig, typename TSpec>
inline typename Member<JournaledStringTree<TSequence, TConfig, TSpec>, JstBaseSequenceMember >::Type &
baseSeq(JournaledStringTree<TSequence, TConfig, TSpec> & jst)
{
    return jst._baseSeq;
}

template <typename TSequence, TConfig, typename TSpec>
inline typename Member<JournaledStringTree<TSequence, TConfig, TSpec> const, JstBaseSequenceMember>::Type &
baseSeq(JournaledStringTree<TSequence, TConfig, TSpec> const & jst)
{
    return jst._baseSeq;
}

// ----------------------------------------------------------------------------
// Function setBaseSeq()
// ----------------------------------------------------------------------------

/*!
 * @fn JournaledStringTree#setBaseSeq
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Sets the new base sequence.
 *
 * @signature void setBaseSeq(jst, seq);
 *
 * @param[in,out] jst   The Journaled-String-Tree to set the new base sequence for.
 * @param[in]     seq   The new sequence to set as base sequence.
 *
 * Setting a new base sequence resets all the node information.
 *
 * @see JournaledStringTree#baseSeq
 */

template <typename TSequence, typename TConfig, typename TSpec, typename TSequence2>
inline void
setBaseSeq(JournaledStringTree<TSequence, TConfig, TSpec> & jst, TSequence2 const & newBaseSeq)
{
    jst._baseSeq = newBaseSeq;
    clear(jst._deltaMap);
}

// ----------------------------------------------------------------------------
// Function depth()
// ----------------------------------------------------------------------------

/*!
 * @fn JournaledStringTree#depth
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Returns the number of sequences represented by the Journaled-String-Tree.
 *
 * @signature TSize depth(jst);
 *
 * @param[in] jst The Journal String Tree to query the depth for.
 *
 * @return TSize The number of sequences represented; of type @link JournaledStringTree#Size @endlink.
 *
 * @see JournaledStringTree#setDepth
 */

template <typename TSequence, typename TConfig, typename TSpec>
inline typename Size<JournaledStringTree<TSequence, TConfig, TSpec>>::Type
depth(JournaledStringTree<TSequence, TConfig, TSpec> & jst)
{
    return getCoverageSize(jst._deltaMap);
}

template <typename TSequence, TConfig, typename TSpec>
inline typename Size<JournaledStringTree<TSequence, TConfig, TSpec> const>::Type
depth(JournaledStringTree<TSequence, TConfig, TSpec> const & jst)
{
    return getCoverageSize(jst._deltaMap);
}

// ----------------------------------------------------------------------------
// Function setDepth()
// ----------------------------------------------------------------------------

/*!
 * @fn JournaledStringTree#setDepth
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Sets the new number of sequences represented by this Journaled-String-Tree.
 *
 * @signature void setDepth(jst, d);
 *
 * @param[in,out] jst   The Journaled-String-Tree to set the new depth for.
 * @param[in]     d     The new depth.
 *
 * Internally parses all nodes of the tree and adapts their coverage accordingly. If the new depth is greater than the
 * old depth the new represented sequences do not cover any of the exising nodes.
 *
 * @see JournaledStringTree#depth
 */

template <typename TSequence, typename TConfig, typename TSpec, typename TSize>
inline void
setDepth(JournaledStringTree<TSequence, TConfig, TSpec> & jst, TSize newDepth)
{
    setCoverageSize(jst._deltaMap, newDepth);
}

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

/*!
 * @fn JournaledStringTree#clear
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Clears the Journaled-String-Tree.
 *
 * @signature void clear(jst);
 *
 * @param[in,out] jst  The Journal String Tree to be cleared.
 */

template <typename TSequence, typename TConfig, typename TSpec>
inline void
clear(JournaledStringTree<TSequence, TConfig, TSpec> & jst)
{
    clear(jst._baseSeq);
    clear(jst._deltaMap);
}

// ----------------------------------------------------------------------------
// Function addNode()
// ----------------------------------------------------------------------------

// Public interface to add a new node.
template <typename TSequence, typename TConfig, typename TSpec, typename TPos, typename TValue, typename TIds,
          typename TDeltaType>
inline SEQAN_ENABLE_IF(Is<SequenceConcept<TIds> >, void)
addNode(JournaledStringTree<TSequence, TConfig, TSpec> & jst,
        TPos basePos,
        TValue const & deltaVal,
        TIds const & ids,
        TDeltaType /*deltaType*/)
{
    typedef typename Member<JournaledStringTree, JstDeltaMapMember>::Type TDeltaMap;
    typedef typename DeltaCoverage<TDeltaMap>::Type TCoverage;

    // Transform the ids into coverage value.
    TCoverage coverage(coverageSize(jst._deltaMap), false, Exact());

    std::for_each(begin(ids, Standard()), end(ids, Standard()),
                  [&](auto const & it)
                  {
                      SEQAN_ASSERT_LT(*it, depth(jst));  // Check that id is valid.
                      coverage[*it] = true;
                  })

    insert(jst._deltaMap, basePos, deltaValue, coverage, TDeltaType());
}

// ----------------------------------------------------------------------------
// Function eraseNode()
// ----------------------------------------------------------------------------

// TODO(rrahn): Implement eraseNode;

// ----------------------------------------------------------------------------
// Function updateNode()
// ----------------------------------------------------------------------------

// TODO(rrahn): implement updateNode() -> updates the coverage of a node.

}

#endif // EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_JOURNALED_STRING_TREE_IMPL_H_
