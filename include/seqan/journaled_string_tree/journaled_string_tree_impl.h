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
 * @brief A virtual string tree over a set of sequences sharing the same base sequence.
 *
 * @signature template <typename TSequence[, typename TConfig][, TSpec]>
 *            class JournaledStringTree<TDeltaStore, TConfig, TSpec>;
 *
 * @tparam TSequence Type of underlying base sequence.
 * @tparam TConfig   A configuration object. Defaults to @link DefaultJstConfig @endlink.
 * @tparam TSpec     The specialization tag for the Journaled-String-Tree. Defaults to @link Default @endlink.
 *
 * This data structure stores delta values between a set of sequences and a common base sequence. A
 * @link DeltaMap @endlink is used to store these information efficiently.
 */

template <typename TSequence, typename TConfig>
class JournaledStringTree<TSequence, TConfig, Default>
{
public:

    typedef JournaledStringTree<TSequence, TConfig, Default>    TJst;
    typedef typename Member<TJst, JstDeltaMapMember>::Type      TDeltaMap;
    typedef typename Member<TJst, JstBaseSequenceMember>::Type  TBaseSeq;
    typedef typename Size<JournaledStringTree>::Type            TSize;

    TSize      _depth;
    TBaseSeq   _baseSeq;
    TDeltaMap  _deltaMap;

    /*!
     * @fn JournaledStringTree::JournaledStringTree
     * @brief Constructor.
     * @headerfile <seqan/journaled_string_tree.h>
     *
     * @signature JournaledStringTree(depth);
     * @signature JournaledStringTree(depth, seq);
     *
     * @param depth  The number of sequences represented by the Journaled-String-Tree.
     * @param seq    The underlying base sequence.
     */
    template <typename TSize>
    JournaledStringTree(TSize newDepth) : _depth(newDepth)
    {}

    template <typename TSize, typename TSeq>
    JournaledStringTree(TSize newDepth, TSeq & seq) : _depth(newDepth), _baseSeq(seq)
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
struct Member<JournaledStringTree<TSequence, TConfig, TSpec>, JstDeltaMapMember>
{
    typedef DeltaMap<TConfig> Type;
};

template <typename TSequence, typename TConfig, typename TSpec>
struct Member<JournaledStringTree<TSequence, TConfig, TSpec> const, JstDeltaMapMember>
{
    typedef DeltaMap<TConfig> const Type;
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
struct Member<JournaledStringTree<TSequence, TConfig, TSpec>, JstBaseSequenceMember>
{
    typedef TSequence Type;
};

template <typename TSequence, typename TConfig, typename TSpec>
struct Member<JournaledStringTree<TSequence, TConfig, TSpec> const, JstBaseSequenceMember>
{
    typedef TSequence const Type;
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

template <typename TSequence, typename TConfig, typename TSpec>
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
 * @remark Setting a new base sequence resets all the node information.
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
 * @return TSize The number of sequences represented by the <tt>jst<\tt>. Of type @link JournaledStringTree#Size @endlink.
 */

template <typename TSequence, typename TConfig, typename TSpec>
inline typename Size<JournaledStringTree<TSequence, TConfig, TSpec> >::Type
depth(JournaledStringTree<TSequence, TConfig, TSpec> & jst)
{
    return jst._depth;
}

template <typename TSequence, typename TConfig, typename TSpec>
inline typename Size<JournaledStringTree<TSequence, TConfig, TSpec> const>::Type
depth(JournaledStringTree<TSequence, TConfig, TSpec> const & jst)
{
    return jst._depth;
}

// ----------------------------------------------------------------------------
// Function maxPos()
// ----------------------------------------------------------------------------

/*!
 * @fn JournaledStringTree#maxPos
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Returns the greatest position allowed within the Journaled-String-Tree.
 *
 * @signature TPos maxPos(jst);
 *
 * @param[in] jst The Journal String Tree to query the maximal position for.
 *
 * @return TPos The maximal position. Of type @link JournaledStringTree#Position @endlink.
 */

template <typename TSequence, typename TConfig, typename TSpec>
inline typename Position<JournaledStringTree<TSequence, TConfig, TSpec> >::Type
maxPos(JournaledStringTree<TSequence, TConfig, TSpec> const /*jst*/)
{
    return MaxValue<typename TConfig::TDeltaPos>::VALUE;
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
    jst._depth = 0;
    clear(jst._baseSeq);
    clear(jst._deltaMap);
}

// ----------------------------------------------------------------------------
// Function addNode()
// ----------------------------------------------------------------------------

/*!
 * @fn JournaledStringTree#addNode
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Adds a new delta node to the Journaled-String-Tree.
 *
 * @signature bool addNode(jst, pos, val, ids, type);
 *
 * @param[in,out]   jst     The Journaled-String-Tree to insert the new node to.
 * @param[in]       pos     The position wihtin the base sequence of the <tt>jst<\tt> at which the new delta node is placed.
 * @param[in]       val     The value represented by the delta node.
 * @param[in]       ids     An object fulfilling the @link StringConcept @endlink representing the ids of the sequences that cover this delta node.
 * @param[in]       type    A tag representing the type of the delta node. One of @link DeltaTypeTags @endlink.
 *
 * Note the value given by <tt>val<\tt> depends on the <tt>type<\tt>, which can be a single character, a string, 
 * an integer or a pair of string and integer, representing a SNP, an insertion, a deletion or 
 * a structural variant respectively.
 *
 * @return bool <tt>true<\tt> if the node could be inserted, or <tt>false<\tt>, in case there exists already a node
 * with at the same <tt>pos<\tt> and with the same <tt>type<\tt>.
 *
 * @remark The insertion time is linear in the number of nodes.
 */

template <typename TSequence, typename TConfig, typename TSpec, typename TPos, typename TValue, typename TIds,
          typename TDeltaType>
inline SEQAN_FUNC_ENABLE_IF(Is<StringConcept<TIds> >, bool)
addNode(JournaledStringTree<TSequence, TConfig, TSpec> & jst,
        TPos basePos,
        TValue const & deltaVal,
        TIds const & ids,
        TDeltaType /*deltaType*/)
{
    typedef JournaledStringTree<TSequence, TConfig, TSpec>  TJst;
    typedef typename Member<TJst, JstDeltaMapMember>::Type  TDeltaMap;
    typedef typename DeltaCoverage<TDeltaMap>::Type         TCoverage;
    typedef typename Value<TIds>::Type                      TID;
    typedef typename MakeUnsigned<TPos>::Type               TMaxPos SEQAN_TYPEDEF_FOR_DEBUG;

    SEQAN_ASSERT_LT(static_cast<TMaxPos>(basePos), maxPos(jst));

    // Transform the ids into coverage value.
    TCoverage coverage;
    resize(coverage, depth(jst), false, Exact());

    forEach(ids, [&jst, &coverage](TID seqId)
                    {
                      SEQAN_ASSERT_LT(seqId, depth(jst));  // Check that id is valid.
                      coverage[seqId] = true;
                    }, Serial());

    return insert(jst._deltaMap, basePos, deltaVal, coverage, TDeltaType());
}

// ----------------------------------------------------------------------------
// Function eraseNode()
// ----------------------------------------------------------------------------

/*!
 * @fn JournaledStringTree#eraseNode
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Erases a delta node in the Journaled-String-Tree.
 *
 * @signature bool eraseNode(jst, pos, type);
 *
 * @param[in,out]   jst     The Journaled-String-Tree to erase the node from.
 * @param[in]       pos     The position wihtin the base sequence of the <tt>jst<\tt> of the delta node to be erased.
 * @param[in]       type    A tag representing the type of the delta node. One of @link DeltaTypeTags @endlink. 
 *                          Note there can be multiple nodes at the same <tt>pos<\tt> with different <tt>type<\tt>.
 *
 * @return bool <tt>true<\tt> if the node could be erased, or <tt>false<\tt>, in case there was no node at the 
 * <tt>pos<\tt> with this <tt>type<\tt>.
 *
 * @remark The deletion time is linear in the number of nodes.
 */

template <typename TSequence, typename TConfig, typename TSpec, typename TPos, typename TDeltaType>
inline bool
eraseNode(JournaledStringTree<TSequence, TConfig, TSpec> & jst,
          TPos basePos,
          TDeltaType /*deltaType*/)
{
    return erase(jst._deltaMap, basePos, TDeltaType());
}

}

#endif // EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_JOURNALED_STRING_TREE_IMPL_H_
