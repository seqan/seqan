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

struct JstBufferMember_;
typedef Tag<JstBufferMember_> JstBufferMember;

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
    typedef typename Container<TJst>::Type                      TDeltaMap;
    typedef typename Source<TJst>::Type                         TSource;
    typedef typename Size<TJst>::Type                           TSize;
//    typedef typename Member<TJst, JstBufferMember>::Type        TBuffer;

    TSize               _historySize;       // The history size which is used to prune the branches after the historySize is reached.
    TSize               _dimension;         // Number of sequences represented by the JST.
    TSource             _source;            // A journaled String representin the baseSequence.
    Holder<TDeltaMap>   _deltaMapHolder;    // The wrapped container.
//    TBuffer         _jsBuffer;          // The buffer representing the current sequence context to iterate through.

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

    JournaledStringTree() = delete;

    // Custom constructor.
    template <typename TSeq, typename TSize>
    JournaledStringTree(TSeq const & source, TSize historySize, TSize dimension) :
        _historySize(historySize),
        _dimension(dimension)
    {
        setHost(_source, source);
        create(_deltaMapHolder);  // Makes delta map holder owner.
    }

    template <typename TSeq, typename TSize>
    JournaledStringTree(TSeq const &source, TSize historySize, TDeltaMap & deltaMap) :
        _historySize(historySize)
    {
        SEQAN_ASSERT(!empty(deltaMap));  // TODO(rrahn): Use exception instead.

        setHost(_source, source);
        setValue(_deltaMapHolder, deltaMap);  // Makes delta map dependent.
        // Receive the dimension from the coverage assuming same coverage size for all mapped delta operations.
        _dimension = length(deltaCoverage(begin(container(*this))));
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

// ----------------------------------------------------------------------------
// Metafunction Container
// ----------------------------------------------------------------------------

/*!
 * @mfn JournaledStringTree#Container
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief The container type.
 *
 * @signature Container<TJst>::Type;
 *
 * @tparam TJst The type of the journal string tree to get the container type for.
 *
 * @return TContainer The type of the container, which is a @link DeltaMap @endlink.
 */
template <typename TSeq, typename TConfig, typename TSpec>
struct Container<JournaledStringTree<TSeq, TConfig, TSpec> >
{
    typedef DeltaMap<TConfig> Type;
};

// ----------------------------------------------------------------------------
// Metafunction Source
// ----------------------------------------------------------------------------

/*!
 * @mfn JournaledStringTree#Source
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief The source sequence type.
 *
 * @signature Source<TJst>::Type;
 *
 * @tparam TJst The type of the journal string tree to get the  source type for.
 *
 * @return TSource The type of the underlying source sequence, which is represented as a @link JournaledString @endlink.
 */
template <typename TSeq, typename TConfig, typename TSpec>
struct Source<JournaledStringTree<TSeq, TConfig, TSpec> >
{
    typedef typename Spec<TSeq>::Type           TSpec_;
    typedef typename Value<TSeq>::Type          TValue_;
    typedef String<TValue_, Journaled<TSpec_> > Type;
};

template <typename TSeq, typename TConfig, typename TSpec>
struct Source<JournaledStringTree<TSeq, TConfig, TSpec> const>
{
    typedef JournaledStringTree<TSeq, TConfig, TSpec> TJst;
    typedef typename Source<TJst>::Type const Type;
};

// ----------------------------------------------------------------------------
// Metafunction Member<TJst, JstBufferMember>
// ----------------------------------------------------------------------------

//template <typename TJst>
//struct Member<TJst, JstBufferMember>
//{
//    typedef JstBuffer<TJst, Forw> Type;
//};

// ============================================================================
// Private Functions
// ============================================================================

namespace impl
{

}

// ============================================================================
// Public Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function source()
// ----------------------------------------------------------------------------

/*!
 * @fn JournaledStringTree#baseSeq
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Returns a reference to the base sequence used for the journaled string tree.
 *
 * @signature TSource source(jst);
 *
 * @param[in] jst    The Journal String Tree.
 *
 * @return TSource A reference to the source sequence of type @link JournaledStringTree#Source @endlink.
 *
 * Note that the returned value is a @link JournaledString @endlink, which is a hosted type and might depend
 * on another resource. So changes applied to the host of this string might invalidate other clients that depend
 * on the same resource.
 */

template <typename TSequence, typename TConfig, typename TSpec>
inline typename Source<JournaledStringTree<TSequence, TConfig, TSpec> >::Type &
source(JournaledStringTree<TSequence, TConfig, TSpec> & jst)
{
    return jst._source;
}

template <typename TSequence, typename TConfig, typename TSpec>
inline typename Source<JournaledStringTree<TSequence, TConfig, TSpec> const>::Type &
source(JournaledStringTree<TSequence, TConfig, TSpec> const & jst)
{
    return jst._source;
}

// ----------------------------------------------------------------------------
// Function container()
// ----------------------------------------------------------------------------

/*!
 * @fn JournaledStringTree#baseSeq
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Returns a reference to the base sequence used for the journaled string tree.
 *
 * @signature TSource source(jst);
 *
 * @param[in] jst    The Journal String Tree.
 *
 * @return TSource A reference to the source sequence of type @link JournaledStringTree#Source @endlink.
 *
 * Note that the returned value is a @link JournaledString @endlink, which is a hosted type and might depend
 * on another resource. So changes applied to the host of this string might invalidate other clients that depend
 * on the same resource.
 */

template <typename TSequence, typename TConfig, typename TSpec>
inline typename Container<JournaledStringTree<TSequence, TConfig, TSpec> >::Type &
container(JournaledStringTree<TSequence, TConfig, TSpec> & jst)
{
    return value(jst._deltaMapHolder);
}

template <typename TSequence, typename TConfig, typename TSpec>
inline typename Container<JournaledStringTree<TSequence, TConfig, TSpec> const>::Type &
container(JournaledStringTree<TSequence, TConfig, TSpec> const & jst)
{
    return value(jst._deltaMapHolder);
}

// ----------------------------------------------------------------------------
// Function dimension()
// ----------------------------------------------------------------------------

/*!
 * @fn JournaledStringTree#dimension
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Returns the number of sequences represented by the Journaled-String-Tree.
 *
 * @signature TSize dimension(jst);
 *
 * @param[in] jst The Journal String Tree to query the dimension for.
 *
 * @return TSize The number of sequences represented by the <tt>jst<\tt>. Of type @link JournaledStringTree#Size @endlink.
 */

template <typename TSequence, typename TConfig, typename TSpec>
inline typename Size<JournaledStringTree<TSequence, TConfig, TSpec> >::Type
dimension(JournaledStringTree<TSequence, TConfig, TSpec> & jst)
{
    return jst._dimension;
}

template <typename TSequence, typename TConfig, typename TSpec>
inline typename Size<JournaledStringTree<TSequence, TConfig, TSpec> const>::Type
dimension(JournaledStringTree<TSequence, TConfig, TSpec> const & jst)
{
    return jst._dimension;
}

// ----------------------------------------------------------------------------
// Function historySize()
// ----------------------------------------------------------------------------

/*!
 * @fn JournaledStringTree#historySize
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Returns the history size of the Journaled-String-Tree.
 *
 * @signature TSize historySize(jst);
 *
 * @param[in] jst The Journal String Tree to query the history size for.
 *
 * @return TSize The history size used to prune branches in order to reprent a tree. Of type @link JournaledStringTree#Size @endlink.
 */

template <typename TSequence, typename TConfig, typename TSpec>
inline typename Size<JournaledStringTree<TSequence, TConfig, TSpec> >::Type
historySize(JournaledStringTree<TSequence, TConfig, TSpec> & jst)
{
    return jst._historySize;
}

template <typename TSequence, typename TConfig, typename TSpec>
inline typename Size<JournaledStringTree<TSequence, TConfig, TSpec> const>::Type
historySize(JournaledStringTree<TSequence, TConfig, TSpec> const & jst)
{
    return jst._historySize;
}

// ----------------------------------------------------------------------------
// Function setHistorySize()
// ----------------------------------------------------------------------------

/*!
 * @fn JournaledStringTree#setHistorySize
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Sets the history size for the Journaled-String-Tree.
 *
 * @signature void historySize(jst, h);
 *
 * @param[in,out] jst The Journal String Tree to set the new history size for.
 * @param[in]     h   The new history size.
 *
 * @note Any iterators, pointers or reference to this object may be invalidated.
 */

template <typename TSequence, typename TConfig, typename TSpec, typename TSize>
inline void
setHistorySize(JournaledStringTree<TSequence, TConfig, TSpec> & jst,
               TSize const h)
{
    jst._historySize = h;
}

// ----------------------------------------------------------------------------
// Function maxSize()
// ----------------------------------------------------------------------------

/*!
 * @fn JournaledStringTree#maxSize
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Returns the maximal length the the Journaled-String-Tree can reach.
 *
 * @signature TSize maxSize(jst);
 *
 * @param[in] jst The Journal String Tree to query the maximal length for.
 *
 * @return TSize The maximal size. Of type @link JournaledStringTree#Size @endlink.
 */

template <typename TSequence, typename TConfig, typename TSpec>
constexpr typename Size<JournaledStringTree<TSequence, TConfig, TSpec> >::Type
maxSize(JournaledStringTree<TSequence, TConfig, TSpec> const /*jst*/)
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
 *
 * Note this function invalidates all other clients that depend on the container of this
 * object.
 */

template <typename TSequence, typename TConfig, typename TSpec>
inline void
clear(JournaledStringTree<TSequence, TConfig, TSpec> & jst)
{
    jst._dimension = 0;
    jst._historySize = 0;
    reset(jst._source);
    clear(jst._deltaMapHolder);
//    clear(jst._buffer);
}

// ----------------------------------------------------------------------------
// Function insertNode()
// ----------------------------------------------------------------------------

/*!
 * @fn JournaledStringTree#addNode
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Adds a new delta node to the Journaled-String-Tree.
 *
 * @signature bool insertNode(jst, pos, val, ids, type);
 *
 * @param[in,out]   jst     The Journaled-String-Tree to insert the new node to.
 * @param[in]       pos     The mapped source position at which the new delta node is placed. This position must not exceed @link JournaledStringTree#maxSize @endlink.
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
 * @note Inserting a new delta might invalidate other clients that depend on the same container.
 *
 * @remark The insertion time is linear in the number of nodes.
 */

template <typename TSequence, typename TConfig, typename TSpec, typename TPos, typename TValue, typename TIds,
          typename TDeltaType>
inline SEQAN_FUNC_ENABLE_IF(Is<StringConcept<TIds> >, bool)
insertNode(JournaledStringTree<TSequence, TConfig, TSpec> & jst,
           TPos srcPos,
           TValue const & deltaVal,
           TIds const & ids,
           TDeltaType /*deltaType*/)
{
    typedef JournaledStringTree<TSequence, TConfig, TSpec>  TJst;
    typedef typename Container<TJst>::Type                  TDeltaMap;
    typedef typename DeltaCoverage<TDeltaMap>::Type         TCoverage;
    typedef typename Value<TIds>::Type                      TID;
    typedef typename MakeUnsigned<TPos>::Type               TMaxPos SEQAN_TYPEDEF_FOR_DEBUG;

    SEQAN_ASSERT_LT(static_cast<TMaxPos>(srcPos), maxSize(jst));

    // Transform the ids into coverage value.
    TCoverage coverage;
    resize(coverage, dimension(jst), false, Exact());

    forEach(ids, [&jst, &coverage](TID seqId)
                    {
                      SEQAN_ASSERT_LT(seqId, dimension(jst));  // Check that id is valid.
                      coverage[seqId] = true;
                    }, Serial());

    return insert(container(jst), srcPos, deltaVal, coverage, TDeltaType());
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
 * @param[in]       pos     The mapped source position of the delta node to be erased.
 * @param[in]       type    A tag representing the type of the delta node. One of @link DeltaTypeTags @endlink. 
 *                          Note there can be multiple nodes at the same <tt>pos<\tt> with different <tt>type<\tt>.
 *
 * @return bool <tt>true<\tt> if the node could be erased, or <tt>false<\tt>, in case there was no node at the 
 * <tt>pos<\tt> with this <tt>type<\tt>.
 *
 * @note Inserting a new delta might invalidate other clients that depend on the same container.
 *
 * @remark The deletion time is linear in the number of nodes.
 */

template <typename TSequence, typename TConfig, typename TSpec, typename TPos, typename TDeltaType>
inline bool
eraseNode(JournaledStringTree<TSequence, TConfig, TSpec> & jst,
          TPos srcPos,
          TDeltaType /*deltaType*/)
{
    return erase(container(jst), srcPos, TDeltaType());
}

}

#endif // EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_JOURNALED_STRING_TREE_IMPL_H_
