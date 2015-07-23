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

// ----------------------------------------------------------------------------
// Class JournaledStringTree                                [StringTreeDefault]
// ----------------------------------------------------------------------------

/*!
 * @class JournaledStringTree Journaled String Tree
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief This container represents a set of strings efficiently in a referentially compressed way.
 *
 * @signature template <typename TSequence[, typename TConfig][, TSpec]>
 *            class JournaledStringTree<TDeltaStore, TConfig, TSpec>;
 *
 * @tparam TSequence Type of underlying base sequence.
 * @tparam TConfig   A configuration objectm, whose traits can be used to specialize the underlying @link DetlaMap @endlink. Defaults to @link DefaultJstConfig @endlink.
 * @tparam TSpec     The specialization tag for the Journaled-String-Tree. Defaults to @link Default @endlink.
 *
 * This data structure stores a set of sequences in a compressed form. The sequences are described in form of a pointer
 * to a common reference sequence and so called delta events, which are stored in a @link DeltaMap @endlink.
 * In addition this data strcuture can be traversed in a data parallel fashion, while identical segments are parsed only
 * once.
 */


// This data structure represents a collection of strings.
// It facilitates an alignment representation of each sequence to a common base sequence.
template <typename TSequence, typename TConfig = DefaultJstConfig<TSequence>, typename TSpec = Default>
class JournaledStringTree
{
public:

    typedef JournaledStringTree<TSequence, TConfig, Default>                        TJst;
    typedef typename Member<TJst, JstDeltaMapMember>::Type                          TDeltaMap;
    typedef typename Member<TJst, JstSourceMember>::Type                            TSource;
    typedef typename Iterator<TSource, Standard>::Type                              TSrcIter;
    typedef typename Size<TJst>::Type                                               TSize;


    TSize       _dimension = 0;         // Number of sequences represented by the JST.
    TSource     _source;            // A journaled String representing the baseSequence.
    TDeltaMap   _map;

    /*!
     * @fn JournaledStringTree::JournaledStringTree
     * @brief Constructor.
     * @headerfile <seqan/journaled_string_tree.h>
     *
     * @signature JournaledStringTree(source, branchLength, dimension);
     * @signature JournaledStringTree(source, branchLength, deltaMap);
     *
     * @param source        The underlying base sequence.
     * @param dimension     The number of sequences represented in this object.
     * @param deltaMap      A delta map of type @link DeltaMap @endlink.
     */

    JournaledStringTree() = delete;

    // Custom constructor.
    template <typename TDim>
    JournaledStringTree(TDim dimension) : _dimension(dimension), _source()
    {}

    // Custom constructor.
    template <typename TDim, typename TSequence2>
    JournaledStringTree(TSequence2 & source, TDim dimension) : _dimension(dimension)
    {
        setHost(_source, source);
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
    typedef typename Member<JournaledStringTree<TSequence, TConfig, TSpec>, JstDeltaMapMember>::Type    TDeltaMap_;
    typedef typename Size<TDeltaMap_>::Type                                                             TDeltaMapSize_;
    typedef typename Member<JournaledStringTree<TSequence, TConfig, TSpec>, JstSourceMember>::Type      TSource_;
    typedef typename Size<TSource_>::Type                                                               TSrcSize_;

    typedef typename If<Eval<sizeof(TDeltaMapSize_) < sizeof(TSrcSize_)>,
                        TSrcSize_,
                        TDeltaMapSize_>::Type Type;  // Size type is the largest size type of delta map or source.
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
    typedef typename Member<JournaledStringTree<TSequence, TConfig, TSpec>, JstDeltaMapMember>::Type    TDeltaMap_;
    typedef typename Position<TDeltaMap_>::Type Type;
};

// ----------------------------------------------------------------------------
// Metafunction Host
// ----------------------------------------------------------------------------

/*!
 * @mfn JournaledStringTree#Host
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief The host type.
 *
 * @signature Host<TJst>::Type;
 *
 * @tparam TJst The type of the journal string tree to get the host type for.
 *
 * @return THost The type of the underlying data structure, which is a @link DeltaMapWrapper @endlink.
 */
template <typename TSeq, typename TConfig, typename TSpec>
struct Host<JournaledStringTree<TSeq, TConfig, TSpec> >
{
    typedef TSeq Type;
};

template <typename TSeq, typename TConfig, typename TSpec>
struct Host<JournaledStringTree<TSeq, TConfig, TSpec> const >
{
    typedef TSeq const Type;
};

// ----------------------------------------------------------------------------
// Metafunction JstSourceMember
// ----------------------------------------------------------------------------

template <typename TSeq, typename TConfig, typename TSpec>
struct Member<JournaledStringTree<TSeq, TConfig, TSpec>, JstSourceMember>
{
    typedef typename Spec<TSeq>::Type                   TSpec_;
    typedef typename Value<TSeq>::Type                  TValue_;
    typedef String<TValue_, Journaled<TSpec_> >         Type;
};

template <typename TSeq, typename TConfig, typename TSpec>
struct Member<JournaledStringTree<TSeq, TConfig, TSpec> const, JstSourceMember>
{
    typedef typename Member<JournaledStringTree<TSeq, TConfig, TSpec>, JstSourceMember>::Type const Type;
};

// ----------------------------------------------------------------------------
// Metafunction Member<TJst, JstBufferMember>
// ----------------------------------------------------------------------------

template <typename TSeq, typename TConfig, typename TSpec>
struct Member<JournaledStringTree<TSeq, TConfig, TSpec>, JstDeltaMapMember>
{
    typedef DeltaMap<TConfig> Type;
};

// TODO(rrahn): Check if there is a default const Member metafunction.
template <typename TSeq, typename TConfig, typename TSpec>
struct Member<JournaledStringTree<TSeq, TConfig, TSpec> const, JstDeltaMapMember>
{
    typedef DeltaMap<TConfig> const Type;
};

// ============================================================================
// Private Functions
// ============================================================================

namespace impl
{

// ----------------------------------------------------------------------------
// Function impl::member();                                   [JstSourceMember]
// ----------------------------------------------------------------------------

template <typename TSequence, typename TConfig, typename TSpec>
inline typename Member<JournaledStringTree<TSequence, TConfig, TSpec>, JstSourceMember>::Type &
member(JournaledStringTree<TSequence, TConfig, TSpec> & jst,
       JstSourceMember const & /*tag*/)
{
    return jst._source;
}

template <typename TSequence, typename TConfig, typename TSpec>
inline typename Member<JournaledStringTree<TSequence, TConfig, TSpec> const, JstSourceMember>::Type &
member(JournaledStringTree<TSequence, TConfig, TSpec> const & jst,
       JstSourceMember const & /*tag*/)
{
    return jst._source;
}

// ----------------------------------------------------------------------------
// Function impl::member();                                 [JstDeltaMapMember]
// ----------------------------------------------------------------------------

template <typename TSequence, typename TConfig, typename TSpec>
inline typename Member<JournaledStringTree<TSequence, TConfig, TSpec>, JstDeltaMapMember>::Type &
member(JournaledStringTree<TSequence, TConfig, TSpec> & jst,
       JstDeltaMapMember const & /*tag*/)
{
    return jst._map;
}

template <typename TSequence, typename TConfig, typename TSpec>
inline typename Member<JournaledStringTree<TSequence, TConfig, TSpec> const, JstDeltaMapMember>::Type &
member(JournaledStringTree<TSequence, TConfig, TSpec> const & jst,
       JstDeltaMapMember const & /*tag*/)
{
    return jst._map;
}

}

// ============================================================================
// Public Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function host()
// ----------------------------------------------------------------------------

/*!
 * @fn JournaledStringTree#host
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Returns the set base sequence.
 *
 * @signature THost host(jst);
 *
 * @param[in] jst The Journal String Tree.
 *
 * @return THost A reference to the host of type @link JournaledStringTree#Host @endlink.
 */

template <typename TSequence, typename TConfig, typename TSpec>
inline typename Host<JournaledStringTree<TSequence, TConfig, TSpec> >::Type &
host(JournaledStringTree<TSequence, TConfig, TSpec> & jst)
{
    return host(jst._source);
}

template <typename TSequence, typename TConfig, typename TSpec>
inline typename Host<JournaledStringTree<TSequence, TConfig, TSpec> const>::Type &
host(JournaledStringTree<TSequence, TConfig, TSpec> const & jst)
{
    return host(jst._source);
}

// ----------------------------------------------------------------------------
// Function setHost()
// ----------------------------------------------------------------------------

/*!
 * @fn JournaledStringTree#setHost
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Sets the base sequence of the Journaled-String-Tree.
 *
 * @signature void setHost(jst, host);
 *
 * @param[in,out] jst   The Journal String Tree.
 * @param[in]     host  The new host to set.
 *
 * When setting a new host, the <tt>jst</tt> will be cleared before.
 */

template <typename TSequence, typename TConfig, typename TSpec, typename THost>
inline void
setHost(JournaledStringTree<TSequence, TConfig, TSpec> & jst,
        THost & host)
{
    clear(jst);
    setHost(jst._source, host);
}

// ----------------------------------------------------------------------------
// Function length()
// ----------------------------------------------------------------------------

/*!
 * @fn JournaledStringTree#length
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Returns the number of sequences represented by the Journaled-String-Tree.
 *
 * @signature TSize length(jst);
 *
 * @param[in] jst The Journal-String-Tree to query the length for.
 *
 * @return TSize The number of sequences represented by the <tt>jst<\tt>. Of type @link JournaledStringTree#Size @endlink.
 */

template <typename TSequence, typename TConfig, typename TSpec>
inline typename Size<JournaledStringTree<TSequence, TConfig, TSpec> const>::Type
length(JournaledStringTree<TSequence, TConfig, TSpec> const & jst)
{
    return jst._dimension;
}

// ----------------------------------------------------------------------------
// Function maxSize()
// ----------------------------------------------------------------------------

/*
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
maxSize(JournaledStringTree<TSequence, TConfig, TSpec> const & /*jst*/)
{
    return maxSize(typename Member<JournaledStringTree<TSequence, TConfig, TSpec>, JstDeltaMapMember>::Type());
}

// ----------------------------------------------------------------------------
// Function size()
// ----------------------------------------------------------------------------

/*
 * @fn JournaledStringTree#size
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Returns the number of delta operations stored in the Journaled-String-Tree.
 *
 * @signature TSize maxSize(jst);
 *
 * @param[in] jst The Journal String Tree to query the maximal length for.
 *
 * @return TSize The maximal size. Of type @link JournaledStringTree#Size @endlink.
 */

template <typename TSequence, typename TConfig, typename TSpec>
inline typename Size<JournaledStringTree<TSequence, TConfig, TSpec> >::Type
size(JournaledStringTree<TSequence, TConfig, TSpec> const jst)
{
    return size(jst._map);
}

// ----------------------------------------------------------------------------
// Function empty()
// ----------------------------------------------------------------------------

/*
 * @fn JournaledStringTree#size
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Returns the delta operations stored in the Journaled-String-Tree.
 *
 * @signature TSize maxSize(jst);
 *
 * @param[in] jst The Journal String Tree to query the maximal length for.
 *
 * @return TSize The maximal size. Of type @link JournaledStringTree#Size @endlink.
 */

template <typename TSequence, typename TConfig, typename TSpec>
inline bool
empty(JournaledStringTree<TSequence, TConfig, TSpec> const jst)
{
    return empty(jst._map);
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
 * Note this function invalidates all other clients that depend on the host of this
 * object.
 */

template <typename TSequence, typename TConfig, typename TSpec>
inline void
clear(JournaledStringTree<TSequence, TConfig, TSpec> & jst)
{
    jst._dimension = 0;
    reset(jst._source);
    clear(jst._map);
}

// ----------------------------------------------------------------------------
// Function insert();
// ----------------------------------------------------------------------------

/*!
 * @fn JournaledStringTree#insertNode
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
 * @note Inserting a new delta might invalidate other clients that depend on the same host.
 *
 * @remark The insertion time is linear in the number of nodes.
 * @see JournaledStringTree#eraseNode
 */

template <typename TSequence, typename TConfig, typename TSpec, typename TPos, typename TValue, typename TIds,
          typename TDeltaType>
inline SEQAN_FUNC_ENABLE_IF(Is<StringConcept<TIds> >, void)
insert(JournaledStringTree<TSequence, TConfig, TSpec> & jst,
       TPos srcPos,
       TValue const & deltaVal,
       TIds const & ids,
       TDeltaType /*deltaType*/)
{
    typedef JournaledStringTree<TSequence, TConfig, TSpec>  TJst;
    typedef typename Member<TJst, JstDeltaMapMember>::Type  TDeltaMap;
    typedef typename DeltaCoverage<TDeltaMap>::Type         TCoverage;
    typedef typename Value<TIds>::Type                      TID;
    typedef typename Size<TJst>::Type                       TSize   SEQAN_TYPEDEF_FOR_DEBUG;

    if (IsSameType<TDeltaType, DeltaTypeIns>::VALUE)
        SEQAN_ASSERT_LEQ(static_cast<TSize>(srcPos), length(impl::member(jst, JstSourceMember())));  // Make sure the delta position does not exceed the source.
    else
        SEQAN_ASSERT_LT(static_cast<TSize>(srcPos), length(impl::member(jst, JstSourceMember())));
    // Transform the ids into coverage value.
    TCoverage coverage;
    resize(coverage, length(jst), false, Exact());

    forEach(ids,[&jst, &coverage](TID seqId)
    {
        SEQAN_ASSERT_LT(static_cast<TSize>(seqId), length(jst));  // Check that id is valid.
        coverage[seqId] = true;
    });

    insert(impl::member(jst, JstDeltaMapMember()), srcPos, deltaVal, coverage, TDeltaType());
}

// ----------------------------------------------------------------------------
// Function erase();
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
 * @note Inserting a new delta might invalidate other clients that depend on the same host.
 *
 * @remark The deletion time is linear in the number of nodes.
 * @see JournaledStringTree#insertNode
 */

template <typename TSequence, typename TConfig, typename TSpec, typename TPos, typename TDeltaType>
inline auto
erase(JournaledStringTree<TSequence, TConfig, TSpec> & jst,
      TPos srcPos,
      TDeltaType /*deltaType*/) -> decltype(erase(impl::member(jst, JstDeltaMapMember()), srcPos, TDeltaType()))
{
    return erase(impl::member(jst, JstDeltaMapMember()), srcPos, TDeltaType());
}

// TODO(rrahn): Implement emplace when needed.
// TODO(rrahn): Implement emplace_hint when needed.
// TODO(rrahn): Implement swap when needed.

}

#endif // EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_JOURNALED_STRING_TREE_IMPL_H_
