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
    typedef typename Host<TJst>::Type                                               TDeltaMap;
    typedef typename Member<TJst, JstSourceMember>::Type                            TSource;
    typedef typename Iterator<TSource, Standard>::Type                              TSrcIter;
    typedef typename Size<TJst>::Type                                               TSize;
    typedef typename Member<TJst, JstBufferMember>::Type                            TBuffer;


    TSize                       _branchLength;      // The history size which is used to prune the branches after the historySize is reached.
    TSize                       _dimension;         // Number of sequences represented by the JST.
    TSource                     _source;            // A journaled String representing the baseSequence.

    Holder<TDeltaMap>           _mapHolder;
    TBuffer                     _buffer;            // The buffer represenging the journaled stings.

    /*!
     * @fn JournaledStringTree::JournaledStringTree
     * @brief Constructor.
     * @headerfile <seqan/journaled_string_tree.h>
     *
     * @signature JournaledStringTree(source, branchLength, dimension);
     * @signature JournaledStringTree(source, branchLength, deltaMap);
     *
     * @param source        The underlying base sequence.
     * @param branchLength  The maximal length of a branch. (Prunes the branch, i.e. defines the maximum size of the sequence context).
     * @param dimension     The number of sequences represented in this object.
     * @param deltaMap      A delta map of type @link DeltaMap @endlink.
     */

    // Custom constructor.
    template <typename TSeq, typename TDim>
    JournaledStringTree(TSeq const & source, TSize branchLength, TDim dimension) :
        _branchLength(branchLength),
        _dimension(dimension),
        _buffer()
    {
        setHost(_source, source);
        create(_mapHolder);
        setDeltaMap(_buffer, value(_mapHolder));
        setSourceBegin(_buffer, begin(_source, Standard()));
        setSourceEnd(_buffer, end(_source, Standard()));
    }

    template <typename TSeq, typename TSize>
    JournaledStringTree(TSeq const & source, TSize branchLength, TDeltaMap & deltaMap) :
        _branchLength(branchLength),
        _buffer()
    {
        SEQAN_ASSERT(!empty(deltaMap));  // TODO(rrahn): Use exception instead.

        setHost(_source, source);
        setValue(_mapHolder, deltaMap);
        setDeltaMap(_buffer, value(_mapHolder));
        setSourceBegin(_buffer, begin(_source, Standard()));
        setSourceEnd(_buffer, end(_source, Standard()));

        _dimension = length(getDeltaCoverage(*begin(host(*this))));
        setDeltaMap(_buffer, host(*this));
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
    typedef DeltaMap<TConfig> Type;
};

template <typename TSeq, typename TConfig, typename TSpec>
struct Host<JournaledStringTree<TSeq, TConfig, TSpec> const >
{
    typedef DeltaMap<TConfig> const Type;
};

// ----------------------------------------------------------------------------
// Metafunction Source
// ----------------------------------------------------------------------------

/*
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
struct Member<JournaledStringTree<TSeq, TConfig, TSpec>, JstSourceMember>
{
    typedef typename Spec<TSeq>::Type           TSpec_;
    typedef typename Value<TSeq>::Type          TValue_;
    typedef String<TValue_, Journaled<TSpec_> > Type;
};

template <typename TSeq, typename TConfig, typename TSpec>
struct Member<JournaledStringTree<TSeq, TConfig, TSpec> const, JstSourceMember>
{
    typedef JournaledStringTree<TSeq, TConfig, TSpec> TJst;
    typedef typename Source<TJst>::Type const Type;
};

// ----------------------------------------------------------------------------
// Metafunction Member<TJst, JstBufferMember>
// ----------------------------------------------------------------------------

template <typename TJst>
struct Member<TJst, JstBufferMember>
{
    typedef JstBuffer_<TJst> Type;
};

// ============================================================================
// Private Functions
// ============================================================================

namespace impl
{

// ----------------------------------------------------------------------------
// Function impl::source();
// ----------------------------------------------------------------------------

template <typename TSequence, typename TConfig, typename TSpec>
inline typename Member<JournaledStringTree<TSequence, TConfig, TSpec>, JstSourceMember>::Type &
source(JournaledStringTree<TSequence, TConfig, TSpec> & jst)
{
    return jst._source;
}

template <typename TSequence, typename TConfig, typename TSpec>
inline typename Member<JournaledStringTree<TSequence, TConfig, TSpec> const, JstSourceMember>::Type &
source(JournaledStringTree<TSequence, TConfig, TSpec> const & jst)
{
    return jst._source;
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
 * @brief Returns a reference to the underlying host holding the delta information of the journaled string tree.
 *
 * @signature THost host(jst);
 *
 * @param[in] jst    The Journal String Tree.
 *
 * @return THost A reference to the host of type @link JournaledStringTree#Host @endlink.
 *
 * Note that other objects might depent on the returned host. Any iterator, pointer or reference that refer to this
 * host may get invalidated if it is changed.
 */

template <typename TSequence, typename TConfig, typename TSpec>
inline typename Host<JournaledStringTree<TSequence, TConfig, TSpec> >::Type &
host(JournaledStringTree<TSequence, TConfig, TSpec> & jst)
{
    return value(jst._mapHolder);
}

template <typename TSequence, typename TConfig, typename TSpec>
inline typename Host<JournaledStringTree<TSequence, TConfig, TSpec> const>::Type &
host(JournaledStringTree<TSequence, TConfig, TSpec> const & jst)
{
    return value(jst._mapHolder);
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
inline typename Size<JournaledStringTree<TSequence, TConfig, TSpec> const>::Type
dimension(JournaledStringTree<TSequence, TConfig, TSpec> const & jst)
{
    return jst._dimension;
}

// ----------------------------------------------------------------------------
// Function branchLength()
// ----------------------------------------------------------------------------

/*!
 * @fn JournaledStringTree#branchLength
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Returns the branch length of the Journaled-String-Tree.
 *
 * @signature TSize branchLength(jst);
 *
 * @param[in] jst The object to query the branch length for.
 *
 * The branch length parameter indicates the length of the spawned branch in the <tt>jst</tt>, which restricts the
 * maximal sequence context size to this value.
 *
 * @return TSize The branch length of type @link JournaledStringTree#Size @endlink.
 */

template <typename TSequence, typename TConfig, typename TSpec>
inline typename Size<JournaledStringTree<TSequence, TConfig, TSpec> const>::Type
branchLength(JournaledStringTree<TSequence, TConfig, TSpec> const & jst)
{
    return jst._branchLength;
}

// ----------------------------------------------------------------------------
// Function setBranchLength()
// ----------------------------------------------------------------------------

/*!
 * @fn JournaledStringTree#setBranchLength
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Sets the branch length for the Journaled-String-Tree.
 *
 * @signature void setBranchLength(jst, l);
 *
 * @param[in,out] jst The object to set the new branch length for.
 * @param[in]     l   The new branch length. Must implement the @link IntegerConcept @endlink.
 *
 * @note After changing the branch length, the <tt>jst</jst> needs to be recreated with @link JournaledStringTree#create 
 * @endlink, otherwise access to the object may result in undefined behaviour.
 */

template <typename TSequence, typename TConfig, typename TSpec, typename TSize>
inline void
setBranchLength(JournaledStringTree<TSequence, TConfig, TSpec> & jst,
                TSize const l)
{
    jst._branchLength = l;
    markModified(jst._buffer);
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
 * Note this function invalidates all other clients that depend on the host of this
 * object.
 */

template <typename TSequence, typename TConfig, typename TSpec>
inline void
clear(JournaledStringTree<TSequence, TConfig, TSpec> & jst)
{
    jst._dimension = 0;
    jst._branchLength = 0;
    reset(jst._source);
    clear(jst._mapHolder);
    clear(jst._buffer);
}

// ----------------------------------------------------------------------------
// Function insertNode()
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
inline SEQAN_FUNC_ENABLE_IF(Is<StringConcept<TIds> >, bool)
insertNode(JournaledStringTree<TSequence, TConfig, TSpec> & jst,
           TPos srcPos,
           TValue const & deltaVal,
           TIds const & ids,
           TDeltaType /*deltaType*/)
{
    typedef JournaledStringTree<TSequence, TConfig, TSpec>  TJst;
    typedef typename Host<TJst>::Type                       TDeltaMap;
    typedef typename DeltaCoverage<TDeltaMap>::Type         TCoverage;
    typedef typename Value<TIds>::Type                      TID;
    typedef typename MakeUnsigned<TPos>::Type               TMaxPos SEQAN_TYPEDEF_FOR_DEBUG;
    typedef typename Size<TJst>::Type                       TSize   SEQAN_TYPEDEF_FOR_DEBUG;

    SEQAN_ASSERT_LT(static_cast<TMaxPos>(srcPos), maxSize(jst));

    // Transform the ids into coverage value.
    TCoverage coverage;
    resize(coverage, dimension(jst), false, Exact());

    forEach(ids,[&jst, &coverage](TID seqId)
    {
        SEQAN_ASSERT_LT(static_cast<TSize>(seqId), dimension(jst));  // Check that id is valid.
        coverage[seqId] = true;
    });
    markModified(jst._buffer);
    return insert(host(jst), srcPos, deltaVal, coverage, TDeltaType());
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
 * @note Inserting a new delta might invalidate other clients that depend on the same host.
 *
 * @remark The deletion time is linear in the number of nodes.
 * @see JournaledStringTree#insertNode
 */

template <typename TSequence, typename TConfig, typename TSpec, typename TPos, typename TDeltaType>
inline bool
eraseNode(JournaledStringTree<TSequence, TConfig, TSpec> & jst,
          TPos srcPos,
          TDeltaType /*deltaType*/)
{
    markModified(jst._buffer);
    return erase(host(jst), srcPos, TDeltaType());
}

// ----------------------------------------------------------------------------
// Function create()
// ----------------------------------------------------------------------------

/*!
 * @fn JournaledStringTree#create
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Creates the sequence context.
 *
 * @signature bool create(jst);
 *
 * @param[in,out]   jst     The object to be created.
 *
 * After modifying the <tt>jst</tt> it needs to be recreated before it can be used.
 *
 * @return bool <tt>true<\tt> if the tree could be created, <tt>false</tt> otherwise.
 */

template <typename TSequence, typename TConfig, typename TSpec>
inline bool
create(JournaledStringTree<TSequence, TConfig, TSpec> & jst)
{
    return create(jst._buffer);
}

// TODO(rrahn): Implement me!
//// ----------------------------------------------------------------------------
//// Function setSourceBegin()
//// ----------------------------------------------------------------------------
//
//template <typename TSequence, typename TConfig, typename TSpec,
//          typename TPosition>
//inline void
//setSourceBegin(JournaledStringTree<TSequence, TConfig, TSpec> & jst,
//               TPosition const beginPos)
//{
//    setSourceBegin(jst._buffer, iter(jst._source, beginPos, Standard()));
//}
//
//// ----------------------------------------------------------------------------
//// Function setSourceEnd()
//// ----------------------------------------------------------------------------
//
//template <typename TSequence, typename TConfig, typename TSpec,
//          typename TPosition>
//inline void
//setSourceEnd(JournaledStringTree<TSequence, TConfig, TSpec> & jst,
//             TPosition const endPos)
//{
//    setSourceEnd(jst._buffer, iter(jst._source, endPos, Standard()));
//}

//// ----------------------------------------------------------------------------
//// Function setSourceBegin()
//// ----------------------------------------------------------------------------
//
//// TODO(rrahn): Consider only const based version.
//template <typename TSequence, typename TConfig, typename TSpec>
//inline typename Iterator<typename Source<JournaledStringTree<TSequence, TConfig, TSpec> >::Type, Standard>::Type &
//sourceBegin(JournaledStringTree<TSequence, TConfig, TSpec> & jst)
//{
//    return jst._sourceBegin;
//}
//
//template <typename TSequence, typename TConfig, typename TSpec>
//inline typename Iterator<typename Source<JournaledStringTree<TSequence, TConfig, TSpec> const>::Type, Standard>::Type &
//sourceBegin(JournaledStringTree<TSequence, TConfig, TSpec> const & jst)
//{
//    return jst._sourceBegin;
//}

//// ----------------------------------------------------------------------------
//// Function sourceEnd()
//// ----------------------------------------------------------------------------
//
//template <typename TSequence, typename TConfig, typename TSpec>
//inline typename Iterator<typename Source<JournaledStringTree<TSequence, TConfig, TSpec> >::Type, Standard>::Type &
//sourceEnd(JournaledStringTree<TSequence, TConfig, TSpec> & jst)
//{
//    return jst._sourceEnd;
//}
//
//template <typename TSequence, typename TConfig, typename TSpec>
//inline typename Iterator<typename Source<JournaledStringTree<TSequence, TConfig, TSpec> const>::Type, Standard>::Type &
//sourceEnd(JournaledStringTree<TSequence, TConfig, TSpec> const & jst)
//{
//    return jst._sourceEnd;
//}
//
//// ----------------------------------------------------------------------------
//// Function sourceEndPosition()
//// ----------------------------------------------------------------------------
//
//template <typename TSequence, typename TConfig, typename TSpec>
//inline typename Position<JournaledStringTree<TSequence, TConfig, TSpec> >::Type
//sourceBegin(JournaledStringTree<TSequence, TConfig, TSpec> const & jst)
//{
//    return position(sourceBegin(jst._buffer));
//}
//
//// ----------------------------------------------------------------------------
//// Function sourceEndPosition()
//// ----------------------------------------------------------------------------
//
//template <typename TSequence, typename TConfig, typename TSpec>
//inline typename Position<JournaledStringTree<TSequence, TConfig, TSpec> >::Type
//sourceEnd(JournaledStringTree<TSequence, TConfig, TSpec> const & jst)
//{
//    return position(sourceEnd(jst._buffer));
//}

// TODO(rrahn): Why do we need a traverser function? We can simply create one. We could create an Iterator, which
// wrapps an traverser with context size 1.

//// ----------------------------------------------------------------------------
//// Function traverser()
//// ----------------------------------------------------------------------------
//
//template <typename TSequence, typename TConfig, typename TSpec,
//          typename TSize,
//          typename TObserver>
//inline typename Traverser<JournaledStringTree<TSequence, TConfig, TSpec>, TObserver>::Type
//traverser(JournaledStringTree<TSequence, TConfig, TSpec> & jst,
//          TSize const contextSize,
//          TObserver & observer)
//{
//    SEQAN_ASSERT_GEQ(contextSize, static_cast<TSize>(1));
//    return typename Traverser<JournaledStringTree<TSequence, TConfig, TSpec>, TObserver>::Type(jst, contextSize, observer);
//}
//
//template <typename TSequence, typename TConfig, typename TSpec,
//          typename TSize,
//          typename TObserver>
//inline typename Traverser<JournaledStringTree<TSequence, TConfig, TSpec> const, TObserver>::Type
//traverser(JournaledStringTree<TSequence, TConfig, TSpec> const & jst,
//          TSize const contextSize,
//          TObserver & observer)
//{
//    SEQAN_ASSERT_GEQ(contextSize, static_cast<TSize>(1));
//    return typename Traverser<JournaledStringTree<TSequence, TConfig, TSpec> const, TObserver>::Type(jst, contextSize, observer);
//}
//
//template <typename TSequence, typename TConfig, typename TSpec,
//          typename TSize>
//inline typename Traverser<JournaledStringTree<TSequence, TConfig, TSpec>, void>::Type
//traverser(JournaledStringTree<TSequence, TConfig, TSpec> & jst,
//          TSize const contextSize)
//{
//    return typename Traverser<JournaledStringTree<TSequence, TConfig, TSpec>, void>::Type(jst);
//}
//
//template <typename TSequence, typename TConfig, typename TSpec>
//inline typename Traverser<JournaledStringTree<TSequence, TConfig, TSpec> const>::Type
//traverser(JournaledStringTree<TSequence, TConfig, TSpec> const & jst)
//{
//    return typename Traverser<JournaledStringTree<TSequence, TConfig, TSpec> const>::Type(jst);
//}

}

#endif // EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_JOURNALED_STRING_TREE_IMPL_H_
