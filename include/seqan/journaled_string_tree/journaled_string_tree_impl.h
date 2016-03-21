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
 * @brief This data structure represents a virtual set of strings using a veritcal compression scheme.
 *
 * @signature template <typename TSequence[, typename TConfig][, TSpec]>
 *            class JournaledStringTree<TDeltaStore, TConfig, TSpec>;
 *
 * @tparam TSequence Type of underlying base sequence.
 * @tparam TConfig   A configuration object, whose traits can be used to specialize the underlying @link DeltaMap @endlink. Defaults to @link DefaultJstConfig @endlink.
 * @tparam TSpec     The specialization tag for the journaled string tree. Defaults to @link Default @endlink.
 *
 * This data structure stores a virtual set of sequences in a compressed form. The sequences are described in form of a pointer
 * to a common reference sequence and in addition delta events, which are stored in a @link DeltaMap @endlink.
 * The interface is a hybrid between a container interface and a associative container interface. So it uses typical 
 * string set intrefaces to access the virtual strings. On the other hand it implements an associative container to
 * store the delta events given their reference position as a key.
 * This data structure can then be traversed efficiently, while exploring only those sequence contexts, that are unique
 * to the set of sequences. This means, that redundant parts, that are shared by one or many sequences is traversed only
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

    TSize       _dimension = 0;     // Number of sequences represented by the JST.
    TSource     _source;            // A journaled String representing the baseSequence.
    TDeltaMap   _map;

    /*!
     * @fn JournaledStringTree::JournaledStringTree
     * @brief Constructor.
     * @headerfile <seqan/journaled_string_tree.h>
     *
     * @signature JournaledStringTree();
     * @signature JournaledStringTree(size);
     * @signature JournaledStringTree(source, size);
     *
     * @param size     The number of sequences represented by the JST.
     * @param source   The underlying base sequence to which the sequences are compressed.
     */

    // Def c'tor.
    JournaledStringTree()
    {}

    // Custom constructor.
    template <typename TDim>
    JournaledStringTree(TDim const dimension) : _dimension(dimension), _source()
    {}

    // Custom constructor.
    template <typename TSequence2, typename TDim>
    JournaledStringTree(TSequence2 && source, TDim const dimension) : _dimension(dimension)
    {
        setHost(_source, std::forward<TSequence2>(source));
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
 * @tparam TJst The JST type to query the specialization type for.
 *
 * @return TSpec The specialization type.
 */

template <typename TSequence, typename TConfig, typename TSpec>
struct Spec<JournaledStringTree<TSequence, TConfig, TSpec> >
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
 * @tparam TJst The JST type to get the size type for.
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
 * @tparam TJst The JST type to get the position type for.
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
 * @tparam TJst The JST type to get the host type for.
 *
 * @return THost The type of the underlying data structure.
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

//// ----------------------------------------------------------------------------
//// Function impl::openJstFromFile();
//// ----------------------------------------------------------------------------
//
//template <typename TSequence, typename TConfig, typename TSpec>
//inline bool openJstFromFile(JournaledStringTree<TSequence, TConfig, TSpec> &, const char *,FileOpenMode const, Vcf);
//
//template <typename TSequence, typename TConfig, typename TSpec>
//inline bool
//openJstFromFile(JournaledStringTree<TSequence, TConfig, TSpec> & /*jst*/,
//                const char * /*fileName*/,
//                FileOpenMode const /*mode*/,
//                TagSelector<> const & /*format*/)
//{
//    return false;
//}
//
//template <typename TSequence, typename TConfig, typename TSpec, typename TTagList>
//inline bool
//openJstFromFile(JournaledStringTree<TSequence, TConfig, TSpec> & jst,
//                const char * fileName,
//                FileOpenMode const mode,
//                TagSelector<TTagList> const & format)
//{
//    typedef typename TTagList::Type TFormat;
//
//    if (isEqual(format, TFormat()))
//        return impl::openJstFromFile(jst, fileName, mode, TFormat());
//    else
//        return openJstFromFile(jst, fileName, mode, static_cast<typename TagSelector<TTagList>::Base const &>(format));
//}

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
 * @param[in] jst The journaled string tree.
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
 * @brief Sets the base sequence of the journaled string tree.
 *
 * @signature void setHost(jst, host);
 *
 * @param[in,out] jst   The journaled string tree.
 * @param[in]     host  The new host to set.
 *
 * When setting a new host, the <tt>jst</tt> will be cleared all previous information are lost.
 */

template <typename TSequence, typename TConfig, typename TSpec, typename THost>
inline void
setHost(JournaledStringTree<TSequence, TConfig, TSpec> & jst,
        THost && host)
{
    clear(jst);
    setHost(jst._source, std::forward<THost>(host));
}

// ----------------------------------------------------------------------------
// Function length()
// ----------------------------------------------------------------------------

/*!
 * @fn JournaledStringTree#length
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Returns the number of sequences represented by the journaled string tree.
 *
 * @signature TSize length(jst);
 *
 * @param[in] jst The journaled string tree to query the length for.
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
// Function resize()
// ----------------------------------------------------------------------------

/*!
 * @fn JournaledStringTree#resize
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Resizes the number of sequences represented by the journaled string tree.
 *
 * @signature void resize(jst, size);
 *
 * @param[in,out] jst   The journaled string tree to resize.
 * @param[in]     size  The new size of the <tt>jst</tt>.
 * 
 * @see JournaledStringTree#length
 *
 * Note, when invoking this method all stored delta events are parsed and their coverage is resized to the given
 * <tt>size</tt>. If <tt>size</tt> is smaller than the previous size this reduces the actual size of the stored coverages.
 * If <tt>size</tt> is bigger than the previous size, the coverages are filled up with 0s, i.e. they represent a copy of
 * the underlying base sequence.
 * Hence the performance is linear in the number of stored delta events times the linear time for resizing the coverages.
 */

template <typename TSequence, typename TConfig, typename TSpec, typename TSize, typename TOverflowSpec>
inline auto
resize(JournaledStringTree<TSequence, TConfig, TSpec> & jst,
       TSize const & newSize,
       Tag<TOverflowSpec> const & /*overflow tag*/) -> decltype(jst._dimension)
{
    auto& deltaMap = impl::member(jst, JstDeltaMapMember());

    for (auto & event : deltaMap)
        resize(getDeltaCoverage(event), newSize, false);

    jst._dimension = newSize;
    return jst._dimension;
}

// ----------------------------------------------------------------------------
// Function maxSize()
// ----------------------------------------------------------------------------

/*
 * @fn JournaledStringTree#maxSize
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Returns the maximal number of delta events that can be stored in the journaled string tree.
 *
 * @signature TSize maxSize(jst);
 *
 * @param[in] jst The journaled string tree to query the maximal length for.
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
 * @brief Returns the number of delta events stored in the journaled string tree.
 *
 * @signature TSize size(jst);
 *
 * @param[in] jst The journaled string tree to query the size for.
 *
 * @return TSize The current number of stored delta events. Of type @link JournaledStringTree#Size @endlink.
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
 * @fn JournaledStringTree#empty
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Checks whether the underlying @link DeltaMap @endlink is empty.
 *
 * @signature bool empty(jst);
 *
 * @param[in] jst The journaled string tree to check for.
 *
 * @return bool <tt>true</tt> if the underlying @link DeltaMap @endlink is empty, otherwise <tt>false</tt>.
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
 * @brief Clears the journaled string tree.
 *
 * @signature void clear(jst);
 *
 * @param[in,out] jst  The journaled string tree to be cleared.
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
 * @fn JournaledStringTree#insert
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Adds a new delta ecent to the journaled string tree.
 *
 * @signature bool insert(jst, pos, val, ids, type);
 *
 * @param[in,out]   jst     The journaled string tree to be modified.
 * @param[in]       pos     The position of the host to place the new delta event at. This position must not exceed length of the @link JournaledStringTree#host @endlink.
 * @param[in]       val     The value of the delta event.
 * @param[in]       ids     An object of @link ContainerConcept @endlink storing the ids of the sequences covering this delta event. The ids must be in the range [0, @link JournaledStringTree#length @endlink).
 * @param[in]       type    A tag representing the type of the delta event. One of @link DeltaTypeTags @endlink.
 *
 * Note the value given by <tt>val<\tt> depends on the <tt>type<\tt>, which can be a single character, a string, 
 * an integer or a pair of integer and string, representing a SNP, an insertion, a deletion or
 * a structural variant respectively.
 *
 * When inserting a deletion of size > 1 (the same holds for the structural variant), the event is split into a left begin and a right begin and two events are inserted into the underlying @link DeltaMap @endlink.
 * Hence the number of inserted events can differ from the actual number of calls to the insert function.
 *
 * @remark The insertion time depends on the implementation of the underlying @link DeltaMap @endlink. See @link DeltaMap#insert @endlink for more details.
 *
 * @see JournaledStringTree#erase
 */

template <typename TSequence, typename TConfig, typename TSpec, typename TPos, typename TValue, typename TIds,
          typename TDeltaType>
inline SEQAN_FUNC_ENABLE_IF(Is<ContainerConcept<TIds> >, void)
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
 * @fn JournaledStringTree#erase
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Erases a delta event in the journaled string tree.
 *
 * @signature bool erase(jst, pos, type);
 *
 * @param[in,out]   jst     The journaled string tree to be modified.
 * @param[in]       pos     The host source position of the delta event.
 * @param[in]       type    A tag representing the type of the delta event. One of @link DeltaTypeTags @endlink.
 *                          Note there can be multiple nodes at the same <tt>pos<\tt> with different <tt>type<\tt>.
 *
 * A delta event is described by it's host position and the type. If a left end of a deletion is deleted its
 * accompanied right end is deleted as well.
 *
 * @remark The deletion time is depends on the implementation of the underlying @link DeltaMap @endlink. 
 * See @link DeltaMap#erase @endlink for detailed information.
 * @see JournaledStringTree#insert
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

//template <typename TSequence, typename TConfig, typename TSpec>
//inline bool
//open(JournaledStringTree<TSequence, TConfig, TSpec> & jst,
//     const char * fileName,
//     FileOpenMode const mode)
//{
//    // Figure out the end of the file name.
//    JstFileFormatSelector selector;
//    guessFormatFromFilename(fileName, selector);
//
//    return impl::openJstFromFile(jst, fileName, mode, selector);
//}

}

#endif // EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_JOURNALED_STRING_TREE_IMPL_H_
