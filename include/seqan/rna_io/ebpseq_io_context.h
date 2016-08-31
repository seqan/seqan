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
// Author: Gianvito Urgese <gianvito.urgese@polito.it>
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_EBSEQ_IO_EBSEQ_IO_CONTEXT_H_
#define SEQAN_INCLUDE_SEQAN_EBSEQ_IO_EBSEQ_IO_CONTEXT_H_

namespace seqan {

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class EbseqIOContext
// ----------------------------------------------------------------------------

/*!
 * @class EbseqIOContext
 * @headerfile <seqan/Ebseq_io.h>
 * @brief The I/O context to use for Ebseq I/O.
 *
 * @signature template <typename TNameStore[, typename TNameStoreCache]>
 *            class EbseqIOContext;
 *
 * @tparam TNameStore      The type used to represent the names.
 * @tparam TNameStoreCache The type used to cache the names. Defaults to @link NameStoreCache @endlink &lt;TNameStore&gtl;.
 */

/*!
 * @fn EbseqIOContext::EbseqIOContext
 * @brief Constructor.
 *
 * @signature EbseqIOContext::EbseqIOContext();
 * @signature EbseqIOContext::EbseqIOContext(sequenceNames, fixedStructureNames);
 *
 * Default constructor or construction with references to sequence and fixedStructure names.
 */

template <typename TNameStore_        = StringSet<CharString>,
          typename TNameStoreCache_   = NameStoreCache<TNameStore_>,
          typename TStorageSpec       = Owner<> >
class EbseqIOContext
{
public:
    typedef TNameStore_ TNameStore;
    typedef TNameStoreCache_ TNameStoreCache;

    typedef typename StorageSwitch<TNameStore, TStorageSpec>::Type      TNameStoreMember;
    typedef typename StorageSwitch<TNameStoreCache, TStorageSpec>::Type TNameStoreCacheMember;

    // Cache for the sequence name lookup.
    TNameStoreMember        _sequenceNames;
    TNameStoreCacheMember   _sequenceNamesCache;

    // Cache for the fixedStructure name lookup.
    TNameStoreMember        _fixedStructureNames;
    TNameStoreCacheMember   _fixedStructureNamesCache;

    // Cache for the baseProbability name lookup.
    TNameStoreMember        _fixedStructureNames;
    TNameStoreCacheMember   _fixedStructureNamesCache;

    CharString              buffer;

    EbseqIOContext() :
        _sequenceNames(TNameStoreMember()),
        _sequenceNamesCache(ifSwitch(typename IsPointer<TNameStoreCacheMember>::Type(),
                                 (TNameStoreCache*)NULL,
                                 _sequenceNames)),
        _fixedStructureNames(TNameStoreMember()),
        _fixedStructureNamesCache(ifSwitch(typename IsPointer<TNameStoreCacheMember>::Type(),
                                 (TNameStoreCache*)NULL,
                                 _fixedStructureNames)),
        _baseProbabilityNames(TNameStoreMember()),
        _baseProbabilityNamesCache(ifSwitch(typename IsPointer<TNameStoreCacheMember>::Type(),
                                 (TNameStoreCache*)NULL,
                                 _fixedStructureNames))
    {}

    EbseqIOContext(TNameStore & nameStore_, TNameStoreCache & nameStoreCache_) :
        _sequenceNames(_referenceCast<typename Parameter_<TNameStoreMember>::Type>(nameStore_)),
        _sequenceNamesCache(ifSwitch(typename IsPointer<TNameStoreCacheMember>::Type(),
                                 &nameStoreCache_,
                                 _sequenceNames)),
        _fixedStructureNames(TNameStoreMember()),
        _fixedStructureNamesCache(ifSwitch(typename IsPointer<TNameStoreCacheMember>::Type(),
                                 (TNameStoreCache*)NULL,
                                 _fixedStructureNames)),
        _baseProbabilityNames(TNameStoreMember()),
        _baseProbabilityNamesCache(ifSwitch(typename IsPointer<TNameStoreCacheMember>::Type(),
                                 (TNameStoreCache*)NULL,
                                 _fixedStructureNames))
    {}

    template <typename TOtherStorageSpec>
    EbseqIOContext(EbseqIOContext<TNameStore, TNameStoreCache, TOtherStorageSpec> & other) :
        _sequenceNames(_referenceCast<typename Parameter_<TNameStoreMember>::Type>(sequenceNames(other))),
        _sequenceNamesCache(ifSwitch(typename IsPointer<TNameStoreCacheMember>::Type(),
                                 &sequenceNamesCache(other),
                                 _sequenceNames)),
        _fixedStructureNames(_referenceCast<typename Parameter_<TNameStoreMember>::Type>(fixedStructureNames(other))),
        _fixedStructureNamesCache(ifSwitch(typename IsPointer<TNameStoreCacheMember>::Type(),
                                 &fixedStructureNamesCache(other),
                                 _fixedStructureNames)),
        _baseProbabilityNames(_referenceCast<typename Parameter_<TNameStoreMember>::Type>(baseProbabilityNames(other))),
        _baseProbabilityNamesCache(ifSwitch(typename IsPointer<TNameStoreCacheMember>::Type(),
                                 &baseProbabilityNamesCache(other),
                                 _baseProbabilityNames))
    {}
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function sequenceNames()
// ----------------------------------------------------------------------------

/*!
 * @fn EbseqIOContext#sequenceNames
 * @brief Return reference to the sequence names from @link EbseqIOContext @endlink.
 *
 * @signature TNameStoreRef sequenceNames(context);
 *
 * @param[in] context The @link EbseqIOContext @endlink to query.
 *
 * @return TNameStoreRef A reference to the <tt>TNameStore</tt> of the context.
 */

template <typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
inline TNameStore &
sequenceNames(EbseqIOContext<TNameStore, TNameStoreCache, TStorageSpec> & context)
{
    return _referenceCast<TNameStore &>(context._sequenceNames);
}

template <typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
inline TNameStore const &
sequenceNames(EbseqIOContext<TNameStore, TNameStoreCache, TStorageSpec> const & context)
{
    return _referenceCast<TNameStore &>(context._sequenceNames);
}

/*!
 * @fn EbseqIOContext#sequenceNamesCache
 * @brief Return reference to sequence names cache from @link EbseqIOContext @endlink.
 *
 * @signature TNameStoreCacheRef sequenceNamesCache(context);
 *
 * @param[in] context The @link BamIOContext @endlink to query.
 *
 * @return TNameStoreCacheRef A reference to the <tt>TNameStoreCache</tt> of the context.
 */

template <typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
inline TNameStoreCache &
sequenceNamesCache(EbseqIOContext<TNameStore, TNameStoreCache, TStorageSpec> & context)
{
    return _referenceCast<TNameStoreCache &>(context._sequenceNamesCache);
}

template <typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
inline TNameStoreCache const &
sequenceNamesCache(EbseqIOContext<TNameStore, TNameStoreCache, TStorageSpec> const & context)
{
    return _referenceCast<TNameStoreCache &>(context._sequenceNamesCache);
}

/*!
 * @fn EbseqIOContext#fixedStructureNames
 * @brief Return reference to the fixedStructure names from @link EbseqIOContext @endlink.
 *
 * @signature TNameStoreRef fixedStructureNames(context);
 *
 * @param[in] context The @link EbseqIOContext @endlink to query.
 *
 * @return TNameStoreRef A reference to the <tt>TNameStore</tt> of the context.
 */

template <typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
inline TNameStore &
fixedStructureNames(EbseqIOContext<TNameStore, TNameStoreCache, TStorageSpec> & context)
{
    return _referenceCast<TNameStore &>(context._fixedStructureNames);
}

template <typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
inline TNameStore const &
fixedStructureNames(EbseqIOContext<TNameStore, TNameStoreCache, TStorageSpec> const & context)
{
    return _referenceCast<TNameStore &>(context._fixedStructureNames);
}

template <typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
inline TNameStoreCache &
fixedStructureNamesCache(EbseqIOContext<TNameStore, TNameStoreCache, TStorageSpec> & context)
{
    return _referenceCast<TNameStoreCache &>(context._fixedStructureNamesCache);
}

template <typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
inline TNameStoreCache const &
fixedStructureNamesCache(EbseqIOContext<TNameStore, TNameStoreCache, TStorageSpec> const & context)
{
    return _referenceCast<TNameStoreCache &>(context._fixedStructureNamesCache);
}

/*!
 * @fn EbseqIOContext#baseProbabilityNames
 * @brief Return reference to the baseProbability names from @link EbseqIOContext @endlink.
 *
 * @signature TNameStoreRef baseProbabilityNames(context);
 *
 * @param[in] context The @link EbseqIOContext @endlink to query.
 *
 * @return TNameStoreRef A reference to the <tt>TNameStore</tt> of the context.
 */

template <typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
inline TNameStore &
baseProbabilityNames(EbseqIOContext<TNameStore, TNameStoreCache, TStorageSpec> & context)
{
    return _referenceCast<TNameStore &>(context._baseProbabilityNames);
}

template <typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
inline TNameStore const &
baseProbabilityNames(EbseqIOContext<TNameStore, TNameStoreCache, TStorageSpec> const & context)
{
    return _referenceCast<TNameStore &>(context._baseProbabilityNames);
}

template <typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
inline TNameStoreCache &
baseProbabilityNamesCache(EbseqIOContext<TNameStore, TNameStoreCache, TStorageSpec> & context)
{
    return _referenceCast<TNameStoreCache &>(context._baseProbabilityNamesCache);
}

template <typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
inline TNameStoreCache const &
baseProbabilityNamesCache(EbseqIOContext<TNameStore, TNameStoreCache, TStorageSpec> const & context)
{
    return _referenceCast<TNameStoreCache &>(context._baseProbabilityNamesCache);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_EBSEQ_IO_EBSEQ_IO_CONTEXT_H_