// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_EXTRAS_INCLUDE_SEQAN_VCF_IO_VCF_IO_CONTEXT_H_
#define SEQAN_EXTRAS_INCLUDE_SEQAN_VCF_IO_VCF_IO_CONTEXT_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class VcfHeader
// ----------------------------------------------------------------------------

/*!
 * @class VcfIOContext
 * @headerfile <seqan/vcf_io.h>
 * @brief The I/O context to use for VCF I/O.
 * 
 * @signature class VcfIOContext;
 * 
 * VcfIOContext objects store the names of (and provide a cache for) reference and sample names.  StringSet of
 * CharString are used for the name stores.
 * 
 * @var TNamesPtr VcfIOContext::sequenceNames;
 * @brief Names of the reference sequences, pointer to @link StringSet @endlink of @link CharString @endlink.
 * 
 * @var TNameStoreCache VcfIOContext::sequenceNamesCache;
 * @brief Name store cache for of the reference names, @link NameStoreCache @endlink of @link CharString @endlink.
 * 
 * @var TNamesPtr VcfIOContext::sampleNames;
 * @brief Names of the samples, pointer to @link StringSet @endlink of @link CharString @endlink.
 * 
 * @var TNameStoreCache VcfIOContext::sampleNamesCache;
 * @brief Name store cache for the sample names, @link NameStoreCache @endlink of @link CharString @endlink.
 */

/*!
 * @fn VcfIOContext::VcfIOContext
 * @brief Constructor.
 * 
 * @signature VcfIOContext::VcfIOContext();
 * @signature VcfIOContext::VcfIOContext(sequenceNames, sampleNames);
 *
 * Default constructor or construction with references to sequence and sample names.
 */

/**
.Class.VcfIOContext
..cat:VCF I/O
..signature:class VcfIOContext
..summary:The I/O context to use for VCF I/O.
..description:
VcfIOContext objects store the names of (and provide a cache for) reference and sample names.
@Class.StringSet@ of @Shortcut.CharString@ are used for the name stores.
..include:seqan/vcf_io.h

.Memfunc.VcfIOContext#VcfIOContext
..class:Class.VcfIOContext
..signature:VcfIOContext()
..signature:VcfIOContext(sequenceNames, sampleNames)
..summary:Constructor.
..remarks:Default constructor or construction with references to sequence and sample names.

.Memvar.VcfIOContext#sequenceNames
..class:Class.VcfIOContext
..summary:Names of the reference sequences.

.Memvar.VcfIOContext#sequenceNamesCache
..class:Class.VcfIOContext
..summary:Name store cache for of the reference names.

.Memvar.VcfIOContext#sampleNames
..class:Class.VcfIOContext
..summary:Names of the samples.

.Memvar.VcfIOContext#sampleNamesCache
..class:Class.VcfIOContext
..summary:Name store cache for the sample names.
*/

template <
    typename TNameStore_ = StringSet<CharString>,
    typename TNameStoreCache_ = NameStoreCache<TNameStore_>,
    typename TStorageSpec = Owner<> >
class VcfIOContext
{
public:
    typedef TNameStore_ TNameStore;
    typedef TNameStoreCache_ TNameStoreCache;

    typedef typename StorageSwitch<TNameStore, TStorageSpec>::Type      TNameStoreMember;
    typedef typename StorageSwitch<TNameStoreCache, TStorageSpec>::Type TNameStoreCacheMember;

    // Cache for the sequence name lookup.
    TNameStoreMember        _contigNames;
    TNameStoreCacheMember   _contigNamesCache;

    // Cache for the sample name lookup.
    TNameStoreMember        _sampleNames;
    TNameStoreCacheMember   _sampleNamesCache;

    CharString              buffer;

    VcfIOContext() :
        _contigNames(TNameStoreMember()),
        _contigNamesCache(ifSwitch(typename IsPointer<TNameStoreCacheMember>::Type(),
                                 (TNameStoreCache*)NULL,
                                 _contigNames)),
        _sampleNames(TNameStoreMember()),
        _sampleNamesCache(ifSwitch(typename IsPointer<TNameStoreCacheMember>::Type(),
                                 (TNameStoreCache*)NULL,
                                 _sampleNames))
    {}

    VcfIOContext(TNameStore & nameStore_, TNameStoreCache & nameStoreCache_) :
        _contigNames(_referenceCast<typename Parameter_<TNameStoreMember>::Type>(nameStore_)),
        _contigNamesCache(ifSwitch(typename IsPointer<TNameStoreCacheMember>::Type(),
                                 &nameStoreCache_,
                                 _contigNames)),
        _sampleNames(TNameStoreMember()),
        _sampleNamesCache(ifSwitch(typename IsPointer<TNameStoreCacheMember>::Type(),
                                 (TNameStoreCache*)NULL,
                                 _sampleNames))
    {}

    template <typename TOtherStorageSpec>
    VcfIOContext(VcfIOContext<TNameStore, TNameStoreCache, TOtherStorageSpec> & other) :
        _contigNames(_referenceCast<typename Parameter_<TNameStoreMember>::Type>(contigNames(other))),
        _contigNamesCache(ifSwitch(typename IsPointer<TNameStoreCacheMember>::Type(),
                                 &contigNamesCache(other),
                                 _contigNames)),
        _sampleNames(_referenceCast<typename Parameter_<TNameStoreMember>::Type>(sampleNames(other))),
        _sampleNamesCache(ifSwitch(typename IsPointer<TNameStoreCacheMember>::Type(),
                                 &sampleNamesCache(other),
                                 _sampleNames))
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

template <typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
inline TNameStore &
contigNames(VcfIOContext<TNameStore, TNameStoreCache, TStorageSpec> & context)
{
    return _referenceCast<TNameStore &>(context._contigNames);
}

template <typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
inline TNameStoreCache &
contigNamesCache(VcfIOContext<TNameStore, TNameStoreCache, TStorageSpec> & context)
{
    return _referenceCast<TNameStoreCache &>(context._contigNamesCache);
}

template <typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
inline TNameStore &
sampleNames(VcfIOContext<TNameStore, TNameStoreCache, TStorageSpec> & context)
{
    return _referenceCast<TNameStore &>(context._sampleNames);
}

template <typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
inline TNameStoreCache &
sampleNamesCache(VcfIOContext<TNameStore, TNameStoreCache, TStorageSpec> & context)
{
    return _referenceCast<TNameStoreCache &>(context._sampleNamesCache);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_EXTRAS_INCLUDE_SEQAN_VCF_IO_VCF_IO_CONTEXT_H_
