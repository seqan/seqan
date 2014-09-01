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
// Class GffIOContext, accessor functions.
// ==========================================================================

#ifndef CORE_INCLUDE_SEQAN_GFF_IO_GFF_IO_CONTEXT_H_
#define CORE_INCLUDE_SEQAN_GFF_IO_GFF_IO_CONTEXT_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/*!
 * @class GffIOContext
 * @headerfile <seqan/bam_io.h>
 * @brief The I/O context to use for GFF I/O.
 * 
 * @signature template <typename TNameStore[, typename TNameStoreCache]>
 *            class GffIOContext;
 * 
 * @tparam TNameStore      The name store class.
 * @tparam TNameStoreCache The name store cache class, default: NameStoreCache&lt;TNameStore&gt;.
 * 
 * @section Examples
 * 
 * Creating a GffIOContext for a raw StringSet of CharString.
 * 
 * @code{.cpp}
 * StringSet<CharString> nameStore;
 * NameStoreCache<StringSet<CharString> > nameStoreCache(nameStore);
 * GffIOContext<StringSet<CharString> > bamIOContext(nameStore, nameStoreCache);
 * // ...
 * @endcode
 *
 * Using a GffIOContext with a FragmentStore.
 * 
 * @code{.cpp}
 * typedef FragmentStore<>::TContigNameStore         TNameStore;
 * typedef NameStoreCache<TNameStore>                TNameStoreCache;
 * FragmentStore<> store;
 * // Optionally, do something with store.
 * typedef GffIOContext<TNameStore, TNameStoreCache> TGffIOContext;
 * TGffIOContext bamIOContext(store.contigNameStore, store.contigNameStoreCache);
 * // ...
 * @endcode
 *
 *
 * @fn GffIOContext::GffIOContext
 * @brief Constructor.
 * @signature GffIOContext::GffIOContext();
 * 
 * @section Remarks
 * 
 * Only the default constructor is provided.
 * 
 *
 * @typedef GffIOContext::TNameStore
 * @brief The name store class.
 *
 * @signature typedef (...) TNameStore;
 * 
 * @see GffIOContext#nameStore
 * 
 *
 * @typedef GffIOContext::TNameStoreCache
 * @brief The name store cache class.
 *
 * @signature typedef (...) TNameStoreCache;
 * 
 * @see GffIOContext#nameStoreCache
 */

/**
.Class.GffIOContext
..cat:GFF I/O
..signature:GffIOContext<TNameStore[, TNameStoreCache]>
..summary:The I/O context to use for GFF I/O.
..param.TNameStore:The name store class.
..param.TNameStoreCache:The name store cache class.
...default:@Class.NameStoreCache@<TNameStore>
..include:bam_io.h
..example.text:Creating a @Class.GffIOContext@ for a raw @Class.StringSet@ of @Shortcut.CharString@.
..example.code:
StringSet<CharString> nameStore;
NameStoreCache<StringSet<CharString> > nameStoreCache(nameStore);
GffIOContext<StringSet<CharString> > bamIOContext(nameStore, nameStoreCache);
// ...
..example.text:Using a @Class.GffIOContext@ with a @Class.FragmentStore@.
..example.code:
typedef FragmentStore<>::TContigNameStore         TNameStore;
typedef NameStoreCache<TNameStore>                TNameStoreCache;
FragmentStore<> store;
// Optionally, do something with store.
typedef GffIOContext<TNameStore, TNameStoreCache> TGffIOContext;
TGffIOContext bamIOContext(store.contigNameStore, store.contigNameStoreCache);
// ...

.Memfunc.GffIOContext#GffIOContext
..class:Class.GffIOContext
..signature:GffIOContext()
..summary:Constructor.
..remarks:Only the default constructor is provided.

.Typedef.GffIOContext#TNameStore
..class:Class.GffIOContext
..summary:The name store class.

.Typedef.GffIOContext#TNameStoreCache
..class:Class.GffIOContext
..summary:The name store cache class.
*/

template <typename TNameStore_, typename TNameStoreCache_ = NameStoreCache<TNameStore_>, typename TStorageSpec = Owner<> >
class GffIOContext
{
public:
    typedef TNameStore_ TNameStore;
    typedef TNameStoreCache_ TNameStoreCache;

    typedef typename StorageSwitch<TNameStore, TStorageSpec>::Type      TNameStoreMember;
    typedef typename StorageSwitch<TNameStoreCache, TStorageSpec>::Type TNameStoreCacheMember;

    TNameStoreMember        _nameStore;
    TNameStoreCacheMember   _nameStoreCache;
    CharString              buffer;

    GffIOContext() :
        _nameStore(TNameStoreMember()),
        _nameStoreCache(TNameStoreCacheMember())
    {}

    GffIOContext(TNameStore & nameStore, TNameStoreCache & nameStoreCache) :
        _nameStore(_referenceCast<typename Parameter_<TNameStoreMember>::Type>(nameStore)),
        _nameStoreCache(_referenceCast<typename Parameter_<TNameStoreCacheMember>::Type>(nameStoreCache))
    {}

    template <typename TOtherStorageSpec>
    GffIOContext(GffIOContext<TNameStore, TNameStoreCache, TOtherStorageSpec> & other) :
        _nameStore(_referenceCast<typename Parameter_<TNameStoreMember>::Type>(nameStore(other))),
        _nameStoreCache(_referenceCast<typename Parameter_<TNameStoreCacheMember>::Type>(nameStoreCache(other)))
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function nameStore()
// ----------------------------------------------------------------------------

/*!
 * @fn GffIOContext#nameStore
 * @brief Return reference to name store from @link GffIOContext @endlink.
 * 
 * @signature TNameStore nameStore(context);
 * 
 * @param[in] context The GffIOContext to query.
 * 
 * @return TNameStore Reference to the name store of the context (type <tt>TNameStore</tt>).
 *
 * @see GffIOContext::TNameStore
 * @see GffIOContext#nameStoreCache
 */

/**
.Function.GffIOContext#nameStore
..class:Class.GffIOContext
..cat:GFF I/O
..summary:Return reference to name store from @Class.GffIOContext@.
..signature:nameStore(context)
..param.context:The @Class.GffIOContext@ to query.
...type:Class.GffIOContext
..see:Typedef.GffIOContext#TNameStore
..include:seqan/bam_io.h
*/

// TODO(holtgrew): Rename to referenceNameStore
template <typename TNameStore, typename TNameStoreCache>
TNameStore &
nameStore(GffIOContext<TNameStore, TNameStoreCache> & context)
{
    return _referenceCast<TNameStore &>(context._nameStore);
}

template <typename TNameStore, typename TNameStoreCache>
TNameStore const &
nameStore(GffIOContext<TNameStore, TNameStoreCache> const & context)
{
    return _referenceCast<TNameStore const &>(context._nameStore);
}

// ----------------------------------------------------------------------------
// Function nameStoreCache()
// ----------------------------------------------------------------------------

/*!
 * @fn GffIOContext#nameStoreCache
 * @brief Return reference to name store cache from @link GffIOContext @endlink.
 * 
 * @signature TNameStoreCache nameStoreCache(context);
 * 
 * @param[in] context The GffIOContext to query.
 *
 * @return TNameStoreCache A reference to the NameStoreCache of the context.
 * 
 * @see GffIOContext::TNameStoreCache
 * @see GffIOContext#nameStore
 */

/**
.Function.GffIOContext#nameStoreCache
..class:Class.GffIOContext
..cat:GFF I/O
..summary:Return reference to name store cache from @Class.GffIOContext@.
..signature:nameStoreCache(context)
..param.context:The @Class.GffIOContext@ to query.
...type:Class.GffIOContext
..see:Typedef.GffIOContext#TNameStoreCache
..include:seqan/bam_io.h
..see:Function.GffIOContext#nameStore
*/

// TODO(holtgrew): Rename to referenceNameStoreCache
template <typename TNameStore, typename TNameStoreCache>
TNameStoreCache &
nameStoreCache(GffIOContext<TNameStore, TNameStoreCache> & context)
{
    return _referenceCast<TNameStoreCache &>(context._nameStoreCache);
}

template <typename TNameStore, typename TNameStoreCache>
TNameStoreCache const &
nameStoreCache(GffIOContext<TNameStore, TNameStoreCache> const & context)
{
    return _referenceCast<TNameStoreCache const &>(context._nameStoreCache);
}

}  // namespace seqan

#endif  // #ifndef CORE_INCLUDE_SEQAN_GFF_IO_GFF_IO_CONTEXT_H_

