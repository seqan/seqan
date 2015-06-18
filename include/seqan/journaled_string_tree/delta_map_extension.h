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

#ifndef INCLUDE_SEQAN_JOURNALED_STRING_TREE_DELTA_MAP_EXTENSION_H_
#define INCLUDE_SEQAN_JOURNALED_STRING_TREE_DELTA_MAP_EXTENSION_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

enum class ExtensionInfo : __uint8
{
    IS_BEGIN,
    IS_END
};

struct ExtensionMapMember_;
typedef Tag<ExtensionMapMember_> ExtensionMapMember;

struct DeltaMapExtensionIterSpec_;
typedef Tag<DeltaMapExtensionIterSpec_> DeltaMapExtensionIterSpec;

template <typename TDeltaMap>
struct ExtendedDeltaEntry_
{
    typedef typename RemoveConst<TDeltaMap>::Type const     TConstMap;
    typedef typename Iterator<TConstMap, Standard>::Type    TMapIter;
    typedef typename Position<TDeltaMap>::Type              TPosition;

    TMapIter        hostIter;
    TPosition       deltaPos;
    ExtensionInfo   info;

    ExtendedDeltaEntry_() = default;

    template <typename TIter, typename TPos>
    ExtendedDeltaEntry_(TIter const & it, TPos const pos, ExtensionInfo inf) :
        hostIter(it),
        deltaPos(pos),
        info(inf)
    {}

    inline bool operator==(ExtendedDeltaEntry_ const & rhs)
    {
        return (hostIter == rhs.hostIter) && (deltaPos == rhs.deltaPos) && (info == rhs.info);
    }

    inline bool operator!=(ExtendedDeltaEntry_ const & rhs)
    {
        return !(*this == rhs);
    }
};

template <typename TDeltaMap>
class DeltaMapExtension_
{
public:
    typedef typename RemoveConst<TDeltaMap>::Type       TNonConstDeltaMap;
    typedef typename Value<DeltaMapExtension_>::Type    TExtendedEntry;
    typedef String<TExtendedEntry>                      TExtensionTable;

    TNonConstDeltaMap const* _mapPtr;       // Not modifiable from within this class.
    TExtensionTable          _extTable;

    DeltaMapExtension_() : _mapPtr(nullptr)
    {}

    DeltaMapExtension_(TDeltaMap const & host) : _mapPtr(&host)
    {}
};

struct DeltaExtensionCompareLessPos_
{

    template <typename TDeltaMap, typename TPosition>
    inline bool
    operator()(ExtendedDeltaEntry_<TDeltaMap> const & lhs, TPosition const & pos)
    {
        return lhs.deltaPos < pos;
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

template <typename TObject>
struct HostIter_;

template <typename TDeltaMap>
struct HostIter_<DeltaMapExtension_<TDeltaMap> >
{
    typedef typename Value<DeltaMapExtension_<TDeltaMap> >::Type    TValue_;
    typedef typename TValue_::TMapIter                              Type;
};

template <typename TDeltaMap>
struct HostIter_<DeltaMapExtension_<TDeltaMap> const>
{
    typedef typename Value<DeltaMapExtension_<TDeltaMap> >::Type    TValue_;
    typedef typename TValue_::TMapIter const                        Type;
};

template <typename TDeltaMap>
struct Host<DeltaMapExtension_<TDeltaMap> >
{
    typedef TDeltaMap const Type;
};

template <typename TDeltaMap>
struct Host<DeltaMapExtension_<TDeltaMap> const>
    : Host<DeltaMapExtension_<TDeltaMap> >{};

template <typename TDeltaMap>
struct Member<DeltaMapExtension_<TDeltaMap>, ExtensionMapMember>
{
    typedef typename Value<DeltaMapExtension_<TDeltaMap> >::Type    TValue_;
    typedef String<TValue_>                                         Type;
};

template <typename TDeltaMap>
struct Value<DeltaMapExtension_<TDeltaMap> >
{
    typedef typename Host<DeltaMapExtension_<TDeltaMap> >::Type THost_;
    typedef ExtendedDeltaEntry_<THost_>                         Type;
};

template <typename  TDeltaMap>
struct Iterator<DeltaMapExtension_<TDeltaMap>, Standard>
{
    typedef Iter<DeltaMapExtension_<TDeltaMap>, DeltaMapExtensionIterSpec> Type;
};

template <typename  TDeltaMap>
struct Iterator<DeltaMapExtension_<TDeltaMap> const, Standard>
{
    typedef Iter<DeltaMapExtension_<TDeltaMap> const, DeltaMapExtensionIterSpec> Type;
};

template <typename  TDeltaMap>
struct Size<DeltaMapExtension_<TDeltaMap> >
{
    typedef typename Size<TDeltaMap>::Type Type;
};

template <typename TDeltaMap>
struct Position<DeltaMapExtension_<TDeltaMap> >
{
    typedef typename Position<TDeltaMap>::Type Type;
};

// ============================================================================
// Private Functions
// ============================================================================

namespace impl
{

template <typename TDeltaMap>
inline typename Size<TDeltaMap>::Type
getDeletionsCount(DeltaMapExtension_<TDeltaMap> & extMap)
{
    SEQAN_ASSERT_NOT(empty(host(extMap)));
    return length(getDeltaStore(host(extMap)._deltaStore, DeltaTypeDel())) +
           length(getDeltaStore(host(extMap)._deltaStore, DeltaTypeSV()));
}

template <typename TDeltaMap, typename TMapIter, typename TTag>
inline void
recordDeletion(DeltaMapExtension_<TDeltaMap> & extMap, TMapIter const & mapIt, TTag const & /*tag*/)
{
    typedef typename Value<DeltaMapExtension_<TDeltaMap> >::Type TEntry;

    TEntry tmp;
    tmp.hostIter = mapIt;
    tmp.deltaPos = getDeltaPosition(*mapIt) + deletionSize(host(extMap)._deltaStore, getStorePosition(*mapIt), TTag());
    tmp.info = ExtensionInfo::IS_END;
    appendValue(extMap._extTable, SEQAN_MOVE(tmp));
}

template <typename TMapExtension>
struct DeltaEndPositionRecorder
{
    typedef typename Host<TMapExtension>::Type              TDeltaMap;
    typedef typename Iterator<TDeltaMap, Standard>::Type    TMapIter;

    TMapExtension & _ext;
    TMapIter        mapIt;

    DeltaEndPositionRecorder(TMapExtension & ext) : _ext(ext)
    {}

    template <typename TTag>
    inline void
    operator()(TTag const & /*deltaType*/)
    {
        // no-op
    }

    inline void
    operator()(DeltaTypeDel const & /*deltaType*/)
    {
        impl::recordDeletion(_ext, mapIt, DeltaTypeDel());
    }

    inline void
    operator()(DeltaTypeSV const & /*deltaType*/)
    {
        impl::recordDeletion(_ext, mapIt, DeltaTypeSV());
    }
};

struct DeltaExtensionEntryLessThanComparator_
{

    template <typename TDeltaMap>
    inline bool
    operator()(ExtendedDeltaEntry_<TDeltaMap> const & lhs, ExtendedDeltaEntry_<TDeltaMap> const & rhs)
    {
        return lhs.deltaPos < rhs.deltaPos;
    }
};
}  // namespace impl

// ============================================================================
// Functions
// ============================================================================

template <typename TDeltaMap>
inline void
sync(DeltaMapExtension_<TDeltaMap> & ext)
{
    SEQAN_ASSERT(ext._mapPtr != nullptr);
    // First clear old values.
    clear(ext._extTable);

    if (empty(host(ext)))
        return;
    // Get the number stored deletions and SVs
    reserve(ext._extTable, impl::getDeletionsCount(ext), Exact());

    impl::DeltaEndPositionRecorder<DeltaMapExtension_<TDeltaMap> > f(ext);
    f.mapIt = begin(host(ext), Standard());
    auto itEnd = end(host(ext), Standard());
    for (; f.mapIt != itEnd; ++f.mapIt)
    {
        DeltaTypeSelector selector;
        applyOnDelta(f, getDeltaType(*f.mapIt), selector);
    }
    std::sort(begin(ext._extTable, Standard()), end(ext._extTable, Standard()),
              impl::DeltaExtensionEntryLessThanComparator_());
    shrinkToFit(ext._extTable);
}

template <typename TDeltaMap>
inline void
setHost(DeltaMapExtension_<TDeltaMap> & ext,
        TDeltaMap const & host)
{
    clear(ext);
    ext._mapPtr = &host;
    sync(ext);
}

// Returns reference to const delta map.
template <typename TDeltaMap>
inline typename Host<DeltaMapExtension_<TDeltaMap> >::Type &
host(DeltaMapExtension_<TDeltaMap> & ext)
{
    return *ext._mapPtr;
}

template <typename TDeltaMap>
inline typename Host<DeltaMapExtension_<TDeltaMap> const>::Type &
host(DeltaMapExtension_<TDeltaMap> const & ext)
{
    return *ext._mapPtr;
}

template <typename TDeltaMap>
inline void
clear(DeltaMapExtension_<TDeltaMap> & ext)
{
    ext._mapPtr = nullptr;
    clear(ext._extTable);
}

// ----------------------------------------------------------------------------
// Function begin()
// ----------------------------------------------------------------------------

template <typename TDeltaMap>
inline typename Iterator<DeltaMapExtension_<TDeltaMap>, Standard>::Type
begin(DeltaMapExtension_<TDeltaMap> & ext, Standard const & /*standard*/)
{
    return typename Iterator<DeltaMapExtension_<TDeltaMap>, Standard>::Type(ext);
}

template <typename TDeltaMap>
inline typename Iterator<DeltaMapExtension_<TDeltaMap> const, Standard>::Type
begin(DeltaMapExtension_<TDeltaMap> const & ext, Standard const & /*standard*/)
{
    return typename Iterator<DeltaMapExtension_<TDeltaMap> const, Standard>::Type(ext);
}

// ----------------------------------------------------------------------------
// Function end()
// ----------------------------------------------------------------------------

template <typename TDeltaMap>
inline typename Iterator<DeltaMapExtension_<TDeltaMap>, Standard>::Type
end(DeltaMapExtension_<TDeltaMap> & ext, Standard const & /*standard*/)
{
    typename Iterator<DeltaMapExtension_<TDeltaMap>, Standard>::Type tmp(ext);
    tmp._hostMapIter = end(host(ext), Standard());
    tmp._extTableIter = end(ext._extTable, Standard());
    return tmp;
}

template <typename TDeltaMap>
inline typename Iterator<DeltaMapExtension_<TDeltaMap> const, Standard>::Type
end(DeltaMapExtension_<TDeltaMap> const & ext, Standard const & /*standard*/)
{
    typename Iterator<DeltaMapExtension_<TDeltaMap> const, Standard>::Type tmp(ext);
    tmp._hostIter = end(host(ext), Standard());
    tmp._endPointIter = end(ext._extTable, Standard());
    return tmp;
}

// ----------------------------------------------------------------------------
// Function operator<<()
// ----------------------------------------------------------------------------

template <typename TStream, typename TDeltaMap>
inline TStream &
operator<<(TStream & stream, ExtendedDeltaEntry_<TDeltaMap> const & entry)
{
    if (entry.info == ExtensionInfo::IS_BEGIN)
        stream << "INF: Begin ";
    else
        stream << "INF: End ";
    stream << "POS: " << entry.deltaPos << " NODE-ID: " << position(entry.hostIter);
    return stream;
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_JOURNALED_STRING_TREE_DELTA_MAP_EXTENSION_H_
