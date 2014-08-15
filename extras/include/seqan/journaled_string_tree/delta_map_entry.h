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
// Implements the delta map entry.
// ==========================================================================

#ifndef EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_DELTA_MAP_ENTRY_H_
#define EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_DELTA_MAP_ENTRY_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

template <typename T>
struct DeltaCoverage;

template <typename T>
struct DeltaPosition;

template <typename T>
struct DeltaRecord;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Tag DeltaMapEntryCoverage
// ----------------------------------------------------------------------------

struct DeltaMapEntryCoverage_;
typedef Tag<DeltaMapEntryCoverage_> DeltaMapEntryCoverage;

// ----------------------------------------------------------------------------
// Tag DeltaMapEntryDeltaStoreValue
// ----------------------------------------------------------------------------

struct DeltaMapEntryDeltaPosition_;
typedef Tag<DeltaMapEntryDeltaPosition_> DeltaMapEntryDeltaPosition;

// ----------------------------------------------------------------------------
// Class DeltaMapEntry
// ----------------------------------------------------------------------------

template <typename TRefPos, typename TStorePos>
class DeltaMapEntry
{
public:

    typedef typename DeltaCoverage<DeltaMapEntry>::Type TCoverage;
    typedef typename DeltaPosition<DeltaMapEntry>::Type TDeltaPos;
    typedef typename DeltaRecord<DeltaMapEntry>::Type TDeltaRecord;

    TDeltaPos    deltaPosition;
    TDeltaRecord deltaRecord;
    TCoverage    deltaCoverage;

    // Default C'tor.
    DeltaMapEntry()
    {}

    // Custom C'tor.
    DeltaMapEntry(TDeltaPos _deltaPos, TDeltaRecord _deltaRecord, TCoverage const & _coverage) :
        deltaPosition(_deltaPos),
        deltaRecord(_deltaRecord),
        deltaCoverage(_coverage)
    {}
};

// ----------------------------------------------------------------------------
// Struct DeltaMapEntryCompareLessByDeltaPosition_
// ----------------------------------------------------------------------------

struct DeltaMapEntryCompareLessByDeltaPosition_
{
    template <typename TRefPos, typename TStorePos>
    inline bool operator()(DeltaMapEntry<TRefPos, TStorePos> const & lhs, DeltaMapEntry<TRefPos, TStorePos> const & rhs)
    {
        return lhs.deltaPosition < rhs.deltaPosition;
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction DeltaPosition
// ----------------------------------------------------------------------------

template <typename TRefPos, typename TStorePos>
struct DeltaPosition<DeltaMapEntry<TRefPos, TStorePos> >
{
    typedef TRefPos Type;
};

template <typename TRefPos, typename TStorePos>
struct DeltaPosition<DeltaMapEntry<TRefPos, TStorePos> const>
{
    typedef TRefPos const Type;
};

// ----------------------------------------------------------------------------
// Metafunction DeltaCoverage
// ----------------------------------------------------------------------------

template <typename TRefPos, typename TStorePos>
struct DeltaCoverage<DeltaMapEntry<TRefPos, TStorePos> >
{
    typedef String<bool, Packed<> > Type;
};

template <typename TRefPos, typename TStorePos>
struct DeltaCoverage<DeltaMapEntry<TRefPos, TStorePos> const>
{
    typedef typename DeltaCoverage<DeltaMapEntry<TRefPos, TStorePos> >::Type const Type;
};

// ----------------------------------------------------------------------------
// Metafunction DeltaRecord
// ----------------------------------------------------------------------------

template <typename TRefPos, typename TStorePos>
struct DeltaRecord<DeltaMapEntry<TRefPos, TStorePos> >
{
    static const unsigned REMAINING_BITS = BitsPerValue<TStorePos>::VALUE - BitsPerValue<DeltaType>::VALUE;
    typedef Pair<DeltaType, TStorePos, BitPacked<BitsPerValue<DeltaType>::VALUE, REMAINING_BITS> > Type;
};

template <typename TRefPos, typename TStorePos>
struct DeltaRecord<DeltaMapEntry<TRefPos, TStorePos> const>
{
    typedef typename DeltaRecord<DeltaMapEntry<TRefPos, TStorePos> >::Type const Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function setRefPosition()
// ----------------------------------------------------------------------------

template <typename TRefPos, typename TStorePos, typename TPosition>
inline void
setDeltaPosition(DeltaMapEntry<TRefPos, TStorePos> & deltaEntry, TPosition newRefPos)
{
    deltaEntry.refPos = newRefPos;
}

// ----------------------------------------------------------------------------
// Function getRefPosition()
// ----------------------------------------------------------------------------

template <typename TRefPos, typename TStorePos>
inline typename DeltaPosition<DeltaMapEntry<TRefPos, TStorePos> >::Type &
getDeltaPosition(DeltaMapEntry<TRefPos, TStorePos> & deltaEntry)
{
    return deltaEntry.deltaPosition;
}

template <typename TRefPos, typename TStorePos>
inline typename DeltaPosition<DeltaMapEntry<TRefPos, TStorePos> const>::Type &
getDeltaPosition(DeltaMapEntry<TRefPos, TStorePos> const & deltaEntry)
{
    return deltaEntry.deltaPosition;
}

// ----------------------------------------------------------------------------
// Function setDeltaCoverage()
// ----------------------------------------------------------------------------

template <typename TRefPos, typename TStorePos, typename TCoverage>
inline void
setDeltaCoverage(DeltaMapEntry<TRefPos, TStorePos> & deltaEntry, TCoverage newDeltaCoverage)
{
    deltaEntry.deltaCoverage = newDeltaCoverage;
}

// ----------------------------------------------------------------------------
// Function getDeltaCoverage()
// ----------------------------------------------------------------------------

template <typename TRefPos, typename TStorePos>
inline typename DeltaCoverage<DeltaMapEntry<TRefPos, TStorePos> >::Type &
getDeltaCoverage(DeltaMapEntry<TRefPos, TStorePos> & deltaEntry)
{
    return deltaEntry.deltaCoverage;
}

template <typename TRefPos, typename TStorePos>
inline typename DeltaCoverage<DeltaMapEntry<TRefPos, TStorePos> const>::Type &
getDeltaCoverage(DeltaMapEntry<TRefPos, TStorePos> const & deltaEntry)
{
    return deltaEntry.deltaCoverage;
}

// ----------------------------------------------------------------------------
// Function setDeltaRecord()
// ----------------------------------------------------------------------------

template <typename TRefPos, typename TStorePos, typename TRecord>
inline void
setDeltaRecord(DeltaMapEntry<TRefPos, TStorePos> & deltaEntry, TRecord newDeltaRecord)
{
    deltaEntry.deltaRecord = newDeltaRecord;
}

// ----------------------------------------------------------------------------
// Function getDeltaCoverage()
// ----------------------------------------------------------------------------

template <typename TRefPos, typename TStorePos>
inline typename DeltaRecord<DeltaMapEntry<TRefPos, TStorePos> >::Type &
getDeltaRecord(DeltaMapEntry<TRefPos, TStorePos> & deltaEntry)
{
    return deltaEntry.deltaRecord;
}

template <typename TRefPos, typename TStorePos>
inline typename DeltaRecord<DeltaMapEntry<TRefPos, TStorePos> const>::Type &
getDeltaRecord(DeltaMapEntry<TRefPos, TStorePos> const & deltaEntry)
{
    return deltaEntry.deltaRecord;
}

// ----------------------------------------------------------------------------
// Function operator==()
// ----------------------------------------------------------------------------

template <typename TRefPos, typename TStorePos>
inline bool operator==(DeltaMapEntry<TRefPos, TStorePos> const & lhs, DeltaMapEntry<TRefPos, TStorePos> const & rhs)
{
    return (lhs.deltaPosition == rhs.deltaPosition) && (lhs.deltaRecord == rhs.deltaRecord) &&
           (lhs.deltaCoverage == rhs.deltaCoverage);
}

// ----------------------------------------------------------------------------
// Function operator!=()
// ----------------------------------------------------------------------------

template <typename TRefPos, typename TStorePos>
inline bool operator!=(DeltaMapEntry<TRefPos, TStorePos> const & lhs, DeltaMapEntry<TRefPos, TStorePos> const & rhs)
{
    return (lhs.deltaPosition != rhs.deltaPosition) && (lhs.deltaRecord != rhs.deltaRecord) &&
           (lhs.deltaCoverage != rhs.deltaCoverage);
}

// ----------------------------------------------------------------------------
// Function operator<<()
// ----------------------------------------------------------------------------

template <typename TStream, typename TRefPos, typename TStorePos>
inline TStream & operator<<(TStream & stream, DeltaMapEntry<TRefPos, TStorePos> const & entry)
{
    stream << "<" << entry.deltaPosition << ", ";
    switch (entry.deltaRecord.i1)
    {
        case DELTA_TYPE_SNP: stream << "SNP ("; break;
        case DELTA_TYPE_DEL: stream << "DEL ("; break;
        case DELTA_TYPE_INS: stream << "INS ("; break;
        case DELTA_TYPE_SV: stream << "SV ("; break;
        default: stream << "NA ("; break;
    }
    stream << entry.deltaRecord.i2 << "), " << entry.deltaCoverage << ">\n";
    return stream;
}

}

#endif // EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_DELTA_MAP_ENTRY_H_
