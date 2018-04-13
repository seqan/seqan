// ==========================================================================
//                      Yara - Yet Another Read Aligner
// ==========================================================================
// Copyright (c) 2011-2018, Enrico Siragusa, FU Berlin
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
//     * Neither the name of Enrico Siragusa or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ENRICO SIRAGUSA OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================

#ifndef APP_YARA_BITS_SEEDS_H_
#define APP_YARA_BITS_SEEDS_H_

using namespace seqan;

// ============================================================================
// Functions
// ============================================================================

// --------------------------------------------------------------------------
// Function getHostPos()
// --------------------------------------------------------------------------

template <typename THost, typename TSpec, typename TPos>
inline typename StringSetPosition<THost>::Type
getHostPos(StringSet<THost, Segment<TSpec> > const & me, TPos pos)
{
    return me.positions[pos];
}

// ----------------------------------------------------------------------------
// Function getSeedIds()
// ----------------------------------------------------------------------------
// getSeqNosByHostSeqNo(seqNo)

template <typename THost, typename TSpec, typename TReadId>
inline Pair<typename Id<StringSet<THost, Segment<TSpec> > const>::Type>
getSeedIds(StringSet<THost, Segment<TSpec> > const & seeds, TReadId readId)
{
    typedef StringSet<THost, Segment<TSpec> > const         TStringSet;
    typedef typename Id<TStringSet>::Type                   TId;
    typedef typename StringSetPositions<THost>::Type const  TPositions;
    typedef typename Iterator<TPositions, Standard>::Type   TPositionsIter;
    typedef typename Value<TPositions>::Type                TPos;

    TPositionsIter posBegin = begin(seeds.positions, Standard());
    TPositionsIter posEnd = end(seeds.positions, Standard());
    TPos pos(readId, 0);

    TPositionsIter seedsBegin = std::lower_bound(posBegin, posEnd, TPos(readId, 0));
    TPositionsIter seedsEnd = std::lower_bound(posBegin, posEnd, TPos(readId + 1, 0));

    return Pair<TId>(position(seedsBegin, seeds.positions), position(seedsEnd, seeds.positions));
}

// ----------------------------------------------------------------------------
// Function getReadSeqId()
// ----------------------------------------------------------------------------

template <typename THost, typename TSpec, typename TSeedId>
inline typename Id<THost>::Type
getReadSeqId(StringSet<THost, Segment<TSpec> > const & seeds, TSeedId seedId)
{
    return getSeqNo(getHostPos(seeds, seedId));
}

// ----------------------------------------------------------------------------
// Function getPosInRead()
// ----------------------------------------------------------------------------

template <typename THost, typename TSpec, typename TSeedId>
inline Pair<typename Position<typename Value<THost>::Type>::Type>
getPosInRead(StringSet<THost, Segment<TSpec> > const & seeds, TSeedId seedId)
{
    typedef typename Position<typename Value<THost>::Type>::Type TPos;

    TPos seedBegin = getSeqOffset(getHostPos(seeds, seedId));
    TPos seedEnd = seedBegin + length(value(seeds, seedId));

    return Pair<TPos>(seedBegin, seedEnd);
}

#endif  // #ifndef APP_YARA_BITS_SEEDS_H_
