// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2013, Knut Reinert, FU Berlin
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
// Author: Hannes Hauswedell <hauswedell@mi.fu-berlin.de>
// ==========================================================================
// Murphy10 reduction tables
// ==========================================================================

#ifndef SEQAN_EXTRAS_REDUCED_AMINOACID_MURPHY10_TABLES_H_
#define SEQAN_EXTRAS_REDUCED_AMINOACID_MURPHY10_TABLES_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// ---------------------------------- N = 10 ------------------------------

template <typename TSpec>
struct TranslateTableRedAAToAscii_<Murphy10, TSpec>
{
    typedef typename ValueSize<ReducedAminoAcid<Murphy10>>::Type Type;
    static constexpr Type
    VALUE[ValueSize<ReducedAminoAcid<Murphy10>>::VALUE]
    {
        'A', // A
        'R', // R K
        'N', // N D Q E
        'C', // C
        'G', // G
        'H', // H
        'I', // I L M V
        'F', // F W Y
        'P', // P
        'S'  // S T
    };
};

template <typename TSpec>
constexpr typename ValueSize<ReducedAminoAcid<Murphy10>>::Type
TranslateTableRedAAToAscii_<Murphy10, TSpec>::VALUE
[ValueSize<ReducedAminoAcid<Murphy10>>::VALUE];

// --
template <typename TSpec>
struct TranslateTableAsciiToRedAA_<Murphy10, TSpec>
{
    typedef typename ValueSize<ReducedAminoAcid<Murphy10>>::Type Type;
    static constexpr Type
    VALUE[256]
    {
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        0,  0,  0,  0,  0,  0,  0,  3,  2,  2,  7,  4,  5,  6,  0,
        1,  6,  6,  2,  0,  8,  2,  1,  9,  9,  0,  6,  7,  0,  7,
        0,  0,  0,  0,  0,  0,  0,  0,  0,  3,  2,  2,  7,  4,  5,
        6,  0,  1,  6,  6,  2,  0,  8,  2,  1,  9,  9,  0,  6,  7,
        0,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        0
    };
};

template <typename TSpec>
constexpr typename ValueSize<ReducedAminoAcid<Murphy10>>::Type
TranslateTableAsciiToRedAA_<Murphy10, TSpec>::VALUE[256];

// --
template <typename TSpec>
struct TranslateTableAAToRedAA_<Murphy10, TSpec>
{
    typedef typename ValueSize<ReducedAminoAcid<Murphy10>>::Type Type;
    static constexpr Type
    VALUE[24]
    {
        0,  1,  2,  2,  3,  2,  2,  4,  5,  6,  6,  1,
        6,  7,  8,  9,  9,  7,  7,  6,  0,  0,  0,  0
    };
};

template <typename TSpec>
constexpr typename ValueSize<ReducedAminoAcid<Murphy10>>::Type
TranslateTableAAToRedAA_<Murphy10, TSpec>::VALUE[24];

// --
template <typename TSpec>
struct TranslateTableByteToRedAA_<Murphy10, TSpec>
{
    typedef typename ValueSize<ReducedAminoAcid<Murphy10>>::Type Type;
    static constexpr Type
    VALUE[256]
    {
        0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  0,  0,  0,  0,  0,
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        0
    };
};

template <typename TSpec>
constexpr typename ValueSize<ReducedAminoAcid<Murphy10>>::Type
TranslateTableByteToRedAA_<Murphy10, TSpec>::VALUE[256];

// ============================================================================
// Functions
// ============================================================================


} // namespace

#endif // SEQAN_EXTRAS_REDUCED_AMINOACID_MURPHY10_TABLES_H_