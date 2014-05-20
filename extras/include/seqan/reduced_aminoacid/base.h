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
// Reduced Versions of the 24-letter amino acid alphabet
// ==========================================================================

#ifndef SEQAN_EXTRAS_REDUCED_AMINOACID_BASE_H_
#define SEQAN_EXTRAS_REDUCED_AMINOACID_BASE_H_

namespace seqan {


// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// -----------------------------------------------------------------------
// Tag ReducedAminoAcid
// -----------------------------------------------------------------------

/*
 * @class ReducedAminoAcid
 * @extends SimpleType
 * @brief Reduced versions of the amino acid alphabet.
 *
 * @signature template <typename TReductionSpec>
 * using ReducedAminoAcid = SimpleType<unsigned char, ReducedAminoAcid_<TReductionSpec> >;
 *
 * @tparam TReductionSpec Either @link Murphy10 @endlink or a specialization of
 * @link ClusterReduction @endlink
 *
 * @section Remarks
 *
 * This module is only available in C++11 and when SEQAN_CXX11_STANDARD is
 * defined.
 *
 * @see ClusterReduction
 *
 * @see Murphy10
 *
 * @headerfile seqan/reduced_aminoacid.h
 */



template <typename TRedSpec>
struct ReducedAminoAcid_ {};

template <typename TRedSpec>
using ReducedAminoAcid = SimpleType<unsigned char, ReducedAminoAcid_<TRedSpec> >;


// ============================================================================
// Metafunctions
// ============================================================================

// -----------------------------------------------------------------------
// Metafunction CompareType
// -----------------------------------------------------------------------

template <typename TRedSpec>
struct CompareType<ReducedAminoAcid<TRedSpec>, __uint8>
{
    typedef ReducedAminoAcid<TRedSpec> Type;
};

template <typename TRedSpec>
struct CompareType<ReducedAminoAcid<TRedSpec>, char>
{
    typedef ReducedAminoAcid<TRedSpec> Type;
};

template <typename TRedSpec>
struct CompareType<ReducedAminoAcid<TRedSpec>, AminoAcid>
{
    typedef ReducedAminoAcid<TRedSpec> Type;
};

template <typename TRedSpec>
struct CompareType<ReducedAminoAcid<TRedSpec>, Unicode>
{
    typedef ReducedAminoAcid<TRedSpec> Type;
};

// -----------------------------------------------------------------------
// Translation Tables (implementations see extra files)
// -----------------------------------------------------------------------

template <typename TRedSpec>
struct TranslateTableAsciiToRedAA_
{
    static char const VALUE[256];
};

template <typename TRedSpec>
struct TranslateTableAAToRedAA_
{
    static char const VALUE[24];
};

template <typename TRedSpec>
struct TranslateTableByteToRedAA_
{
    static char const VALUE[256];
};

// needs to be overwritten
template <typename TRedSpec>
struct TranslateTableRedAAToAscii_
{
    static char const VALUE[24];
};

// ============================================================================
// Functions
// ============================================================================


}
#endif // def SEQAN_EXTRAS_REDUCED_AMINOACID_BASE_H_
