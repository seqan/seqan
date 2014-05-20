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
// Murphy10 reduction of AminoAcid alphabet
// ==========================================================================

#ifndef SEQAN_EXTRAS_REDUCED_AMINOACID_MURPHY10_BASE_H_
#define SEQAN_EXTRAS_REDUCED_AMINOACID_MURPHY10_BASE_H_

namespace seqan {


// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// -----------------------------------------------------------------------
// Tag Murphy10
// -----------------------------------------------------------------------

/*
 * @tag Murphy10
 * @brief Specialization for @link ReducedAminoAcid @endlink
 *
 * @signature typedef Murphy10 Tag<Murphy10_>;
 *
 * @section Remarks
 * @subsection Background
 *
 * This is the 10-character reduction defined by Murphy et al,
 * 2000, <a href="http://www.ncbi.nlm.nih.gov/pubmed/10775656">http://www.ncbi.nlm.nih.gov/pubmed/10775656</a>
 *
 * This alphabet is used by many tools, e.g. Rapsearch2.
 *
 * It looks like this :
 *
 * A, (R, K), (N, D, Q, E), C, G, H, (I, L, M, V), (F, W, Y), P, (S, T)
 *
 * @headerfile seqan/reduced_aminoacid.h
 */


struct Murphy10_ {};

typedef Tag<Murphy10_> Murphy10;

// ============================================================================
// Metafunctions
// ============================================================================

// -----------------------------------------------------------------------
// Metafunction ValueSize
// -----------------------------------------------------------------------

template <>
struct ValueSize<ReducedAminoAcid<Murphy10> >
{
    typedef __uint8 Type;
    static const Type VALUE = 10;
};

// -----------------------------------------------------------------------
// Metafunction BitPerValue
// -----------------------------------------------------------------------

template <>
struct BitsPerValue<ReducedAminoAcid<Murphy10> >
{
    typedef __uint8 Type;
    static const Type VALUE = 4;
};

// -----------------------------------------------------------------------
// Translation Tables (implementations see extra files)
// -----------------------------------------------------------------------

template <>
struct TranslateTableRedAAToAscii_<Murphy10 >
{
    static char const VALUE[10];
};

// ============================================================================
// Functions
// ============================================================================


}
#endif // def SEQAN_EXTRAS_REDUCED_AMINOACID_MURPHY10_BASE_H_
