// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2017, Knut Reinert, FU Berlin
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
// Author: Kristin Knorr <kristin.knorr@fu-berlin.de>
// ==========================================================================
// Bacardit5cn (Contact Number) reduction of AminoAcid alphabet
// ==========================================================================

#ifndef SEQAN_REDUCED_AMINOACID_BACARDIT5CN_BASE_H_
#define SEQAN_REDUCED_AMINOACID_BACARDIT5CN_BASE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// -----------------------------------------------------------------------
// Tag Bacardit5cn
// -----------------------------------------------------------------------

/*!
 * @tag Bacardit5cn
 * @brief Specialization for @link ReducedAminoAcid @endlink#
 * @headerfile seqan/reduced_aminoacid.h
 *
 * @signature typedef Bacardit5cn Tag<Bacardit5cn_>;
 *
 * This is the 5-character reduction depending on Contact Number defined by Bacardit et al,
 * 2009, <a href="https://www.ncbi.nlm.nih.gov/pubmed/19126227">https://www.ncbi.nlm.nih.gov/pubmed/19126227</a>
 *
 * Since it was created from the 20-letter alphabet the clusters in SeqAn are
 * not identical (they contain more symbols). This is the clustering:
 * @code{.txt}
 *   'A', // A I S X
 *   'B', // B D E P Q Y Z
 *   'C', // C H J L U V
 *   'F', // F G M W *
 *   'K', // K N O R T
 * @endcode
 *
 */

struct Bacardit5cn_ {};

typedef Tag<Bacardit5cn_> Bacardit5cn;

// ============================================================================
// Metafunctions
// ============================================================================

// -----------------------------------------------------------------------
// Metafunction ValueSize
// -----------------------------------------------------------------------

template <>
struct ValueSize<SimpleType<unsigned char, ReducedAminoAcid_<Bacardit5cn> > >
{
    typedef uint8_t Type;
    static const Type VALUE = 5;
};

// -----------------------------------------------------------------------
// Metafunction BitPerValue
// -----------------------------------------------------------------------

template <>
struct BitsPerValue<SimpleType<unsigned char, ReducedAminoAcid_<Bacardit5cn> > >
{
    typedef uint8_t Type;
    static const Type VALUE = 3;
};

// -----------------------------------------------------------------------
// Translation Tables (implementations see extra files)
// -----------------------------------------------------------------------

// ============================================================================
// Functions
// ============================================================================

}
#endif // def SEQAN_REDUCED_AMINOACID_BACARDIT5CN_BASE_H_
