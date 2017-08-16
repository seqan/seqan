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
// Bacardit5sa (Solvent Accessibility) reduction of AminoAcid alphabet
// ==========================================================================

#ifndef SEQAN_REDUCED_AMINOACID_BACARDIT5SA_BASE_H_
#define SEQAN_REDUCED_AMINOACID_BACARDIT5SA_BASE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// -----------------------------------------------------------------------
// Tag Bacardit5sa
// -----------------------------------------------------------------------

/*!
 * @tag Bacardit5sa
 * @brief Specialization for @link ReducedAminoAcid @endlink#
 * @headerfile seqan/reduced_aminoacid.h
 *
 * @signature typedef Bacardit5sa Tag<Bacardit5sa_>;
 *
 * This is the 5-character reduction depending on solvent accessibility defined by Bacardit et al,
 * 2009, <a href="https://www.ncbi.nlm.nih.gov/pubmed/19126227">https://www.ncbi.nlm.nih.gov/pubmed/19126227</a>
 *
 * Since it was created from the 20-letter alphabet the clusters in SeqAn are
 * not identical (they contain more symbols). This is the clustering:
 * @code{.txt}
 *   'A', // A C I U X Y
 *   'B', // B D E S W Z
 *   'F', // F J L M V *
 *   'G', // G H N T
 *   'K', // K O P Q R
 * @endcode
 *
 */

struct Bacardit5sa_ {};

typedef Tag<Bacardit5sa_> Bacardit5sa;

// ============================================================================
// Metafunctions
// ============================================================================

// -----------------------------------------------------------------------
// Metafunction ValueSize
// -----------------------------------------------------------------------

template <>
struct ValueSize<SimpleType<unsigned char, ReducedAminoAcid_<Bacardit5sa> > >
{
    typedef uint8_t Type;
    static const Type VALUE = 5;
};

// -----------------------------------------------------------------------
// Metafunction BitPerValue
// -----------------------------------------------------------------------

template <>
struct BitsPerValue<SimpleType<unsigned char, ReducedAminoAcid_<Bacardit5sa> > >
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
#endif // def SEQAN_REDUCED_AMINOACID_BACARDIT5SA_BASE_H_
