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
// Bacardit5sa (Solvent Accessibility) reduction tables
// ==========================================================================

#ifndef SEQAN_REDUCED_AMINOACID_BACARDIT5SA_TABLES_H_
#define SEQAN_REDUCED_AMINOACID_BACARDIT5SA_TABLES_H_

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

template <typename TSpec>
struct TranslateTableRedAAToChar_<Bacardit5sa, TSpec>
{
    static const char VALUE[ValueSize<SimpleType<unsigned char, ReducedAminoAcid_<Bacardit5sa> > >::VALUE];
};

template <typename TSpec>
struct TranslateTableCharToRedAA_<Bacardit5sa, TSpec>
{
    static const char VALUE[256];
};

template <typename TSpec>
struct TranslateTableAAToRedAA_<Bacardit5sa, TSpec>
{
    static const char VALUE[27];
};

template <typename TSpec>
struct TranslateTableByteToRedAA_<Bacardit5sa, TSpec>
{
    static const char VALUE[256];
};

// ---------------------------------- N = 5 ------------------------------

template <typename TVoidSpec>
char const TranslateTableRedAAToChar_<Bacardit5sa, TVoidSpec>::VALUE[5] =
{
    'A', // A C I U X Y
    'B', // B D E S W Z
    'F', // F J L M V *
    'G', // G H N T
    'K', // K O P Q R
};

template <typename TVoidSpec>
char const TranslateTableCharToRedAA_<Bacardit5sa, TVoidSpec>::VALUE[256] =
{
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  1,  0,  1,  1,  2,  3,  3,  0,  2,
     4,  2,  2,  3,  4,  4,  4,  4,  1,  3,  0,  2,  1,  0,  0,
     1,  0,  0,  0,  0,  0,  0,  0,  1,  0,  1,  1,  2,  3,  3,
     0,  2,  4,  2,  2,  3,  4,  4,  4,  4,  1,  3,  0,  2,  1,
     0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
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

template <typename TVoidSpec>
char const TranslateTableAAToRedAA_<Bacardit5sa, TVoidSpec>::VALUE[27] =
{
     0,  1,  0,  1,  1,  2,  3,  3,  0,  2,  4,  2,  2,
     3,  4,  4,  4,  4,  1,  3,  0,  2,  1,  0,  1,  0,  2
};

template <typename TVoidSpec>
char const TranslateTableByteToRedAA_<Bacardit5sa, TVoidSpec>::VALUE[256] =
{
     0,  1,  2,  3,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
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

// ============================================================================
// Functions
// ============================================================================

} // namespace

#endif // SEQAN_REDUCED_AMINOACID_BACARDIT5SA_TABLES_H_
