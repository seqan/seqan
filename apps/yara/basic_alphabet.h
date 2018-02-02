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
// This file contains Dna5 specializations to deal with uncalled bases (N).
// ==========================================================================

#ifndef APP_YARA_BASIC_ALPHABET_H_
#define APP_YARA_BASIC_ALPHABET_H_

namespace seqan {

// ============================================================================
// Operators
// ============================================================================

static unsigned char __MASK_DNA5_EQ[]  = {1, 2, 4, 8, 0};
static unsigned char __MASK_DNA5_LT[]  = {0, 1, 2, 3, 4};
static unsigned char __MASK_DNA5Q_LT[] = {0, 1, 2, 3, 5};

// ----------------------------------------------------------------------------
// Operators ==, !=, <, >, <=, >=                [Dna5 vs Dna5Q, Dna5Q vs Dna5]
// ----------------------------------------------------------------------------

template <>
inline bool operator==(Dna5 const & left_, Dna5Q const & right_)
{
    return __MASK_DNA5_EQ[ordValue(left_)] & __MASK_DNA5_EQ[ordValue(right_)];
}

template <>
inline bool operator==(Dna5Q const & left_, Dna5 const & right_)
{
    return __MASK_DNA5_EQ[ordValue(left_)] & __MASK_DNA5_EQ[ordValue(right_)];
}

template <>
inline bool operator!=(Dna5 const & left_, Dna5Q const & right_)
{
    return !(__MASK_DNA5_EQ[ordValue(left_)] & __MASK_DNA5_EQ[ordValue(right_)]);
}

template <>
inline bool operator!=(Dna5Q const & left_, Dna5 const & right_)
{
    return !(__MASK_DNA5_EQ[ordValue(left_)] & __MASK_DNA5_EQ[ordValue(right_)]);
}

template <>
inline bool operator<(Dna5 const & left_, Dna5Q const & right_)
{
    return __MASK_DNA5_LT[ordValue(left_)] < __MASK_DNA5Q_LT[ordValue(right_)];
}

template <>
inline bool operator<(Dna5Q const & left_, Dna5 const & right_)
{
    return __MASK_DNA5Q_LT[ordValue(left_)] < __MASK_DNA5_LT[ordValue(right_)];
}

template <>
inline bool operator>(Dna5 const & left_, Dna5Q const & right_)
{
    return __MASK_DNA5_LT[ordValue(left_)] > __MASK_DNA5Q_LT[ordValue(right_)];
}

template <>
inline bool operator>(Dna5Q const & left_, Dna5 const & right_)
{
    return __MASK_DNA5Q_LT[ordValue(left_)] > __MASK_DNA5_LT[ordValue(right_)];
}

template <>
inline bool operator<=(Dna5 const & left_, Dna5Q const & right_)
{
    return __MASK_DNA5_LT[ordValue(left_)] <= __MASK_DNA5Q_LT[ordValue(right_)];
}

template <>
inline bool operator<=(Dna5Q const & left_, Dna5 const & right_)
{
    return __MASK_DNA5Q_LT[ordValue(left_)] <= __MASK_DNA5_LT[ordValue(right_)];
}

template <>
inline bool operator>=(Dna5 const & left_, Dna5Q const & right_)
{
    return __MASK_DNA5_LT[ordValue(left_)] >= __MASK_DNA5Q_LT[ordValue(right_)];
}

template <>
inline bool operator>=(Dna5Q const & left_, Dna5 const & right_)
{
    return __MASK_DNA5Q_LT[ordValue(left_)] >= __MASK_DNA5_LT[ordValue(right_)];
}

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Functions ordLess/Equal/Greater()                             [Dna5 vs Dna5]
// ----------------------------------------------------------------------------
    
template <>
inline bool ordLess(Dna5 const & left_, Dna5 const & right_)
{
    return __MASK_DNA5_LT[ordValue(left_)] < __MASK_DNA5Q_LT[ordValue(right_)];
}

template <>
inline bool ordEqual(Dna5 const & left_, Dna5 const & right_)
{
    return __MASK_DNA5_EQ[ordValue(left_)] & __MASK_DNA5_EQ[ordValue(right_)];
}

template <>
inline bool ordGreater(Dna5 const & left_, Dna5 const & right_)
{
    return __MASK_DNA5_LT[ordValue(left_)] > __MASK_DNA5Q_LT[ordValue(right_)];
}

// ----------------------------------------------------------------------------
// Function ordEqual()
// ----------------------------------------------------------------------------
// This function is overloaded to avoid casting TValue2 to Dna.

template <typename TValue2>
inline bool ordEqual(Dna const & left, TValue2 const & right)
{
    return ordValue(left) == ordValue(right);
}

}

#endif  // #ifndef APP_YARA_BASIC_ALPHABET_H_

