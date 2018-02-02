// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Basic definitions for the module random.
// ==========================================================================

#ifndef SEQAN_RANDOM_RANDOM_BASE_H_
#define SEQAN_RANDOM_RANDOM_BASE_H_

namespace seqan {

// ===========================================================================
// Forwards, Tags.
// ===========================================================================

// ===========================================================================
// Classes
// ===========================================================================

// ===========================================================================
// Metafunctions
// ===========================================================================

/*!
 * @mfn GetDefaultRng
 * @headerfile <seqan/random.h>
 * @brief Return the default Rng to use in a given class, specialization, or algorithm.
 *
 * @signature GetDefaultRng<T>::Type;
 *
 * @tparam T The type or tag to get the default Rng for.
 *
 * @return Type The Rng type to use. Defaults to <tt> std::mt19937 </tt>.
 *
 * @see defaultRng
 */

template <typename T>
struct GetDefaultRng
{
    typedef std::mt19937 Type;
};

// ===========================================================================
// Functions
// ===========================================================================

/*!
 * @fn defaultRng
 * @headerfile <seqan/random.h>
 * @brief Return the default random number generator object of a given type.
 *
 * @signature TRng defaultRng<TRng>();
 *
 * @tparam TRng The random number generator type to construct.
 *
 * @return TRng Default random number generator of the given type TRng.
 *
 * @section Remarks
 *
 * The random number generator will be default constructed, i.e. with the default seed.
 *
 * This function is NOT thread-safe!  Also, data structures and functions using defaultRng are not thread-safe.  Data
 * structures using global random number generator state should use pointers.  This way, the random number generator
 * state to be used can be set to be thread-local.
 *
 * @see GetDefaultRng
 */

template <typename TRng>
inline TRng &
defaultRng()
{
    static TRng x;
    return x;
}


template <typename T>
inline typename GetDefaultRng<T>::Type &
defaultRng(T const &)
{
    typedef typename GetDefaultRng<T>::Type TRng;
    return defaultRng<TRng>();
}

}  // namespace seqan

#endif  // SEQAN_RANDOM_RANDOM_BASE_H_
