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
// Author: Rene Rahn <rene.rahn@fu-berlin.de>
// ==========================================================================
// Shuffling and conversion
// ==========================================================================

#ifndef INCLUDE_SEQAN_RANDOM_RANDOM_UTIL_H_
#define INCLUDE_SEQAN_RANDOM_RANDOM_UTIL_H_

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

// ===========================================================================
// Functions
// ===========================================================================

/*!
 * @fn shuffle
 * @headerfile <seqan/random.h>
 * @brief Shuffle the given container.
 *
 * @signature void shuffle(container, rng);
 *
 * @param[in,out] container The container to shuffle.
 * @param[in,out] rng       The random number generator to use.
 */

template <typename TContainer, typename TRNG>
void shuffle(TContainer & container, TRNG & rng)
{
    typedef typename Position<TContainer>::Type TPosition;
    typedef typename Value<TContainer>::Type TValue;

    TValue tmp;
    for (TPosition i = 0, iend = length(container); i < iend; ++i) {
        std::uniform_int_distribution<TPosition> uniformDist(i, iend - 1);
        TPosition j = uniformDist(rng);
        // swap
        move(tmp, container[i]);
        move(container[i], container[j]);
        move(container[j], tmp);
    }
}

/*!
 * @fn calcBetaDistParam
 * @headerfile <seqan/random.h>
 * @brief Computes alpha and beta paramters for the beta distribution from mean and stddev of the underlying distribution.
 *
 * @signature param cvtBetaDistParam(mean, stddev);
 *
 * @param[in] mean   The mean of the underlying distribution.
 * @param[in] stddev The standard deviation of the underlying distribution.
 *
 * @return param The distribution specific paramter of type <tt> typename beta_distribution::param_type </tt>.
 */

template<typename TRealType>
auto cvtBetaDistParam(TRealType const & mean, TRealType const & stddev)
    -> typename BetaDistribution<TRealType>::param_type
{
    using TParamType = typename BetaDistribution<TRealType>::param_type;
    return TParamType((((1 - mean) / stddev / stddev - 1 / mean) * mean * mean), ((((1 - mean) / stddev / stddev - 1 / mean) * mean * mean) * (1 / mean - 1)));
}

/*!
 * @fn calcLogNormalDistParam
 * @headerfile <seqan/random.h>
 * @brief Computes <tt>m</tt> and <tt>s</tt> paramters for the lognormal distribution from mean and stddev of the underlying distribution.
 *
 * @signature param cvtLogNormalDistParam(mean, stddev);
 *
 * @param[in] mean   The mean of the underlying distribution.
 * @param[in] stddev The standard deviation of the underlying distribution.
 *
 * @return param The distribution specific paramter of type <tt> typename lognormal_distribution::param_type </tt>.
 *
 */
// compute paramters m and s for lognorm from mean and stddev.
template<typename TRealType>
auto cvtLogNormalDistParam(TRealType const & mean, TRealType const & stddev)
    -> typename std::lognormal_distribution<TRealType>::param_type
{
    using TParamType = typename std::lognormal_distribution<TRealType>::param_type;
    return TParamType(static_cast<TRealType>(std::log(mean) - 0.5 * std::log(1.0 + stddev * stddev / mean / mean)),
                      static_cast<TRealType>(std::sqrt(std::log(1.0 + stddev * stddev / mean / mean))));
}

}  // namespace seqan

#endif  // INCLUDE_SEQAN_RANDOM_RANDOM_UTIL_H_
