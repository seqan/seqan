// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
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
// Author: Rene Rahn <rene.rahn@fu-berlin.de>
// ==========================================================================
// beta distribution for random module
// build by two std::gamma_distribution's
// source: http://stackoverflow.com/questions/15165202/random-number-generator-with-beta-distribution.
// ==========================================================================

#ifndef INCLUDE_SEQAN_MISC_RANDOM_UTIL_H_
#define INCLUDE_SEQAN_MISC_RANDOM_UTIL_H_

#include <iostream>
#include <sstream>
#include <string>
#include <seqan/basic.h>
#include <random>

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================



/*!
 * @class BetaDistribution
 * @headerfile <seqan/misc/random_util.h>
 * @brief Wrapper for beta distribution which is not part of the stl random module.
 *
 * @signature class BetaDistribution<TRealType>;
 * @tparam TRealType A floating-point type. Aliased as member type result_type. By default, this is double.
 *
 * This distribution produces random numbers that are distributed according to the beta distribution.
 *
 * @fn BetaDistribution::BetaDistribution
 * @brief Constructor for the beta distribution.
 *
 * @signature explicit BetaDistribution(a, b)
 * @signature explicit BetaDistribution(param)
 *
 * @param[in] a Tha alpha value of the beta distribution. Defaults to 2.0.
 * @param[in] b Tha beta value of the beta distribution. Defaults to 2.0.
 * @param[in] param Tha param type of the beta distribution.
 *
 */

template <typename TRealType = double>
class BetaDistribution
{
public:
    typedef TRealType result_type;

    // __ ParamType. ___________________________________________________________

    class param_type
    {
    public:
        typedef BetaDistribution distribution_type;

        explicit param_type(TRealType a = 2.0, TRealType b = 2.0) :
            a_param(a), b_param(b) { }

        TRealType a() const { return a_param; }
        TRealType b() const { return b_param; }

        bool operator==(const param_type& other) const
        {
            return (a_param == other.a_param && b_param == other.b_param);
        }

        bool operator!=(const param_type& other) const
        {
            return !(*this == other);
        }

    private:
        TRealType a_param, b_param;
    };

    // __ Constructor. ___________________________________________________________

    explicit BetaDistribution(TRealType a = 2.0, TRealType b = 2.0) : a_gamma(a), b_gamma(b)
    {}

    explicit BetaDistribution(const param_type& param) : a_gamma(param.a()), b_gamma(param.b())
    {}

    // ___ Member functions. _____________________________________________________

    void reset() {}

    /*!
     * @fn BetaDistribution#param
     * @brief Returns the parameters of this distribution.
     *
     * @signature TParamType param() const;
     *
     * @return TParamType The parameter of this distribution. Of type <tt> typename BetaDistribution\<TRealType\>::param_type </tt>.
     * @see BetaDistribution#setParam
     */

    param_type param() const
    {
        return param_type(alpha(), beta());
    }

    /*!
     * @fn BetaDistribution#setParam
     * @brief Sets the given parameters to this distribution.
     *
     * @signature void setParam(param);
     *
     * @param[in] param The parameters to be set to this distribution.
     * @see BetaDistribution#param
     */

    void setParam(const param_type& param)
    {
        a_gamma = gamma_dist_type(param.a());
        b_gamma = gamma_dist_type(param.b());
    }


    /*!
     * @fn BetaDistribution#operator()
     * @brief Generate random number.
     *
     * @signature TRealType operator(rng);
     *            TRealType operator(rng, param);
     * @param[in,out] rng The random number generator engine.
     * @param[in] param The parameter to used for the distribution.
     *
     * @return TRealType The generated random number.
     */
    template <typename URNG>
    result_type operator()(URNG& engine)
    {
        return generate(engine, a_gamma, b_gamma);
    }

    template <typename URNG>
    result_type operator()(URNG& engine, const param_type& param)
    {
        gamma_dist_type a_param_gamma(param.a()), b_param_gamma(param.b());
        return generate(engine, a_param_gamma, b_param_gamma);
    }

    /*!
     * @fn BetaDistribution#min
     * @brief Minimum value.
     *
     * @signature TRealType min();
     *
     * @return TRealType Returns the greatest lower bound of the range of values potentially returned by member operator(), which for <tt>BetaDistribution</tt> is always zero.
     */
    result_type min() const { return 0.0; }

    /*!
     * @fn BetaDistribution#max
     * @brief Maximum value.
     *
     * @signature TRealType max();
     *
     * @return TRealType Returns the least upper bound of the range of values potentially returned by member operator().
     */
    result_type max() const { return 1.0; }

    /*!
     * @fn BetaDistribution#alpha
     * @brief Returns the alpha parameter associated with the beta distribution.
     *
     * @signature TRealType alpha();
     *
     * @return TRealType Returns the alpha parameter.
     */
    result_type alpha() const { return a_gamma.alpha(); }

    /*!
     * @fn BetaDistribution#beta
     * @brief Returns the beta parameter associated with the beta distribution.
     *
     * @signature TRealType beta();
     *
     * @return TRealType Returns the beta parameter.
     */
    result_type beta() const { return b_gamma.alpha(); }

    bool operator==(const BetaDistribution<result_type>& other) const
    {
        return (param() == other.param() &&
                a_gamma == other.a_gamma &&
                b_gamma == other.b_gamma);
    }

    bool operator!=(const BetaDistribution<result_type>& other) const
    {
        return !(*this == other);
    }

private:
    typedef std::gamma_distribution<result_type> gamma_dist_type;

    gamma_dist_type a_gamma, b_gamma;

    template <typename URNG>
    result_type
    generate(URNG& engine,
             gamma_dist_type& x_gamma,
             gamma_dist_type& y_gamma)
    {
        result_type x = x_gamma(engine);
        return x / (x + y_gamma(engine));
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

template <typename CharT, typename TRealType>
std::basic_ostream<CharT>&
operator<<(std::basic_ostream<CharT>& os,
           const BetaDistribution<TRealType>& beta)
{
    os << "~Beta(" << beta.a() << "," << beta.b() << ")";
    return os;
}

template <typename CharT, typename TRealType>
std::basic_istream<CharT>&
operator>>(std::basic_istream<CharT>& is,
           BetaDistribution<TRealType>& beta)
{
    std::string str;
    TRealType a, b;
    if (std::getline(is, str, '(') && str == "~Beta" &&
        is >> a && is.get() == ',' && is >> b && is.get() == ')') {
        beta = BetaDistribution<TRealType>(a, b);
    } else {
        is.setstate(std::ios::failbit);
    }
    return is;
}

/*!
 * @fn calcBetaDistParam
 * @headerfile <seqan/misc/random_util.h>
 * @brief Computes paramters alpha and beta for the beta distribution from mean and stddev of the underlying distribution.
 *
 * @signature TParamType calcBetaDistParam(mean, stddev);
 *
 * @param[in] mean   The mean of the underlying distribution.
 * @param[in] stddev The standard deviation of the underlying distribution.
 *
 * @return TParamType The return type of type <tt> typename beta_distribution::param_type </tt>.
 */

template<typename TRealType>
auto calcBetaDistParam(TRealType const & mean, TRealType const & stddev)
    -> typename BetaDistribution<TRealType>::param_type
{
    using TParamType = typename BetaDistribution<TRealType>::param_type;
    return TParamType((((1 - mean) / stddev / stddev - 1 / mean) * mean * mean), ((((1 - mean) / stddev / stddev - 1 / mean) * mean * mean) * (1 / mean - 1)));
}

/*!
 * @fn calcLogNormalDistParam
 * @headerfile <seqan/misc/random_util.h>
 * @brief Computes paramters <tt>m</tt> and <tt>s</tt> for the lognormal distribution from mean and stddev of the underlying distribution.
 *
 * @signature TParamType calcLogNormalDistParam(mean, stddev);
 *
 * @param[in] mean   The mean of the underlying distribution.
 * @param[in] stddev The standard deviation of the underlying distribution.
 *
 * @return TParamType The return type of type <tt> typename lognormal_distribution::param_type </tt>.
 *
 */
// compute paramters m and s for lognorm from mean and stddev.
template<typename TRealType>
auto calcLogNormalDistParam(TRealType const & mean, TRealType const & stddev)
    -> typename std::lognormal_distribution<TRealType>::param_type
{
    using TParamType = typename std::lognormal_distribution<TRealType>::param_type;
    return TParamType(static_cast<TRealType>(std::log(mean) - 0.5 * std::log(1.0 + stddev * stddev / mean / mean)),
                      static_cast<TRealType>(std::sqrt(std::log(1.0 + stddev * stddev / mean / mean))));
}

} // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_MISC_RANDOM_UTIL_H_
