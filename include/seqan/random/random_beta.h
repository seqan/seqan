
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
// Author: Rene Rahn <rene.rahn@fu-berlin.de>
// ==========================================================================
// beta distribution for random module
// build by two std::gamma_distribution's
// source: http://stackoverflow.com/questions/15165202/random-number-generator-with-beta-distribution.
// ==========================================================================

#ifndef INCLUDE_SEQAN_RANDOM_RANDOM_BETA_DIST_H_
#define INCLUDE_SEQAN_RANDOM_RANDOM_BETA_DIST_H_

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
 * @headerfile <seqan/random.h>
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
    typedef TRealType result_type;  // Conform the stl standard.

    // __ ParamType. ___________________________________________________________

    class param_type
    {
    public:
        typedef BetaDistribution distribution_type;  // Conform the stl standard.

        explicit param_type(TRealType a = 2.0, TRealType b = 2.0) : aParam(a), bParam(b)
        {}

        TRealType alpha() const { return aParam; }
        TRealType beta() const { return bParam; }

        bool operator==(const param_type& other) const
        {
            return (aParam == other.aParam && bParam == other.bParam);
        }

        bool operator!=(const param_type& other) const
        {
            return !(*this == other);
        }
    private:
        TRealType aParam, bParam;
    };

    // __ Constructor. ___________________________________________________________

    explicit BetaDistribution(TRealType a = 2.0, TRealType b = 2.0) : aGamma(a), bGamma(b)
    {}

    explicit BetaDistribution(const param_type& param) : aGamma(param.alpha()), bGamma(param.beta())
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
        aGamma = TGammaDist(param.alpha());
        bGamma = TGammaDist(param.beta());
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
    template <typename TRng>
    result_type operator()(TRng& engine)
    {
        return generate(engine, aGamma, bGamma);
    }

    template <typename TRng>
    result_type operator()(TRng& engine, param_type const& param)
    {
        TGammaDist aParamGamma(param.alpha()), bParamGamma(param.beta());
        return generate(engine, aParamGamma, bParamGamma);
    }

    /*!
     * @fn BetaDistribution#min
     * @brief Minimum value.
     *
     * @signature TRealType min();
     *
     * @return TRealType Returns the greatest lower bound of the range of values potentially returned by member operator(), which for <tt>BetaDistribution</tt> is always zero.
     */
    constexpr result_type min() const { return 0.0; }

    /*!
     * @fn BetaDistribution#max
     * @brief Maximum value.
     *
     * @signature TRealType max();
     *
     * @return TRealType Returns the least upper bound of the range of values potentially returned by member operator().
     */
    constexpr result_type max() const { return 1.0; }

    /*!
     * @fn BetaDistribution#alpha
     * @brief Returns the alpha parameter associated with the beta distribution.
     *
     * @signature TRealType alpha();
     *
     * @return TRealType Returns the alpha parameter.
     */
    result_type alpha() const { return aGamma.alpha(); }

    /*!
     * @fn BetaDistribution#beta
     * @brief Returns the beta parameter associated with the beta distribution.
     *
     * @signature TRealType beta();
     *
     * @return TRealType Returns the beta parameter.
     */
    result_type beta() const { return bGamma.alpha(); }

    bool operator==(const BetaDistribution<result_type>& other) const
    {
        return (aGamma == other.aGamma &&
                bGamma == other.bGamma);
    }

    bool operator!=(const BetaDistribution<result_type>& other) const
    {
        return !(*this == other);
    }

private:
    typedef std::gamma_distribution<result_type> TGammaDist;

    TGammaDist aGamma, bGamma;

    template <typename TRng>
    result_type
    generate(TRng& engine,
             TGammaDist& x_gamma,
             TGammaDist& y_gamma)
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

// ----------------------------------------------------------------------------
// Function operator<<
// ----------------------------------------------------------------------------

/*!
 * @fn BetaDistribution#operator<<
 * @brief Writes to output stream.
 *
 * @signature stream operator<<(stream, beta);
 *
 * Writes a textual representation of the distribution parameters and any other internal data kept by the object, 
 * in such a way that if this same text is inserted using <tt>operator>></tt> into an object of the same type, 
 * the same sequence of random numbers would be generated by equivalent invocations of <tt>operator()</tt>.
 *
 * @param[in,out] stream The stream to write to.
 * @param[in]     beta   An instance of the beta distribution.
 *
 * @return stream The updated stream to write too.
 * @see BetaDistribution#operator>>
 */

template <typename TStream, typename TRealType>
inline TStream &
operator<<(TStream& target,
           BetaDistribution<TRealType> const & beta)
{
    typename DirectionIterator<TStream, Output>::Type it = directionIterator(target, Output());
    write(it, "~Beta(");
    write(it, beta.alpha());
    write(it, ",");
    write(it, beta.beta());
    write(it, ")");
    return target;
}

// ----------------------------------------------------------------------------
// Function operator>>
// ----------------------------------------------------------------------------

/*!
 * @fn BetaDistribution#operator>>
 * @brief Reads from input stream.
 *
 * @signature stream operator>>(stream, beta);
 *
 * Restores the distribution parameters and any other internal data into distr from the textual 
 * representation provided by is.
 *
 * <tt>beta</tt> will generate the same sequence of random numbers as if equivalent invocations 
 * of <tt>operator()</tt> were performed on the object of the same type from which the provided 
 * text was obtained.
 *
 * @param[in,out] stream The stream to read from.
 * @param[in,out] beta   An instance of the beta distribution to read from input.
 *
 * @return stream The updated stream to write too.
 * @see BetaDistribution#operator>>
 */

template <typename TStream, typename TRealType>
inline bool
_streamBetaDistHelper(TStream & source,
                      BetaDistribution<TRealType> & beta)
{
    using TParam = typename BetaDistribution<TRealType>::param_type;
    typename DirectionIterator<TStream, Input>::Type it = directionIterator(source, Input());
    CharString buf;
    TRealType a, b;
    readUntil(buf, it, OrFunctor<EqualsChar<'('>, IsNewline>());
    if (buf != "~Beta" && *it != '(')
        return false;
    skipOne(it);
    clear(buf);
    readUntil(buf, it, OrFunctor<EqualsChar<','>, IsNewline>());
    if (*it != ',')
        return false;
    skipOne(it);
    lexicalCast(a, buf);
    clear(buf);
    readUntil(buf, it, OrFunctor<EqualsChar<')'>, IsNewline>());
    if (*it != ')')
        return false;
    skipOne(it);
    lexicalCast(b, buf);

    beta.setParam(TParam(a, b));
    return true;
}

template <typename TStream, typename TRealType>
inline TStream&
operator>>(TStream & source,
           BetaDistribution<TRealType> & beta)
{
    if (!_streamBetaDistHelper(source, beta))
        source.setstate(std::ios::failbit);
    return source;
}
    
} // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_RANDOM_RANDOM_BETA_DIST_H_
