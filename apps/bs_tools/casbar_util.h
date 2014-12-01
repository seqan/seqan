/*==========================================================================

Functors and functions for likelihood function evaluation

==========================================================================*/

#ifndef APPS_BS_TOOLS_CASBAR_UTIL_H__
#define APPS_BS_TOOLS_CASBAR_UTIL_H__

#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>

using namespace std;
using namespace seqan;

template<typename TDna>
bool
isTransition(TDna base1, TDna base2)  // why not error with & ?
{
    if (base1 == 'A' && base2 == 'G') return true;
    else if (base1 == 'G' && base2 == 'A') return true;
    else if (base1 == 'C' && base2 == 'T') return true;
    else if (base1 == 'T' && base2 == 'C') return true;
    else return false;
}

template<typename TMethOptions, typename TOptions>
void
computeGenotypePriors(TMethOptions &methOptions, TOptions &options)
{
    resize(methOptions.genPriors, 4*4*4, -1.0, Exact());
    long double r_j;
    long double r_k;

    if(options._debugLevel > 1)
    {
        std::cout << "Genotype prior probabilities: "<< std::endl;
        std::cout << "Ref." << '\t' << "Allel1" << '\t' << "Allele2" << std::endl;
    }

    if (!methOptions.uniformGenPriors)
    {
        for (unsigned i = 0; i < 4; ++i)
        {
            for (unsigned j = 0; j < 4; ++j)
            {
                for (unsigned k = 0; k < 4; ++k)
                {   // transitions are 4 x as frequent as transversions, but keep sum the same
                    r_j = (isTransition((Dna)i, (Dna)j)? (double)(12.0/6.0):(double)(3.0/6.0));
                    r_k = (isTransition((Dna)i, (Dna)k)? (double)(12.0/6.0):(double)(3.0/6.0));

                    if (j == i && k == i)
                        methOptions.genPriors[i<<4| j<<2| k] = (1.0-options.pHetSnp -options.pHomoSnp - options.pHetSnp*options.pHomoSnp);
                    else if (j == i)
                        methOptions.genPriors[i<<4| j<<2| k] = options.pHetSnp*r_k*(3.0/16.0);
                    else if (k == i)
                        methOptions.genPriors[i<<4| j<<2| k] = options.pHetSnp*r_j*(3.0/16.0);
                    else if (j == k)
                        methOptions.genPriors[i<<4| j<<2| k] = options.pHomoSnp*r_j*(3.0/16.0);    // + pError ?
                    else
                        methOptions.genPriors[i<<4| j<<2| k] = options.pHetSnp * options.pHomoSnp * (6.0/16.0);  // pHetSnp*(4/6)*0.0005*(4/6) + pRefError*pHetSnp*(4/6)

                    if(options._debugLevel > 1)
                    {
                        std::cout << (Dna)i << '\t' << (Dna)j << '\t' << (Dna)k << "\t\t" << (i<<4| j<<2| k) << '\t' << std::setprecision (25) << methOptions.genPriors[i<<4| j<<2| k] << std::endl;
                    }
                }
            }
        }
    }
    else
    {
        for (unsigned i = 0; i < 4; ++i)
        {
            for (unsigned j = 0; j < 4; ++j)
            {
                for (unsigned k = 0; k < 4; ++k)
                {
                    if (j == i && k == i)
                        methOptions.genPriors[i<<4| j<<2| k] = (1.0-options.pHetSnp -options.pHomoSnp - options.pHetSnp*options.pHomoSnp);
                    else if (j == i || k == i)
                        methOptions.genPriors[i<<4| j<<2| k] = options.pHetSnp*(6.0/16.0);
                    else if (j == k)
                        methOptions.genPriors[i<<4| j<<2| k] = options.pHomoSnp*(3.0/16.0);
                    else
                        methOptions.genPriors[i<<4| j<<2| k] = options.pHetSnp * options.pHomoSnp * (6.0/16.0);
               }
            }
        }
    }
}


// Prob. function
struct NaiveMult_;
typedef Tag<NaiveMult_> NaiveMult;

struct LogFunction_;
typedef Tag<LogFunction_> LogFunction;


// Optimize
struct Sampling_;
typedef Tag<Sampling_> Sampling;

struct Newton_;
typedef Tag<Newton_> Newton;


// NSpace
// Functor for  naive evaliuation
// Returns values for f(x)
template <typename TValue>
struct FctNaive_0N
{
    FctNaive_0N(String<String<TValue> > const& constants) : constants(constants)
    { // Constructor
    }
    TValue operator()(TValue const&beta)
    { // beta is estimate so far.
        TValue f = 1.0;
        for (unsigned i = 0; i < length(constants[0]); ++i)
        {
            f *= (constants[0][i] + beta*constants[1][i]);
        }

        return f;
    }
private:
    String<String<TValue> > constants;
};


// Functor for log function
// Returns values for f(x)
template <typename TValue>
struct FctLog_02N
{
    FctLog_02N(String<String<TValue> > const& constants) : constants(constants)
    { // Constructor
    }
    ::boost::math::tuple<TValue, TValue> operator()(TValue const&beta)
    {
        TValue f_0 = 0;
        TValue f_2 = 0;
        for (unsigned i = 0; i < length(constants[0]); ++i)
        {
            f_0 += std::log10(constants[0][i] + constants[1][i]*beta);
            f_2 += -pow(constants[1][i],2)/( pow(constants[0][i] , 2) + 2*constants[0][i]*constants[1][i]*beta + pow(constants[1][i], 2)*pow(beta, 2) );
        }
        f_2 *= 1.0/(std::log((double)10));

        return  boost::math::make_tuple(f_0, f_2);
    }
private:
    String<String<TValue> > constants;
};

// Functor for log function
// Returns values for f'(x) and f''(x)
template <typename TValue>
struct FctLog_12N
{
    FctLog_12N(String<String<TValue> > const& constants) : constants(constants)
    { // Constructor
    }
    ::boost::math::tuple<TValue, TValue> operator()(TValue const&beta)
    {
        TValue f_1 = 0;
        TValue f_2 = 0;
        for (unsigned i = 0; i < length(constants[0]); ++i)
        {
            f_1 += constants[1][i]/(constants[0][i] + constants[1][i]*beta);
            f_2 += -pow(constants[1][i],2)/( pow(constants[0][i], 2) + 2*constants[0][i]*constants[1][i]*beta + pow(constants[1][i], 2)*pow(beta, 2) );
        }

        f_1 *= 1.0/(std::log((double)10)) ;  // 1/(ln*b)
        f_2 *= 1.0/(std::log((double)10)) ;

        return boost::math::make_tuple(f_1, f_2);
    }
private:
    String<String<TValue> > constants;
};


#endif
