// ==========================================================================
//                                   ANISE
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
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

#include "separating_columns.h"

#include "rep_sep/local_variation_store.h"

namespace rep_sep {

namespace {  // anonymous namespace

// Numerical Recipes in C, Second Edition (1992)
inline double _gammln(float xx)
{
    double x,y,tmp,ser;
    static double cof[6] = { 76.18009172947146, -86.50532032941677, 24.01409824083091,
                             -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5 };
    y = x = xx;
    tmp = x+ 5.5;
    tmp -= (x + 0.5) * std::log(tmp);
    ser = 1.000000000190015;
    for (int j = 0; j <= 5; ++j)
        ser += cof[j] / ++y;

    return -tmp + log(2.5066282746310005 * ser / x);
}

// Numerical Recipes in C, Second Edition (1992)
// Returns: ln(n!)
inline double _factln(int n)
{
    static float a[101];

    SEQAN_ASSERT_GEQ_MSG(n, 0, "Negative factorial requested");
    if (n <= 1)
        return 0.0;
    else if (n <= 100)
        return a[n] ? a[n] : (a[n] = _gammln(n + 1.0));
    else
        return _gammln(n + 1.0);
}

#if 0  // unused, was used for Kececegliou method

// Numerical Recipes in C, Second Edition (1992)
// Numerical implementation of binomial coefficient
inline double binomialCoefficient(int n, int k)
{
    if (k > n)
        return 0;
    // The floor function cleans up roundoff errors
    return std::floor(0.5 + std::exp(_factln(n)-_factln(k)-_factln(n-k)));
}

// The cumulative binomial distribution function; numerical recipes.
inline double cumulativeBinomialDistribution(double p, int m, int n)
{
    double res = 0;
    for (int i = m; i <= n; ++i)
    {
        double x = std::pow(p, i);
        double y = std::pow((1 - p), (n - i));
        if (x != 0 && y != 0)
            res += binomialCoefficient(n, i) * x * y;
    }
    return res;
}
#endif  // #if 0

// Poisson probability function.
inline double pPoisson(double lambda, double k)
{
    return std::pow(lambda,k) / std::exp(_factln(k)) * std::exp(-lambda);
}

/*
 * Numerical Recipes in C, Second Edition (1992)
 * Numerical implementation of binomial coefficient
 */
inline double _nchoose(int n, int k)
{
	if (k > n)
        return 0;
	// The floor function cleans up roundoff errors
	return std::floor(0.5 + std::exp(_factln(n) - _factln(k) - _factln(n - k)));
}

// Hypergeometric probability function.
inline double pHyper(unsigned N, unsigned m, unsigned n, unsigned k)
{
    if (k > m)
        return 0;
    return (_nchoose(m, k) * _nchoose(N - m, n - k)) / _nchoose(N, n);
}

// From PHRED quality to error probability.
inline double errorProbability(char quality)
{
    double frac = (quality - 'I');
    frac *= -1;
    frac /= 10;
    double result = std::pow(10.0,frac);
    return result;
}

// ----------------------------------------------------------------------------
// Class SeparatingColumnEnumerator
// ----------------------------------------------------------------------------

// Enumerate separating columngs for a LocalVariationStore.
//
// N's are not counted as deviations from the consensus.

class SeparatingColumnEnumerator
{
public:
    // The configuration for the column separation.
    ReadSeparatorOptions options;
    // The LocalVariationStore to enumerate the columns for.
    LocalVariationStore const * store;

    SeparatingColumnEnumerator() : store()
    {}

    SeparatingColumnEnumerator(LocalVariationStore const & store) : store(&store)
    {}

    // Write out the (sorted) indices of the separating columns of store into columnIds);
    inline void run(seqan::String<unsigned> & columnIds);

    // Test columns for being a separating pairs using Tammi's algorithm, depending on the configuration in options.
    //
    // readIds are the reads aligning to both columns u and v.
    inline bool testColumns(unsigned u, unsigned v, std::vector<unsigned> const & readIds);
};

bool SeparatingColumnEnumerator::testColumns(
        unsigned u, unsigned v, std::vector<unsigned> const & readIds)
{
    using namespace seqan;

    // The number of deviating values in columns u and v.
    unsigned nU = 0, nV = 0;
    // The number of correlating deviations (not necessarily to the same non-consensus base) in columns u and v.
    unsigned cObs = 0;

    // Compute nU, nV, and cObs.
    typedef LocalVariationStore::TColumnEntry TColumnEntry;
    typedef StringSet<String<TColumnEntry> > TColumns;
    typedef seqan::Value<TColumns>::Type TColumn;
    typedef seqan::Iterator<TColumn const, seqan::Standard>::Type TColumnIter;
    if (options.verbosity >= 3)
        std::cerr << "u == " << u << "\tv ==" << v << "\n";

    // Build columns, reduced to the reads in readIds.  The columns are sorted by read ID (member i1).
    std::set<unsigned> readIDSet(readIds.begin(), readIds.end());
    TColumn colU, colV;
    for (TColumnIter itU = begin(store->columns[u], seqan::Standard()); itU != end(store->columns[u], seqan::Standard()); ++itU)
        if (readIDSet.count(itU->i1))
            appendValue(colU, *itU);
    for (TColumnIter itV = begin(store->columns[v], seqan::Standard()); itV != end(store->columns[v], seqan::Standard()); ++itV)
        if (readIDSet.count(itV->i1))
            appendValue(colV, *itV);

    if (options.verbosity >= 3)
    {
        std::cerr << "column " << u << "\n";
        for (unsigned i = 0; i < length(colU); ++i)
            std::cerr << colU[i].i1 << "\t" << colU[i].i2 << "\t" << colU[i].i3 << "\n";
        std::cerr << "column " << v << "\n";
        for (unsigned i = 0; i < length(colV); ++i)
            std::cerr << colV[i].i1 << "\t" << colV[i].i2 << "\t" << colV[i].i3 << "\n";
    }

#if SEQAN_ENABLE_DEBUG
    // Check that the columns have the same size and are for the same reads.
    SEQAN_ASSERT_EQ(length(colU), length(colV));
    for (unsigned i = 0; i < length(colU); ++i)
        SEQAN_ASSERT_EQ(colU[i].i1, colV[i].i1);
#endif  // #if SEQAN_ENABLE_DEBUG

    // Re-compute consensus for both columns, might differ from global ones.
    LocalVariationStore::TProfileChar profileU, profileV;
    for (unsigned i = 0; i < length(colU); ++i)
    {
        SEQAN_ASSERT(store->compressionMap.count(colU[i].i1));
        unsigned numReads = length(store->compressionMap.find(colU[i].i1)->second);

        profileU.count[ordValue(colU[i].i2)] += numReads;
        profileV.count[ordValue(colV[i].i2)] += numReads;
    }
    LocalVariationStore::TConsensusAlphabet consensusU(_getMaxIndex(profileU)), consensusV(_getMaxIndex(profileV));

    // Count number of deviations.
    for (unsigned i = 0; i < length(colU); ++i)
    {
        SEQAN_ASSERT(store->compressionMap.count(colU[i].i1));
        unsigned numReads = length(store->compressionMap.find(colU[i].i1)->second);

        if (options.verbosity >= 3)
            std::cerr << "colU[i].i2 == " << colU[i].i2 << ", store->consensus[u] == " << store->consensus[u]
                      << ", colV[i].i2 == " << colV[i].i2 << ", store->consensus[v] == " << store->consensus[v]
                      << ", numReads == " << numReads << "\n";
        
        bool devU = (colU[i].i2 != 'N' && colU[i].i2 != consensusU);
        bool devV = (colV[i].i2 != 'N' && colV[i].i2 != consensusV);
        nU += devU * numReads;
        nV += devV * numReads;
        cObs += (devU && devV) * numReads;
    }

    if (options.verbosity >= 3)
        std::cerr << "u == " << u << ", v == " << v << "\n"
                  << "nU == " << nU << ", nV == " << nV << ", cObs == " << cObs << "\n\n";

    // We will compute the values of lambdaU, lambdaV, and pCorr depending on th method selected.
    double lambdaU = 0, lambdaV = 0, pCorr = 0;

    if (options.tammiMethod == ReadSeparatorOptions::TAMMI_SIMPLE)
    {
        // The number of reads in both columns.
        unsigned overlap = 0;
        for (unsigned i = 0; i < readIds.size(); ++i)
        {
            SEQAN_ASSERT(store->compressionMap.count(readIds[i]));
            overlap += length(store->compressionMap.find(readIds[i])->second);
        }

        // Compute lambdaU and lambdaV for constant per-position error probability.
        lambdaU = overlap * options.pErr;
        lambdaV = overlap * options.pErr;

        if (options.verbosity >= 3)
            std::cerr << "overlap == " << overlap << ", lambdaU == " << lambdaU << ", lambdaV == " << lambdaV << "\n";

        // Compute p^{corr}.
        pCorr = 1.0;
        for (unsigned i = 0; i < cObs; ++i)
            pCorr -= pHyper(overlap, nV, nU, i);
        pCorr = std::max(pCorr, 0.0);  // against underflow.

        if (options.verbosity >= 3)
            std::cerr << "pCorr == " << pCorr << "\n";
    }
    else if (options.tammiMethod == ReadSeparatorOptions::TAMMI_PHRED)
    {
        // Build set of read ids to quickly test for being contained in overlap.
        // TODO(holtgrew): Use scanning algorithm as above to get rid of log(n) factor?
        std::set<unsigned> readIdSet;
        std::copy(readIds.begin(), readIds.end(), std::inserter(readIdSet, readIdSet.end()));
        std::vector<double> probsU, probsV;

        // Compute lambdaU and lambdaV.
        for (unsigned i = 0; i < length(colU); ++i)
            if (readIdSet.count(colU[i].i1))
            {
                for (unsigned j = 0; j < length(store->compressionMap.find(colU[i].i1)->second); ++j)
                {
                    probsU.push_back(errorProbability(colU[i].i3));
                    lambdaU += probsU.back();
                }
            }
        for (unsigned i = 0; i < length(colV); ++i)
            if (readIdSet.count(colV[i].i1))
            {
                for (unsigned j = 0; j < length(store->compressionMap.find(colU[i].i1)->second); ++j)
                {
                    probsV.push_back(errorProbability(colV[i].i3));
                    lambdaV += probsV.back();
                }
            }
        SEQAN_CHECK(probsU.size() == probsV.size(), "Invalid columns (this probably is a bug)!");

        // Compute E.
        double E = 0;
        for (unsigned i = 0; i < probsU.size(); ++i)
        {
            double puj = probsU[i];
            double pvj = probsV[i];
            E += ( (cObs * puj) / ( (cObs * puj) + ( (lambdaU - puj) * (1 - puj) ) ) )
                    * ( (cObs * pvj) / ( (cObs * pvj) + ( (lambdaV - pvj) * (1 - pvj) ) ) );
        }

        // Compute p^{corr}
        pCorr = 1.0;
        for (unsigned i = 0; i < cObs; ++i)
            pCorr -= pPoisson(E, i);
        pCorr = std::max(pCorr, 0.0);  // against underflow
    }
    else
    {
        SEQAN_FAIL("Only Tammi's phred and simple method are implemented.");
    }

    // Compute p_u^{col} and p_v^{col}.
    double pColU = 1.0;
    for (unsigned i = 0; i < nU; ++i)
        pColU -= pPoisson(lambdaU, i);
    pColU = std::max(pColU, 0.0);  // against underflow.
    double pColV = 1.0;
    for (unsigned i = 0; i < nV; ++i)
        pColV -= pPoisson(lambdaV, i);
    pColV = std::max(pColV, 0.0);  // against underflow.

    // Compute final p^{tot}.
    double pCol = pColU + pColV - pColU * pColV;
    double pTot = pCol * pCorr;

    if (options.verbosity >= 3)
        std::cerr << "pColU == " << pColU << ", pColV == " << pColV
                  << ", pCol == " << pCol << ", pCorr == " << pCorr
                  << ", pTot == " << pTot << "\n";

    if (options.verbosity >= 3)
        std::cerr << "SEPARATING\t" << u << "\t" << v << "? => " << (pTot <= options.maxRandomCorrelation) << "\n";
    return (pTot <= options.maxRandomCorrelation);
}

void SeparatingColumnEnumerator::run(seqan::String<unsigned> & columnIds)
{
    using namespace seqan;

    std::set<unsigned> result;

    std::vector<unsigned> tmp;
    for (unsigned i = 0; i < length(store->coveringReads); ++i)
    {
        for (unsigned j = i + 1; j < length(store->coveringReads); ++j)
        {
            tmp.clear();
            std::set_intersection(
                    begin(store->coveringReads[i], seqan::Standard()),
                    end(store->coveringReads[i], seqan::Standard()),
                    begin(store->coveringReads[j], seqan::Standard()),
                    end(store->coveringReads[j], seqan::Standard()),
                    std::back_inserter(tmp));
            if (tmp.size() < options.minCommonReads)
                break;  // columns right of j cannot overlap sufficiently with i
            if (result.count(i) && result.count(j))
                continue;  // skip, both columns already in result

            if (options.verbosity >= 3)
                std::cerr << "Testing i=" << i << " (" << store->positions[i].first << ", " << store->positions[i].second << ") "
                          << "j=" << j << "(" << store->positions[j].first << ", " << store->positions[j].second << ")\n";

            if (testColumns(i, j, tmp))
            {
                result.insert(i);
                result.insert(j);
            }
        }
    }

    // Write out result from vector into string.
    resize(columnIds, result.size(), 0);
    std::copy(result.begin(), result.end(), begin(columnIds, seqan::Standard()));
}

}  // anonymous namespace


// ----------------------------------------------------------------------------
// Function enumerateSeparatingColumns()
// ----------------------------------------------------------------------------

void enumerateSeparatingColumns(seqan::String<unsigned> & result,
                                LocalVariationStore const & store)
{
    SeparatingColumnEnumerator sepEnum(store);
    sepEnum.run(result);
}

void enumerateSeparatingColumns(seqan::String<unsigned> & result,
                                LocalVariationStore const & store,
                                ReadSeparatorOptions const & options)
{
    SeparatingColumnEnumerator sepEnum(store);
    sepEnum.options = options;
    sepEnum.run(result);
}

}  // namespace rep_sep
