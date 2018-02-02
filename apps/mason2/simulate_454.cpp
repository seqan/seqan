// ==========================================================================
//                         Mason - A Read Simulator
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

#include "sequencing.h"

// Maximal homopolymer length we will observe.
const unsigned MAX_HOMOPOLYMER_LEN = 40;

// ===========================================================================
// Class ThresholdMatrix
// ===========================================================================

// 454 Model

class ThresholdMatrix
{
public:
    // The scaling parameter k.
    double _k;
    // Whether or not to use the sqrt for the std deviation computation.
    bool _useSqrt;
    // Mean of the log normally distributed noise.
    double _noiseMu;
    // Standard deviation of the log normally distributed noise.
    double _noiseSigma;
    // The edge length of the matrix.
    mutable unsigned _size;
    // The data of the matrix.
    mutable seqan::String<double> _data;

    ThresholdMatrix()
            : _k(0), _useSqrt(false), _noiseMu(0), _noiseSigma(0), _size(0)
    {}

    ThresholdMatrix(double k, bool useSqrt, double noiseMu, double noiseSigma)
            : _k(k), _useSqrt(useSqrt), _noiseMu(noiseMu), _noiseSigma(noiseSigma), _size(0)
    {}

    inline double
    computeThreshold(unsigned r1, unsigned r2) const
    {
        if (r1 > r2)
            return computeThreshold(r2, r1);
        // The epsilon we use for convergence detection.
        const double EPSILON = 0.00001;

        // In i, we will count the number of iterations so we can limit the maximal
        // number of iterations.
        unsigned i = 0;

        // f1 is the density function for r1 and f2 the density function for r2.

        // Pick left such that f1(left) > f2(left).
        double left = r1;
        if (left == 0) left = 0.23;
        while (dispatchDensityFunction(r1, left) <= dispatchDensityFunction(r2, left))
            left /= 2.0;
        // And pick right such that f1(right) < f2(right).
        double right = r2;
        if (right == 0) right = 0.5;
        while (dispatchDensityFunction(r1, right) >= dispatchDensityFunction(r2, right))
            right *= 2.;

        // Now, search for the intersection point.
        while (true)
        {
            SEQAN_ASSERT_LT_MSG(i, 1000u, "Too many iterations (%u)! r1 = %u, r2 = %u.", i, r1, r2);
            i += 1;

            double center = (left + right) / 2;
            double fCenter1 = dispatchDensityFunction(r1, center);
            double fCenter2 = dispatchDensityFunction(r2, center);
            double delta = fabs(fCenter1 - fCenter2);
            if (delta < EPSILON)
                return center;

            if (fCenter1 < fCenter2)
                right = center;
            else
                left = center;
        }
    }

    inline void
    extendThresholds(unsigned dim) const
    {
        // Allocate new data array for matrix.  Then compute values or copy
        // over existing ones.
        seqan::String<double> newData;
        resize(newData, dim * dim);
        for (unsigned i = 0; i < dim; ++i) {
            for (unsigned j = 0; j < dim; ++j) {
                if (i == j)
                    continue;
                if (i < _size && j < _size)
                    newData[i * dim + j] = _data[i * _size + j];
                else
                    newData[i * dim + j] = computeThreshold(i, j);
            }
        }
        // Update matrix.
        assign(_data, newData);
        _size = dim;
    }

    inline double
    getThreshold(unsigned r1, unsigned r2) const
    {
        if (_size <= r1 || _size <= r2)
            extendThresholds(std::max(r1, r2) + 1);
        return _data[r1 * _size + r2];
    }

    inline void
    setK(double k)
    {
        _k = k;
    }

    inline void
    setUseSqrt(bool useSqrt)
    {
        _useSqrt = useSqrt;
    }

    inline void
    setNoiseMu(double mu)
    {
        _noiseMu = mu;
    }

    inline void
    setNoiseSigma(double sigma)
    {
        _noiseSigma = sigma;
    }

    inline void
    setNoiseMeanStdDev(double mean, double stdDev)
    {
        auto tmp = seqan::cvtLogNormalDistParam(mean, stdDev);
        _noiseMu = tmp.m();
        _noiseSigma = tmp.s();
    }

    inline double
    normalDensityF(double x, double mu, double sigma) const
    {
        const double PI = 3.14159265;
        double sigma2 = sigma * sigma;
        return exp(- (x - mu) * (x - mu) / (2 * sigma2)) / sqrt(2 * PI * sigma2);
    }

    inline double
    lognormalDensityF(double x, double mu, double sigma) const
    {
        if (x <= 0)
            return 0;
        const double PI = 3.14159265;
        double sigma2 = sigma * sigma;
        double log_mu2 = (log(x) - mu) * (log(x) - mu);
        return exp(-log_mu2 / (2 * sigma2)) / (x * sigma * sqrt(2 * PI));
    }

    inline double
    dispatchDensityFunction(unsigned r, double x) const
    {
        if (r == 0) {
            return lognormalDensityF(x, _noiseMu, _noiseSigma);
        } else {
            double rd = static_cast<double>(r);
            return normalDensityF(x, rd, (_useSqrt ? sqrt(rd) : rd));
        }
    }
};

// ===========================================================================
// Class Roche454Model
// ===========================================================================

// Stores the threshold matrix.

class Roche454Model
{
public:
    ThresholdMatrix thresholdMatrix;
};

// ===========================================================================
// Class Roche454SequencingSimulator
// ===========================================================================

// ---------------------------------------------------------------------------
// Constructor Roche454SequencingSimulator::Roche454SequencingSimulator()
// ---------------------------------------------------------------------------

Roche454SequencingSimulator::Roche454SequencingSimulator(
        TRng & rng,
        TRng & methRng,
        SequencingOptions const & seqOptions,
        Roche454SequencingOptions const & roche454Options) :
        SequencingSimulator(rng, methRng, seqOptions), roche454Options(roche454Options), model(new Roche454Model())
{
    _initModel();
}

// ---------------------------------------------------------------------------
// Function Roche454SequencingSimulator::_initModel()
// ---------------------------------------------------------------------------

// Initialize the threshold matrix.
void Roche454SequencingSimulator::_initModel()
{
    model->thresholdMatrix.setK(roche454Options.k);
    model->thresholdMatrix.setUseSqrt(roche454Options.sqrtInStdDev);
    model->thresholdMatrix.setNoiseMeanStdDev(roche454Options.backgroundNoiseMean, roche454Options.backgroundNoiseStdDev);
}

// ---------------------------------------------------------------------------
// Function Roche454SequencingSimulator::readLength()
// ---------------------------------------------------------------------------

unsigned Roche454SequencingSimulator::readLength()
{
    if (roche454Options.lengthModel == Roche454SequencingOptions::UNIFORM)
    {
        // Pick uniformly.
        double minLen = roche454Options.minReadLength;
        double maxLen = roche454Options.maxReadLength;
        std::uniform_real_distribution<double> dist(minLen, maxLen);
        double len = dist(rng);
        return static_cast<unsigned>(round(len));
    }
    else
    {
        // Pick normally distributed.
        std::normal_distribution<double> dist(roche454Options.meanReadLength,
                                              roche454Options.stdDevReadLength);
        double len = dist(rng);
        return static_cast<unsigned>(round(len));
    }
}

// ---------------------------------------------------------------------------
// Function Roche454SequencingSimulator::simulateRead()
// ---------------------------------------------------------------------------

void Roche454SequencingSimulator::simulateRead(
        TRead & seq, TQualities & quals, SequencingSimulationInfo & info,
        TFragment const & frag, Direction dir, Strand strand)
{
    clear(seq);
    clear(quals);

    // Compute read length and check whether it fits in fragment.
    unsigned sampleLength = this->readLength();
    if (sampleLength > length(frag))
    {
        throw std::runtime_error("454 read is too long, increase fragment length");
    }

    // Get a copy of the to be sequenced base stretch.
    TRead haplotypeInfix;
    if (dir == LEFT)
        haplotypeInfix = prefix(frag, sampleLength);
    else
        haplotypeInfix = suffix(frag, length(frag) - sampleLength);
    if (strand == REVERSE)
        reverseComplement(haplotypeInfix);

    // In the flow cell simulation, we will simulate light intensities which will be stored in observedIntensities.
    seqan::String<double> observedIntensities;
    reserve(observedIntensities, 4 * sampleLength);
    seqan::Dna5String observedBases;
    // We also store the real homopolymer length.
    seqan::String<unsigned> realBaseCount;

    // Probability density function to use for the background noise.
    std::lognormal_distribution<double> distNoise(seqan::cvtLogNormalDistParam(roche454Options.backgroundNoiseMean,
                                                                               roche454Options.backgroundNoiseStdDev));

    // Initialize information about the current homopolymer length.
    unsigned homopolymerLength = 0;
    seqan::Dna homopolymerType = haplotypeInfix[0];
    while (homopolymerLength < length(haplotypeInfix) && haplotypeInfix[homopolymerLength] == homopolymerType)
        ++homopolymerLength;

    // Simulate flowcell.
    for (unsigned i = 0, j = 0; i < sampleLength; ++j, j = j % 4)  // i indicates first pos of current homopolymer, j indicates flow phase
    {
        if (ordValue(homopolymerType) == j)
        {
            // Simulate positive flow observation.
            double l = homopolymerLength;
            double sigma = roche454Options.k * (roche454Options.sqrtInStdDev ? sqrt(l) : l);
            std::normal_distribution<double> distIntensity(homopolymerLength, sigma);
            double intensity = distIntensity(rng);
            intensity += distNoise(rng);  // Add noise.
            appendValue(observedIntensities, intensity);
            appendValue(realBaseCount, homopolymerLength);
            // Get begin pos and length of next homopolymer.
            i += homopolymerLength;
            if (i < length(haplotypeInfix))
            {
                homopolymerType = haplotypeInfix[i];
                homopolymerLength = 0;
                while (((i + homopolymerLength) < length(haplotypeInfix)) && haplotypeInfix[i + homopolymerLength] == homopolymerType)
                    ++homopolymerLength;
            }
        }
        else
        {
            // Simulate negative flow observation.
            //
            // Constants taken from MetaSim paper which have it from the
            // original 454 publication.
            double intensity = std::max(0.0, distNoise(rng));
            appendValue(observedIntensities, intensity);
            appendValue(realBaseCount, 0);
        }
    }

    seqan::String<seqan::CigarElement<> > cigar;

    // Call bases, from this build the edit string and maybe qualities.  We only support the "inter" base calling
    // method which was published by the MetaSim authors in the PLOS paper.
    typedef seqan::Iterator<seqan::String<double>, seqan::Standard>::Type IntensitiesIterator;
    int i = 0;  // Flow round, Dna(i % 4) gives base.
    for (IntensitiesIterator it = begin(observedIntensities); it != end(observedIntensities); ++it, ++i)
    {
        double threshold = model->thresholdMatrix.getThreshold(static_cast<unsigned>(floor(*it)), static_cast<unsigned>(ceil(*it)));
        unsigned calledBaseCount = static_cast<unsigned>(*it < threshold ? floor(*it) : ceil(*it));
        // Add any matches.
        unsigned j = 0;
        for (; j < std::min(calledBaseCount, realBaseCount[i]); ++j)
        {
            appendOperation(cigar, 'M');
            appendValue(seq, seqan::Dna(i % 4));
        }
        // Add insertions, if any.
        for (; j < calledBaseCount; ++j)
        {
            appendOperation(cigar, 'I');
            appendValue(seq, seqan::Dna(i % 4));
        }
        // Add deletions, if any.
        for (; j < realBaseCount[i]; ++j)
            appendOperation(cigar, 'D');
        // Simulate qualities if configured to do so.
        if (seqOptions->simulateQualities)
        {
            // Compute likelihood for calling the bases, given this intensity and the Phred score from this.
            double densitySum = 0;
            for (unsigned j = 0; j <= std::max(4u, 2 * MAX_HOMOPOLYMER_LEN); ++j)  // Anecdotally through plot in maple: Enough to sum up to 4 or 2 times the maximal homopolymer length.
                densitySum += model->thresholdMatrix.dispatchDensityFunction(j, *it);
            double x = 0;  // Probability of seeing < (j+1) bases.
            for (unsigned j = 0; j < calledBaseCount; ++j) {
                x += model->thresholdMatrix.dispatchDensityFunction(j, *it);
                int q = -static_cast<int>(10 * ::std::log10(x / densitySum));
                q = std::max(0, std::min(40, q));
                appendValue(quals, (char)('!' + q));
            }
        }
    }

    // Write out extended sequencing information info if configured to do so.  We always write out the sample position
    // and alignment information.
    info.cigar = cigar;
    unsigned len = 0;
    _getLengthInRef(len, cigar);
    info.beginPos = (dir == LEFT) ? beginPosition(frag) : (beginPosition(frag) + length(frag) - len);
    info.isForward = (strand == FORWARD);

    if (seqOptions->embedReadInfo)
    {
        if (dir == LEFT)
            info.sampleSequence = prefix(frag, len);
        else
            info.sampleSequence = suffix(frag, length(frag) - len);
        if (strand == REVERSE)
            reverseComplement(info.sampleSequence);
    }
}
