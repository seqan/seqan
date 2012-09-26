// ==========================================================================
//                          Mason - A Read Simulator
// ==========================================================================
// Copyright (C) 2010 Manuel Holtgrewe, FU Berlin
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License along
// with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// ==========================================================================
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Base-calling support code for 454 reads.
// ==========================================================================

#ifndef SIMULATE_454_BASE_CALLING_H_
#define SIMULATE_454_BASE_CALLING_H_

#include <cmath>

using namespace seqan;

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
    mutable String<double> _data;

    ThresholdMatrix()
            : _k(0), _useSqrt(false), _noiseMu(0), _noiseSigma(0), _size(0)
    {}
    
    ThresholdMatrix(double k, bool useSqrt, double noiseMu, double noiseSigma)
            : _k(k), _useSqrt(useSqrt), _noiseMu(noiseMu), _noiseSigma(noiseSigma), _size(0)
    {}

    ThresholdMatrix(double k, bool useSqrt, double noiseMean, double noiseStdDev, MeanStdDev const &)
            : _k(k), _useSqrt(useSqrt),
              _noiseMu(::std::log(noiseMean) - 0.5 * ::std::log(1.0 + noiseStdDev * noiseStdDev / noiseMean / noiseMean)),
              _noiseSigma(::std::sqrt(::std::log(1.0 + noiseStdDev * noiseStdDev / noiseMean / noiseMean))),
              _size(0)
    {}
};

inline double
normalDensityF(double x, double mu, double sigma)
{
    const double PI = 3.14159265;
    double sigma2 = sigma * sigma;
    return exp(- (x - mu) * (x - mu) / (2 * sigma2)) / sqrt(2 * PI * sigma2);
}

inline double
lognormalDensityF(double x, double mu, double sigma)
{
    if (x <= 0)
        return 0;
    const double PI = 3.14159265;
    double sigma2 = sigma * sigma;
    double log_mu2 = (log(x) - mu) * (log(x) - mu);
    return exp(-log_mu2 / (2 * sigma2)) / (x * sigma * sqrt(2 * PI));
}

inline double
dispatchDensityFunction(ThresholdMatrix const & matrix, unsigned r, double x)
{
    if (r == 0) {
        return lognormalDensityF(x, matrix._noiseMu, matrix._noiseSigma);
    } else {
        double rd = static_cast<double>(r);
        return normalDensityF(x, rd, (matrix._useSqrt ? sqrt(rd) : rd));
    }
}

inline double
computeThreshold(ThresholdMatrix const & matrix, unsigned r1, unsigned r2)
{
    if (r1 > r2)
        return computeThreshold(matrix, r2, r1);
    // The epsilon we use for convergence detection.
    const double EPSILON = 0.00001;

    // In i, we will count the number of iterations so we can limit the maximal
    // number of iterations.
    unsigned i = 0;

    // f1 is the density function for r1 and f2 the density function for r2.

    // Pick left such that f1(left) > f2(left).
    double left = r1;
    if (left == 0) left = 0.23;
    while (dispatchDensityFunction(matrix, r1, left) <= dispatchDensityFunction(matrix, r2, left))
        left /= 2.0;
    // And pick right such that f1(right) < f2(right).
    double right = r2;
    if (right == 0) right = 0.5;
    while (dispatchDensityFunction(matrix, r1, right) >= dispatchDensityFunction(matrix, r2, right))
        right *= 2.;

    // Now, search for the intersection point.
    while (true) {
        SEQAN_ASSERT_LT_MSG(i, 1000u, "Too many iterations (%u)! r1 = %u, r2 = %u.", i, r1, r2);
        i += 1;

        double center = (left + right) / 2;
        double fCenter1 = dispatchDensityFunction(matrix, r1, center);
        double fCenter2 = dispatchDensityFunction(matrix, r2, center);
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
extendThresholds(ThresholdMatrix const & matrix, unsigned dim)
{
    // Allocate new data array for matrix.  Then compute values or copy
    // over existing ones.
    String<double> newData;
    resize(newData, dim * dim, Exact());
    for (unsigned i = 0; i < dim; ++i) {
        for (unsigned j = 0; j < dim; ++j) {
            if (i == j)
                continue;
            if (i < matrix._size && j < matrix._size)
                newData[i * dim + j] = matrix._data[i * matrix._size + j];
            else
                newData[i * dim + j] = computeThreshold(matrix, i, j);
        }
    }
    // Update matrix.
    assign(matrix._data, newData);
    matrix._size = dim;
}

inline double
getThreshold(ThresholdMatrix const & matrix, unsigned r1, unsigned r2)
{
    if (matrix._size <= r1 || matrix._size <= r2)
        extendThresholds(matrix, _max(r1, r2) + 1);
    return matrix._data[r1 * matrix._size + r2];
}

inline void
setK(ThresholdMatrix & matrix, double k)
{
    matrix._k = k;
}

inline void
setUseSqrt(ThresholdMatrix & matrix, bool useSqrt)
{
    matrix._useSqrt = useSqrt;
}

inline void
setNoiseMu(ThresholdMatrix & matrix, double mu)
{
    matrix._noiseMu = mu;
}

inline void
setNoiseSigma(ThresholdMatrix & matrix, double sigma)
{
    matrix._noiseSigma = sigma;
}

inline void
setNoiseMeanStdDev(ThresholdMatrix & matrix, double mean, double stdDev)
{
    matrix._noiseMu = ::std::log(mean) - 0.5 * ::std::log(1.0 + stdDev * stdDev / mean / mean);
    matrix._noiseSigma = ::std::sqrt(::std::log(1.0 + stdDev * stdDev / mean / mean));
}

#endif  // SIMULATE_454_BASE_CALLING_H_
