// ==========================================================================
//                         Mason - A Read Simulator
// ==========================================================================
// Copyright (c) 2006-2024, Knut Reinert, FU Berlin
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
// Simulator for the sequencing process.
// ==========================================================================

#ifndef APPS_MASON2_SEQUENCING_H_
#define APPS_MASON2_SEQUENCING_H_

#include <stdexcept>
#include <random>

#include <seqan/bam_io.h>
#include <seqan/seq_io.h>
#include <seqan/modifier.h>
#include <seqan/sequence.h>

#include "mason_options.h"
#include "methylation_levels.h"

// ===========================================================================
// Class IlluminaSequencingOptions
// ===========================================================================

class IlluminaModel
{
public:
    // Probabilities for a mismatch at a given position.
    seqan2::String<double> mismatchProbabilities;

    // Standard deviations for the normal distributions of base qualities for the mismatch case.
    seqan2::String<double> mismatchQualityMeans;
    // Standard deviations for the normal distributions of base qualities for the mismatch case.
    seqan2::String<double> mismatchQualityStdDevs;

    // Standard deviations for the normal distributions of base qualities for the non-mismatch case.
    seqan2::String<double> qualityMeans;
    // Standard deviations for the normal distributions of base qualities for the non-mismatch case.
    seqan2::String<double> qualityStdDevs;

    IlluminaModel()
    {}
};

// ===========================================================================
// Class ThresholdMatrix
// ===========================================================================

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
    mutable seqan2::String<double> _data;

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
        [[maybe_unused]] unsigned i = 0;

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
        seqan2::String<double> newData;
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
        auto tmp = seqan2::cvtLogNormalDistParam(mean, stdDev);
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

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

typedef seqan2::Dna5String TRead;
typedef seqan2::CharString TQualities;
typedef std::mt19937 TRng;
typedef seqan2::Infix<seqan2::Dna5String const>::Type TFragment;
typedef seqan2::String<seqan2::CigarElement<> > TCigarString;

// ----------------------------------------------------------------------------
// Class SequencingSimulationInfo
// ----------------------------------------------------------------------------

// Composition of verbose information for sequencing simulation.
//
// Objects of this type store information such as the CIGAR string and the original sampled sequence for debug and
// evaluation purposes.

struct SequencingSimulationInfo
{
    // Originally sampled sequence, that together with the errors introduced by sequencing gives the read sequence.
    TRead sampleSequence;

    // The CIGAR string, MXID for matches, mismatches, insertions, deletions with respect to the reference.
    TCigarString cigar;

    // Whether or not this comes from the forward strand.
    bool isForward;

    // The contig and haplotype ID.
    int rID, hID;

    // The begin position of the sequence.
    int beginPos;

    // Number of bases covering SNPs and overlapping with indels.
    int snpCount, indelCount;

    SequencingSimulationInfo() : isForward(false), rID(-1), hID(-1), beginPos(-1), snpCount(0), indelCount(0)
    {}

    // Returns number in reference covered by this read.
    int lengthInRef() const
    {
        int result = 0;
        for (unsigned i = 0; i < length(cigar); ++i)
            if (cigar[i].operation == 'M' || cigar[i].operation == 'M' || cigar[i].operation == 'D')
                result += cigar[i].count;
        return result;
    }

    template <typename TStream>
    void serialize(TStream & stream) const
    {
        stream << "SEQUENCE=" << rID << " HAPLOTYPE=" << hID << " BEGIN_POS=" << beginPos
               << " SAMPLE_SEQUENCE=" << sampleSequence << " CIGAR=";
        for (unsigned i = 0; i < length(cigar); ++i)
            stream << cigar[i].count << cigar[i].operation;
        stream << " STRAND=" << (isForward ? 'F' : 'R');
        stream << " NUM_SNPS=" << snpCount << " NUM_INDELS=" << indelCount;
    }
};

// ----------------------------------------------------------------------------
// Class SequencingSimulator
// ----------------------------------------------------------------------------

// Responsible for the simulation of sequencing.
//
// We pick the read length independently of the error since the sequencing simulator might have context dependent
// errors.

class SequencingSimulator
{
public:
    // Sequencing direction from left or right side.
    enum Direction
    {
        LEFT,
        RIGHT
    };

    static Direction direction(Direction s)
    {
        return (s == LEFT) ? RIGHT : LEFT;
    }

    // Strand for sequencing.
    enum Strand
    {
        FORWARD,
        REVERSE
    };

    static Strand toggle(Strand s)
    {
        return (s == FORWARD) ? REVERSE : FORWARD;
    }

    // The random number generators for sequences and methylation/BS-treatment simulation.
    TRng & rng;
    TRng & methRng;
    // Overall sequencing options.
    SequencingOptions const * seqOptions;

    // Buffer for the materialization of BS-seq treated fragments.
    seqan2::Dna5String methFrag;

    SequencingSimulator(TRng & rng, TRng & methRng, SequencingOptions const & _options) :
            rng(rng), methRng(methRng), seqOptions(&_options)
    {}

    virtual ~SequencingSimulator() = 0;

    // Pick read length for the fragment to be simulated.
    virtual unsigned readLength() = 0;

    // Simulate paired-end sequencing from a fragment.
    //
    // If BS-seq is enabled in seqOptions->bsSeqOptions then levels must be != 0.
    void simulatePairedEnd(TRead & seqL, TQualities & qualsL, SequencingSimulationInfo & infoL,
                           TRead & seqR, TQualities & qualsR, SequencingSimulationInfo & infoR,
                           TFragment const & frag,
                           MethylationLevels const * levels = 0);

    // Simulate single-end sequencing from a fragment.
    //
    // If BS-seq is enabled in seqOptions->bsSeqOptions then levels must be != 0.
    void simulateSingleEnd(TRead & seq, TQualities & quals, SequencingSimulationInfo & info,
                           TFragment const & frag,
                           MethylationLevels const * levels = 0);

    // Actually simulate read and qualities from fragment and direction forward/reverse strand.
    //
    // seq -- target sequence of the read to simulate
    // quals -- target qualities of the read to simulate
    // frag -- source fragment
    // strand -- the strand of the fragment, coordinates are relative to forward and will be affected by this flag
    // dir -- whether this is the left or right read
    virtual void simulateRead(TRead & seq, TQualities & quals, SequencingSimulationInfo & info,
                              TFragment const & frag, Direction dir, Strand strand) = 0;

private:
    // Simulate BS-seq treatment on forward/reverse strand of frag with the given methylation levels.
    //
    // The result is a DNA string with the translations.
    void _simulateBSTreatment(seqan2::Dna5String & methFragment,
                              TFragment const & frag,
                              MethylationLevels const & levels,
                              bool reverse);

    // Implementation of the single-end and paired-end sequencing.  The functions without an underscore forward here.
    void _simulatePairedEnd(TRead & seqL, TQualities & qualsL, SequencingSimulationInfo & infoL,
                            TRead & seqR, TQualities & qualsR, SequencingSimulationInfo & infoR,
                            TFragment const & frag, bool isForward);

    // Simulate single-end sequencing from a fragment.
    //
    // If BS-seq is enabled in seqOptions->bsSeqOptions then levels must be != 0.
    void _simulateSingleEnd(TRead & seq, TQualities & quals, SequencingSimulationInfo & info,
                            TFragment const & frag, bool isForward);
};

inline SequencingSimulator::~SequencingSimulator() {}

// ----------------------------------------------------------------------------
// Class IlluminaSequencingSimulator
// ----------------------------------------------------------------------------

// Illumina read simulation.

class IlluminaSequencingSimulator : public SequencingSimulator
{
public:
    // Configuration for Illumina sequencing.
    IlluminaSequencingOptions illuminaOptions;

    // Storage for the Illumina simulation.
    std::unique_ptr<IlluminaModel> model;

    IlluminaSequencingSimulator(TRng & rng, TRng & methRng, SequencingOptions const & seqOptions,
                                IlluminaSequencingOptions const & illuminaOptions);

    // Pick read length for the fragment to be simulated.
    virtual unsigned readLength()
    {
        return illuminaOptions.readLength;
    }

    // Actually simulate read and qualities from fragment and direction forward/reverse strand.
    virtual void simulateRead(TRead & seq, TQualities & quals, SequencingSimulationInfo & info,
                              TFragment const & frag, Direction dir, Strand strand);

private:
    // Initialize the model.
    void _initModel();

    // Simulate PHRED qualities from the CIGAR string.
    void _simulateQualities(TQualities & quals, TCigarString const & cigar);

    // Simulate CIGAR string.  We can do this with position specific parameters only and thus independent of any
    // context.
    void _simulateCigar(TCigarString & cigar);
};

// ----------------------------------------------------------------------------
// Class Roche454SequencingSimulator
// ----------------------------------------------------------------------------

// 454 read simulation.

class Roche454SequencingSimulator : public SequencingSimulator
{
public:
    // Configuration for Roche454 sequencing.
    Roche454SequencingOptions roche454Options;

    // Precomputed model data for 454 Sequencing.
    std::unique_ptr<Roche454Model> model;

    Roche454SequencingSimulator(TRng & rng, TRng & methRng,
                                SequencingOptions const & seqOptions,
                                Roche454SequencingOptions const & roche454Options);

    // Initialize the threshold matrix.
    void _initModel();

    // Pick read length for the sequence to be sampled from fragments.
    virtual unsigned readLength();

    // Actually simulate read and qualities from fragment and direction forward/reverse strand.
    virtual void simulateRead(TRead & seq, TQualities & quals, SequencingSimulationInfo & info,
                              TFragment const & frag, Direction dir, Strand strand);
};

// ----------------------------------------------------------------------------
// Class SangerSequencingSimulator
// ----------------------------------------------------------------------------

// Sanger read simulation.

class SangerSequencingSimulator : public SequencingSimulator
{
public:
    // Configuration for Sanger sequencing.
    SangerSequencingOptions sangerOptions;

    SangerSequencingSimulator(TRng & rng, TRng & methRng,
                              SequencingOptions const & seqOptions,
                              SangerSequencingOptions const & sangerOptions) :
            SequencingSimulator(rng, methRng, seqOptions), sangerOptions(sangerOptions)
    {}

    // Pick read length for the sequence to be sampled from fragments.
    virtual unsigned readLength();

    // Actually simulate read and qualities from fragment and direction forward/reverse strand.
    virtual void simulateRead(TRead & seq, TQualities & quals, SequencingSimulationInfo & info,
                              TFragment const & frag, Direction dir, Strand strand);

    // Simulate CIGAR string.  We can do this with position specific parameters only and thus independent of any
    // context.
    void _simulateCigar(TCigarString & cigar, unsigned sampleLength);

    void _simulateQualities(TQualities & quals, TCigarString const & cigar, unsigned sampleLength);
};

// ----------------------------------------------------------------------------
// Class SequencingSimulatorFactory
// ----------------------------------------------------------------------------

// Factory for SequencingSimulator objects.op

class SequencingSimulatorFactory
{
public:
    TRng & rng;
    TRng & methRng;

    SequencingOptions const & seqOptions;
    IlluminaSequencingOptions const & illuminaOptions;
    Roche454SequencingOptions const & roche454Options;
    SangerSequencingOptions const & sangerOptions;

    SequencingSimulatorFactory(TRng & rng, TRng & methRng,
                               SequencingOptions const & seqOptions,
                               IlluminaSequencingOptions const & illuminaOptions,
                               Roche454SequencingOptions const & roche454Options,
                               SangerSequencingOptions const & sangerOptions) :
            rng(rng), methRng(methRng), seqOptions(seqOptions), illuminaOptions(illuminaOptions),
            roche454Options(roche454Options), sangerOptions(sangerOptions)
    {}

    std::unique_ptr<SequencingSimulator> make();
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function appendOrientation()
// ----------------------------------------------------------------------------

// Returns (a, b) where a is the difference in resulting read length and b is the difference in used up input/reference
// sequence length.

inline std::pair<int, int> appendOperation(TCigarString & cigar, char op)
{
    // Canceling out of events.  This can only happen if the CIGAR string is not empty.
    if (!empty(cigar) && ((back(cigar).operation == 'I' && op == 'D') ||
                          (back(cigar).operation == 'D' && op == 'I')))
    {
        if (back(cigar).count > 1)
            back(cigar).count -= 1;
        else
            eraseBack(cigar);
        return (op == 'D') ? std::make_pair(-1, 0) : std::make_pair(0, -1);
    }

    // No canceling out of events.  The read length increases by one if the operation is no deletion and one base of
    // input sequence is used up if the operation is not an insertion.
    if (!empty(cigar) && back(cigar).operation == op)
        back(cigar).count += 1;
    else
        appendValue(cigar, seqan2::CigarElement<>(op, 1));
    return std::make_pair((op != 'D'), (op != 'I'));
}

// ---------------------------------------------------------------------------
// Function _simulateSequence()
// ---------------------------------------------------------------------------

// Simulate the characters that polymorphisms turn into and inserted characters.
//
// Through the usage of ModifiedString, we will always go from the left to the right end.
template <typename TFrag>
void _simulateSequence(TRead & read, TRng & rng, TFrag const & frag,
                       TCigarString const & cigar)
{
    clear(read);

    typedef typename seqan2::Iterator<TFrag>::Type TFragIter;
    TFragIter it = begin(frag, seqan2::Standard());

    for (unsigned i = 0; i < length(cigar); ++i)
    {
        //unsigned numSimulate = 0;
        if (cigar[i].operation == 'M')
        {
            for (unsigned j = 0; j < cigar[i].count; ++j, ++it)
                appendValue(read, *it);
            continue;
        }
        else if (cigar[i].operation == 'D')
        {
            it += cigar[i].count;
            continue;
        }

        // Otherwise, we have insertions or mismatches.
        for (unsigned j = 0; j < cigar[i].count; ++j)
        {
            // Pick a value between 0 and 1.
            std::uniform_real_distribution<double> dist(0, 1);
            double x = 1.0;
            while (x == 1.0)
                x = dist(rng);
            int num = static_cast<int>(x / 0.25);

            // NOTE: We can only insert CGAT, but we can have a polymorphism to N.

            if (cigar[i].operation == 'I')
                appendValue(read, seqan2::Dna5(num));
            else
                appendValue(read, seqan2::Dna5(num + (num == ordValue(*it))));
        }

        if (cigar[i].operation == 'X')
            it += cigar[i].count;
    }
}

#endif  // #ifndef APPS_MASON2_SEQUENCING_H_
