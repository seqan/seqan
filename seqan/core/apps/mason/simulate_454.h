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
// Code specific to the 454 read simulation.
// ==========================================================================

#ifndef SIMULATE_454_H_
#define SIMULATE_454_H_

#include "simulate_454_base_calling.h"

// Maximal homopolymer length we will observe.
const unsigned MAX_HOMOPOLYMER_LEN = 40;

// ============================================================================
// Enums, Tags, Classes.
// ============================================================================

struct _LS454Reads;
typedef Tag<_LS454Reads> LS454Reads;

template <>
struct Options<LS454Reads> : public Options<Global>
{
    // Read Length Parameters.

    // Iff true, read lengths follow a uniform distribution, otherwise a
    // standard distribution will be used.
    bool readLengthIsUniform;

    // Average read length.
    double readLengthMean;

    // For standard distributed read lengths, this is the standard deviation,
    // for uniform read length the interval around the average to use for
    // picking the read lengths.
    double readLengthError;

    // Base Calling Error Model Parameters.

    // If set, $\sigma = k * \sqrt(r)$, otherwise $\sigma = k * r$ is used.
    bool sqrtInStdDev;

    // Proportionality factor for calculating standard deviation proportional
    // to sqrt(homopolymer length).
    double k;

    // Noise parameters.  We take the default values 0.23 and 0.15 from Metasim.
    
    // The mean of the lognormal distribution for the noise.
    double backgroundNoiseMean;

    // The standard deviation of the lognormal distribution for the noise.
    double backgroundNoiseStdDev;

    Options()
            : readLengthIsUniform(false),
              readLengthMean(400),
              readLengthError(40),
              sqrtInStdDev(true),
              k(0.15),
              backgroundNoiseMean(0.23),
              backgroundNoiseStdDev(0.15)
    {}
};

template <>
struct ReadSimulationInstruction<LS454Reads> : public ReadSimulationInstruction<Global> {
    // For each insertion in the edit string, this string provides the
    // nucleotide types to be inserted at this point.
    String<Dna5> insertionNucleotides;
};

template<>
struct ModelParameters<LS454Reads> : public ModelParameters<Global>
{
    ThresholdMatrix thresholdMatrix;
};

// ============================================================================
// Metafunctions.
// ============================================================================

// ============================================================================
// Functions.
// ============================================================================

template <typename TStream>
TStream & operator<<(TStream & stream, Options<LS454Reads> const & options) {
    stream << static_cast<Options<Global> >(options);
    stream << "454-options {" << std::endl
           << "  readLengthIsUniform:   " << options.readLengthIsUniform << std::endl
           << "  readLengthMean:        " << options.readLengthMean << std::endl
           << "  readLengthError:       " << options.readLengthError << std::endl
           << "  sqrtInStdDev:          " << options.sqrtInStdDev << std::endl
           << "  k:                     " << options.k << std::endl
           << "  backgroundNoiseMean:   " << options.backgroundNoiseMean << std::endl
           << "  backgroundNoiseStdDev: " << options.backgroundNoiseStdDev << std::endl
           << "}" << std::endl;
    return stream;
}

void setUpArgumentParser(ArgumentParser & parser,
                         LS454Reads const &)
{
    setUpArgumentParser(parser);
    addUsageLine(parser, "454 [\\fIOPTIONS\\fP] \\fBSEQUENCE\\fP");

    addSection(parser, "454 Read Length Parameters");

    addOption(parser, ArgParseOption("nu",  "read-length-uniform", "If set, the read lengths are simulated with a uniform distribution, with standard distribution otherwise."));
    addOption(parser, ArgParseOption("nm",  "read-length-mean", "The mean of the read lengths.", ArgParseOption::DOUBLE, "REAL"));
    setDefaultValue(parser, "read-length-mean", "400");
    addOption(parser, ArgParseOption("ne",  "read-length-error", "The standard deviation (for standard distribution) and interval length (for uniform distribution) for the read length.", ArgParseOption::DOUBLE, "REAL"));
    setDefaultValue(parser, "read-length-error", "40");

    addSection(parser, "454 Error Model Parameters");

    addOption(parser, ArgParseOption("nsq",  "no-sqrt-in-std-dev", "If set, no square root is used in error calculation."));
    addOption(parser, ArgParseOption("k",    "proportionality-factor", "Proportionality factor for calculating standard deviation proportional to sqrt(homopolymer length).", ArgParseOption::DOUBLE));
    setDefaultValue(parser, "proportionality-factor", "0.15");
    addOption(parser, ArgParseOption("bm",  "background-noise-mean", "Background noise mean.", ArgParseOption::DOUBLE));
    setDefaultValue(parser, "background-noise-mean", "0.23");
    addOption(parser, ArgParseOption("bs",  "background-noise-stddev", "Background noise std dev.", ArgParseOption::DOUBLE));
    setDefaultValue(parser, "background-noise-stddev", "0.15");
}

seqan::ArgumentParser::ParseResult
parseArgumentsAndCheckModelSpecific(Options<LS454Reads> & options,
                                    ArgumentParser & parser)
{
    options.readLengthIsUniform = isSet(parser, "read-length-uniform");
    getOptionValue(options.readLengthMean, parser, "read-length-mean");
    getOptionValue(options.readLengthError, parser, "read-length-error");

    options.sqrtInStdDev = !isSet(parser, "no-sqrt-in-std-dev");
    getOptionValue(options.k, parser, "proportionality-factor");
    getOptionValue(options.backgroundNoiseMean, parser, "background-noise-mean");
    getOptionValue(options.backgroundNoiseStdDev, parser, "background-noise-stddev");

    return seqan::ArgumentParser::PARSE_OK;
}

// For 454 reads, we do not need model specific data (yet?).
int simulateReadsSetupModelSpecificData(ModelParameters<LS454Reads> & parameters,
                                        Options<LS454Reads> const & options)
{
    setK(parameters.thresholdMatrix, options.k);
    setUseSqrt(parameters.thresholdMatrix, options.sqrtInStdDev);
    setNoiseMeanStdDev(parameters.thresholdMatrix, options.backgroundNoiseMean, options.backgroundNoiseStdDev);
    return 0;
}

template <typename TRNG>
inline
unsigned pickReadLength(TRNG & rng, Options<LS454Reads> const & options)
{
    if (options.readLengthIsUniform) {
        // Pick uniformly.
        double minLen = options.readLengthMean - options.readLengthError;
        double maxLen = options.readLengthMean + options.readLengthError;
        double len = pickRandomNumber(rng, Pdf<Uniform<double> >(minLen, maxLen));
        return static_cast<unsigned>(round(len));
    } else {
        // Pick normally distributed.
        double len = pickRandomNumber(rng, Pdf<Normal>(options.readLengthMean, options.readLengthError));
        return static_cast<unsigned>(round(len));
    }
}


template <typename TRNG, typename TContig>
void buildSimulationInstructions(ReadSimulationInstruction<LS454Reads> & inst, TRNG & rng, unsigned readLength, TContig const & contig, ModelParameters<LS454Reads> const & parameters, Options<LS454Reads> const & options)
{
    typedef Iterator<String<Dna5>, Standard>::Type TIterator;
    
    if (inst.endPos == inst.beginPos)
        return;

    //
    // Perform Flowcell Simulation.
    //
    reserve(inst.editString, readLength, Generous());
    clear(inst.editString);
    if (options.simulateQualities) {
        reserve(inst.qualities, readLength, Generous());
        clear(inst.qualities);
    }

    // Get a copy of the haplotype region we are considering.
    String<Dna5> haplotypeInfix = infix(contig, inst.beginPos, inst.endPos);
	SEQAN_ASSERT_LEQ(inst.endPos, length(contig));
    SEQAN_ASSERT_EQ(readLength, inst.endPos - inst.beginPos);

    // In the flow cell simulation, we will simulate light intensities which
    // will be stored in observedIntensities.
    String<double> observedIntensities;
    reserve(observedIntensities, 4 * readLength);
    String<Dna5> observedBases;
    // We also store the real homopolymer length.
    String<unsigned> realBaseCount;

    // Probability density function to use for the background noise.
    Pdf<LogNormal> noisePdf(options.backgroundNoiseMean, options.backgroundNoiseStdDev, MeanStdDev());

    // Initialize information about the current homopolymer length.
    unsigned homopolymerLength = 0;
    Dna homopolymerType = haplotypeInfix[0];
    while (haplotypeInfix[homopolymerLength] == homopolymerType)
        ++homopolymerLength;

    // Simulate flowcell.
    for (unsigned i = 0, j = 0; i < readLength; ++j, j = j % 4) {  // i indicates first pos of current homopolymer, j indicates flow phase
        if (ordValue(homopolymerType) == j) {
            // Simulate positive flow observation.
            double l = homopolymerLength;
            double sigma = options.k * (options.sqrtInStdDev ? sqrt(l) : l);
            double intensity = pickRandomNumber(rng, Pdf<Normal>(homopolymerLength, sigma));
            intensity += pickRandomNumber(rng, noisePdf);  // Add noise.
            appendValue(observedIntensities, intensity);
            appendValue(realBaseCount, homopolymerLength);
            // Get begin pos and length of next homopolymer.
            i += homopolymerLength;
			if (i < length(haplotypeInfix)) {
				homopolymerType = haplotypeInfix[i];
				homopolymerLength = 0;
				while (((i + homopolymerLength) < length(haplotypeInfix)) && haplotypeInfix[i + homopolymerLength] == homopolymerType)
					++homopolymerLength;
			}
        } else {
            // Simulate negative flow observation.
            //
            // Constants taken from MetaSim paper which have it from the
            // original 454 publication.
            double intensity = _max(0.0, pickRandomNumber(rng, noisePdf));
            appendValue(observedIntensities, intensity);
            appendValue(realBaseCount, 0);
        }
    }

    inst.mismatchCount = 0;
    inst.insCount = 0;
    inst.delCount = 0;
	clear(inst.insertionNucleotides);

    // Call bases, from this build the edit string and maybe qualities.  We
    // only support the "inter" base calling method which was published by
    // the MetaSim authors in the PLOS paper.
    typedef Iterator<String<double>, Standard>::Type IntensitiesIterator;
    int i = 0;  // Flow round, Dna(i % 4) gives base.
    for (IntensitiesIterator it = begin(observedIntensities); it != end(observedIntensities); ++it, ++i) {
        double threshold = getThreshold(parameters.thresholdMatrix, static_cast<unsigned>(floor(*it)), static_cast<unsigned>(ceil(*it)));
        unsigned calledBaseCount = static_cast<unsigned>(*it < threshold ? floor(*it) : ceil(*it));
        // Add any matches.
        unsigned j = 0;
        for (; j < _min(calledBaseCount, realBaseCount[i]); ++j)
            appendValue(inst.editString, ERROR_TYPE_MATCH);
        // Add insertions, if any.
        for (; j < calledBaseCount; ++j) {
            appendValue(inst.insertionNucleotides, Dna(i % 4));
            appendValue(inst.editString, ERROR_TYPE_INSERT);
            ++inst.insCount;
        }
        // Add deletions, if any.
        for (; j < realBaseCount[i]; ++j) {
            appendValue(inst.editString, ERROR_TYPE_DELETE);
            ++inst.delCount;
        }
        // Simulate qualities if configured to do so.
        if (options.simulateQualities) {
            // Compute likelihood for calling the bases, given this intensity and the Phred score from this.
            double densitySum = 0;
            for (unsigned j = 0; j <= _max(4u, 2 * MAX_HOMOPOLYMER_LEN); ++j)  // Anecdotally through plot in maple: Enough to sum up to 4 or 2 times the maximal homopolymer length.
                densitySum += dispatchDensityFunction(parameters.thresholdMatrix, j, *it);
            double x = 0;  // Probability of seeing < (j+1) bases.
            for (unsigned j = 0; j < calledBaseCount; ++j) {
                x += dispatchDensityFunction(parameters.thresholdMatrix, j, *it);
                unsigned phredScore = -static_cast<int>(10 * ::std::log10(x / densitySum));
                appendValue(inst.qualities, phredScore);
            }
        }
    }
	
    SEQAN_ASSERT_EQ(length(inst.insertionNucleotides), inst.insCount);
    if (options.simulateQualities)
        SEQAN_ASSERT_EQ(length(inst.qualities) + inst.delCount, length(inst.editString));
}

template <typename TRNG, typename TString>
void applySimulationInstructions(TString & read, TRNG & /*rng*/, ReadSimulationInstruction<LS454Reads> const & inst, Options<LS454Reads> const & options)
{
    typedef typename Value<TString>::Type TAlphabet;

    SEQAN_ASSERT_EQ(length(inst.insertionNucleotides), inst.insCount);
    if (options.simulateQualities)
        SEQAN_ASSERT_EQ(length(inst.qualities) + inst.delCount, length(inst.editString));
    SEQAN_ASSERT_EQ(length(read), length(inst.editString) - inst.insCount);
    
    TString tmp;
    reserve(tmp, length(read) + inst.insCount - inst.delCount);
    unsigned j = 0;  // Index in read
    unsigned k = 0;  // Index in inst.insertionNucleotides
    unsigned l = 0;  // Index in inst.qualities
    for (unsigned i = 0; i < length(inst.editString); ++i) {
        SEQAN_ASSERT_LEQ(j, i);

        TAlphabet c;
        switch (inst.editString[i]) {
            case ERROR_TYPE_MATCH:
                SEQAN_ASSERT_LT_MSG(j, length(read), "i = %u", i);
                appendValue(tmp, read[j]);
                if (options.simulateQualities)
                    assignQualityValue(back(tmp), inst.qualities[l++]);
                j += 1;
                break;
            case ERROR_TYPE_MISMATCH:
                SEQAN_ASSERT_FAIL("No mismatches should occur for 454 reads!");
                break;
            case ERROR_TYPE_INSERT:
                appendValue(tmp, inst.insertionNucleotides[k++]);
                if (options.simulateQualities)
                    assignQualityValue(back(tmp), inst.qualities[l++]);
                break;
            case ERROR_TYPE_DELETE:
                j += 1;
                break;
            default:
                SEQAN_ASSERT_FAIL("Invalid error type.");
        }
    }
    SEQAN_ASSERT_EQ(j, length(read));

    move(read, tmp);
}

#endif  // SIMULATE_454_H_
