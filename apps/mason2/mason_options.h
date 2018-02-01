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
// We define the options of all Mason components in a central place.
//
//
// There is one *Options struct for each simulation core component.  Each
// program has its own Mason*Options struct for storing its configuration.
//
// Each struct has a function addOptions() and getOptions() function for
// adding options to an ArgumentParser and getting the values back from the
// parser after parsing.
//
// Also, some helper functions are provided for converting booleans and enums
// to human readable strings.
// ==========================================================================

#ifndef APPS_MASON2_MASON_OPTIONS_H_
#define APPS_MASON2_MASON_OPTIONS_H_

#include <seqan/arg_parse.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/bam_io.h>
#include <seqan/gff_io.h>
#include <seqan/vcf_io.h>

#include "mason_types.h"

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class MasonSimulateGenomeOptions
// ----------------------------------------------------------------------------

// Configuration for simulating genomes.
//
// This is separate from the application options of the simulate_genome program to allow simulation from other parts of
// the program when necessary.

struct MasonSimulateGenomeOptions
{
    // A list of the lengths of the contigs to simulate.
    seqan::String<int> contigLengths;
    // The seed to use for random number generation.
    int seed;

    MasonSimulateGenomeOptions() : seed(0)
    {}
};

// ----------------------------------------------------------------------------
// Class BSSeqOptions
// ----------------------------------------------------------------------------

// Configuration for BS-Seq treatment simulation.

struct BSSeqOptions
{
    // Verbosity: 0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    int verbosity;

    // The BS simulation type to use.
    enum BSProtocol
    {
        DIRECTIONAL,
        UNDIRECTIONAL
    };

    // Whether or not BS-treatment is enabled or not.
    bool bsSimEnabled;
    // The rate that unmethylated Cs to become Ts.
    double bsConversionRate;
    // The protocol to use for the simulation.
    BSProtocol bsProtocol;

    BSSeqOptions() : verbosity(1), bsSimEnabled(false), bsConversionRate(1.0), bsProtocol(DIRECTIONAL)
    {}

    // Add options to the argument parser.
    void addOptions(seqan::ArgumentParser & parser) const;

    // Add possible text sections to the argument parser.
    void addTextSections(seqan::ArgumentParser & parser) const;

    // Get option values from the argument parser.
    void getOptionValues(seqan::ArgumentParser const & parser);

    // Print settings to out.
    void print(std::ostream & out) const;
};

// ----------------------------------------------------------------------------
// Class MethylationLevelSimulatorOptions
// ----------------------------------------------------------------------------

// Configuration for methylation level computation.

struct MethylationLevelSimulatorOptions
{
    // Verbosity: 0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    int verbosity;

    // Whether or not to simulate methylation levels.
    bool simulateMethylationLevels;
    // Median and standard deviation for picking methylation level for CpGs.
    double methMuCG, methSigmaCG;
    // Median and standard deviation for picking methylation level for CHGs.
    double methMuCHG, methSigmaCHG;
    // Median and standard deviation for picking methylation level for CHHs.
    double methMuCHH, methSigmaCHH;

    MethylationLevelSimulatorOptions() :
            verbosity(1), simulateMethylationLevels(false), methMuCG(0), methSigmaCG(0),
            methMuCHG(0), methSigmaCHG(0), methMuCHH(0), methSigmaCHH(0)
    {}

    // Add options to the argument parser.
    void addOptions(seqan::ArgumentParser & parser) const;

    // Add possible text sections to the argument parser.
    void addTextSections(seqan::ArgumentParser & parser) const;

    // Get option values from the argument parser.
    void getOptionValues(seqan::ArgumentParser const & parser);

    // Print settings to out.
    void print(std::ostream & out) const;
};


// ----------------------------------------------------------------------------
// Class MaterializerOptions
// ----------------------------------------------------------------------------

// Configuration for the contig- and haplotype-wise materialization of variants from a VCF file to a reference sequence
// FASTA file.

struct MaterializerOptions
{
    // Verbosity: 0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    int verbosity;

    // Path to reference file.  Required.
    seqan::CharString fastaFileName;
    // Path to VCF file.  No variation is applied if empty.
    seqan::CharString vcfFileName;

    // TODO(holtgrew): Add options for methylation levels FASTA input here?

    MaterializerOptions() : verbosity(1)
    {}

    // Add options to the argument parser.
    void addOptions(seqan::ArgumentParser & parser) const;

    // Add possible text sections to the argument parser.
    void addTextSections(seqan::ArgumentParser & parser) const;

    // Get option values from the argument parser.
    void getOptionValues(seqan::ArgumentParser const & parser);

    // Print settings to out.
    void print(std::ostream & out) const;
};

// ----------------------------------------------------------------------------
// Class FragmentSamplerOptions
// ----------------------------------------------------------------------------

// Configuration for the sampling of uniformly or normally distributed fragments from a sequence.

struct FragmentSamplerOptions
{
    // Verbosity: 0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    int verbosity;

    // Type for selecting the fragment size distribution: uniformly or normally distributed.
    enum FragmentSizeModel
    {
        UNIFORM,
        NORMAL
    };

    /*
    // The BS simulation type to use.
    enum BSProtocol
    {
        DIRECTIONAL,
        UNDIRECTIONAL
    };
    */

    // Lower bound on fragment size to simulate, used for requiring minimal size when simulating Illumina reads.  A
    // value of 0 indicates no bound.
    int fragSizeLowerBound;

    // Smallest fragment size if uniformly distributed.
    int minFragmentSize;
    // Maximal fragment size if uniformly distributed.
    int maxFragmentSize;

    // Mean fragment size if normally distributed.
    int meanFragmentSize;
    // Standard deviation of fragment size if normally distributed.
    int stdDevFragmentSize;

    // The model to use for the fragment size.
    FragmentSizeModel model;
    /*
    // Wether or not to simulate BS-seq.
    bool bsSimEnabled;
    // The rate that unmethylated Cs to become Ts.
    double bsConversionRate;
    // The protocol to use for the simulation.
    BSProtocol bsProtocol;
    */

    FragmentSamplerOptions() :
            verbosity(1), fragSizeLowerBound(0), minFragmentSize(0), maxFragmentSize(0), meanFragmentSize(0),
            stdDevFragmentSize(0), model(UNIFORM)
    {}

    // Add options to the argument parser.
    void addOptions(seqan::ArgumentParser & parser) const;

    // Add possible text sections to the argument parser.
    void addTextSections(seqan::ArgumentParser & parser) const;

    // Get option values from the argument parser.
    void getOptionValues(seqan::ArgumentParser const & parser);

    // Print settings to out.
    void print(std::ostream & out) const;
};

// ----------------------------------------------------------------------------
// Class SequencingOptions
// ----------------------------------------------------------------------------

// Configuration for generic read simulation.
//
// Generic configuration for read simulation, such as mate orientation, enabling or disabling quality and paired
// simulation, and selecting specific strands.

struct SequencingOptions
{
    // Verbosity: 0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    int verbosity;

    // Enum for selecting the relative read orientation.
    enum MateOrientation
    {
        FORWARD_REVERSE = 0,  // R1 --> <-- R2
        REVERSE_FORWARD = 1,  // R1 <-- --> R2
        FORWARD_FORWARD = 2,  // R1 --> --> R2
        FORWARD_FORWARD2 = 3  // R2 --> --> R1
    };

    // Enum for selecting forward and/or reverse strand.
    enum SourceStrands
    {
        BOTH,
        FORWARD,
        REVERSE
    };

    // Enum for selecting the sequencing technology to simulate.
    enum SequencingTechnology
    {
        ILLUMINA,
        ROCHE_454,
        SANGER
    };

    // Options for BS-Seq simulation.
    BSSeqOptions bsSeqOptions;

    // Prefix to give all reads.
    seqan::CharString readNamePrefix;
    // Whether or not to simulate qualities.
    bool simulateQualities;
    // Whether or not to simulate mate pairs.
    bool simulateMatePairs;
    // Whether or not to gather read information.
    bool embedReadInfo;
    // Mate orientation.
    MateOrientation mateOrientation;
    // Whether to simulate from forward/reverse strand or both.
    SourceStrands strands;
    // The sequencing technology.
    SequencingTechnology sequencingTechnology;

    SequencingOptions() :
            verbosity(1), simulateQualities(false), simulateMatePairs(false), embedReadInfo(false),
            mateOrientation(FORWARD_REVERSE), strands(BOTH), sequencingTechnology(ILLUMINA)
    {}

    // Add options to the argument parser.
    void addOptions(seqan::ArgumentParser & parser) const;

    // Add possible text sections to the argument parser.
    void addTextSections(seqan::ArgumentParser & parser) const;

    // Get option values from the argument parser.
    void getOptionValues(seqan::ArgumentParser const & parser);

    // Print settings to out.
    void print(std::ostream & out) const;
};

// ----------------------------------------------------------------------------
// Class IlluminaSequencingOptions
// ----------------------------------------------------------------------------

// Configuration for Illumina read simulation.
//
// Configuration specific to the Illumina sequencing model.

// TODO(holtgrew): Allow for giving a FASTQ file as the input for qualities and N-patterns.

struct IlluminaSequencingOptions
{
    // Verbosity: 0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    int verbosity;

    // Length of the reads to simulate.
    unsigned readLength;

    // Path to file with positional error probabilities.
    seqan::CharString probabilityMismatchFile;

    // The default orientation for Illumina paired-end reads.
    SequencingOptions::MateOrientation defaultOrientation;

    // -----------------------------------------------------------------------
    // Base Calling Error Model Parameters.
    // -----------------------------------------------------------------------

    // Probability of an insertion.
    double probabilityInsert;
    // Probability of a deletion.
    double probabilityDelete;

    // TODO(holtgrew): Load probabilities from file or use FASTQ file as quality/N template.

    // Scale factor to apply to mismatch probabilities,
    double probabilityMismatchScale;

    // Probability of a mismatch (single-base polymorphism).
    double probabilityMismatch;
    // Probability for a mismatch in the first base.
    double probabilityMismatchBegin;
    // Probability for a mismatch in the last base.
    double probabilityMismatchEnd;
    // Relative position in the read between 0 and 1 where the steeper curve begins.
    double positionRaise;

    // // If set then no Ns will be introduced into the read.
    // bool illuminaNoN;

    // Paths to left/right template FASTQ files.  The qualities will be used to compute positional qualities, patterns
    // of Ns will be applied to the simulated reads.  If set, this will be used instead of the built-in model.
    seqan::CharString leftTemplateFastq, rightTemplateFastq;

    // -----------------------------------------------------------------------
    // Base Calling Quality Model Parameters.
    // -----------------------------------------------------------------------

    // Mean quality for non-mismatches at the first base.
    double meanQualityBegin;
    // Mean quality for non-mismatches at last base.
    double meanQualityEnd;
    // Standard deviation quality for non-mismatches at the first base.
    double stdDevQualityBegin;
    // Standard deviation quality for non-mismatches at the last base.
    double stdDevQualityEnd;

    // Mean quality for mismatches at the first base.
    double meanMismatchQualityBegin;
    // Mean quality for mismatches at last base.
    double meanMismatchQualityEnd;
    // Standard deviation quality for mismatches at the first base.
    double stdDevMismatchQualityBegin;
    // Standard deviation quality for mismatches at the last base.
    double stdDevMismatchQualityEnd;

    IlluminaSequencingOptions() :
            verbosity(1),
            readLength(0),
            defaultOrientation(SequencingOptions::FORWARD_REVERSE),
            // Base Calling Error Model Parameters
            probabilityInsert(0.001),
            probabilityDelete(0.001),
            probabilityMismatchScale(1.0),
            probabilityMismatch(0.004),
            probabilityMismatchBegin(0.002),
            probabilityMismatchEnd(0.012),
            positionRaise(0.66),
            // illuminaNoN(false),
            // Base Calling Quality Model Parameters
            meanQualityBegin(40),
            meanQualityEnd(39.5),
            stdDevQualityBegin(0.05),
            stdDevQualityEnd(10),
            meanMismatchQualityBegin(39.5),
            meanMismatchQualityEnd(30),
            stdDevMismatchQualityBegin(3),
            stdDevMismatchQualityEnd(15)
    {}

    // Add options to the argument parser.
    void addOptions(seqan::ArgumentParser & parser) const;

    // Add possible text sections to the argument parser.
    void addTextSections(seqan::ArgumentParser & parser) const;

    // Get option values from the argument parser.
    void getOptionValues(seqan::ArgumentParser const & parser);

    // Print settings to out.
    void print(std::ostream & out) const;
};

// ----------------------------------------------------------------------------
// Class SangerSequencingOptions
// ----------------------------------------------------------------------------

// Configuration for Sanger read simulation.
//
// Configuration specific to the Sanger sequencing model.

struct SangerSequencingOptions
{
    // Verbosity: 0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    int verbosity;

    // The default orientation for Illumina paired-end reads.
    SequencingOptions::MateOrientation defaultOrientation;

    // Read Length Parameters.

    // Iff true, read lengths follow a uniform distribution, otherwise a
    // standard distribution will be used.
    bool readLengthIsUniform;

    // Average read length for normal distribution.
    double readLengthMean;

    // For standard distributed read lengths, this is the standard deviation.
    double readLengthError;

    // Minimal and maximal read lenght in case of uniform distribution.
    double readLengthMin;
    double readLengthMax;

    // Base Calling Error Model Parameters.

    // Mismatch probability ramp.
    double probabilityMismatchBegin;
    double probabilityMismatchEnd;

    // Insert probability ramp.
    double probabilityInsertBegin;
    double probabilityInsertEnd;

    // Delete probability ramp.
    double probabilityDeleteBegin;
    double probabilityDeleteEnd;

    // Quality mean and standard deviation ramp specification for matches.
    double qualityMatchStartMean;
    double qualityMatchEndMean;
    double qualityMatchStartStdDev;
    double qualityMatchEndStdDev;

    // Quality mean and standard deviation ramp specification for errors.
    double qualityErrorStartMean;
    double qualityErrorEndMean;
    double qualityErrorStartStdDev;
    double qualityErrorEndStdDev;

    SangerSequencingOptions() :
            verbosity(1),
            defaultOrientation(SequencingOptions::FORWARD_REVERSE),
            readLengthIsUniform(false),
            readLengthMean(400),
            readLengthError(40),
            readLengthMin(100),
            readLengthMax(200),
            probabilityMismatchBegin(0.005),
            probabilityMismatchEnd(0.01),
            probabilityInsertBegin(0.0025),
            probabilityInsertEnd(0.005),
            probabilityDeleteBegin(0.0025),
            probabilityDeleteEnd(0.005),
            qualityMatchStartMean(40),
            qualityMatchEndMean(39),
            qualityMatchStartStdDev(0.1),
            qualityMatchEndStdDev(2),
            qualityErrorStartMean(30),
            qualityErrorEndMean(20),
            qualityErrorStartStdDev(2),
            qualityErrorEndStdDev(5)
    {}

    // Add options to the argument parser.
    void addOptions(seqan::ArgumentParser & parser) const;

    // Add possible text sections to the argument parser.
    void addTextSections(seqan::ArgumentParser & parser) const;

    // Get option values from the argument parser.
    void getOptionValues(seqan::ArgumentParser const & parser);

    // Print settings to out.
    void print(std::ostream & out) const;
};

// ----------------------------------------------------------------------------
// Class Roche454SequencingOptions
// ----------------------------------------------------------------------------

// Configuration for Roche 454 read simulation.
//
// Configuration specific to the Roche 454 sequencing model.

struct Roche454SequencingOptions
{
    // Verbosity: 0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    int verbosity;

    // The default orientation for Roche 454 paired-end reads.
    SequencingOptions::MateOrientation defaultOrientation;

    // Enum for selecting the read length model.
    enum ReadLengthModel
    {
        UNIFORM,
        NORMAL
    };

    // Read Length Parameters

    // The model for the read length.
    ReadLengthModel lengthModel;
    // The minimal and maximal read length if uniformly distributed.
    int minReadLength, maxReadLength;
    // The mean read length if normally distributed.
    double meanReadLength;
    // The read length standard deviation if normally distributed.
    double stdDevReadLength;

    // Base-Calling Error Model Parameters

    // If true then $\sigma = k * \sqrt(r)$, otherwise $\sigma = k * r$ is used.
    bool sqrtInStdDev;

    // Proportionality factor for calculating standard deviation proportional
    // to sqrt(homopolymer length).
    double k;

    // Noise parameters.  We take the default values 0.23 and 0.15 from Metasim.

    // The mean of the lognormal distribution for the noise.
    double backgroundNoiseMean;

    // The standard deviation of the lognormal distribution for the noise.
    double backgroundNoiseStdDev;

    Roche454SequencingOptions() :
            verbosity(1),
            defaultOrientation(SequencingOptions::FORWARD_FORWARD2),
            lengthModel(UNIFORM),
            minReadLength(0),
            maxReadLength(0),
            meanReadLength(0),
            stdDevReadLength(0),
            sqrtInStdDev(true),
            k(0),
            backgroundNoiseMean(0),
            backgroundNoiseStdDev(0)
    {}

    // Add options to the argument parser.
    void addOptions(seqan::ArgumentParser & parser) const;

    // Add possible text sections to the argument parser.
    void addTextSections(seqan::ArgumentParser & parser) const;

    // Get option values from the argument parser.
    void getOptionValues(seqan::ArgumentParser const & parser);

    // Print settings to out.
    void print(std::ostream & out) const;
};

// ----------------------------------------------------------------------------
// Class MasonSimulatorOptions
// ----------------------------------------------------------------------------

// Configuration for the program mason_simulator.
//
// mason_simulator reads in a VCF file and a genome file and then materializes the haplotype using the core components
// of mason_materializer.  Then, fragments are sampled from the materialized file and reads are sequencing is simulated
// to yield the reads and alignmetns.  For this, the core components of mason_fragments and mason_seq_fragments are
// used.
//
// Thus, we store only a few more settings besides re-using the configuration.

struct MasonSimulatorOptions
{
    // Verbosity: 0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    int verbosity;
    // The seed for the random number generator.
    int seed;
    // The seed for the random number generator for methylation simulation.
    int methSeed;
    // The spacing of the see when using multi-threading.  Thread i (beginning with 0) will get (seed + i) as its
    // initial seed.  While not crytographically safe, this should be OK for read simulation.
    int seedSpacing;
    // The number of threads to use for the simulation.
    int numThreads;
    // Number of reads/pairs to simulate in one chunk
    int chunkSize;

    // Number of reads/pairs to simulate.
    int numFragments;

    // Whether to force single-end simulation although !empty(outFileNameRight).
    bool forceSingleEnd;

    // Path to input methylation FASTA file.
    seqan::CharString methFastaInFile;
    // Path to output sequence files for left (and single end) and right reads.
    seqan::CharString outFileNameLeft, outFileNameRight;
    // Path to output SAM file.
    seqan::CharString outFileNameSam;

    // Configuration for the reading of the reference and application of the variants from the VCF file.
    MaterializerOptions matOptions;
    // Configuration for the methylation simulation.  Required for repairing methylation levels after variation.
    MethylationLevelSimulatorOptions methOptions;

    // Configuration for sampling references.
    FragmentSamplerOptions fragSamplerOptions;

    // Generic sequencing configuration.
    SequencingOptions seqOptions;
    // Configuration of the Illumina read simulation.
    IlluminaSequencingOptions illuminaOptions;
    // Configuration of the Sanger read simulation.
    SangerSequencingOptions sangerOptions;
    // Configuration of the Roche 454 read simulation.
    Roche454SequencingOptions rocheOptions;

    MasonSimulatorOptions() :
            verbosity(1), seed(0), methSeed(0), seedSpacing(2048), numThreads(1), chunkSize(64*1024), numFragments(0),
            forceSingleEnd(false)
    {}

    // Add options to the argument parser.  Calls addOptions() on the nested *Options objects.
    void addOptions(seqan::ArgumentParser & parser) const;

    // Add possible text sections to the argument parser, calls addTextSections() on *Options objects.
    void addTextSections(seqan::ArgumentParser & parser) const;

    // Get option values from the argument parser.  Calls getOptionValues() on the nested *Option objects.
    void getOptionValues(seqan::ArgumentParser const & parser);

    // Print settings to out.
    void print(std::ostream & out) const;
};

// ----------------------------------------------------------------------------
// Class MasonMaterializerOptions
// ----------------------------------------------------------------------------

// Configuration for the program mason_materializer.

struct MasonMaterializerOptions
{
    // Verbosity: 0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    int verbosity;

    // Seed to use in RNG.
    int seed;
    // Seed to use in RNG for methylation simulation.
    int methSeed;

    // Options for the materializer.
    MaterializerOptions matOptions;
    // Options for the methylation simulation.
    MethylationLevelSimulatorOptions methOptions;

    // Path to output file.
    seqan::CharString outputFileName;
    // Path to TSV file to write the resulting breakpoints in variant genomes to.
    seqan::CharString outputBreakpointFile;
    // Separator between contig names and haplotype number.
    seqan::CharString haplotypeNameSep;
    // FASTA file to load the methylation levels from.
    seqan::CharString methFastaInFile;
    // FASTA file to write the methylation levels to.
    seqan::CharString methFastaOutFile;

    MasonMaterializerOptions() : verbosity(1), seed(0), methSeed(0)
    {}

    // Add options to the argument parser.
    void addOptions(seqan::ArgumentParser & parser) const;

    // Add possible text sections to the argument parser.
    void addTextSections(seqan::ArgumentParser & parser) const;

    // Get option values from the argument parser.
    void getOptionValues(seqan::ArgumentParser const & parser);

    // Print settings to out.
    void print(std::ostream & out) const;
};

// ----------------------------------------------------------------------------
// Class MasonSplicingOptions
// ----------------------------------------------------------------------------

// Configuration for the program mason_materializer.

struct MasonSplicingOptions
{
    // Verbosity: 0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    int verbosity;

    // Seed to use in RNG.
    int seed;

    // Options for the materializer.
    MaterializerOptions matOptions;

    // Path to input GFF/GTF file.
    seqan::CharString inputGffFile;
    // Type of the annotations to splice.
    seqan::CharString gffType;
    // Name of the group-by column.
    seqan::CharString gffGroupBy;

    // Path to output file.
    seqan::CharString outputFileName;
    // Separator between contig names and haplotype number.
    seqan::CharString haplotypeNameSep;

    MasonSplicingOptions() : verbosity(1), seed(0)
    {}

    // Add options to the argument parser.
    void addOptions(seqan::ArgumentParser & parser) const;

    // Add possible text sections to the argument parser.
    void addTextSections(seqan::ArgumentParser & parser) const;

    // Get option values from the argument parser.
    void getOptionValues(seqan::ArgumentParser const & parser);

    // Print settings to out.
    void print(std::ostream & out) const;
};

// ----------------------------------------------------------------------------
// Class MasonFragmentSequencingOptions
// ----------------------------------------------------------------------------

// Configuration for the program mason_frag_simulator.

struct MasonFragmentSequencingOptions
{
    // Verbosity: 0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    int verbosity;

    // The seed to use for random number generation.
    int seed;

    // Whether to force single-end simulation although !empty(outFileNameRight).
    bool forceSingleEnd;

    // Path to input file.
    seqan::CharString inputFileName;

    // Path to output sequence files for left (and single end) and right reads.
    seqan::CharString outFileNameLeft, outFileNameRight;

    // Generic sequencing configuration.
    SequencingOptions seqOptions;
    // Configuration of the Illumina read simulation.
    IlluminaSequencingOptions illuminaOptions;
    // Configuration of the Sanger read simulation.
    SangerSequencingOptions sangerOptions;
    // Configuration of the Roche 454 read simulation.
    Roche454SequencingOptions rocheOptions;

    MasonFragmentSequencingOptions() : verbosity(1), seed(0)
    {}

    // Add options to the argument parser.  Calls addOptions() on the nested *Options objects.
    void addOptions(seqan::ArgumentParser & parser) const;

    // Add possible text sections to the argument parser, calls addTextSections() on *Options objects.
    void addTextSections(seqan::ArgumentParser & parser) const;

    // Get option values from the argument parser.  Calls getOptionValues() on the nested *Option objects.
    void getOptionValues(seqan::ArgumentParser const & parser);

    // Print settings to out.
    void print(std::ostream & out) const;
};

// --------------------------------------------------------------------------
// Class MasonMethylationOptions
// --------------------------------------------------------------------------

// This struct stores the options from the command line.
//
// You might want to rename this to reflect the name of your app.

struct MasonMethylationOptions
{
    // Verbosity level.  0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    int verbosity;

    // The seed to use for the RNG.
    int seed;

    // Methylation simulation options.
    MethylationLevelSimulatorOptions methOptions;

    // FASTA file to import.
    seqan::CharString fastaInFile;
    // FASTA file to write the methylation levels to.
    seqan::CharString methFastaOutFile;

    MasonMethylationOptions() : verbosity(1), seed(0)
    {}

    // Add options to the argument parser.  Calls addOptions() on the nested *Options objects.
    void addOptions(seqan::ArgumentParser & parser) const;

    // Add possible text sections to the argument parser, calls addTextSections() on *Options objects.
    void addTextSections(seqan::ArgumentParser & parser) const;

    // Get option values from the argument parser.  Calls getOptionValues() on the nested *Option objects.
    void getOptionValues(seqan::ArgumentParser const & parser);

    // Print settings to out.
    void print(std::ostream & out) const;
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

char const * getYesNoStr(bool b);
char const * getVerbosityStr(int verbosity);
char const * getFragmentSizeModelStr(FragmentSamplerOptions::FragmentSizeModel model);
char const * getMateOrientationStr(SequencingOptions::MateOrientation orientation);
char const * getSourceStrandsStr(SequencingOptions::SourceStrands strands);
char const * getSequencingTechnologyStr(SequencingOptions::SequencingTechnology technology);
char const * getFragmentSizeModelStr(Roche454SequencingOptions::ReadLengthModel model);

// ----------------------------------------------------------------------------
// Function setDateAndVersion()
// ----------------------------------------------------------------------------

inline
void setDateAndVersion(seqan::ArgumentParser & parser)
{
#ifdef SEQAN_APP_VERSION
    #ifdef SEQAN_REVISION
        setVersion(parser, SEQAN_APP_VERSION " [" SEQAN_REVISION "]");
    #else
        setVersion(parser, SEQAN_APP_VERSION);
    #endif
#endif
#ifdef SEQAN_DATE
    setDate(parser, SEQAN_DATE);
#endif
}

#endif  // #ifndef APPS_MASON2_MASON_OPTIONS_H_
