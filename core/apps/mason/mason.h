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
// Globally shared code for the Mason read simulator.
// ==========================================================================

#ifndef MASON_H_
#define MASON_H_

#include <numeric>
#include <map>

#include <seqan/random.h>
#include <seqan/sequence_journaled.h>
#include <seqan/parallel.h>
#include <seqan/arg_parse.h>

#include "store_config.h"
#include "util.h"

using namespace seqan;

// ============================================================================
// Enums, Tags, Classes.
// ============================================================================

// Enum describing the read type to be simulated.
enum ReadsType
{
    READS_TYPE_ILLUMINA,
    READS_TYPE_454,
    READS_TYPE_SANGER
};

// Naming scheme of paired-end reads.
enum ReadNaming
{
    READ_NAMING_NOSUFFIX,
    READ_NAMING_SLASH_SUFFIX0, // /0, /1
    READ_NAMING_SLASH_SUFFIX1  // /1, /2
};

// Tag for global options.
typedef void Global;

// Class options.  We will use template inheritance for specializing
// the options class and C++ inheritance to prevent redundant code.
template <typename TTag>
struct Options;

// Program-wide options.
template <>
struct Options<Global>
{
    // true iff verbosity is enabled.
    bool verbose;
    // true if very verbose is enabled.
    bool veryVerbose;

    // The type of the reads to be simulated.
    ReadsType readsType;

    // Basic Read Simulation Parameters.

    // If set to true then reads are also sampled from regions in the
    // genome that have Ns in them.
    bool allowNFromGenome;
    // Seed to use for the random number generator.
    unsigned seed;
    // Number of reads (pairs) to simulate.
    unsigned numReads;
    // true iff a random reference sequence is to be used.
    bool useRandomSequence;
    // Probability for A, C, G in the random simulated sequence.
    double sourceProbabilityA;
    double sourceProbabilityC;
    double sourceProbabilityG;
    // Length of random sequence to be simulated.
    unsigned randomSourceLength;
    // true iff only the forward strand is to be simulated.
    bool onlyForward;
    // true iff only the reverse strand is to be simulated.
    bool onlyReverse;
    // The output file.  Defaults to REFERENCE-FILE.reads.fastq,
    // possibly with a ".1" or ".2" before the ".fastq" if mate pairs
    // are simulated.
    CharString outputFile;
    // Path to the Sam file to generate.  Defaults to fastq file name
    // with suffix ".sam"
    CharString samFile;
    // Path to the VCF file to generate the haplotype from
    CharString vcfFile;
    // true iff qualities are to be simulated.
    bool simulateQualities;
    // true iff additional information is to be included in the reads file.
    bool includeReadInformation;
    // Path to file with sample counts from each contig.  The file is a
    // TSV file containing pairs mapping the contig id (non-whitespace
    // sequence in FASTA header) to a read count.  Overrides parameter
    // numReads if given.
    CharString sampleCountsFilename;

    // Mate-Pair Related Options.

    // true iff generating mate pairs is enabled.
    bool generateMatePairs;
    // true iff mate pair library sizes are to be uniformly distributed,
    // otherwise standard distribution is used.
    bool libraryLengthIsUniform;
    // Mate-pair library mean length.
    double libraryLengthMean;
    // Mate-pair library length error.  Standard deviation for normally
    // distributed library lengths, interval length around mean for uniform
    // distribution.
    double libraryLengthError;
    // Set how to name the reads in the FASTA/FASTQ files.
    ReadNaming readNaming;
    // The prefix for read names.  The default is the file name.
    CharString readNamePrefix;

    // Haplotype parameters.

    // Number of haplotypes to generated.  All are generated from the input genome.
    unsigned numHaplotypes;
    // SNP rate.
    double haplotypeSnpRate;
    // Indel rate.
    double haplotypeIndelRate;
    // Smallest number of indels.
    unsigned haplotypeIndelRangeMin;
    // Largest number of indels.
    unsigned haplotypeIndelRangeMax;
    // If true then no Ns are substituted or inserted into the haplotype.
    bool haplotypeNoN;

    Options()
            : verbose(false),
              veryVerbose(false),
              allowNFromGenome(false),
              seed(0),
              numReads(1000),
              useRandomSequence(false),
              sourceProbabilityA(0.25),
              sourceProbabilityC(0.25),
              sourceProbabilityG(0.25),
              randomSourceLength(1000*1000),
              onlyForward(false),
              onlyReverse(false),
              outputFile(""),
              samFile(""),
              vcfFile(""),
              simulateQualities(false),
              includeReadInformation(false),
              generateMatePairs(false),
              libraryLengthIsUniform(false),
              libraryLengthMean(1000),
              libraryLengthError(100),
              readNaming(READ_NAMING_NOSUFFIX),
              numHaplotypes(1),
              haplotypeSnpRate(0.001),
              haplotypeIndelRate(0.001),
              haplotypeIndelRangeMin(1),
              haplotypeIndelRangeMax(6),
              haplotypeNoN(false)
    {}
};

// Use this container for model specific parameters generated before the actual
// simulation.  Setup in void simulateReadsSetupModelSpecificData(...).
template <typename TTag>
struct ModelParameters;

// Global model parameters.
template <>
struct ModelParameters<Global>
{
    // If non-empty, sampleCounts[i] gives the number of reads to
    // sample from contig i (the i-th sequence in the FASTA input
    // file).
    String<size_t> sampleCounts;
};

// Enum describing the type of an error.
enum ErrorType
{
    ERROR_TYPE_MATCH    = 0,
    ERROR_TYPE_MISMATCH = 1,
    ERROR_TYPE_INSERT   = 2,
    ERROR_TYPE_DELETE   = 3
};

// Class for storing the read simulation instructions.  Will be
// specialized for each technology to simulate.
template <typename TReadTypeTag>
struct ReadSimulationInstruction;

// Read simulation instructions used by all simulated technologies.
template <>
struct ReadSimulationInstruction<Global>
{
    // Index of the contig to sample the read from.
    unsigned contigId;
    // Begin and end position of the infix to sample the read from.
    size_t beginPos;
    size_t endPos;
    // Index of the haplotype to sample the read from.
    unsigned char haplotype;
    // Whether or not the read is to be sampled from the forward strand.
    bool isForward;
    // Number of characters added/removed to the string by indels.
    unsigned short mismatchCount;
    unsigned short delCount;
    unsigned short insCount;
    // Edit string of the read.
    String<ErrorType> editString;
    // String of qualities for the bases written out.
    String<int> qualities;

    ReadSimulationInstruction() : mismatchCount(0), delCount(0), insCount(0) {}
};

// Class to represent variations on a contig
struct Snp
{
    unsigned int    virtualPos;
    unsigned short  length;
//    unsigned char   contigId; // SNPs are grouped by contigId, hence superfluous
    unsigned char   type;
};

// Functor to sort variations
struct SnpLess
{
    inline bool operator() (Snp const &a, Snp const &b) const
    {
//        if (a.contigId < b.contigId) return true;
//        if (a.contigId > b.contigId) return false;
        if (a.virtualPos < b.virtualPos) return true;
        if (a.virtualPos > b.virtualPos) return false;
        return false;
    }
};



// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// Prints the global options to stream.
template <typename TStream>
TStream & operator<<(TStream & stream, Options<Global> const & options) {
    char const * READ_NAMINGS[] = {"no-suffix", "suffix-zero-based", "suffix-one-based"};
    stream << "global-options {" << std::endl
           << "  allowNFromGenome:       " << (options.allowNFromGenome ? "true" : "false") << std::endl
           << "  seed:                   " << options.seed << std::endl
           << "  numReads:               " << options.numReads << std::endl
           << "  useRandomSequence:      " << (options.useRandomSequence ? "true" : "false") << std::endl
           << "  sourceProbabilityA:     " << options.sourceProbabilityA << std::endl
           << "  sourceProbabilityC:     " << options.sourceProbabilityC << std::endl
           << "  sourceProbabilityG:     " << options.sourceProbabilityG << std::endl
           << "  useRandomSequence:      " << (options.useRandomSequence ? "true" : "false") << std::endl
           << "  randomSourceLength:     " << options.randomSourceLength << std::endl
           << "  onlyForward:            " << (options.onlyForward ? "true" : "false") << std::endl
           << "  onlyReverse:            " << (options.onlyReverse ? "true" : "false") << std::endl
           << "  outputFile:             \"" << options.outputFile << "\"" << std::endl
           << "  samFile:                \"" << options.samFile << "\"" << std::endl
           << "  vcfFile:                \"" << options.vcfFile << "\"" << std::endl
           << "  simulateQualities:      " << (options.simulateQualities ? "true" : "false") << std::endl
           << "  includeReadInformation: " << options.includeReadInformation << std::endl
           << "  generateMatePairs:      " << (options.generateMatePairs ? "true" : "false") << std::endl
           << "  libraryLengthMean:      " << options.libraryLengthMean << std::endl
           << "  libraryLengthError:     " << options.libraryLengthError << std::endl
           << "  libraryLengthIsUniform: " << options.libraryLengthIsUniform << std::endl
           << "  readNaming:             " << READ_NAMINGS[(int)options.readNaming] << std::endl
           << "  numHaplotypes:          " << options.numHaplotypes << std::endl
           << "  haplotypeSnpRate:       " << options.haplotypeSnpRate << std::endl
           << "  haplotypeIndelRate:     " << options.haplotypeIndelRate << std::endl
           << "  haplotypeIndelRangeMin: " << options.haplotypeIndelRangeMin << std::endl
           << "  haplotypeIndelRangeMax: " << options.haplotypeIndelRangeMax << std::endl
           << "  haplotypeNoN:           " << options.haplotypeNoN << std::endl
           << "  sampleCountsFilename:   " << options.sampleCountsFilename << std::endl
           << "  readNamePrefix          " << options.readNamePrefix << std::endl
           << "}" << std::endl;
    return stream;
}

// Stream operator for simulation instructions.
template <typename TStream>
TStream & operator<<(TStream & stream, ReadSimulationInstruction<Global> const & inst) {
    stream << "(haplotype=" << inst.haplotype << ", contigId=" << inst.contigId << ", isForward=" << inst.isForward << ", beginPos=" << inst.beginPos << ", endPos=" << inst.endPos << ", insCount=" << inst.insCount << ", delCount=" << inst.delCount << ", editString=";
    for (unsigned i = 0; i < length(inst.editString); ++i) {
        stream << "MEID"[inst.editString[i]];
    }
    stream << ", qualities=[";
    for (unsigned i = 0; i < length(inst.qualities); ++i) {
        if (i != 0)
            stream << ", ";
        stream << inst.qualities[i];
    }
    stream << "])";
    return stream;
}

// Initialize the command line parser for the global options.
void setUpArgumentParser(ArgumentParser & parser)
{
    setShortDescription(parser, "A Read Simulator");
    setVersion(parser, "0.1.1");
    setDate(parser, "July 27, 2012");
    setCategory(parser, "Simulation");

    addDescription(parser, "Use 'random.fasta' for the SEQUENCE file name to generate it randomly.");

    addSection(parser, "Main Options");
    
    addOption(parser, ArgParseOption("aNg",  "allow-N-from-genome", "Allow N from genome.  Default: false."));
    addOption(parser, ArgParseOption("s", "seed", "The seed for Rng.", ArgParseOption::INTEGER, "INT"));
    setDefaultValue(parser, "seed", "0");
    addOption(parser, ArgParseOption("N",  "num-reads", "Number of reads (or mate pairs) to simulate.", ArgParseOption::INTEGER));
    setDefaultValue(parser, "num-reads", "1000");
    addOption(parser, ArgParseOption("sn", "source-length", "Length of random source sequence.", ArgParseOption::INTEGER));
    setDefaultValue(parser, "source-length", "1000000");
    addOption(parser, ArgParseOption("spA", "source-probability-A", "Propabilibty for A in randomly generated sequence.", ArgParseOption::DOUBLE));
    setDefaultValue(parser, "source-probability-A", "0.25");
    addOption(parser, ArgParseOption("spC", "source-probability-C", "Propabilibty for C in randomly generated sequence.", ArgParseOption::DOUBLE));
    setDefaultValue(parser, "source-probability-C", "0.25");
    addOption(parser, ArgParseOption("spG", "source-probability-G", "Propabilibty for G in randomly generated sequence.", ArgParseOption::DOUBLE));
    setDefaultValue(parser, "source-probability-G", "0.25");
    addOption(parser, ArgParseOption("snN", "source-no-N", "If set then no Ns are generated in the random source sequence."));
    addOption(parser, ArgParseOption("f",  "forward-only", "Simulate from forward strand only."));
    addOption(parser, ArgParseOption("r",  "reverse-only", "Simulate from reverse strand only."));
    addOption(parser, ArgParseOption("o",  "output-file", "Write results to PARAM.fasta file instead of SEQUENCE.reads.fasta.", ArgParseOption::OUTPUTFILE));
    setValidValues(parser, "output-file", "fasta fa fastq fq");
    addOption(parser, ArgParseOption("sq", "simulate-qualities", "Simulate qualities, generate FASTQ instead of FASTA."));
    addOption(parser, ArgParseOption("i", "include-read-information", "Include additional read information in reads file."));
    addOption(parser, ArgParseOption("v", "verbose", "Verbosity mode."));
    addOption(parser, ArgParseOption("vv", "very-verbose", "High verbosity mode, implies verbosity mode."));
    addOption(parser, ArgParseOption("scf", "sample-counts-file", "Path to TSV file that maps contig ids to read counts.", ArgParseOption::STRING));
    addOption(parser, ArgParseOption("vcf", "vcf-allele-file", "Path to VCF file to synthesize the haplotype from.", ArgParseOption::STRING));

    addSection(parser, "Mate-Pair Options");

    addOption(parser, ArgParseOption("ll", "library-length-mean", "Mate-pair mean library length.", ArgParseOption::DOUBLE));
    setDefaultValue(parser, "library-length-mean", "1000");
    addOption(parser, ArgParseOption("le", "library-length-error", "Mate-pair library tolerance.", ArgParseOption::DOUBLE));
    setDefaultValue(parser, "library-length-error", "100");
    addOption(parser, ArgParseOption("lu", "library-length-uniform", "Mate-pair library length is uniform."));
    addOption(parser, ArgParseOption("mp", "mate-pairs", "Enable mate pair simulation."));
    addOption(parser, ArgParseOption("rn", "read-naming", "Read naming scheme in FASTQ/FASTA files. See section on read naming below.", ArgParseOption::INTEGER));
    setMinValue(parser, "read-naming", "0");
    setMaxValue(parser, "read-naming", "2");
    setDefaultValue(parser, "read-naming", "0");
    addOption(parser, ArgParseOption("rnp", "read-name-prefix", "Read name prefix. Default is output file name.", ArgParseOption::STRING));

    addSection(parser, "Haplotype Options");

    addOption(parser, ArgParseOption("hn", "num-haplotypes", "Number of haplotypes to simulate.", ArgParseOption::INTEGER));
    setDefaultValue(parser, "num-haplotypes", "1");
    addOption(parser, ArgParseOption("hs", "haplotype-snp-rate", "Haplotype SNP rate.", ArgParseOption::DOUBLE));
    setDefaultValue(parser, "haplotype-snp-rate", "0.001");
    addOption(parser, ArgParseOption("hi", "haplotype-indel-rate", "Haplotype indel rate.", ArgParseOption::DOUBLE));
    setDefaultValue(parser, "haplotype-indel-rate", "0.001");
    addOption(parser, ArgParseOption("hm", "haplotype-indel-range-min", "Haplotype indel size min.", ArgParseOption::INTEGER));
    setDefaultValue(parser, "haplotype-indel-range-min", "1");
    addOption(parser, ArgParseOption("hM", "haplotype-indel-range-max", "Haplotype indel size max.", ArgParseOption::INTEGER));
    setDefaultValue(parser, "haplotype-indel-range-max", "6");
    addOption(parser, ArgParseOption("hnN", "haplotype-no-N", "Do not allow Ns to be substituted or inserted into N.  Default: Is allowed."));

    addTextSection(parser, "Read Naming Scheme");
	addText(parser, "Note that the suffixes will not appear in the SAM file.  In the SAM format, the FLAG field is used to mark the leftmost/rightmost fragment/read.");
	addListItem(parser, "0", "No suffix, left-end and right-end read will be named 'example.fastq.0000', for example.");
	addListItem(parser, "1", "Add zero-based slash-suffix, i.e. names will end on '/0' and '/1'");
	addListItem(parser, "2", "Add one-based slash-suffix, i.e. names will end on '/1' and '/2'");

    // Need reads type and {SEQUENCE.fasta, random.fasta}.
    addArgument(parser, ArgParseArgument(ArgParseArgument::STRING));
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE));
    setValidValues(parser, 1, "fasta fa");
}

// Parse command line for global options and perform some checks.
template <typename TOptions>
ArgumentParser::ParseResult
parseArgumentsAndCheck(TOptions & options,
                       CharString & referenceFilename,
                       ArgumentParser & parser,
                       const int argc,
                       const char * argv[]) {
    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    // If parsing was not successful then exit with code 1 if there were errors.  Otherwise, exit with code 0 (e.g. help
    // was printed).
    if (res != ArgumentParser::PARSE_OK)
        return res;

    options.allowNFromGenome = isSet(parser, "allow-N-from-genome");
    getOptionValue(options.seed, parser, "seed");
    getOptionValue(options.numReads, parser, "num-reads");
    getOptionValue(options.randomSourceLength, parser, "source-length");
    getOptionValue(options.sourceProbabilityA, parser, "source-probability-A");
    getOptionValue(options.sourceProbabilityC, parser, "source-probability-C");
    getOptionValue(options.sourceProbabilityG, parser, "source-probability-G");
    options.onlyForward = isSet(parser, "forward-only");
    options.onlyReverse = isSet(parser, "reverse-only");
    getOptionValue(options.outputFile, parser, "output-file");
    options.simulateQualities = isSet(parser, "simulate-qualities");
    options.includeReadInformation = isSet(parser, "include-read-information");
    options.verbose = isSet(parser, "verbose") || isSet(parser, "very-verbose");
    options.veryVerbose = isSet(parser, "very-verbose");
    
    getOptionValue(options.libraryLengthMean, parser, "library-length-mean");
    getOptionValue(options.libraryLengthError, parser, "library-length-error");
    options.libraryLengthIsUniform = isSet(parser, "library-length-uniform");
    options.generateMatePairs = isSet(parser, "mate-pairs");
    if (isSet(parser, "read-naming")) {
        int mode = 0;
        getOptionValue(mode, parser, "read-naming");
        switch (mode) {
            case 1:
                options.readNaming = READ_NAMING_SLASH_SUFFIX0;
                break;
            case 2:
                options.readNaming = READ_NAMING_SLASH_SUFFIX1;
                break;
            default:
                options.readNaming = READ_NAMING_NOSUFFIX;
        }
    }
    getOptionValue(options.readNamePrefix, parser, "read-name-prefix");
    getOptionValue(options.numHaplotypes, parser, "num-haplotypes");
    getOptionValue(options.haplotypeSnpRate, parser, "haplotype-snp-rate");
    getOptionValue(options.haplotypeIndelRate, parser, "haplotype-indel-rate");
    getOptionValue(options.haplotypeIndelRangeMin, parser, "haplotype-indel-range-min");
    getOptionValue(options.haplotypeIndelRangeMax, parser, "haplotype-indel-range-max");
    options.haplotypeNoN = isSet(parser, "haplotype-no-N");

    getOptionValue(options.sampleCountsFilename, parser, "sample-counts-file");
    getOptionValue(options.vcfFile, parser, "vcf-allele-file");

    // First argument is "illumina", second one name of reference file.
    getArgumentValue(referenceFilename, parser, 1);

    if (referenceFilename == "random.fasta")
        options.useRandomSequence = true;

    return parseArgumentsAndCheckModelSpecific(options, parser);
}

// Write a random DNA sequence of the given length to the file with the given name.
template <typename TRNG, typename TOptions>
int writeRandomSequence(TRNG & rng, size_t length, CharString const & fileName, TOptions const & options)
{
//IOREV uses custom IO, should adapt generic IO
    DnaString randomSequence;
    reserve(randomSequence, length);

    // Simulate the random sequence with the given background probabilities.
    double xA = options.sourceProbabilityA;
    double xC = xA + options.sourceProbabilityC;
    double xG = xC + options.sourceProbabilityG;

    Pdf<Uniform<double> > pdf(0, 1);
    for (size_t i = 0; i < length; ++i) {
        double x = pickRandomNumber(rng, pdf);
        if (x < xA)
            appendValue(randomSequence, Dna('A'));
        else if (x < xC)
            appendValue(randomSequence, Dna('C'));
        else if (x < xG)
            appendValue(randomSequence, Dna('G'));
        else
            appendValue(randomSequence, Dna('T'));
    }

    std::ofstream file;
    file.open(toCString(fileName), std::ios_base::out | std::ios_base::trunc);
    if (!file.is_open()) {
        std::cerr << "Failed to write random sequence to " << fileName << std::endl;
        return 1;
    }
    write(file, randomSequence, "random_sequence", Fasta());
    file.close();
    return 0;
}

// Load sample counts as integers from a file.
template <typename TFragmentStore, typename TSpec, typename TOptions>
int loadSampleCounts(ModelParameters<TSpec> & modelParameters, TFragmentStore /*const*/ & fragmentStore, TOptions const & options)
{
    std::ifstream file;
    file.open(toCString(options.sampleCountsFilename), std::ios_base::in);
    if (!file.is_open()) {
        std::cerr << "Failed to open sample counts file " << options.sampleCountsFilename << std::endl;
        return 1;
    }

    resize(modelParameters.sampleCounts, length(fragmentStore.contigNameStore), 0);

    char c = _streamGet(file);
    CharString contigName;
    size_t sampleCount;

    refresh(fragmentStore.contigNameStoreCache);

    while (!_streamEOF(file)) {
        clear(contigName);

        _parseReadSamIdentifier(file, contigName, c);
        _parseSkipUntilChar(file, '\t', c);
        sampleCount = _parseReadNumber(file, c);
        _parseSkipLine(file, c);

        size_t contigId = 0;
        if (!getIdByName(fragmentStore.contigNameStore, contigName, contigId, fragmentStore.contigNameStoreCache)) {
            std::cerr << "ERROR: Could not find contig with name \"" << contigName << "\" (from read counts file) in contigs." << std::endl;
            return 1;
        }

        if (options.veryVerbose)
            std::cout << "Sample count for contig " << contigName << " (id=" << contigId << ") is " << sampleCount << std::endl;
        modelParameters.sampleCounts[contigId] = sampleCount;
    }

    return 0;
}

// Top level dispatched function that mainly contains the input/output
// logic and dispatches to simulateReadsMain() for the actual
// simulation steps.
template <typename TOptions, typename TReadsTypeTag>
int simulateReads(TOptions options, CharString referenceFilename, TReadsTypeTag const &) {
    // Print options.
    std::cerr << options;
    std::cerr << "reference file: " << referenceFilename << std::endl;

    // Create a Rng with the given seed.
    Rng<MersenneTwister> rng(options.seed);

    // Load or randomly generate the reference sequence.
    FragmentStore<MyFragmentStoreConfig> fragmentStore;
    if (options.useRandomSequence) {
        referenceFilename = "random.fasta";
        std::cerr << "Generating random sequence of length " << options.randomSourceLength
                  << " to file \"" << referenceFilename << "\"." << std::endl;
        int ret = writeRandomSequence(rng, options.randomSourceLength, referenceFilename, options);
        if (ret != 0)
            return ret;
    }
    // Generate output file name if necessary.
    if (options.outputFile == "") {
        options.outputFile = referenceFilename;
        if (options.simulateQualities)
            append(options.outputFile, ".fastq");
        else
            append(options.outputFile, ".fasta");
    }
    // Set read name prefix to file name if necessary.
    if (empty(options.readNamePrefix))
        options.readNamePrefix = options.outputFile;
    // Generate SAM file if necessary.
    if (options.samFile == "") {
        options.samFile = options.outputFile;
        append(options.samFile, ".sam");
    }
    std::cerr << "Loading reference sequence from \"" << referenceFilename << "\"" << std::endl;
    if (!loadContigs(fragmentStore, referenceFilename)) {
        std::cerr << "Could not load reference sequence from " << referenceFilename << std::endl;
        return 1;
    }

    // Load and/or build the model specific parameters for the simulation.
    ModelParameters<TReadsTypeTag> modelParameters;
    int ret = simulateReadsSetupModelSpecificData(modelParameters, options);
    if (ret != 0)
        return ret;

    // Load sample counts file if given.
    if (length(options.sampleCountsFilename) > 0) {
        ret = loadSampleCounts(modelParameters, fragmentStore, options);
        if (ret != 0)
            return ret;
    }
    
    // Kick off the read generation.
    ret = simulateReadsMain(fragmentStore, rng, options, modelParameters);
    if (ret != 0)
        return ret;

    // Write out results.
    if (options.generateMatePairs) {
        // Find first space in read names.  Is the same for all.
        unsigned pos = 0;
        if (options.readNaming != READ_NAMING_NOSUFFIX && !empty(fragmentStore.readNameStore)) {
            while (pos < length(fragmentStore.readNameStore[0]) && fragmentStore.readNameStore[0][pos] != ' ')
                pos += 1;
        }
        // Build filename with '.1.' infix.
        CharString mateFilename = options.outputFile;
        size_t dotPos = 0;
        for (size_t i = 0; i < length(mateFilename); ++i)
            if (mateFilename[i] == '.')
                dotPos = i;
        replace(mateFilename, dotPos, dotPos + 1, "_1.");
        // Write out first mates.
        std::cerr << "Writing resulting reads to \"" << mateFilename << "\" mates/1" << std::endl;
        StringSet<String<Dna5Q>, Dependent<> > reads;
        StringSet<CharString, Owner<> > readNames;
        for (size_t i = 0; i < length(fragmentStore.matePairStore); ++i) {
            size_t readId = fragmentStore.matePairStore[i].readId[0];
            appendValue(readNames, fragmentStore.readNameStore[readId]);
            if (options.readNaming != READ_NAMING_NOSUFFIX) {
                if (options.readNaming == READ_NAMING_SLASH_SUFFIX0)
                    insert(back(readNames), pos, "/0");
                else if (options.readNaming == READ_NAMING_SLASH_SUFFIX1)
                    insert(back(readNames), pos, "/1");
            }
            appendValue(reads, fragmentStore.readSeqStore[readId]);
        }
        {
            std::fstream fstrm(toCString(mateFilename), std::ios_base::out);
            if (!fstrm.is_open()) {
                std::cerr << "Could not open out file \"" << mateFilename << "\"" << std::endl;
                return 1;
            }
            if (options.simulateQualities)
                write(fstrm, readNames, reads, Fastq());
            else
                write(fstrm, readNames, reads, Fasta());
        }
        // Build filename with '_2.' infix.
        replace(mateFilename, dotPos + 1, dotPos + 2, "2");
        // Write out second mates.
        std::cerr << "Writing resulting reads to \"" << mateFilename << "\" mates/2" << std::endl;
        clear(reads);
        clear(readNames);
        for (size_t i = 0; i < length(fragmentStore.matePairStore); ++i) {
            size_t readId = fragmentStore.matePairStore[i].readId[1];
            appendValue(readNames, fragmentStore.readNameStore[readId]);
            if (options.readNaming != READ_NAMING_NOSUFFIX) {
                if (options.readNaming == READ_NAMING_SLASH_SUFFIX0)
                    insert(back(readNames), pos, "/1");
                else if (options.readNaming == READ_NAMING_SLASH_SUFFIX1)
                    insert(back(readNames), pos, "/2");
            }
            appendValue(reads, fragmentStore.readSeqStore[readId]);
        }
        {
            std::fstream fstrm(toCString(mateFilename), std::ios_base::out);
            if (!fstrm.is_open()) {
                std::cerr << "Could not open out file \"" << mateFilename << "\"" << std::endl;
                return 1;
            }
            if (options.simulateQualities)
                write(fstrm, readNames, reads, Fastq());
            else
                write(fstrm, readNames, reads, Fasta());
        }
    } else {
        std::cerr << "Writing resulting reads to \"" << options.outputFile << "\"" << std::endl;
        std::fstream fstrm(toCString(options.outputFile), std::ios_base::out);
        if (!fstrm.is_open()) {
            std::cerr << "Could not open out file \"" << options.outputFile << "\"" << std::endl;
            return 1;
        }
        if (options.simulateQualities)
            write(fstrm, fragmentStore.readNameStore, fragmentStore.readSeqStore, Fastq());
        else
            write(fstrm, fragmentStore.readNameStore, fragmentStore.readSeqStore, Fasta());
    }
    std::cerr << "Writing Sam file \"" << options.samFile << "\"" << std::endl;
    {
        std::fstream fstrm(toCString(options.samFile), std::ios_base::out);
        if (!fstrm.is_open()) {
            std::cerr << "Could not open Sam file \"" << options.samFile << "\"" << std::endl;
            return 1;
        }
        write(fstrm, fragmentStore, Sam());
    }
    return 0;
}

// Load a haplotype from a VCF file.
void loadHaplotype(StringSet<String<Dna5, Journaled<Alloc<> > > > & haplotype,
                    FragmentStore<MyFragmentStoreConfig> & fragmentStore,
                    String<String<Snp> > & snpSet,
                    Options<Global> const & options)
{
    refresh(fragmentStore.contigNameStoreCache);

    clear(snpSet);
    resize(snpSet, length(fragmentStore.contigStore), Exact());

    clear(haplotype);
    resize(haplotype, length(fragmentStore.contigStore), Exact());
    for (unsigned i = 0; i < length(fragmentStore.contigStore); ++i)
        setHost(haplotype[i], fragmentStore.contigStore[i].seq);
    
    std::ifstream file(toCString(options.vcfFile), std::ios_base::in | std::ios_base::binary);
    if (!file.is_open())
    {
        std::cerr << "Couldn't read VCF file " << options.vcfFile << std::endl;
        return;
    }
    
    char c;
    CharString contigName, oldContigName, dummy;
    CharString ref, altStr, formatStr, valuesStr;
    StringSet<CharString> alt, format, values, genotypes;
    int delta = 0; // stores the difference between genomic and allele position
    int oldContigId = -2;
    __int64 ignoreUpTo = -1;
    Snp snp;
    
    // statistics
    unsigned numInconsistencies = 0;
    unsigned numSNPs = 0;
    unsigned numIndels = 0;
    unsigned indelLenSum = 0;
    unsigned unknownContigs = 0;
    
    while (!_streamEOF(file)) 
    {
        // ignore comments
        _parseSkipWhitespace(file, c);
        
        if (c == '#')
        {
            _parseSkipLine(file, c);
            continue;
        }
        
        // read column 1: contig name
        clear(contigName);
        _parseCharsUntilWhitespace(file, c, contigName);
        _parseSkipWhitespace(file, c);
        
        // read column 2: contig position (fwd strand, 1-based counting)
        __int64 pos = _parseReadNumber(file, c);
        --pos;
        _parseSkipWhitespace(file, c);
        
        // skip column 3: id
        clear(dummy);
        _parseCharsUntilWhitespace(file, c, dummy);
        _parseSkipWhitespace(file, c);
        
        // read column 4: ref sequence
        clear(ref);
        _parseCharsUntilWhitespace(file, c, ref);
        _parseSkipWhitespace(file, c);
        
        // read column 5: alt sequence
        clear(altStr);
        _parseCharsUntilWhitespace(file, c, altStr);
        _parseSkipWhitespace(file, c);

        clear(alt);
        strSplit(alt, altStr, ',');
        
        // skip column 6: qual
        clear(dummy);
        _parseCharsUntilWhitespace(file, c, dummy);
        _parseSkipWhitespace(file, c);
        
        // skip column 7: filter
        clear(dummy);
        _parseCharsUntilWhitespace(file, c, dummy);
        _parseSkipWhitespace(file, c);
        
        // skip column 8: info
        clear(dummy);
        _parseCharsUntilWhitespace(file, c, dummy);
        _parseSkipWhitespace(file, c);
        
        // skip column 9: format
        clear(formatStr);
        _parseCharsUntilWhitespace(file, c, formatStr);
        _parseSkipWhitespace(file, c);
        
        // skip column 10: values
        clear(valuesStr);
        _parseCharsUntilWhitespace(file, c, valuesStr);
        _parseSkipLine(file, c);
        
        clear(format);
        clear(values);
        strSplit(format, formatStr, ':');
        strSplit(values, valuesStr, ':');
        unsigned idx;
        for (idx = 0; idx < length(format); ++idx)
            if (format[idx] == "GT")
                break;
        
        // extract genotype information
        if (idx >= length(values))
            continue;
        
        clear(genotypes);
        strSplit(genotypes, values[idx], ':');
        unsigned gt;
        for (gt = 0; gt < length(genotypes); ++gt)
            if (!empty(genotypes[gt]) && genotypes[gt][0] != '0')
                break;
        
        if (gt == length(genotypes))
            continue;

        SEQAN_ASSERT(!empty(genotypes[gt]));
        unsigned allele = genotypes[gt][0] - '1';

        // trim common left characters
        while (!empty(ref))
        {
            char leftChar = ref[0];
            unsigned i;
//            for (i = 0; i < length(alt); ++i)
//            {
//                if (empty(alt[i]) || leftChar != alt[i][0])
//                    break;
//            }
            // cut first characters if they are all equal
//            if (i == length(alt))
            if (!empty(alt[allele]) && leftChar == alt[allele][0])
            {
                ++pos;
                erase(ref, 0);
                i = allele;
//                for (i = 0; i < length(alt); ++i)
                    erase(alt[i], 0);
            }
            else
                break;
        }

        // trim common right characters
        while (!empty(ref))
        {
            char rightChar = back(ref);
            unsigned i;
//            for (i = 0; i < length(alt); ++i)
//            {
//                if (empty(alt[i]) || rightChar != back(alt[i]))
//                    break;
//            }
            // cut last characters if they are all equal
//            if (i == length(alt))
            if (!empty(alt[allele]) && rightChar == back(alt[allele]))
            {
                eraseBack(ref);
                i = allele;
//                for (i = 0; i < length(alt); ++i)
                    eraseBack(alt[i]);
            }
            else
                break;
        }
        
        int contigId;
        
        // cached lookup of contig name 
        if (oldContigName == contigName)
            contigId = oldContigId;
        else
        {
            oldContigName = contigName;
            if (!getIdByName(fragmentStore.contigNameStore, contigName, contigId, fragmentStore.contigNameStoreCache)) 
            {
                dummy = "chr";
                append(dummy, contigName);
                if (!getIdByName(fragmentStore.contigNameStore, dummy, contigId, fragmentStore.contigNameStoreCache)) 
                {
                    contigId = -1;
                    if (oldContigId == -1)
                        oldContigId = -2;
                }
            }
        }
        
        if (oldContigId != contigId)
        {
            if (oldContigId >= 0)
            {
                std::cerr << "    " << fragmentStore.contigNameStore[oldContigId] << "\tSNPs:" << numSNPs << "\tindels:" << numIndels << "\tindel len sum:\t" << indelLenSum << "\tvrate:" << (numSNPs+indelLenSum)/(double)length(fragmentStore.contigStore[oldContigId].seq);
                if (numInconsistencies > 0)
                    std::cerr << "\tinconsistencies:" << numInconsistencies;
                std::cerr << std::endl;
            }
            if (contigId < 0)
            {
                ++unknownContigs;
                std::cerr << "    WARNING: contig " << contigName << " not found" << std::endl;
            }
            delta = 0;
            numIndels = 0;
            indelLenSum = 0;
            numSNPs = 0;
            numInconsistencies = 0;
            ignoreUpTo = -1;
            oldContigId = contigId;
        }
        
        if (contigId < 0)
            continue;

        if (!empty(ref))
        {
            if (pos < ignoreUpTo) continue;
            if (infix(fragmentStore.contigStore[contigId].seq, pos, pos + length(ref)) != ref)
            {
                ++numInconsistencies;
//                std::cerr << "DIFFERENCE: " << contigName << '\t' << pos << "\t." << infix(fragmentStore.contigStore[contigId].seq, pos, pos + length(ref)) << ".\t!=\t." << ref << ".\t" << delta <<std::endl;
            } else
            {
//                std::cerr << "OK: " << contigName << '\t' << pos << "\t." << infix(fragmentStore.contigStore[contigId].seq, pos, pos + length(ref)) << ".\t=>\t." << alt[allele] << ".\t" << delta << std::endl;
            }
            erase(haplotype[contigId], pos + delta, pos + delta + length(ref));
            ignoreUpTo = pos + length(ref);
        }
        if (!empty(alt[allele]))
        {
            insert(haplotype[contigId], pos + delta, alt[allele]);
//            if (empty(ref))
//                std::cerr << "OK: " << contigName << '\t' << pos << "\t..\t=>\t." << alt[allele] << ".\t" << delta << std::endl;
        }        
        
//        snp.contigId = contigId;
        snp.virtualPos = pos + delta;
        if (length(ref) == length(alt[allele]))
        {
            // HAPTYPE  ----GATTACA---- 
            // REF      ----CACACAC----
            snp.type = ERROR_TYPE_MISMATCH;
            snp.length = length(ref);
            indelLenSum += snp.length;
            ++numSNPs;
        }
        else if (length(ref) < length(alt[allele]))
        {
            // HAPTYPE  ----GATTACA---- 
            // REF      ----       ----
            snp.type = ERROR_TYPE_INSERT;
            snp.length = length(alt[allele]);
            indelLenSum += snp.length;
            ++numIndels;
        }
        else
        {
            // HAPTYPE  ----       ----
            // REF      ----GATTACA---- 
            snp.type = ERROR_TYPE_DELETE;
            snp.length = length(ref);
            ++numIndels;
        }
        appendValue(snpSet[contigId], snp);

        delta += length(alt[allele]) - length(ref);
    }
    if (oldContigId >= 0)
    {
        std::cerr << "    " << fragmentStore.contigNameStore[oldContigId] << "\tSNPs:" << numSNPs << "\tindels:" << numIndels << "\tindel len sum:\t" << indelLenSum << "\tvrate:" << (numSNPs+indelLenSum)/(double)length(fragmentStore.contigStore[oldContigId].seq);
        if (numInconsistencies > 0)
            std::cerr << "\tinconsistencies:" << numInconsistencies;
        std::cerr << std::endl;
    }
    if (unknownContigs > 0)
        std::cerr << "    " << unknownContigs << " contigs not found" << std::endl;

    SEQAN_OMP_PRAGMA(parallel for schedule(static,1))
    for (int i = 0; i < (int)length(snpSet); ++i)
        std::sort(begin(snpSet[i], Standard()), end(snpSet[i], Standard()), SnpLess());
}



// Build a haplotype, based on the contigs from the given fragment store.
template <typename TRNG>
void buildHaplotype(StringSet<String<Dna5, Journaled<Alloc<> > > > & haplotype,
                    FragmentStore<MyFragmentStoreConfig> & fragmentStore,
                    String<String<Snp> > & snpSet,
                    TRNG & rng,
                    Options<Global> const & options)
{
    resize(haplotype, length(fragmentStore.contigStore), Exact());
    String<Dna5> buffer;
    reserve(buffer, options.haplotypeIndelRangeMax);
    Snp snp;

    for (unsigned i = 0; i < length(fragmentStore.contigStore); ++i)
    {
        // statistics
        unsigned numSNPs = 0;
        unsigned numIndels = 0;
        unsigned indelLenSum = 0;

        std::cout << "    contig # " << i+1 << "/" << length(fragmentStore.contigStore) << '\t' << fragmentStore.contigNameStore[i] << std::flush;
        clear(haplotype[i]);
        setHost(haplotype[i], fragmentStore.contigStore[i].seq);
        String<Dna5> const & contig = fragmentStore.contigStore[i].seq;
        String<Dna5, Journaled<Alloc<> > > & haplotypeContig = haplotype[i];

        // Only generate Ns in the haplotype if allowed.
        int maxOrdValue = options.haplotypeNoN ? 3 : 4;

        // j is position in original sequence, k is position in haplotype
        for (size_t j = 0, k = 0; j < length(contig);)
        {
            double x = pickRandomNumber(rng, Pdf<Uniform<double> >(0, 1));
            snp.virtualPos = k;
            if (x < options.haplotypeSnpRate) {
                // SNP
                ++numSNPs;
                snp.length = 1;

                Dna5 c = Dna5(pickRandomNumber(rng, Pdf<Uniform<int> >(0, maxOrdValue - 1)));
                if (c == contig[j])
                    c = Dna5(ordValue(c) + 1);
                if (options.haplotypeNoN)
                    SEQAN_ASSERT(c != Dna5('N'));
                assignValue(haplotypeContig, k, c);
                ++j;
                ++k;

                // HAPTYPE  ----GATTACA---- 
                // REF      ----CACACAC----
                snp.type = ERROR_TYPE_MISMATCH;
                appendValue(snpSet[i], snp);
            }
            else if (x < options.haplotypeSnpRate + options.haplotypeIndelRate) 
            {
                // Indel of random length.
                ++numIndels;
                unsigned rangeLen = options.haplotypeIndelRangeMax - options.haplotypeIndelRangeMin;
                unsigned indelLen = options.haplotypeIndelRangeMin + static_cast<unsigned>(pickRandomNumber(rng, Pdf<Uniform<double> >(0, 1)) * rangeLen);
                snp.length = indelLen;
                indelLenSum += indelLen;
                if (pickRandomNumber(rng, Pdf<Uniform<double> >(0, 1)) < 0.5) 
                {
                    // Insertion.
                    clear(buffer);
                    for (unsigned ii = 0; ii < indelLen; ++ii)
                        appendValue(buffer, Dna5(pickRandomNumber(rng, Pdf<Uniform<int> >(0, maxOrdValue))));
                    insert(haplotypeContig, k, buffer);
                    k += indelLen;

                    // HAPTYPE  ----GATTACA---- 
                    // REF      ----       ----
                    snp.type = ERROR_TYPE_INSERT;
                }
                else
                {
                    // Deletion.
                    indelLen = _min(indelLen, length(haplotypeContig) - k);
                    erase(haplotypeContig, k, k + indelLen);
                    j += indelLen;

                    // HAPTYPE  ----       ----
                    // REF      ----GATTACA---- 
                    snp.type = ERROR_TYPE_DELETE;
                }
                appendValue(snpSet[i], snp);
            }
            else 
            {
                // Match.
                j += 1;
                k += 1;
            }
        }
        std::cout << "\tSNPs:" << numSNPs << "\tindels:" << numIndels << "\tindel len sum:\t" << indelLenSum << "\tvrate:" << (numSNPs+indelLenSum)/(double)length(contig) << std::endl;
    }
}

// Build read simulation instructions for a haplotype.
//
// pick a contig, probability is proportional to the length
// pick a start position, end position = start position + read length
// pick whether to match on the forward or reverse strand
// simulate edit string
// build quality values
// possibly adjust mate if left read has insert at the beginning or right read has insert at the right
//
// Note that we set the isForward flag of the simulation instructions here.  Notably, the technology specific simulation
// parts do not interpret this flag but only provide simulation instructions on the forward strand.  Later, when actually
// cutting out the sequences and modifying the sampled read, we have to take this in mind.  This means, there we have to
// reverse the edit string etc. there.
template <typename TReadsTag, typename TRNG, typename THaplotypeSpec>
int buildReadSimulationInstruction(
        String<ReadSimulationInstruction<TReadsTag> > & instructions,
        TRNG & rng,
        unsigned const & haplotypeId,
        StringSet<String<Dna5, THaplotypeSpec> > const & haplotype,
        String<double> const & relativeContigLengths,
        size_t const & contigId,
        bool fixedContigId,
        ModelParameters<TReadsTag> const & parameters,
        Options<TReadsTag> const & options)
{
    ReadSimulationInstruction<TReadsTag> inst;
    inst.haplotype = haplotypeId;

    // We have to retry simulation if the mate pair did not fit in.
    bool invalid = false;
    do {
        clear(instructions);
        invalid = false;  // By default, we do not want to repeat.
        if (fixedContigId) {
            // Use precomputed contig id.
            inst.contigId = contigId;
        } else {
            // Pick contig id, probability is proportional to the length.
            double x = pickRandomNumber(rng, Pdf<Uniform<double> >(0, 1));
            for (unsigned i = 0; i < length(relativeContigLengths); ++i) {
                if (x < relativeContigLengths[i]) {
                    inst.contigId = i - 1;
                    break;
                }
            }
        }
        // Pick whether on forward or reverse strand.
        if (options.onlyForward)
            inst.isForward = true;
        else if (options.onlyReverse)
            inst.isForward = false;
        else
            inst.isForward = pickRandomNumber(rng, Pdf<Uniform<int> >(0, 1));
        // Pick the length in the haplotype infix the read comes from, possibly randomly.
        unsigned readLength = pickReadLength(rng, options);
        // This cannot work if the haplotype is shorter than the length of the read to simulate.
        if (length(haplotype[inst.contigId]) < readLength) {
            std::cerr << "ERROR: haplotype (== " << length(haplotype[inst.contigId]) << ") < read length!" << std::endl;
            return 1;
        }
        // Pick a start and end position.
        inst.beginPos = pickRandomNumber(rng, Pdf<Uniform<size_t> >(0, length(haplotype[inst.contigId]) - readLength - 1));
        inst.endPos = inst.beginPos + readLength;
        // Simulate the read with these parameters.
        buildSimulationInstructions(inst, rng, readLength, haplotype[inst.contigId], parameters, options);
        // Append read to result list.
        appendValue(instructions, inst);

        ReadSimulationInstruction<TReadsTag> inst2(inst);

        // Maybe create a mate for this read.
        if (options.generateMatePairs) {
            // Pick a read length, possibly randomly.
            unsigned readLength = pickReadLength(rng, options);
            // Pick a library length, according to the options.
            size_t libraryLength = pickLibraryLength(rng, options);
            // Compute start and end position.
            if (inst.isForward)
            {
                inst.endPos = inst.beginPos + libraryLength;
                inst.beginPos = inst.endPos - readLength;
            }
            else
            {
                if (inst.endPos < libraryLength)
                    invalid = true;
                inst.beginPos = inst.endPos - libraryLength;
                inst.endPos = inst.beginPos + readLength;
            }
            // Set orientation of second mate.
            inst.isForward = !back(instructions).isForward;
            // Verify that the mate fits right of the originally simulated read.
            size_t contigLength = length(haplotype[inst.contigId]);
            if ((inst.beginPos > contigLength) || (inst.endPos > contigLength)) {
                // Mate did not fit!  Remove previously added read and set
                // invalid to true so we repeat this simulation.
                SEQAN_ASSERT_GT(length(instructions), 0u);
                eraseBack(instructions);
                invalid = true;
                if (options.verbose) {
                    std::cerr << "INFO: Mate did not fit! Repeating..." << std::endl;
                    std::cerr << "      inst2.beginPos == " << inst2.beginPos << ", inst2.endPos == " << inst2.endPos << std::endl;
                    std::cerr << "       inst.beginPos == " << inst.beginPos << ",  inst.endPos == " << inst.endPos << std::endl;
                }
                continue;
            }
            // Simulate the read with these parameters.
            buildSimulationInstructions(inst, rng, readLength, haplotype[inst.contigId], parameters, options);
            // Append read to result list.
            appendValue(instructions, inst);
        }

        // Check whether there are Ns in the selected areas.
        if (!options.allowNFromGenome) {
            for (unsigned i = 0; i < length(instructions); ++i) {
                String<Dna5, THaplotypeSpec> const & contig = haplotype[instructions[i].contigId];
                typedef typename Position<Dna5String>::Type TPosition;
                TPosition beginPos = instructions[i].beginPos;
                TPosition endPos = instructions[i].endPos;
                if (beginPos > endPos)
                    std::swap(beginPos, endPos);
                for (unsigned i = beginPos; i != endPos; ++i) {
                    if (contig[i] == Dna5('N')) {
                        invalid = true;
                        break;
                    }
                }
            }
        }
    } while (invalid);
	
	if (options.generateMatePairs)
		SEQAN_ASSERT_EQ(length(instructions), 2u);
	else
		SEQAN_ASSERT_EQ(length(instructions), 1u);
	
    return 0;
}

// Pick library length, based on the configuration in options.
template <typename TRNG>
inline
unsigned pickLibraryLength(TRNG & rng, Options<Global> const & options)
{
    if (options.libraryLengthIsUniform) {
        // Pick uniformly.
        double minLen = options.libraryLengthMean - options.libraryLengthError;
        double maxLen = options.libraryLengthMean + options.libraryLengthError;
        double len = pickRandomNumber(rng, Pdf<Uniform<double> >(minLen, maxLen));
        return static_cast<unsigned>(_max(0.0, round(len)));
    } else {
        // Pick normally distributed.
        double len = pickRandomNumber(rng, Pdf<Normal>(options.libraryLengthMean, options.libraryLengthError));
        return static_cast<unsigned>(_max(0.0, round(len)));
    }
}

// Performs the actual read simulation.
template <typename TRNG, typename TReadsTag, typename TOptions>
int simulateReadsMain(FragmentStore<MyFragmentStoreConfig> & fragmentStore,
                      TRNG & rng,
                      TOptions const & options,
                      ModelParameters<TReadsTag> const & parameters) {
    typedef FragmentStore<MyFragmentStoreConfig> TFragmentStore;
    typedef Value<TFragmentStore::TMatePairStore>::Type TMatePairStoreElement;
    typedef typename Value<typename TFragmentStore::TAlignedReadStore>::Type	TAlignedElement;
    typedef typename TAlignedElement::TGapAnchors								TReadGapAnchors;

    if (options.verbose)
        std::cerr << "Simulating reads..." << std::endl;

    typedef Position<CharString>::Type TPos;

    // Number of reads comes from command line by default.  If sample
    // counts are given, we compute it from this instead.
    size_t numReads = options.numReads;
    if (length(parameters.sampleCounts) > 0u) {
        numReads = 0;
        for (unsigned i = 0; i < length(parameters.sampleCounts); ++i)
            numReads += parameters.sampleCounts[i];
    }

    // First, we randomly pick the haplotype for each read to be
    // simulated or read it from the sample counts in parameters.
    String<unsigned> haplotypeIds;
    // Pick random haplotype origin.
    reserve(haplotypeIds, numReads);
    for (size_t i = 0; i < numReads; ++i)
        appendValue(haplotypeIds, pickRandomNumber(rng, Pdf<Uniform<unsigned> >(0, options.numHaplotypes - 1)));

    // Maybe pick contig ids to sample from.
    String<unsigned> contigIds;
    if (length(parameters.sampleCounts) > 0) {
        for (unsigned i = 0; i < length(fragmentStore.contigNameStore); ++i) {
            resize(contigIds, length(contigIds) + parameters.sampleCounts[i], i);
            if (options.veryVerbose)
                std::cerr << parameters.sampleCounts[i] << " reads from haplotype " << i << "..." << std::endl;
        }
        shuffle(contigIds, rng);
    }

    // We do not build all haplotypes at once since this could cost a
    // lot of memory.
    //
    // TODO(holtgrew): Would only have to switch pointers to journals which is possible.
    //
    // for each haplotype id
    //   simulate haplotype
    //   for each simulation instruction for this haplotype:
    //     build simulated read
    reserve(fragmentStore.readSeqStore, numReads, Exact());
    reserve(fragmentStore.readNameStore, numReads, Exact());
    CharString readNameBuf;
    char outFileName[151];
    snprintf(outFileName, 150, "%s", toCString(options.readNamePrefix));
    for (unsigned haplotypeId = 0; haplotypeId < options.numHaplotypes; ++haplotypeId)
    {
        std::cerr << "Simulating for haplotype #" << haplotypeId << "..." << std::endl;
        std::cout << "  Building haplotype..." << std::endl;
        StringSet<String<Dna5, Journaled<Alloc<> > > > haplotypeContigs;

        typedef String<Snp> TContigSnpSet;
        String<TContigSnpSet> snpSet;
        resize(snpSet, length(fragmentStore.contigNameStore));

        double buildStart = sysTime();

        if (empty(options.vcfFile))
            buildHaplotype(haplotypeContigs, fragmentStore, snpSet, rng, options);
        else
            loadHaplotype(haplotypeContigs, fragmentStore, snpSet, options);
        std::cout << "  Finished haplotype creation in " << (sysTime() - buildStart) << 's' << std::endl;

        // TODO(holtgrew): Assigning of string set with compatible string should be possible.
        StringSet<String<Dna5> > haplotypeContigsCopy;
        for (unsigned i = 0; i < length(haplotypeContigs); ++i)
            appendValue(haplotypeContigsCopy, haplotypeContigs[i]);

        // Build partial sums over relative contig lengths so we can pick the contigs later on.
        size_t totalLength = 0;
        for (unsigned i = 0; i < length(fragmentStore.contigStore); ++i)
            totalLength += length(fragmentStore.contigStore[i].seq);
        String<double> relativeContigLengths;
        resize(relativeContigLengths, length(fragmentStore.contigStore) + 1, Exact());
        front(relativeContigLengths) = 0.0;
        for (unsigned i = 0; i < length(fragmentStore.contigStore); ++i) {
            double l = static_cast<double>(length(fragmentStore.contigStore[i].seq));
            relativeContigLengths[i + 1] = l / totalLength;
        }
        std::partial_sum(begin(relativeContigLengths), end(relativeContigLengths), begin(relativeContigLengths));
        back(relativeContigLengths) = 1.0;

        // Simulate the reads...
        std::cerr << "  Simulating reads for haplotype #" << haplotypeId << "..." << std::endl;

//         std::cerr << "Journal: " << haplotypeContigs[0]._journalEntries << std::endl;

        Snp searchSnp;
        char tagBuffer[40];

        for (unsigned j = 0; j < length(haplotypeIds); ++j) {
            if (haplotypeIds[j] != haplotypeId)
                continue;  // Guard against instructions on wrong haplotype.

            // Build simulation instructions.
            String<ReadSimulationInstruction<TReadsTag> > instructions;
            // TODO(holtgrew): Pick contig id outside of instructions.
            size_t contigId = 0;
            bool fixedContigId = false;
            if (length(parameters.sampleCounts) > 0) {
                fixedContigId = true;
                if (options.generateMatePairs)
                    contigId = contigIds[length(fragmentStore.readSeqStore) / 2];
                else
                    contigId = contigIds[length(fragmentStore.readSeqStore)];
            }
            int res = buildReadSimulationInstruction(instructions, rng, haplotypeId, haplotypeContigsCopy, relativeContigLengths, contigId, fixedContigId, parameters, options);
            if (res != 0)
                return res;

            for (unsigned k = 0; k < length(instructions); ++k) {
                ReadSimulationInstruction<TReadsTag> & inst = instructions[k];
                // Reverse edit string and qualities string if this instruction simulates from reverse strand.
                if (!inst.isForward)
                {
                    reverse(inst.editString);
                    reverse(inst.qualities);
                }
                // Apply simulation instructions.
                SEQAN_ASSERT_EQ(length(fragmentStore.readSeqStore), length(fragmentStore.readNameStore));
                // Cut out segment from haplotype.
                String<Dna5Q> read = infix(haplotypeContigsCopy[inst.contigId], inst.beginPos, inst.endPos);
                String<Dna5Q> haplotypeInfix = read;  // Copy for printing later on.
                applySimulationInstructions(read, rng, inst, options);
                if (!inst.isForward)
                {
                    reverseComplement(read);
                    // Reconstruct edit string and qualities in sequencing direction.
                    reverse(inst.editString);
                    reverse(inst.qualities);
                }
                // Append read sequence to read seq store and mate pair to read name store.  This also yields the read
                // id.  We will generate and append the read name below, depending on the read id.
                unsigned readId;
                if (options.generateMatePairs)
                    readId = appendRead(fragmentStore, read, length(fragmentStore.matePairStore));
                else
                    readId = appendRead(fragmentStore, read);

                // Get expected begin/end position in the original sequence.
                TPos origBeginPos = virtualToHostPosition(haplotypeContigs[inst.contigId], inst.beginPos);
                TPos origEndPos = virtualToHostPosition(haplotypeContigs[inst.contigId], inst.endPos);
                
                unsigned numSNPs = 0;
                unsigned numIndels = 0;
                
                if (options.includeReadInformation && inst.contigId < length(snpSet))
                {
                    typedef Iterator<TContigSnpSet, Standard>::Type TIter;                    
                    TIter allBeg = begin(snpSet[inst.contigId], Standard());
                    TIter allEnd = end(snpSet[inst.contigId], Standard());

//                    searchSnp.contigId = inst.contigId;
                    searchSnp.virtualPos = inst.beginPos;
                    TIter snpBeg = std::lower_bound(allBeg, allEnd, searchSnp, SnpLess());
                    searchSnp.virtualPos = inst.endPos;
                    TIter snpEnd = std::upper_bound(allBeg, allEnd, searchSnp, SnpLess());
                    
                    // find first SNP that overlaps the read
                    while (snpBeg != allBeg && (snpBeg - 1)->virtualPos + (snpBeg - 1)->length > inst.beginPos)
                        --snpBeg;
                    
//                    if (readId == 88704)
//                    {
//                        std::cout << "begPos:\t" << inst.beginPos << std::endl;
//                        std::cout << "endPos:\t" << inst.endPos << std::endl;
//                        for (TIter it = snpBeg; it != snpEnd; ++it)
//                            std::cout << it->virtualPos << '\t' << (int)it->type << '\t' << it->length << std::endl;
//                    }

                
                    for (TIter it = snpBeg; it != snpEnd; ++it)
                    {
                        int len = it->length;

                        // cut at the right end of the read
                        if (it->virtualPos + len > inst.endPos)
                            len = inst.endPos - it->virtualPos;

                        // cut at the left end of the read
                        if (it->virtualPos <= inst.beginPos)
                        {
                            // deletions must occur in the read to affect it
                            if (it->type == ERROR_TYPE_DELETE)
                                continue;
                            len -= inst.beginPos - it->virtualPos;
                        }
                        
                        if (len <= 0) continue;
                        
                        if (it->type == ERROR_TYPE_MISMATCH)
                            numSNPs += len;
                        else if (it->type == ERROR_TYPE_INSERT)
                            numIndels += len;
                        else // if (it->type == ERROR_TYPE_DELETE)
                            // at least one base left and right of the deletion is required
                            if (it->virtualPos < inst.endPos)
                                numIndels += it->length;
                    }
                }
                

                // Generate read name.
                // TODO(holtgrew): Remove mateNum, not required?
                if (options.generateMatePairs)
                {
                    // Generate the mate num \in {1, 2}, randomly but consistent so two entries belonging together have
                    // different nums.
                    if (options.includeReadInformation)
                    {
                        resize(readNameBuf, 1024 + length(haplotypeInfix) + length(inst.editString));
                        sprintf(&readNameBuf[0], "%s.%09u contig=%s haplotype=%u length=%lu orig_begin=%lu orig_end=%lu snps=%u indels=%u haplotype_infix=%s edit_string=", outFileName, readId / 2, toCString(fragmentStore.contigNameStore[inst.contigId]), haplotypeId, static_cast<long unsigned>(length(read)), static_cast<long unsigned>(origBeginPos), static_cast<long unsigned>(origEndPos),                             numSNPs, numIndels, toCString(CharString(haplotypeInfix)));
                    }
                    else
                    {
                        resize(readNameBuf, 1024);
                        sprintf(&readNameBuf[0], "%s.%09u", outFileName, readId / 2);
                    }
                } else {
                    if (options.includeReadInformation)
                    {
                        resize(readNameBuf, 1024 + length(haplotypeInfix) + length(inst.editString));
                        sprintf(&readNameBuf[0], "%s.%09u contig=%s haplotype=%u length=%lu orig_begin=%lu orig_end=%lu snps=%u indels=%u haplotype_infix=%s edit_string=", outFileName, readId, toCString(fragmentStore.contigNameStore[inst.contigId]), haplotypeId, static_cast<long unsigned>(length(read)), static_cast<long unsigned>(origBeginPos), static_cast<long unsigned>(origEndPos),                             numSNPs, numIndels, toCString(CharString(haplotypeInfix)));
                    }
                    else
                    {
                        resize(readNameBuf, 1024);
                        sprintf(&readNameBuf[0], "%s.%09u", outFileName, readId);
                    }
                }
                if (options.includeReadInformation) 
                {
                    char buffer[2] = "*";
                    for (unsigned i = 0; i < length(inst.editString); ++i) {
                        buffer[0] = "MEID"[static_cast<int>(inst.editString[i])];
                        strcat(&readNameBuf[0], buffer);
                    }
                }
                
                CharString readName(&readNameBuf[0]);
                if (options.includeReadInformation)
                {
                    if (inst.isForward)
                        append(readName, " strand=forward");
                    else
                        append(readName, " strand=reverse");
                }
                appendValue(fragmentStore.readNameStore, readName);

                // Print info about read and haplotype.
                if (options.veryVerbose) {
                    std::cout << ",-- Read #" << readId << std::endl
                              << "| inst.beginPos    " << inst.beginPos << std::endl
                              << "| inst.endPos      " << inst.endPos << std::endl
                              << "| origBeginPos     " << origBeginPos << std::endl
                              << "| origEndPos       " << origEndPos << std::endl
                              << "| numSNPs          " << numSNPs << std::endl
                              << "| numIndels        " << numIndels << std::endl
                              << "| isgapinhost      " << isGapInHost(haplotypeContigs[inst.contigId], inst.beginPos-1) << std::endl
                              << "| isgapinhost      " << isGapInHost(haplotypeContigs[inst.contigId], inst.beginPos) << std::endl
                              << "| isgapinhost      " << isGapInHost(haplotypeContigs[inst.contigId], inst.beginPos+1) << std::endl
                              << "| name:            " << readName << std::endl
                              << "| original infix:  " << infix(fragmentStore.contigStore[inst.contigId].seq, origBeginPos, origEndPos) << std::endl
                              << "| haplotype infix: " << infix(haplotypeContigsCopy[inst.contigId], inst.beginPos, inst.endPos) << std::endl
                              << "| read:            " << read << std::endl
                              << "`-- " << std::endl;
                }

                // Swap original begin and end position if from reverse strand.
                if (!inst.isForward)
                    std::swap(origBeginPos, origEndPos);

                // Add matches to aligned read store.
                if (options.generateMatePairs)
                    appendAlignedRead(fragmentStore, readId, inst.contigId, origBeginPos, origEndPos, length(fragmentStore.matePairStore));
                else
                    appendAlignedRead(fragmentStore, readId, inst.contigId, origBeginPos, origEndPos);
                
                sprintf(tagBuffer, "XE:i:%u\tXS:i:%u\tXI:i:%u", (inst.mismatchCount + inst.insCount + inst.delCount), numSNPs, numIndels);
                appendValue(fragmentStore.alignedReadTagStore, tagBuffer);

                // Adding mate pair information.
                if (options.generateMatePairs && readId % 2 == 1) {  // Only append mate pair info after simulating second mate.
                    // Append mate pair element to fragment store's mate pair store.
                    TMatePairStoreElement matePair;
                    matePair.readId[0] = readId - 1;
                    matePair.readId[1] = readId;
                    appendValue(fragmentStore.matePairStore, matePair);
                }
            }
			if (options.generateMatePairs)
            {
                // When generating mate pairs, an even number of reads is generated in each step.
 				SEQAN_ASSERT_EQ(length(fragmentStore.alignedReadStore) % 2, 0u);
				SEQAN_ASSERT_EQ(length(fragmentStore.readNameStore) % 2, 0u);
				SEQAN_ASSERT_EQ(length(fragmentStore.readSeqStore) % 2, 0u);
			}
		}
    }

    // Last but not least, convert the matches collected before to a global alignment.
    convertMatchesToGlobalAlignment(fragmentStore, Score<int, EditDistance>(), True());
	
	// AlignedReadLayout layout;
	// layoutAlignment(layout, fragmentStore);
	// printAlignment(std::cout, Raw(), layout, fragmentStore, 0, 0, 300, 0, 100);
    
    if (options.verbose)
        std::cerr << "Simulated " << length(fragmentStore.readSeqStore) << " reads" << std::endl;

    return 0;
}

#endif  // MASON_H_
