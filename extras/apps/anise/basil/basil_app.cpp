// ==========================================================================
//                                 BASIL
// ==========================================================================
// Copyright (c) 2006-2014, Knut Reinert, FU Berlin
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

#include "basil_app.h"

#include <set>
#include <thread>
#include <functional>

#include <seqan/bam_io.h>
#include <seqan/vcf_io.h>

#include "progress_indicator.h"
#include "cluster_matching.h"
#include "thread_safe_queue.h"

namespace  // anonymous
{

// ----------------------------------------------------------------------------
// Function insertOrErase()
// ----------------------------------------------------------------------------

void insertOrErase(std::set<seqan::BamAlignmentRecord *> & set, seqan::BamAlignmentRecord * ptr)
{
    // printf("set.count(%p) == %d\n", ptr, (int)set.count(ptr));
    if (set.count(ptr))
    {
        set.erase(ptr);
        // static int i = 0;
        // fprintf(stderr, "%d deletions\n", ++i);
        // fprintf(stderr, "DEALLOCATED\t%p\n", ptr);
        delete ptr;
    }
    else
    {
        set.insert(ptr);
    }
    // fprintf(stderr, "set.size() == %d\n", (int)set.size());
    // fprintf(stderr, "    set.count(%p) == %d\n", ptr, (int)set.count(ptr));
}

// ----------------------------------------------------------------------------
// Function getVerbosityStr()
// ----------------------------------------------------------------------------

char const * getVerbosityStr(int i)
{
    switch (i)
    {
        case 0:
            return "quiet";
        case 1:
            return "normal";
        case 2:
            return "verbose";
        default:
            return "very-verbose";
    }
}

// ----------------------------------------------------------------------------
// Function getYesNo()
// ----------------------------------------------------------------------------

char const * getYesNo(bool b)
{
    return b ? "YES" : "NO";
}

// ----------------------------------------------------------------------------
// Function getClusterSelection()
// ----------------------------------------------------------------------------

char const * getClusterSelection(OeaClusterSelection s)
{
    switch (s)
    {
        case OeaClusterSelection::CHAINING:
            return "CHAINING";
        case OeaClusterSelection::SET_COVER:
            return "SET COVER";
        default:
            return "<invalid>";
    }
}

}  // anonymous namespace

// ----------------------------------------------------------------------------
// Class BasilOptions
// ----------------------------------------------------------------------------

void BasilOptions::print(std::ostream & stream)
{
    stream << "____OPTIONS___________________________________________________________________\n"
           << "\n"
           << "VERBOSITY                 \t" << getVerbosityStr(verbosity) << "\n"
           << "\n"
           // << "NUM THREADS               \t" << numThreads << "\n"
           // << "CHUNK SIZE                \t" << localAlignmentChunkSize << "\n"
           // << "\n"
           << "IN FILE                   \t" << inputFile << "\n"
           << "REFERENCE FILE            \t" << referenceFile << "\n"
           << "OUTPUT VCF FILE           \t" << outVcfFile << "\n"
           << "OUTPUT DEBUG DIR          \t" << outputDebugDir << "\n"
           << "\n"
           << "AUTO LIBRARY INFO         \t" << getYesNo(autoLibraryInfo) << "\n"
           << "AUTO LIBRARY NUM RECORDS  \t" << autoLibraryNumRecords << "\n"
           << "CMD LINE LIBRARY INFO\n"
           << "  FRAGMENT LEN MEDIAN     \t" << libraryInfo.median << "\n"
           << "  FRAGMENT LEN STDDEV     \t" << libraryInfo.stdDev << "\n"
           // << "  FRAGMENT LEN MAX NORMAL \t" << libraryInfo.maxNormalISize << "\n"
           << "  DEFAULT ORIENTATION     \t" << getOrientationStr(libraryInfo.defaultOrient) << "\n"
           << "FRAGMENT SIZE FACTOR      \t" << fragmentSizeFactor << "\n"
           << "MAX FRAGMENT SIZE         \t" << maxFragmentSize() << "\n"
           << "\n"
           << "REALIGNMENT NUM THREADS   \t" << realignmentNumThreads << "\n"
           << "\n"
           << "MAX COVERAGE              \t" << filterMaxCoverage << "\n"
           << "MIN ALN QUALITY           \t" << filterMinAlignmentQuality << "\n"
           << "\n"
           << "MAX ALIGNMENT LENGTH      \t" << maxAlignmentLength << "\n"
           << "\n"
           << "CLIPPING CLUSTERING\n"
           << "  CLIPPING WINDOW RADIUS  \t" << clippingWindowRadius << "\n"
           << "  CLIPPING MIN LENGTH     \t" << clippingMinLength << "\n"
           << "  CLIPPING MIN COVERAGE   \t" << clippingMinCoverage << "\n"
           << "\n"
           << "BREAKPOINT WINDOW RADIUS\t" << breakpointWindowRadius << "\n";
    if (!breakpointWindowRadius)
        stream << "  => USING FRAGMENT SIZE\n";
    stream << "\n"
           << "OEA CLUSTER SELECTION     \t" << getClusterSelection(oeaClusterSelection) << "\n"
           << "MIN OEA SUPPORT           \t" << oeaMinSupport << "\n"
           << "MIN OEA SIDE SUPPORT      \t" << oeaMinSupportEachSide << "\n";
}

seqan::ArgumentParser::ParseResult BasilOptions::parseCommandLine(int argc, char const ** argv)
{
    // Declare ArgumentParser variable.
    seqan::ArgumentParser parser("basil");

    // Set short description, version, and date.
    setShortDescription(parser, "BASe-resolution Insert Locator");
#ifdef SEQAN_REVISION
        setVersion(parser, "0.2.0-beta.1 [" + std::string(SEQAN_REVISION) + "]");
#else
        setVersion(parser, "0.2.0-beta.1");
#endif
#ifdef SEQAN_DATE
        setDate(parser, SEQAN_DATE);
#endif
    // Set usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fB-ir\\fP \\fIREF\\fP \\fB-im\\fP \\fIMAPPING\\fP \\fB-ov\\fP \\fIOUT.vcf\\fP");
    addDescription(parser, "Scan SAM/BAM file \\fIMAPPING\\fP for signatures of structural variations.  The reference ");
    addDescription(parser, "is given by the FASTA file \\fIREF\\fP.");

    // Define BasilOptions -- General
    addOption(parser, seqan::ArgParseOption("q", "quiet", "Only print on errors."));
    addOption(parser, seqan::ArgParseOption("v", "verbose", "Higher verbosity."));
    addOption(parser, seqan::ArgParseOption("vv", "very-verbose", "Highest verbosity."));

    // // Define Parallelism Options.
    // addOption(parser, seqan::ArgParseOption("nt", "num-threads", "Number of threads to use.", seqan::ArgParseOption::INTEGER, "NUM"));
    // setMinValue(parser, "num-threads", "1");
    // setDefaultValue(parser, "num-threads", "1");
    // addOption(parser, seqan::ArgParseOption("", "local-alignment-chunk-size", "Chunk size for local alignment.", seqan::ArgParseOption::INTEGER, "NUM"));
    // setMinValue(parser, "local-alignment-chunk-size", "1");
    // setDefaultValue(parser, "local-alignment-chunk-size", "65536");

    // Define BasilOptions -- Input / Output Related
    addSection(parser, "Input / Output BasilOptions");
    addOption(parser, seqan::ArgParseOption("ir", "input-reference", "FASTA file with the reference. Required when using", seqan::ArgParseOption::INPUTFILE, "REF.fa"));
    setRequired(parser, "input-reference");
    setValidValues(parser, "input-reference", "fasta fa");
    addOption(parser, seqan::ArgParseOption("im", "input-mapping", "SAM/BAM file to use as the input.", seqan::ArgParseOption::INPUTFILE, "IN"));
    setRequired(parser, "input-mapping");
    setValidValues(parser, "input-mapping", "sam bam");
    addOption(parser, seqan::ArgParseOption("ov", "out-vcf", "VCF file to write variations to. Use \"-\" for stdout.", seqan::ArgParseOption::OUTPUTFILE, "OUT"));
    setDefaultValue(parser, "out-vcf", "-");
    setValidValues(parser, "out-vcf", "vcf");
    addOption(parser, seqan::ArgParseOption("", "output-debug-dir", "Directory for debug output files.  Created if required.", seqan::ArgParseOption::STRING, "PATH"));

    // Define BasilOptions -- Library Info
    addSection(parser, "Library Info");

    addOption(parser, seqan::ArgParseOption("", "fragment-size-factor", "Factor to multiple fragment size stddev with to get allowed error.", seqan::ArgParseOption::DOUBLE, "FACTOR"));
    setMinValue(parser, "fragment-size-factor", "0");
    setDefaultValue(parser, "fragment-size-factor", "8");

    addOption(parser, seqan::ArgParseOption("", "auto-library-num-records", "Number of records to use for automatic library evaluation.  Set to 0 to evaluate all.", seqan::ArgParseOption::INTEGER, "NUM"));
    setMinValue(parser, "auto-library-num-records", "0");
    setDefaultValue(parser, "auto-library-num-records", "100000");

    addOption(parser, seqan::ArgParseOption("", "fragment-size-median", "Median fragment size.", seqan::ArgParseOption::DOUBLE, "SIZE"));
    setMinValue(parser, "fragment-size-median", "0");
    setDefaultValue(parser, "fragment-size-median", "250");

    addOption(parser, seqan::ArgParseOption("", "fragment-size-std-dev", "Fragment size standard deviation.", seqan::ArgParseOption::DOUBLE, "STDDEV"));
    setMinValue(parser, "fragment-size-std-dev", "0");
    setDefaultValue(parser, "fragment-size-std-dev", "3");

    addOption(parser, seqan::ArgParseOption("", "fragment-default-orientation", "Default orientation.", seqan::ArgParseOption::STRING, "FACTOR"));
    setValidValues(parser, "fragment-default-orientation", "F+ F- R+ R-");
    setDefaultValue(parser, "fragment-default-orientation", "R+");

    // Define BasilOptions -- Coverage Filter
    addSection(parser, "Coverage Filter");
    addOption(parser, seqan::ArgParseOption("", "filter-max-coverage", "Filter out calls at locations that have a higher coverage than this number.  Set to 0 to disable filter.", seqan::ArgParseOption::INTEGER, "COVERAGE"));
    setMinValue(parser, "filter-max-coverage", "0");
    setDefaultValue(parser, "filter-max-coverage", "200");
    addOption(parser, seqan::ArgParseOption("", "filter-min-aln-quality", "Ignore alignments with a quality below this value.", seqan::ArgParseOption::INTEGER, "QUAL"));
    setMinValue(parser, "filter-min-aln-quality", "0");
    setDefaultValue(parser, "filter-min-aln-quality", "0");

    // Define options -- Output Data BasilOptions
    addSection(parser, "Output Data BasilOptions");
    addOption(parser, seqan::ArgParseOption("", "out-variation-type", "The types of variants to write out.", seqan::ArgParseOption::STRING, "TYPE", true));
    // TODO(holtgrew): Enable more types.
    // setValidValues(parser, "out-variant-type", "deletion insertion inversion tandem_duplication translocation");
    setValidValues(parser, "out-variation-type", "insertion");
    setDefaultValue(parser, "out-variation-type", "insertion");
    addOption(parser, seqan::ArgParseOption("", "out-individual-name", "The name of the individual in the output.", seqan::ArgParseOption::STRING, "ID"));
    setDefaultValue(parser, "out-individual-name", "individual");
    // addOption(parser, seqan::ArgParseOption("", "out-min-score", "Minimal score that records must have to be written out.", seqan::ArgParseOption::INTEGER, "NUM"));
    // setDefaultValue(parser, "out-min-score", "1");

    // Define BasilOptions -- Soft-Clipping-Clustering
    addSection(parser, "Clipping Position Clustering");
    addOption(parser, seqan::ArgParseOption("", "clipping-window-radius", "Window radius to use for clipping position clustering.", seqan::ArgParseOption::INTEGER, "RADIUS"));
    setDefaultValue(parser, "clipping-window-radius", "20");
    setMinValue(parser, "clipping-window-radius", "0");
    addOption(parser, seqan::ArgParseOption("", "max-alignment-length", "Maximal alignment length.", seqan::ArgParseOption::INTEGER, "LEN"));
    setDefaultValue(parser, "max-alignment-length", "2000");
    setMinValue(parser, "max-alignment-length", "0");
    addOption(parser, seqan::ArgParseOption("", "clipping-min-length", "Smallest number of characters that have to be soft-clipped such that the alignment is not ignored.", seqan::ArgParseOption::INTEGER, "COUNT"));
    setDefaultValue(parser, "clipping-min-length", "15");
    addOption(parser, seqan::ArgParseOption("", "clipping-min-coverage", "The number of clipping positions to find in one window such that the window is not discarded.", seqan::ArgParseOption::INTEGER, "COUNT"));
    setDefaultValue(parser, "clipping-min-coverage", "5");

    // Define BasilOptions -- OEA Clustering.
    addSection(parser, "One-End Anchor Clustering");

    addOption(parser, seqan::ArgParseOption("", "oea-cluster-selection", "Algorithm for OEA cluster selection.", seqan::ArgParseOption::STRING, "KIND"));
    setDefaultValue(parser, "oea-cluster-selection", "chaining");
    setValidValues(parser, "oea-cluster-selection", "chaining set_cover");

    addOption(parser, seqan::ArgParseOption("", "oea-min-support", "Smallest number of EOA reads to support an insertion.", seqan::ArgParseOption::INTEGER, "COUNT"));
    setDefaultValue(parser, "oea-min-support", "4");
    setMinValue(parser, "oea-min-support", "1");
    addOption(parser, seqan::ArgParseOption("", "oea-min-support-each-side", "Smallest number of EOA reads on each side to support an insertion.", seqan::ArgParseOption::INTEGER, "COUNT"));
    setDefaultValue(parser, "oea-min-support-each-side", "2");
    setMinValue(parser, "oea-min-support-each-side", "1");

    // Define BasilOptions -- Realignment Options
    addSection(parser, "Realignment");

    addOption(parser, seqan::ArgParseOption("", "realignment-num-threads", "Number of threads to use for the realignment.", seqan::ArgParseOption::INTEGER, "COUNT"));
    setDefaultValue(parser, "realignment-num-threads", "1");
    setMinValue(parser, "realignment-num-threads", "1");

    // Define BasilOptions -- Breakpoint Labeling
    addSection(parser, "Breakpoint Labeling");

    addOption(parser, seqan::ArgParseOption("", "breakpoint-window-radius", "Radius around breakpoints to extract for remapping.  Set to 0 to use maximal fragment size.", seqan::ArgParseOption::INTEGER, "RADIUS"));
    setDefaultValue(parser, "breakpoint-window-radius", "0");
    setMinValue(parser, "breakpoint-window-radius", "0");

    // Adding section on library property detection.
    addTextSection(parser, "Library Properties");
    addText(parser,
            "The terms insert size, fragment, and template length all denote the length of the physical fragment "
            "that was extracted and is then sequenced from both sides to yield paired reads.");
    addText(parser,
            "Note that if you set \\fB--fragment-size-mean\\fP or \\fB--fragment-size-std-dev\\fP then you have "
            "to set both.");

    // Adding references section.
    addTextSection(parser, "References");
    addText(parser,
            "Hajirasouliha, I., Hormozdiari, F., Alkan, C., Kidd, J.M., Birol, I., Eichler, E.E., Sahinalp, S.C.  "
            "Detection and characterization of sequence insertions using paired-end next-generation sequencing.  "
            "Bioinformatics 2010 May; 15;26(10):1277-83.");

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    // General BasilOptions
    if (isSet(parser, "quiet"))
        verbosity = 0;
    if (isSet(parser, "verbose"))
        verbosity = 2;
    if (isSet(parser, "very-verbose"))
        verbosity = 3;

    // // Parallelism 
    // getOptionValue(numThreads, parser, "num-threads");
    // getOptionValue(localAlignmentChunkSize, parser, "local-alignment-chunk-size");

    // Input / Output BasilOptions
    getOptionValue(referenceFile, parser, "input-reference");
    getOptionValue(inputFile, parser, "input-mapping");
    getOptionValue(outVcfFile, parser, "out-vcf");
    getOptionValue(outputDebugDir, parser, "output-debug-dir");

    // Library Basil
    getOptionValue(fragmentSizeFactor, parser, "fragment-size-factor");
    getOptionValue(autoLibraryNumRecords, parser, "auto-library-num-records");
    autoLibraryInfo = !(isSet(parser, "fragment-size-median") || isSet(parser, "fragment-size-std-dev"));
    bool allSet = !(isSet(parser, "fragment-size-median") || isSet(parser, "fragment-size-std-dev"));
    if (!autoLibraryInfo && !allSet)
    {
        std::cerr << "ERROR: If you set one of --fragment-size-median and --fragment-size-std-dev "
                  << "then you have to set both!\n";
        return seqan::ArgumentParser::PARSE_ERROR;
    }
    getOptionValue(libraryInfo.median, parser, "fragment-size-median");
    getOptionValue(libraryInfo.stdDev, parser, "fragment-size-std-dev");
    {
        seqan::CharString tmp;
        getOptionValue(tmp, parser, "fragment-default-orientation");
        if (tmp == "R-")
            libraryInfo.defaultOrient = BamLibraryInfo::R_MINUS;
        else if (tmp == "R+")
            libraryInfo.defaultOrient = BamLibraryInfo::R_PLUS;
        else if (tmp == "F-")
            libraryInfo.defaultOrient = BamLibraryInfo::F_MINUS;
        else // if (tmp == "F+")
            libraryInfo.defaultOrient = BamLibraryInfo::F_PLUS;
    }

    // Coverage Filter Basil
    getOptionValue(filterMaxCoverage, parser, "filter-max-coverage");
    getOptionValue(filterMinAlignmentQuality, parser, "filter-min-aln-quality");

    // Get individual name.
    getOptionValue(individualName, parser, "out-individual-name");

    // Soft-Clipping-Clustering Configuration
    getOptionValue(maxAlignmentLength, parser, "max-alignment-length");
    getOptionValue(clippingWindowRadius, parser, "clipping-window-radius");
    getOptionValue(clippingMinLength, parser, "clipping-min-length");
    getOptionValue(clippingMinCoverage, parser, "clipping-min-coverage");

    // OEA Clustering Configuration
    seqan::CharString tmp;
    getOptionValue(tmp, parser, "oea-cluster-selection");
    oeaClusterSelection = (tmp == "set_cover") ? OeaClusterSelection::SET_COVER : OeaClusterSelection::CHAINING;
    getOptionValue(oeaMinSupport, parser, "oea-min-support");
    getOptionValue(oeaMinSupportEachSide, parser, "oea-min-support-each-side");

    // Realignment
    getOptionValue(realignmentNumThreads, parser, "realignment-num-threads");

    // Breakpoint window radius.
    getOptionValue(breakpointWindowRadius, parser, "breakpoint-window-radius");

    return seqan::ArgumentParser::PARSE_OK;
}

// ----------------------------------------------------------------------------
// Class BasilAppImpl
// ----------------------------------------------------------------------------

class BasilAppImpl
{
public:
    BasilOptions options;

    // Entry point for the applications.
    int run(int argc, char const ** argv);

private:
    friend class BasilApp;

    // Prepare run.
    void prepare();

    // Execute and fill members for the results.
    void execute(std::vector<OeaClusterRecord> & oeaClusters,
                 std::vector<ClippingClusterRecord> & clippingClusters);
    
    // Write out the results to disk and quit.
    void finalize(std::vector<std::pair<OeaClusterRecord, ClippingClusterRecord> > const & records);

    // Match clusters and build final result.
    void matchClusters(std::vector<std::pair<OeaClusterRecord, ClippingClusterRecord> > & matchedClusters,
                       std::vector<OeaClusterRecord> const & oeaClusters,
                       std::vector<ClippingClusterRecord> const & clippingClusters);

    // Write resulting matched clusters to VCF file.
    void writeResultVcf(std::vector<std::pair<OeaClusterRecord, ClippingClusterRecord> > const & records);
};

int BasilAppImpl::run(int argc, char const ** argv)
{
    // Parse options, quit on error or if a special option is requested.
    seqan::ArgumentParser::ParseResult res = options.parseCommandLine(argc, argv);
    if (res != seqan::ArgumentParser::PARSE_OK)
        return (res == seqan::ArgumentParser::PARSE_ERROR);

    // Prepare the exeuction.
    prepare();

    // Perform the actual process.
    std::vector<OeaClusterRecord> oeaClusters;
    std::vector<ClippingClusterRecord> clippingClusters;
    execute(oeaClusters, clippingClusters);

    // Match clusters and build final result.
    std::vector<std::pair<OeaClusterRecord, ClippingClusterRecord> > matchedClusters;
    matchClusters(matchedClusters, oeaClusters, clippingClusters);

    // Write out the results.
    finalize(matchedClusters);

    return 0;
}

void BasilAppImpl::matchClusters(std::vector<std::pair<OeaClusterRecord, ClippingClusterRecord> > & matchedClusters,
                                 std::vector<OeaClusterRecord> const & oeaClusters,
                                 std::vector<ClippingClusterRecord> const & clippingClusters)
{
    std::cerr << "\n"
              << "____ MATCHING CLUSTERS _______________________________________________________\n"
              << "\n";

    std::cerr << "OEA CLUSTERS # " << oeaClusters.size() << "\n"
              << "CLIPPING CLUSTERS # " << clippingClusters.size() << "\n"
              << "\n";

    std::cerr << "Running cluster matching..." << std::flush;
    ClusterMatching clusterMatching;
    clusterMatching.run(matchedClusters, oeaClusters, clippingClusters);
    std::cerr << " OK\n";
}

void BasilAppImpl::prepare()
{
    std::cerr << "BASIL\n"
              << "=====\n"
              << "\n";

    BamLibraryEstimator estimator(options.autoLibraryNumRecords);
    if (estimator.run(options.libraryInfo, toCString(options.inputFile)) != 0)
        throw BasilAppException("Could not estimate library info!");

    seqan::FaiIndex faiIndex;
    if (read(faiIndex, toCString(options.referenceFile)) != 0)
    {
        if (options.verbosity >= BasilOptions::NORMAL)
            std::cerr << "Building FAI index...\n";
        if (build(faiIndex, toCString(options.referenceFile)) != 0)
            throw BasilAppException("Problem building FAI index.");
        if (write(faiIndex) != 0)
            throw BasilAppException("Could not write FAI index to file.");
    }

    options.print(std::cerr);

    std::cerr << "\n"
              << "____ SEARCHING FOR SV SIGNATURES _____________________________________________\n"
              << "\n";
}

void BasilAppImpl::execute(std::vector<OeaClusterRecord> & oeaClusters,
                           std::vector<ClippingClusterRecord> & clippingClusters)
{
    std::unique_ptr<ProgressBar> progress;

    BamFilterOptions bamFilterOptions;
    bamFilterOptions.maxCoverage = options.filterMaxCoverage;
    bamFilterOptions.minAlignmentQuality = options.filterMinAlignmentQuality;
    bamFilterOptions.clippingMinLength = options.clippingMinLength;
    bamFilterOptions.realignmentNumThreads = options.realignmentNumThreads;
    bamFilterOptions.realignmentChunkSize = 1024;
    bamFilterOptions.maxAlignmentLength = options.maxAlignmentLength;
    bamFilterOptions.libraryInfo = options.libraryInfo;
    bamFilterOptions.fragmentSizeFactor = options.fragmentSizeFactor;
    bamFilterOptions.refFileName = options.referenceFile;

    BamReaderOptions bamReaderOptions;
    bamReaderOptions.bamFileName = options.inputFile;

    // Maximal size of the buffer queues.
    unsigned const BUFFER_QUEUE_SIZE = 1;

    // Build pipeline for clipping scanning.
    //
    // We have one in and one out queue for both the OEA and clipping cluster algorithm.  We pass the filtered records
    // through both algorithms and at the end we collect the pointers and free the records again.

    typedef ThreadSafeQueue<std::vector<seqan::BamAlignmentRecord *>> TQueue;
    TQueue oeaInQueue(BUFFER_QUEUE_SIZE);
    TQueue clippingInQueue(BUFFER_QUEUE_SIZE);
    TQueue oeaOutQueue(BUFFER_QUEUE_SIZE);
    TQueue clippingOutQueue(BUFFER_QUEUE_SIZE);

    OeaClusterOptions oeaClusterOptions;
    oeaClusterOptions.libraryInfo = options.libraryInfo;
    oeaClusterOptions.fragmentSizeFactor = options.fragmentSizeFactor;
    oeaClusterOptions.debug = (options.verbosity >= BasilOptions::VERBOSE);
    oeaClusterOptions.clusterSelection = options.oeaClusterSelection;
    oeaClusterOptions.oeaMinSupport = options.oeaMinSupport;
    oeaClusterOptions.oeaMinSupportEachSide = options.oeaMinSupportEachSide;
    OeaClusterAlgo oeaClusterAlgo(oeaClusterOptions);
    ClippingClusterOptions clippingClusterOptions;
    clippingClusterOptions.maxAlignmentLength = options.maxAlignmentLength;
    clippingClusterOptions.minCoverageWindowLength = options.clippingWindowRadius;;
    clippingClusterOptions.minCoverage = options.clippingMinCoverage;
    clippingClusterOptions.minClippingLength = options.clippingMinLength;
    ClippingClusterAlgo clippingClusterAlgo(clippingClusterOptions);

    // Pipeline of threads.
    std::vector<std::thread> threads;

    // Pass all buffers through OEA cluster algorithm and push them to the oeaOutQueue.  From there, the records will be
    // freed again (if seen in both oeaOutQueue and clippingOutQueue.
    threads.push_back(std::thread(
            [&]() {
                TQueue::QueueResult res;
                std::vector<seqan::BamAlignmentRecord *> buffer;
                while ((res = oeaInQueue.pop(buffer)) == TQueue::OK)
                {
                    oeaClusterAlgo.push(buffer);
                    oeaOutQueue.push(std::move(buffer));
                }
                if (res == TQueue::CLOSED)
                    oeaOutQueue.close();
                else
                    throw std::runtime_error("Error reading from queue in OEA cluster step.");
            }));

    // Pass all buffers through clipping cluster algorithm and push them to the oeaOutQueue.  From there, the records will be
    // freed again (if seen in both oeaOutQueue and clippingOutQueue.
    threads.push_back(std::thread(
            [&]() {
                TQueue::QueueResult res;
                std::vector<seqan::BamAlignmentRecord *> buffer;
                while ((res = clippingInQueue.pop(buffer)) == TQueue::OK)
                {
                    clippingClusterAlgo.push(buffer);
                    clippingOutQueue.push(std::move(buffer));
                }
                if (res == TQueue::CLOSED)
                    clippingOutQueue.close();
                else
                    throw std::runtime_error("Error reading from queue in clipping cluster step.");
            }));
    
    // Build filter pipline that reads records from BAM file, filters them, and then pushes the resulting records to
    // oeaInQueue and clippingInQueue.
    TQueue queue0;  // initial queue
    // Setup BamReader with progress callback.
    BamReader bamReader(bamReaderOptions);
    // Setup progress indicator.
    unsigned long long const MIB = 1024 * 1024;
    progress.reset(new ProgressBar(std::cerr, 0, fileSize(bamReader.bamStream()) / MIB,
                                   (options.verbosity == BasilOptions::NORMAL)));
    progress->setLabel(toCString(bamReaderOptions.bamFileName));
    progress->updateDisplay();
    auto callback = [&bamReader,&progress,MIB] { progress->advanceTo(positionInFile(bamReader.bamStream()) / MIB); };
    bamReader.setProgressCallback(callback);

    // Build the filter pipeline.
    int BUFFER_SIZE = 2;
    BamFilter bamFilter(bamReader.bamIOContext(), bamFilterOptions);
    Pipeline<std::vector<seqan::BamAlignmentRecord *>> pipeline = bamFilter.makePipeline(BUFFER_SIZE);
    pipeline.start();

    // Read records from BAM File and push them to the queue back.
    threads.push_back(std::thread(
            [&]() {
                std::vector<seqan::BamAlignmentRecord *> buffer;
                while (!bamReader.atEnd())
                {
                    bamReader.read(buffer);
                    pipeline.frontQueue().push(std::move(buffer));  // one copy is enough
                }
                // Close the next pipeline steps.
                pipeline.frontQueue().close();
            }));

    // Consume records and copy them to both the clipping and the OEA cluster algorithm.
    TQueue queue6;  // initial queue
    threads.push_back(std::thread(
            [&]() {
                TQueue::QueueResult res;
                std::vector<seqan::BamAlignmentRecord *> buffer;
                while ((res = pipeline.backQueue().pop(buffer)) == TQueue::OK)
                {
                    // Push to both OEA and clipping cluster queue.
                    // fprintf(stderr, "INTO OEA/CLIPPING IN QUEUE: %d\n", (unsigned)in.size());
                    oeaInQueue.push(buffer);
                    clippingInQueue.push(std::move(buffer));
                }
                // Close the next pipeline queues.
                oeaInQueue.close();
                clippingInQueue.close();
            }));

    // Create a thread that collects all records and frees them.
    std::set<seqan::BamAlignmentRecord *> ptrs;
    bool oeaDone = false, clippingDone = false;
    auto popAndHandle = [&](TQueue & queue, bool & doneFlag)
            {
                std::vector<seqan::BamAlignmentRecord *> el;
                // fprintf(stderr, "el.size() == %d\n", (int)el.size());
                TQueue::QueueResult status = queue.pop(el);
                if (status != TQueue::OK)
                    doneFlag = true;
                else
                    for (auto ptr : el)
                        insertOrErase(ptrs, ptr);
            };
    threads.push_back(
        std::thread([&]() {
            while (!oeaDone || !clippingDone)
            {
                popAndHandle(oeaOutQueue, oeaDone);
                popAndHandle(clippingOutQueue, clippingDone);
                // printf("oea done %i, clipping done %i\n", oeaDone, clippingDone);
            }
        }));
    SEQAN_CHECK(ptrs.empty(), "All records must be freed!");

    // Wait for all threads and the pipeline having finished.
    for (auto & thread : threads)
        thread.join();
    pipeline.join();

    // Finalize clustering objects.
    oeaClusterAlgo.finish();
    clippingClusterAlgo.finish();

    // Finalize progress display.
    progress->finish();

    oeaClusters = oeaClusterAlgo.result();
    clippingClusters = clippingClusterAlgo.result();
}

void BasilAppImpl::finalize(std::vector<std::pair<OeaClusterRecord, ClippingClusterRecord> > const & matchedClusters)
{
    std::cerr << "\n"
              << "____ WRITING RESULTS _________________________________________________________\n"
              << "\n";

    writeResultVcf(matchedClusters);

    std::cerr << "\n"
              << "All done!\n"
              << "Have a nice day!\n";
}

void BasilAppImpl::writeResultVcf(std::vector<std::pair<OeaClusterRecord, ClippingClusterRecord> > const & records)
{
    seqan::VcfStream vcfStream;
    open(vcfStream, toCString(options.outVcfFile), seqan::VcfStream::WRITE);
    if (!isGood(vcfStream))
        throw BasilAppException("Could not open VCF file.");

    seqan::FaiIndex faiIndex;
    if (read(faiIndex, toCString(options.referenceFile)) != 0)
        throw BasilAppException("Could not load FAI index from file.");

    // -----------------------------------------------------------------------
    // Write Header
    // -----------------------------------------------------------------------
    
    appendValue(vcfStream.header.headerRecords, seqan::VcfHeaderRecord("fileformat", "VCFv4.1"));
    appendValue(vcfStream.header.headerRecords, seqan::VcfHeaderRecord("source", "BASIL"));
    appendValue(vcfStream.header.headerRecords, seqan::VcfHeaderRecord("reference", options.referenceFile));
    appendValue(vcfStream.header.headerRecords, seqan::VcfHeaderRecord(
            "INFO", "<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">"));
    appendValue(vcfStream.header.headerRecords, seqan::VcfHeaderRecord(
            "INFO", "<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">"));
    appendValue(vcfStream.header.headerRecords, seqan::VcfHeaderRecord(
            "INFO", "<ID=OEA_ONLY,Number=0,Type=Flag,Description=\"Breakpoint support by OEA signature only\">"));

    appendValue(vcfStream.header.headerRecords, seqan::VcfHeaderRecord(
            "ALT", "<ID=INS,Description=\"Insertion of novel sequence\">"));

    // appendValue(vcfStream.header.headerRecords, seqan::VcfHeaderRecord(
    //         "FORMAT", "<ID=GT,Number=1,Type=String,Description=\"Genotype\">"));
    appendValue(vcfStream.header.headerRecords, seqan::VcfHeaderRecord(
            "FORMAT", "<ID=GSCORE,Number=1,Type=String,Description=\"Sum of Geometric score means (see BASIL documentation)\">"));
    appendValue(vcfStream.header.headerRecords, seqan::VcfHeaderRecord(
            "FORMAT", "<ID=CLEFT,Number=1,Type=String,Description=\"Clipped alignments supporting call from left side.\">"));
    appendValue(vcfStream.header.headerRecords, seqan::VcfHeaderRecord(
            "FORMAT", "<ID=CRIGHT,Number=1,Type=String,Description=\"Clipped alignments supporting call from right side.\">"));
    appendValue(vcfStream.header.headerRecords, seqan::VcfHeaderRecord(
            "FORMAT", "<ID=OEALEFT,Number=1,Type=String,Description=\"One-end anchored alignments supporting call from left side.\">"));
    appendValue(vcfStream.header.headerRecords, seqan::VcfHeaderRecord(
            "FORMAT", "<ID=OEARIGHT,Number=1,Type=String,Description=\"One-end anchored alignments supporting call from right side.\">"));

    // Copy over sequence names.
    for (unsigned i = 0; i < numSeqs(faiIndex); ++i)
    {
        seqan::CharString contigStr = "<ID=";
        append(contigStr, sequenceName(faiIndex, i));
        append(contigStr, ",length=");
        std::stringstream ss;
        ss << sequenceLength(faiIndex, i);
        append(contigStr, ss.str());
        append(contigStr, ">");
        appendValue(vcfStream.header.headerRecords, seqan::VcfHeaderRecord("contig", contigStr));
        appendName(*vcfStream._context.sequenceNames,
                   sequenceName(faiIndex, i),
                   vcfStream._context.sequenceNamesCache);
    }

    appendName(*vcfStream._context.sampleNames, options.individualName);

    // -----------------------------------------------------------------------
    // Write Records
    // -----------------------------------------------------------------------

    unsigned siteId = 0;
    typedef std::vector<std::pair<OeaClusterRecord, ClippingClusterRecord> >::const_iterator TIterator;
    for (TIterator it = records.begin(); it != records.end(); ++it)
    {
        seqan::VcfRecord vcfRecord;
        if (it->second.region.seqId == -1)
        {
            vcfRecord.rID = it->first.region.seqId;
            vcfRecord.beginPos = it->first.region.beginPos + (it->first.region.endPos - it->first.region.beginPos) / 2;
        }
        else
        {
            vcfRecord.rID = it->second.region.seqId;
            vcfRecord.beginPos = it->second.region.beginPos + (it->second.region.endPos - it->second.region.beginPos) / 2;
        }
        if (readRegion(vcfRecord.ref, faiIndex, vcfRecord.rID, vcfRecord.beginPos, vcfRecord.beginPos + 1) != 0)
            throw BasilAppException("Could not read reference character!");

        vcfRecord.alt = "<INS>";
        std::stringstream idSS;
        idSS << "site_" << siteId++;
        vcfRecord.id = idSS.str();
        vcfRecord.filter = "PASS";
        vcfRecord.info = "IMPRECISE;SVTYPE=INS";
        if (it->second.region.seqId == -1)
            append(vcfRecord.info, ";OEA_ONLY");
        vcfRecord.format = "GSCORE:CLEFT:CRIGHT:OEALEFT:OEARIGHT";
        
        appendValue(vcfRecord.genotypeInfos, "");
        std::stringstream ss;
        ss << (it->first.score() + it->second.score()) << ":"
           << it->second.leftWeight << ":"
           << it->second.rightWeight << ":"
           << it->first.leftWeight << ":"
           << it->first.rightWeight;
        append(vcfRecord.genotypeInfos[0], ss.str());

        if (writeRecord(vcfStream, vcfRecord) != 0)
            throw BasilAppException("ERROR: Problem writing to VCF file!");
    }
}

// ----------------------------------------------------------------------------
// Class BasilApp
// ----------------------------------------------------------------------------

BasilApp::BasilApp() : impl(new BasilAppImpl())
{}

BasilApp::~BasilApp()
{}

int BasilApp::run(int argc, char const ** argv)
{
    double startTime = seqan::sysTime();
    int result = impl->run(argc, argv);
    std::cerr << "\nTook " << (seqan::sysTime() - startTime) << " s\n";
    return result;
}
