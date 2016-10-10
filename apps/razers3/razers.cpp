/*==========================================================================
             RazerS - Fast Read Mapping with Controlled Loss Rate
                   http://www.seqan.de/projects/razers.html

 ============================================================================
  Copyright (C) 2008 by David Weese

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your options) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ==========================================================================*/

#define SEQAN_PROFILE                   // enable time measuring
//#define SEQAN_DEBUG_SWIFT				// test SWIFT correctness and print bucket parameters
//#define RAZERS_DEBUG					// print verification regions
#define RAZERS_PRUNE_QGRAM_INDEX        // ignore highly abundant k-mers
//#define RAZERS_MEMOPT					// optimize memory consumption
#define RAZERS_MASK_READS               // remove matches with max-hits optimal hits on-the-fly
//#define NO_PARAM_CHOOSER				// disable loss-rate parameter choosing
#define RAZERS_ISLAND_CRITERION         // island match criterion
#define RAZERS_NOOUTERREADGAPS          // enforce the alignment of the first and last base (determines the lakes)

#define RAZERS_OPENADDRESSING           // enables open addressing for the k-mer index as well as the possibility to set the load factor (-lf)
#define RAZERS_BANDED_MYERS             // uses a banded version of Myers bitvector algorithm (analogous to H. Hyyr\"o, 2001)
//#define SEQAN_OPENADDRESSING_COMPACT	// saves some memory for the openaddressing index / faster hash table access (if undefined)s
//#define RAZERS_DEBUG_MATEPAIRS
#define RAZERS_DEFER_COMPACTION         // mask duplicates on the fly and defer compaction
#define RAZERS_EXTERNAL_MATCHES         // use external memory algorithms for managing matches

//#define RAZERS_PROFILE                // Extensive profiling information.
//#define RAZERS_TIMER					// output information on how fast filtration and verification as well as waiting times
//#define RAZERS_WINDOW					// use the findWindownext function on the "normal" index

#define RAZERS_MATEPAIRS                // enable paired-end matching
//#define SEQAN_USE_SSE2_WORDS			// use SSE2 128-bit integers for MyersBitVector

// Warn the user about missing OpenMP.  This can be suppressed by setting the
// CXX flag "SEQAN_IGNORE_MISSING_OPENMP=1".

#ifdef _OPENMP
#include <omp.h>
#else
#if !defined(SEQAN_IGNORE_MISSING_OPENMP) || (SEQAN_IGNORE_MISSING_OPENMP == 0)
#pragma message("OpenMP not found! Shared-memory parallelization will be disabled in RazerS3.")
#endif  // #if !defined(SEQAN_IGNORE_MISSING_OPENMP) || (SEQAN_IGNORE_MISSING_OPENMP == 0)
#endif

#include <iostream>
#include <sstream>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/store.h>
#include <seqan/arg_parse.h>
#include <seqan/parallel.h>

#ifdef STDLIB_VS
#include <process.h>
#endif

#include <memory>
#include <unordered_map>

#include "razers.h"
#include "outputFormat.h"
#include "paramChooser.h"

#ifdef RAZERS_WINDOW
#include "razers_window.h"
#endif

#ifdef RAZERS_MATEPAIRS
#include "razers_matepairs.h"
#ifdef _OPENMP
#include "razers_matepairs_parallel.h"
#endif
#endif

#ifdef _OPENMP
#include "razers_parallel.h"
#endif

#ifdef RAZERS_PROFILE
#include "profile_timeline.h"
#endif  // #ifdef RAZERS_PROFILE

using namespace std;
using namespace seqan;

struct MyFragStoreConfig :
    public FragmentStoreConfig<>
{
    typedef String<Dna5> TContigSeq;
};

//////////////////////////////////////////////////////////////////////////////
// Main read mapper function
template <typename TSpec>
int mapReads(
    StringSet<CharString> & genomeFileNames,
    StringSet<CharString> & readFileNames,  // NULL terminated list of one/two read files (single/mate-pairs)
    RazerSOptions<TSpec> & options)
{
    FragmentStore<MyFragStoreConfig>    store;          // stores all of the tables
    String<String<unsigned short> >     stats;      // needed for mapping quality calculation

    // dump configuration in verbose mode
    if (options._debugLevel >= 1)
    {
        CharString bitmap;
        Shape<Dna, GenericShape> shape;
        stringToShape(shape, options.shape);
        shapeToString(bitmap, shape);

        cerr << "___SETTINGS____________" << endl;
        cerr << "Genome file:                     \t" << genomeFileNames[0] << endl;
        if (length(readFileNames) < 2)
            cerr << "Read file:                       \t" << readFileNames[0] << endl;
        else
        {
            cerr << "Read files:                      \t" << readFileNames[0] << endl;
            for (unsigned i = 1; i < length(readFileNames); ++i)
                cerr << "                                 \t" << readFileNames[i] << endl;
        }
        cerr << "Compute forward matches:         \t";
        if (options.forward)
            cerr << "YES" << endl;
        else
            cerr << "NO" << endl;
        cerr << "Compute reverse matches:         \t";
        if (options.reverse)
            cerr << "YES" << endl;
        else
            cerr << "NO" << endl;
        cerr << "Allow Indels:                    \t";
        if (options.gapMode == RAZERS_GAPPED)
            cerr << "YES" << endl;
        else
            cerr << "NO" << endl;
        cerr << "Error rate:                      \t" << options.errorRate << endl;
        if (options.threshold > 0)
            cerr << "Minimal threshold:               \t" << options.threshold << endl;
        else
            cerr << "Pigeonhole mode with overlap:    \t" << options.overlap << endl;
        cerr << "Shape:                           \t" << bitmap << endl;
        cerr << "Repeat threshold:                \t" << options.repeatLength << endl;
        cerr << "Overabundance threshold:         \t" << options.abundanceCut << endl;
        if (options.threshold > 0)
            cerr << "Taboo length:                    \t" << options.tabooLength << endl;
        if (options._debugLevel >= 1)
        {
#ifdef STDLIB_VS
            int pid = _getpid();
#else // #ifdef STDLIB_VS
            int pid = getpid();
#endif // #ifdef STDLIB_VS
            cerr << "Program PID:                     \t" << pid << endl;
        }
        cerr << endl;
    }

    // circumvent numerical obstacles
    options.errorRate += 0.0000001;

#ifdef RAZERS_PROFILE
    timelineBeginTask(TASK_LOAD);
#endif  // #ifdef RAZERS_PROFILE
    //////////////////////////////////////////////////////////////////////////////
    // Step 1: Load reads
    SEQAN_PROTIMESTART(load_time);

#ifdef RAZERS_MATEPAIRS
    if (length(readFileNames) == 2)
    {
        if (!loadReads(store, options.readFile, toCString(readFileNames[1]), options))
        {
            //if (!loadReads(readSet, readQualities, readNames, readFileNames[0], readFileNames[1], options)) {
            cerr << "Failed to load reads" << endl;
            return RAZERS_READS_FAILED;
        }
    }
    else
#endif
    {
        if (!loadReads(store, options.readFile, options))
        {
            cerr << "Failed to load reads" << endl;
            return RAZERS_READS_FAILED;
        }
    }

    if (options._debugLevel >= 1)
        cerr << lengthSum(store.readSeqStore) << " bps of " << length(store.readSeqStore) << " reads loaded." << endl;
    options.timeLoadFiles = SEQAN_PROTIMEDIFF(load_time);

    if (options._debugLevel >= 1)
        cerr << "Loading reads took               \t" << options.timeLoadFiles << " seconds" << endl;

    #ifdef RAZERS_MEMOPT
    if (length(store.readSeqStore) > 16777216)
    {
        cerr << "more than 2^24 reads. Switch of RAZERS_MEMOPT in razers.cpp or use less." << std::endl;
        return 1;
    }
    #endif

    //////////////////////////////////////////////////////////////////////////////
    // Step 2: Load genomes
    if (length(genomeFileNames) == 1)
    {
        int result = getGenomeFileNameList(genomeFileNames[0], genomeFileNames, options);
        if (result == RAZERS_GENOME_FAILED)
        {
            cerr << "Failed to open genome file " << genomeFileNames[0] << endl;
            return result;
        }
    }
#ifdef RAZERS_PROFILE
    timelineEndTask(TASK_LOAD);
#endif  // #ifdef RAZERS_PROFILE

    //////////////////////////////////////////////////////////////////////////////
    // Step 3: Find matches using SWIFT
    loadContigs(store, genomeFileNames, false); // add filenames to the contig store (they are loaded on-demand)
    int error = _mapReads(store, stats, (RazerSCoreOptions<TSpec>&)options);
    if (error != 0)
    {
        switch (error)
        {
        case RAZERS_GENOME_FAILED:
            cerr << "Failed to load genomes" << endl;
            break;

        case RAZERS_INVALID_SHAPE:
            cerr << "Invalid Shape" << endl;
            break;
        }
        return error;
    }

    //////////////////////////////////////////////////////////////////////////////
    // Step 4: Remove duplicates and output matches
#ifdef RAZERS_PROFILE
    timelineBeginTask(TASK_DUMP_MATCHES);
#endif  // #ifdef RAZERS_PROFILE
    if (!options.spec.DONT_DUMP_RESULTS)
        dumpMatches(store, stats, readFileNames[0], options);
#ifdef RAZERS_PROFILE
    timelineEndTask(TASK_DUMP_MATCHES);
#endif  // #ifdef RAZERS_PROFILE

    return 0;
}

inline void whichMacros()
{
#ifdef RAZERS_OPENADDRESSING
    std::cerr << "Index:    Open addressing" << std::endl;
#else
    std::cerr << "Index:    Normal" << std::endl;
#endif

#ifdef RAZERS_TIMER
    std::cerr << "Timer:    ON" << std::endl;
#else
    std::cerr << "Timer:    OFF" << std::endl;
#endif

#ifdef _OPENMP
    std::cerr << "OpenMP:   ON" << std::endl;
#else
    std::cerr << "OpenMP:   OFF" << std::endl;
#endif

#ifdef RAZERS_BANDED_MYERS
    std::cerr << "Myers:    Banded" << std::endl;
#else
    std::cerr << "Myers:    Unbanded" << std::endl;
#endif

#ifdef RAZERS_PROFILE
    std::cerr << "Timeline: ON" << std::endl;
#else
    std::cerr << "Timeline: OFF" << std::endl;
#endif
    std::cerr << std::endl;
}

void setUpArgumentParser(ArgumentParser & parser, RazerSOptions<> & options, ParamChooserOptions const & pm_options)
{
    setAppName(parser, "razers3");
    setShortDescription(parser, "Faster, fully sensitive read mapping");
    setCategory(parser, "Read Mapping");
    options.version = "3.3";
	setVersion(parser, SEQAN_APP_VERSION " [" SEQAN_REVISION "]");
    setDate(parser, SEQAN_DATE);

    // Need genome and reads (hg18.fa reads.fq)
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE));
    setValidValues(parser, 0, seqan::SeqFileIn::getFileExtensions());
    setHelpText(parser, 0, "A reference genome file.");
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "READS", true));
    setValidValues(parser, 1, seqan::SeqFileIn::getFileExtensions());
    setHelpText(parser, 1, "Either one (single-end) or two (paired-end) read files.");

    addUsageLine(parser, "[\\fIOPTIONS\\fP] <\\fIGENOME FILE\\fP> <\\fIREADS FILE\\fP>");
#ifdef RAZERS_MATEPAIRS
    addUsageLine(parser, "[\\fIOPTIONS\\fP] <\\fIGENOME FILE\\fP> <\\fIPE-READS FILE1\\fP> <\\fIPE-READS FILE2\\fP>");
#endif

    addDescription(parser, "RazerS 3 is a versatile full-sensitive read mapper based on k-mer counting and seeding filters. "
                           "It supports single and paired-end mapping, shared-memory parallelism, and optimally parametrizes "
                           "the filter based on a user-defined minimal sensitivity. "
                           "See \\fIhttp://www.seqan.de/projects/razers\\fP for more information.");

    addDescription(parser, "Input to RazerS 3 is a reference genome file and either one file with single-end reads "
                           "or two files containing left or right mates of paired-end reads. Use - to read single-end "
                           "reads from stdin.");

    addDescription(parser, "(c) Copyright 2009-2014 by David Weese.");

    addSection(parser, "Main Options");
    addOption(parser, ArgParseOption("i", "percent-identity", "Percent identity threshold.", ArgParseOption::DOUBLE));
    setMinValue(parser, "percent-identity", "50");
    setMaxValue(parser, "percent-identity", "100");
    setDefaultValue(parser, "percent-identity", 100 - (100.0 * options.errorRate));
    addOption(parser, ArgParseOption("rr", "recognition-rate", "Percent recognition rate.", ArgParseOption::DOUBLE));
    setMinValue(parser, "recognition-rate", "80");
    setMaxValue(parser, "recognition-rate", "100");
    setDefaultValue(parser, "recognition-rate", 100 - (100.0 * pm_options.optionLossRate));
    addOption(parser, ArgParseOption("ng", "no-gaps", "Allow only mismatches, no indels. Default: allow both."));
    addOption(parser, ArgParseOption("f", "forward", "Map reads only to forward strands."));
    addOption(parser, ArgParseOption("r", "reverse", "Map reads only to reverse strands."));
    addOption(parser, ArgParseOption("m", "max-hits", "Output only <\\fINUM\\fP> of the best hits.", ArgParseOption::INTEGER));
    setMinValue(parser, "max-hits", "1");
    setDefaultValue(parser, "max-hits", options.maxHits);
    addOption(parser, ArgParseOption("", "unique", "Output only unique best matches (-m 1 -dr 0 -pa)."));
    addOption(parser, ArgParseOption("tr", "trim-reads", "Trim reads to given length. Default: off.", ArgParseOption::INTEGER));
    setMinValue(parser, "trim-reads", "14");
    addOption(parser, ArgParseOption("o", "output", "Mapping result filename (use - to dump to stdout in razers format). Default: <\\fIREADS FILE\\fP>.razers.", ArgParseOption::OUTPUT_FILE));
    setValidValues(parser, "output", ".razers .eland .fa .fasta .gff .sam .bam .afg");
    addOption(parser, ArgParseOption("v", "verbose", "Verbose mode."));
    addOption(parser, ArgParseOption("vv", "vverbose", "Very verbose mode."));

#ifdef RAZERS_MATEPAIRS
    addSection(parser, "Paired-end Options");
    addOption(parser, ArgParseOption("ll", "library-length", "Paired-end library length.", ArgParseOption::INTEGER));
    setMinValue(parser, "library-length", "1");
    setDefaultValue(parser, "library-length", options.libraryLength);
    addOption(parser, ArgParseOption("le", "library-error", "Paired-end library length tolerance.", ArgParseOption::INTEGER));
    setMinValue(parser, "library-error", "0");
    setDefaultValue(parser, "library-error", options.libraryError);
#endif

    addSection(parser, "Output Format Options");
    addOption(parser, ArgParseOption("a", "alignment", "Dump the alignment for each match (only \\fIrazer\\fP or \\fIfasta\\fP format)."));
    addOption(parser, ArgParseOption("pa", "purge-ambiguous", "Purge reads with more than <\\fImax-hits\\fP> best matches."));
    addOption(parser, ArgParseOption("dr", "distance-range", "Only consider matches with at most NUM more errors compared to the best. Default: output all.", ArgParseOption::INTEGER));
    addOption(parser, ArgParseOption("gn", "genome-naming", "Select how genomes are named (see Naming section below).", ArgParseOption::INTEGER));
    setMinValue(parser, "genome-naming", "0");
    setMaxValue(parser, "genome-naming", "1");
    setDefaultValue(parser, "genome-naming", options.genomeNaming);
//	addHelpLine(parser, "0 = use Fasta id");
//	addHelpLine(parser, "1 = enumerate beginning with 1");
    addOption(parser, ArgParseOption("rn", "read-naming", "Select how reads are named (see Naming section below).", ArgParseOption::INTEGER));
    setDefaultValue(parser, "read-naming", options.readNaming);
    setMinValue(parser, "read-naming", "0");
    setMaxValue(parser, "read-naming", "3");

//	addHelpLine(parser, "0 = use Fasta id");
//	addHelpLine(parser, "1 = enumerate beginning with 1");
//	addHelpLine(parser, "2 = use the read sequence (only for short reads!)");
//	addHelpLine(parser, "3 = use the Fasta id, do NOT append '/L' or '/R' for mate pairs");
    addOption(parser, ArgParseOption("", "full-readid", "Use the whole read id (don't clip after whitespace)."));
    addOption(parser, ArgParseOption("so", "sort-order", "Select how matches are sorted (see Sorting section below).", ArgParseOption::INTEGER));
    setDefaultValue(parser, "sort-order", options.sortOrder);
    setMinValue(parser, "sort-order", "0");
    setMaxValue(parser, "sort-order", "1");
//	addHelpLine(parser, "0 = 1. read number, 2. genome position");
//	addHelpLine(parser, "1 = 1. genome position, 2. read number");
    addOption(parser, ArgParseOption("pf", "position-format", "Select begin/end position numbering (see Coordinate section below).", ArgParseOption::INTEGER));
    setMinValue(parser, "position-format", "0");
    setMaxValue(parser, "position-format", "1");
    setDefaultValue(parser, "position-format", options.sortOrder);
//	addHelpLine(parser, "0 = gap space");
//	addHelpLine(parser, "1 = position space");
    addOption(parser, ArgParseOption("ds", "dont-shrink-alignments", "Disable alignment shrinking in SAM.  This is required for generating a gold mapping for Rabema."));
    addOption(parser, ArgParseOption("ga", "global-alignment", "Compute global alignment (in SAM output)."));
    hideOption(parser, "global-alignment");

    addSection(parser, "Filtration Options");
    addOption(parser, ArgParseOption("fl", "filter", "Select k-mer filter.", ArgParseOption::STRING));
    setValidValues(parser, "filter", "pigeonhole swift");
    setDefaultValue(parser, "filter", "pigeonhole");
    addOption(parser, ArgParseOption("mr", "mutation-rate", "Set the percent mutation rate (\\fIpigeonhole\\fP).", ArgParseOption::DOUBLE));
    setMinValue(parser, "mutation-rate", "0");
    setMaxValue(parser, "mutation-rate", "20");
    setDefaultValue(parser, "mutation-rate", 100.0 * options.mutationRate);
    addOption(parser, ArgParseOption("ol", "overlap-length", "Manually set the overlap length of adjacent k-mers (\\fIpigeonhole\\fP).", ArgParseOption::INTEGER));
    setMinValue(parser, "overlap-length", "0");
#ifndef NO_PARAM_CHOOSER
    addOption(parser, ArgParseOption("pd", "param-dir", "Read user-computed parameter files in the directory <\\fIDIR\\fP> (\\fIswift\\fP).", ArgParseOption::STRING, "DIR"));
#endif
    addOption(parser, ArgParseOption("t", "threshold", "Manually set minimum k-mer count threshold (\\fIswift\\fP).", ArgParseOption::INTEGER));
    setMinValue(parser, "threshold", "1");
    addOption(parser, ArgParseOption("tl", "taboo-length", "Set taboo length (\\fIswift\\fP).", ArgParseOption::INTEGER));
    setMinValue(parser, "taboo-length", "1");
    setDefaultValue(parser, "taboo-length", options.tabooLength);
    addOption(parser, ArgParseOption("s", "shape", "Manually set k-mer shape.", ArgParseOption::STRING, "BITSTRING"));
//    setDefaultValue(parser, "shape", options.shape);  // <-- doesn't work with KNIME which always sets default values
    addOption(parser, ArgParseOption("oc", "overabundance-cut", "Set k-mer overabundance cut ratio.", ArgParseOption::INTEGER));
    setMinValue(parser, "overabundance-cut", "0");
    setMaxValue(parser, "overabundance-cut", "1");
    setDefaultValue(parser, "overabundance-cut", options.abundanceCut);
    addOption(parser, ArgParseOption("rl", "repeat-length", "Skip simple-repeats of length <\\fINUM\\fP>.", ArgParseOption::INTEGER));
    setMinValue(parser, "repeat-length", "1");
    setDefaultValue(parser, "repeat-length", options.repeatLength);
#ifdef RAZERS_OPENADDRESSING
    addOption(parser, ArgParseOption("lf", "load-factor", "Set the load factor for the open addressing k-mer index.", ArgParseOption::DOUBLE));
    setMinValue(parser, "load-factor", "1");
    setDefaultValue(parser, "load-factor", options.loadFactor);
#endif

    addSection(parser, "Verification Options");
    addOption(parser, ArgParseOption("mN", "match-N", "N matches all other characters. Default: N matches nothing."));
    addOption(parser, ArgParseOption("ed", "error-distr", "Write error distribution to \\fIFILE\\fP.", ArgParseOption::STRING, "FILE"));
    addOption(parser, ArgParseOption("mf", "mismatch-file", "Write mismatch patterns to \\fIFILE\\fP.", ArgParseOption::STRING, "FILE"));

    addSection(parser, "Misc Options");
    addOption(parser, ArgParseOption("cm", "compact-mult", "Multiply compaction threshold by this value after reaching and compacting.", ArgParseOption::DOUBLE));
    setMinValue(parser, "compact-mult", "0");
    setDefaultValue(parser, "compact-mult", options.compactMult);
    addOption(parser, ArgParseOption("ncf", "no-compact-frac", "Don't compact if in this last fraction of genome.", ArgParseOption::DOUBLE));
    setMinValue(parser, "no-compact-frac", "0");
    setMaxValue(parser, "no-compact-frac", "1");
    setDefaultValue(parser, "no-compact-frac", options.noCompactFrac);

#ifdef _OPENMP
    addSection(parser, "Parallelism Options");
#endif  // #ifdef _OPENMP
    addOption(parser, ArgParseOption("tc", "thread-count", "Set the number of threads to use (0 to force sequential mode).", ArgParseOption::INTEGER));
    setMinValue(parser, "thread-count", "0");
#ifndef _OPENMP
    hideOption(parser, "tc");
#endif  // #ifndef _OPENMP
#ifdef _OPENMP
    setDefaultValue(parser, "thread-count", options.threadCount);
    addOption(parser, ArgParseOption("pws", "parallel-window-size", "Collect candidates in windows of this length.", ArgParseOption::INTEGER));
    setMinValue(parser, "parallel-window-size", "1");
    setDefaultValue(parser, "parallel-window-size", options.windowSize);
    addOption(parser, ArgParseOption("pvs", "parallel-verification-size", "Verify candidates in packages of this size.", ArgParseOption::INTEGER));
    setMinValue(parser, "parallel-verification-size", "1");
    setDefaultValue(parser, "parallel-verification-size", options.verificationPackageSize);
    addOption(parser, ArgParseOption("pvmpc", "parallel-verification-max-package-count", "Largest number of packages to create for verification per thread-1.", ArgParseOption::INTEGER));
    setMinValue(parser, "parallel-verification-max-package-count", "1");
    setDefaultValue(parser, "parallel-verification-max-package-count", options.maxVerificationPackageCount);
//    addHelpLine(parser, "Go over package size if this limit is reached.");
    addOption(parser, ArgParseOption("amms", "available-matches-memory-size", "Bytes of main memory available for storing matches.", ArgParseOption::INTEGER));
    setMinValue(parser, "available-matches-memory-size", "-1");
    setDefaultValue(parser, "available-matches-memory-size", options.availableMatchesMemorySize);
//    addHelpLine(parser, "Used to switch to external sorting.");
//    addHelpLine(parser, "-1 = always external");
//    addHelpLine(parser, " 0 = never");
//    addHelpLine(parser, " x = use other value x as threshold.");
    addOption(parser, ArgParseOption("mhst", "match-histo-start-threshold", "When to start histogram.", ArgParseOption::INTEGER));
    setMinValue(parser, "match-histo-start-threshold", "1");
    setDefaultValue(parser, "match-histo-start-threshold", options.matchHistoStartThreshold);
#endif

    addTextSection(parser, "Formats, Naming, Sorting, and Coordinate Schemes");

    addText(parser, "RazerS 3 supports various output formats. The output format is detected automatically from the file name suffix.");
	addListItem(parser, ".razers", "Razer format");
	addListItem(parser, ".fa, .fasta", "Enhanced Fasta format");
	addListItem(parser, ".eland", "Eland format");
	addListItem(parser, ".gff", "GFF format");
	addListItem(parser, ".sam", "SAM format");
	addListItem(parser, ".bam", "BAM format");
	addListItem(parser, ".afg", "Amos AFG format");

    addText(parser, "By default, reads and contigs are referred by their Fasta ids given in the input files. "
                    "With the \\fB-gn\\fP and \\fB-rn\\fP options this behaviour can be changed:");
    addListItem(parser, "0", "Use Fasta id.");
    addListItem(parser, "1", "Enumerate beginning with 1.");
    addListItem(parser, "2", "Use the read sequence (only for short reads!).");
    addListItem(parser, "3", "Use the Fasta id, do NOT append /L or /R for mate pairs.");

    addText(parser, "");
    addText(parser, "The way matches are sorted in the output file can be changed with the \\fB-so\\fP option for the following formats: "
                    "\\fBrazers\\fP, \\fBfasta\\fP, \\fBsam\\fP, and \\fBafg\\fP. Primary and secondary sort keys are:");
    addListItem(parser, "0", "1. read number, 2. genome position");
    addListItem(parser, "1", "1. genome position, 2. read number");

    addText(parser, "");
    addText(parser, "The coordinate space used for begin and end positions can be changed with the "
                    "\\fB-pf\\fP option for the \\fBrazer\\fP and \\fBfasta\\fP formats:");
    addListItem(parser, "0", "Gap space. Gaps between characters are counted from 0.");
    addListItem(parser, "1", "Position space. Characters are counted from 1.");


    addTextSection(parser, "Examples");
    addListItem(parser,
                "\\fBrazers3\\fP \\fB-i\\fP \\fB96\\fP \\fB-tc\\fP \\fB12\\fP \\fB-o\\fP \\fBmapped.razers\\fP \\fBhg18.fa\\fP \\fBreads.fq\\fP",
                "Map single-end reads with 4% error rate using 12 threads.");
    addListItem(parser,
                "\\fBrazers3\\fP \\fB-i\\fP \\fB95\\fP \\fB-no-gaps\\fP \\fB-o\\fP \\fBmapped.razers\\fP \\fBhg18.fa\\fP \\fBreads.fq.gz\\fP",
                "Map single-end gzipped reads with 5% error rate and no indels.");
    addListItem(parser,
                "\\fBrazers3\\fP \\fB-i\\fP \\fB94\\fP \\fB-rr\\fP \\fB95\\fP \\fB-tc\\fP \\fB12\\fP \\fB-ll\\fP \\fB280\\fP \\fB--le\\fP \\fB80\\fP \\fB-o\\fP \\fBmapped.razers\\fP \\fBhg18.fa\\fP \\fBreads_1.fq\\fP \\fBreads_2.fq\\fP",
                "Map paired-end reads with up to 6% errors, 95% sensitivity, 12 threads, and only output aligned pairs with an outer distance of 200-360bp.");
}

ArgumentParser::ParseResult
extractOptions(
    StringSet<CharString> & genomeFileNames, StringSet<CharString> & readFileNames,
    RazerSOptions<> & options, ParamChooserOptions & pm_options,
    ArgumentParser const & parser)
{
    //////////////////////////////////////////////////////////////////////////////
    // Extract options

    bool stop = false;
    getOptionValue(options.forward, parser, "forward");
    getOptionValue(options.reverse, parser, "reverse");
    getOptionValue(options.errorRate, parser, "percent-identity");
    getOptionValue(options.mutationRate, parser, "mutation-rate");
#ifndef NO_PARAM_CHOOSER
    getOptionValue(pm_options.optionLossRate, parser, "recognition-rate");
    getOptionValue(pm_options.paramFolder, parser, "param-dir");
    // append slash/backslash
    if (!empty(pm_options.paramFolder))
    {
        if (back(pm_options.paramFolder) != '/' && back(pm_options.paramFolder) != '\\')
        {
#ifdef STDLIB_VS
            appendValue(pm_options.paramFolder, '\\');
#else
            appendValue(pm_options.paramFolder, '/');
#endif
        }
    }
#endif
    options.gapMode = (isSet(parser, "no-gaps")) ? RAZERS_UNGAPPED : RAZERS_GAPPED;
#ifdef RAZERS_MATEPAIRS
    getOptionValue(options.libraryLength, parser, "library-length");
    getOptionValue(options.libraryError, parser, "library-error");
#endif
    getOptionValue(options.maxHits, parser, "max-hits");
    getOptionValue(options.purgeAmbiguous, parser, "purge-ambiguous");
    getOptionValue(options.scoreDistanceRange, parser, "distance-range");
    if (isSet(parser, "distance-range"))
        options.scoreDistanceRange++;
    getOptionValue(options.dumpAlignment, parser, "alignment");

    getOptionValue(options.sortOrder, parser, "sort-order");
    getOptionValue(options.genomeNaming, parser, "genome-naming");
    getOptionValue(options.readNaming, parser, "read-naming");
    getOptionValue(options.fullFastaId, parser, "full-readid");
    getOptionValue(options.positionFormat, parser, "position-format");
    getOptionValue(options.compactMult, parser, "compact-mult");
    getOptionValue(options.noCompactFrac, parser, "no-compact-frac");
    getOptionValue(options.dontShrinkAlignments, parser, "dont-shrink-alignments");
    getOptionValue(options.computeGlobal, parser, "global-alignment");
    getOptionValue(options.shape, parser, "shape");
    getOptionValue(options.overlap, parser, "overlap-length");
    getOptionValue(options.abundanceCut, parser, "overabundance-cut");
    getOptionValue(options.repeatLength, parser, "repeat-length");

#ifdef _OPENMP
    getOptionValue(options.threadCount, parser, "thread-count");
    getOptionValue(options.windowSize, parser, "parallel-window-size");
    getOptionValue(options.verificationPackageSize, parser, "parallel-verification-size");
    getOptionValue(options.maxVerificationPackageCount, parser, "parallel-verification-max-package-count");
    getOptionValue(options.availableMatchesMemorySize, parser, "available-matches-memory-size");
    getOptionValue(options.matchHistoStartThreshold, parser, "match-histo-start-threshold");
#else
    options.threadCount = 0;
#endif

#ifdef RAZERS_OPENADDRESSING
    getOptionValue(options.loadFactor, parser, "load-factor");
#endif
    getOptionValue(options.trimLength, parser, "trim-reads");
    getOptionValue(options.tabooLength, parser, "taboo-length");
    getOptionValue(options.matchN, parser, "match-N");
    getOptionValue(options.errorPrbFileName, parser, "error-distr");
    getOptionValue(options.mismatchFilename, parser, "mismatch-file");

    if (isSet(parser, "verbose"))
        options._debugLevel = max(options._debugLevel, 1);
    if (isSet(parser, "vverbose"))
        options._debugLevel = max(options._debugLevel, 3);
    if (isSet(parser, "unique"))
    {
        options.maxHits = 1;
        options.scoreDistanceRange = 1;
        options.purgeAmbiguous = true;
    }
    if (!options.forward && !options.reverse)  // enable both per default
    {
        options.forward = true;
        options.reverse = true;
    }

#ifdef RAZERS_MATEPAIRS
    unsigned maxReadFiles = 2;
#else
    unsigned maxReadFiles = 1;
#endif
    resize(genomeFileNames, length(genomeFileNames) + 1);
    getArgumentValue(back(genomeFileNames), parser, 0);
    resize(readFileNames, _min(getArgumentValueCount(parser, 1), maxReadFiles), Exact());
    for (unsigned i = 0; i < length(readFileNames); ++i)
        getArgumentValue(readFileNames[i], parser, 1, i);

    // Get output file name from command line if set.  Otherwise, autogenerate from input file name.
    if (isSet(parser, "output"))
    {
        getOptionValue(options.output, parser, "output");
    }
    else
    {
        options.output = readFileNames[0];
        append(options.output, ".razers");
    }

    // Get lower case of the output file name.  File endings are accepted in both upper and lower case.
    CharString tmp = options.output;
    toLower(tmp);

    if (endsWith(tmp, ".fa") || endsWith(tmp, ".fasta"))
        options.outputFormat = 1;
    else if (endsWith(tmp, ".eland"))
        options.outputFormat = 2;
    else if (endsWith(tmp, ".gff"))
        options.outputFormat = 3;
    else if (endsWith(tmp, ".sam") || endsWith(tmp, ".bam"))
        options.outputFormat = 4;
    else if (endsWith(tmp, ".afg"))
        options.outputFormat = 5;
    else
        options.outputFormat = 0;   // default is ".razers"

    // don't append /L/R in SAM mode
    if (!isSet(parser, "read-naming") && options.outputFormat == 4)
        options.readNaming = 3;

    CharString filter;
    getOptionValue(filter, parser, "filter");
    if (filter == "pigeonhole")
    {
        if (isSet(parser, "threshold") && (stop = true))
            cerr << "k-mer threshold can only be set for the swift filter (-fl swift)" << endl;
        options.threshold = 0;
    }
    else
    {
        getOptionValue(options.threshold, parser, "threshold");
        if (isSet(parser, "overlap-length") && (stop = true))
            cerr << "k-mer overlap length can only be set for the pigeonhole filter (-fl pigeonhole)" << endl;
    }

    if (isSet(parser, "shape"))
    {
        unsigned ones = 0;
        unsigned zeros = 0;
        for (unsigned i = 0; i < length(options.shape); ++i)
            switch (options.shape[i])
            {
            case '0':
                ++zeros;
                break;

            case '1':
                ++ones;
                break;

            default:
                cerr << "Shape must be a binary string" << endl;
                stop = true;
                i = length(options.shape);
            }
        if ((ones == 0 || ones > 31) && !stop)
        {
            cerr << "Invalid Shape" << endl;
            stop = true;
        }

        unsigned maxOnes = 14;
#ifdef RAZERS_OPENADDRESSING
        maxOnes = 31;
#endif
        if ((ones < 7 || ones > maxOnes) && !stop)
            cerr << "Warning: Shape should contain at least 7 and at most " << maxOnes << " '1's" << endl;
        options.delta = ones + zeros;
    }
    if (getArgumentValueCount(parser, 1) == 1)
        options.libraryLength = -1;     // only 1 readset -> disable mate-pair mapping
    if ((getArgumentValueCount(parser, 1) > maxReadFiles) && (stop = true))
        cerr << "More than " << maxReadFiles << " read files specified." << endl;
    if ((getArgumentValueCount(parser, 1) == 0) && (stop = true))
        cerr << "No read files specified." << endl;

    options.compMask[4] = (options.matchN) ? 15 : 0;
    options.errorRate = (100.0 - options.errorRate) / 100.0;
    options.mutationRate = options.mutationRate / 100.0;
    options.lossRate = pm_options.optionLossRate = (100.0 - pm_options.optionLossRate) / 100.0;

    return (stop) ? ArgumentParser::PARSE_ERROR : ArgumentParser::PARSE_OK;
}

//////////////////////////////////////////////////////////////////////////////
// Command line parsing and parameter choosing
int main(int argc, const char * argv[])
{
    //whichMacros();

#ifdef RAZERS_PROFILE
    initTimeline();
    unsigned x = timelineAddTaskType("ON_CONTIG", "Work on contig.");
    (void)x;  // Disable warning if assertions are disable.
    SEQAN_ASSERT_EQ(x, 1u);  // The following will be OK, too.
    timelineAddTaskType("INIT", "Initialization.");
    timelineAddTaskType("REVCOMP", "Reverse-complementing contig.");
    timelineAddTaskType("FILTER", "Filtration using SWIFT.");
    timelineAddTaskType("VERIFY", "Verification of SWIFT hits.");
    timelineAddTaskType("WRITEBACK", "Write back to block-local store.");
    timelineAddTaskType("COMPACT", "Compaction");
    timelineAddTaskType("DUMP_MATCHES", "Dump matches.");
    timelineAddTaskType("LOAD", "Load input.");
    timelineAddTaskType("SORT", "Sorting.");
    timelineAddTaskType("COPY_FINDER", "Copy SWIFT Finder.");
#endif  // #ifndef RAZERS_PROFILE

    RazerSOptions<>         options;
    ParamChooserOptions     pm_options;
    StringSet<CharString>   genomeFileNames;
    StringSet<CharString>   readFileNames;

    // Change defaults
    options.forward = false;
    options.reverse = false;

    // Set up command line parser.
    ArgumentParser argParser;
    setUpArgumentParser(argParser, options, pm_options);

    // Parse command line.
    ArgumentParser::ParseResult res = parse(argParser, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
    {
        if (res == ArgumentParser::PARSE_ERROR)
            cerr << "Exiting ..." << endl;
        return (res == ArgumentParser::PARSE_ERROR) ? RAZERS_INVALID_OPTIONS : 0;
    }
    // Extract options.
    res = extractOptions(genomeFileNames, readFileNames, options, pm_options, argParser);
    if (res != ArgumentParser::PARSE_OK)
    {
        cerr << "Exiting ..." << endl;
        return RAZERS_INVALID_OPTIONS;
    }

#ifdef _OPENMP
    // Set maximal number of threads.
    int oldMaxThreads = omp_get_max_threads();
    omp_set_num_threads(options.threadCount == 0 ? 1 : options.threadCount);
#endif  // #ifdef _OPENMP

	//////////////////////////////////////////////////////////////////////////////
	// open left reads file

    bool success;
    if (!isEqual(readFileNames[0], "-"))
        success = open(options.readFile, toCString(readFileNames[0]));
    else
        success = open(options.readFile, std::cin);

    if (!success)
        return RAZERS_READS_FAILED;

    //////////////////////////////////////////////////////////////////////////////
    // get read length
    int readLength = estimateReadLength(options.readFile);
    if (readLength == RAZERS_READS_FAILED)
    {
        cerr << "Failed to open reads file " << readFileNames[0] << endl;
        cerr << "Exiting ..." << endl;
        return RAZERS_READS_FAILED;
    }
    if (readLength == 0)
    {
        cerr << "Failed to read the first read sequence." << endl;
        cerr << "Exiting ..." << endl;
        return RAZERS_READS_FAILED;
    }

    if (options.trimLength > readLength)
        options.trimLength = readLength;

#ifndef NO_PARAM_CHOOSER
    if (options.threshold != 0 && !(isSet(argParser, "shape") || isSet(argParser, "threshold")))
    {
        if (options.lowMemory)
            pm_options.maxWeight = 13;
        pm_options.verbose = (options._debugLevel >= 1);
        pm_options.optionErrorRate = options.errorRate;
        if (options.gapMode == RAZERS_UNGAPPED)
        {
            pm_options.optionProbINSERT = (ParamChooserOptions::TFloat)0.0;
            pm_options.optionProbDELETE = (ParamChooserOptions::TFloat)0.0;
        }
        else
        {
            pm_options.optionProbINSERT = (ParamChooserOptions::TFloat)0.01;    //this number is basically without meaning, any value > 0 will lead to
            pm_options.optionProbDELETE = (ParamChooserOptions::TFloat)0.01;    //edit distance parameter choosing
        }

        if (options.trimLength > 0)
            readLength = options.trimLength;
        if (readLength > 0)
        {
/*			if(options.maqMapping && readLength != options.artSeedLength)
                pm_options.totalN = options.artSeedLength;
            else*/
            pm_options.totalN = readLength;
            if (options._debugLevel >= 1)
                cerr << "___PARAMETER_CHOOSING__" << endl;
            if (!chooseParams(options, pm_options))
            {
                if (pm_options.verbose)
                    cerr << "Couldn't find preprocessed parameter files. Please configure manually (options --shape and --threshold)." << endl;
                if (options._debugLevel >= 1)
                    cerr << "Using default configurations (shape = " << options.shape << " and k-mer lemma)." << endl;
            }
            if (options._debugLevel >= 1)
                cerr << endl;
        }
        else
        {
            cerr << "Failed to load reads" << endl;
            cerr << "Exiting ..." << endl;
            return RAZERS_READS_FAILED;
        }
    }
#endif

    if (argc > 1)
        options.commandLine = argv[1];
    for (int i = 2; i < argc; ++i)
    {
        appendValue(options.commandLine, ' ');
        append(options.commandLine, argv[i]);
    }

    int result = mapReads(genomeFileNames, readFileNames, options);
    if (result != 0)
        cerr << "Exiting ..." << endl;

#ifdef _OPENMP
    // Restoring number of threads for side-effect freeness.
    omp_set_num_threads(oldMaxThreads);
#endif  // #ifdef _OPENMP

#ifdef RAZERS_PROFILE
    dumpTimeline("razers.profile.txt", true);
#endif  // #ifndef RAZERS_PROFILE

    return result;
}
