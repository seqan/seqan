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

#define SEQAN_PROFILE                    // enable time measuring
//#define SEQAN_DEBUG_SWIFT                // test SWIFT correctness and print bucket parameters
//#define RAZERS_DEBUG                    // print verification regions
#define RAZERS_PRUNE_QGRAM_INDEX        // ignore highly abundant q-grams
#define RAZERS_CONCATREADS                // use <ConcatDirect> StringSet to store reads
#define RAZERS_MEMOPT                    // optimize memory consumption
#define RAZERS_MASK_READS                // remove matches with max-hits optimal hits on-the-fly
//#define NO_PARAM_CHOOSER                // disable loss-rate parameter choosing
#define RAZERS_MATEPAIRS                // enable paired-end matching
//#define RAZERS_DIRECT_MAQ_MAPPING
//#define SEQAN_USE_SSE2_WORDS            // use SSE2 128-bit integers for MyersBitVector
#define RAZERS_OPENADDRESSING
//#define RAZERS_SPLICED

#ifdef RAZERS_SPLICED
  #define RAZERS_MATEPAIRS
#endif

#include <cctype>

#include <seqan/platform.h>
#ifdef STDLIB_VS
    #define SEQAN_DEFAULT_TMPDIR "C:\\TEMP\\"
#else
    #define SEQAN_DEFAULT_TMPDIR "./"
#endif

#include <seqan/arg_parse.h>
#include "razers.h"
#include "outputFormat.h"
#include "paramChooser.h"

#ifdef RAZERS_MATEPAIRS
#include "razers_matepairs.h"
#endif

#ifdef RAZERS_SPLICED
#include "razers_spliced.h"
#endif



#include <iostream>
#include <sstream>

using namespace std;
using namespace seqan;


//////////////////////////////////////////////////////////////////////////////
// Read a list of genome file names
template<typename TSpec>
int getGenomeFileNameList(CharString filename, StringSet<CharString> & genomeFileNames, RazerSOptions<TSpec> &options)
{
    ifstream file;
    file.open(toCString(filename),ios_base::in | ios_base::binary);
    if(!file.is_open())
        return RAZERS_GENOME_FAILED;

    DirectionIterator<std::fstream, Input>::Type reader(file);
    if (!atEnd(reader))
        return 0;

    clear(genomeFileNames);
    if (*reader == '>' && *reader != '@')    //if file does not start with a fasta header --> list of multiple reference genome files
    {
        if(options._debugLevel >=1)
            cout << endl << "Reading multiple genome files:" <<endl;

        unsigned i = 1;
        CharString line;
        while(!atEnd(reader))
        {
            readLine(line, reader);
            cropOuter(line, IsWhitespace());
            appendValue(genomeFileNames, line);
            if(options._debugLevel >=2)
                cout <<"Genome file #"<< i <<": " << back(genomeFileNames) << endl;
            ++i;
        }
        if(options._debugLevel >=1)
            cout << i-1 << " genome files total." <<endl;
    }
    else        //if file starts with a fasta header --> regular one-genome-file input
        appendValue(genomeFileNames, filename, Exact());
    file.close();
    return 0;
}

//////////////////////////////////////////////////////////////////////////////
// Main read mapper function
template <typename TSpec>
int mapReads(
    StringSet<CharString> & genomeFileNames,
    StringSet<CharString> & readFileNames,    // NULL terminated list of one/two read files (single/mate-pairs)
    CharString & errorPrbFileName,
    RazerSOptions<TSpec> &options)
{
    TReadSet                readSet;
    StringSet<CharString>    genomeNames;    // genome names, taken from the Fasta file
    StringSet<CharString>    readNames;        // read names, taken from the Fasta file
    TMatches                matches;        // resulting forward/reverse matches
    String<String<unsigned short> >     stats;        // needed for mapping quality calculation

    // dump configuration in verbose mode
    if (options._debugLevel >= 1)
    {
        CharString bitmap;
        Shape<Dna, GenericShape> shape;
        stringToShape(shape, options.shape);
        shapeToString(bitmap, shape);

        cerr << "___SETTINGS____________" << endl;
        cerr << "Genome file:                     \t" << genomeFileNames[0] << endl;
        if (length(readFileNames) > 1u && empty(readFileNames[1]))
            cerr << "Read file:                       \t" << readFileNames[0] << endl;
        else
        {
            cerr << "Read files:                      \t" << readFileNames[0] << endl;
            for (unsigned i = 1; i < length(readFileNames); ++i)
                cerr << "                                 \t" << readFileNames[i] << endl;
        }
        cerr << "Compute forward matches:         \t";
        if (options.forward)    cerr << "YES" << endl;
        else                cerr << "NO" << endl;
        cerr << "Compute reverse matches:         \t";
        if (options.reverse)        cerr << "YES" << endl;
        else                cerr << "NO" << endl;
        cerr << "Error rate:                      \t" << options.errorRate << endl;
        cerr << "Minimal threshold:               \t" << options.threshold << endl;
        cerr << "Shape:                           \t" << bitmap << endl;
        cerr << "Repeat threshold:                \t" << options.repeatLength << endl;
        cerr << "Overabundance threshold:         \t" << options.abundanceCut << endl;
        cerr << "Taboo length:                    \t" << options.tabooLength << endl;
        cerr << endl;
    }

    // circumvent numerical obstacles
    options.errorRate += 0.0000001;

    //////////////////////////////////////////////////////////////////////////////
    // Step 1: Load fasta files and determine genome file type
    SEQAN_PROTIMESTART(load_time);

#ifdef RAZERS_MATEPAIRS
    if (length(readFileNames) == 2)
    {
        if (!loadReads(readSet, readNames, options.readFile, toCString(readFileNames[1]), options)) {
        //if (!loadReads(readSet, readQualities, readNames, readFileNames[0], readFileNames[1], options)) {
            cerr << "Failed to load reads" << endl;
            return RAZERS_READS_FAILED;
        }
    }
    else
#endif
    {
        if (!loadReads(readSet, readNames, options.readFile, options)) {
        //if (!loadReads(readSet, readQualities, readNames, readFileNames[0], readFileNames[1], options)) {
            cerr << "Failed to load reads" << endl;
            return RAZERS_READS_FAILED;
        }
    }

    if (options._debugLevel >= 1) cerr << lengthSum(readSet) << " bps of " << length(readSet) << " reads loaded." << endl;
    options.timeLoadFiles = SEQAN_PROTIMEDIFF(load_time);

    if (options._debugLevel >= 1)
        cerr << "Loading reads took               \t" << options.timeLoadFiles << " seconds" << endl;

    if (length(genomeFileNames) == 1)
    {
        int result = getGenomeFileNameList(genomeFileNames[0], genomeFileNames, options);
        if (result == RAZERS_GENOME_FAILED)
        {
            cerr << "Failed to open genome file " << genomeFileNames[0] << endl;
            return result;
        }
    }

    //////////////////////////////////////////////////////////////////////////////
    // Step 2: Find matches using SWIFT
#ifdef RAZERS_PARALLEL
    typedef typename RazerSOptions<TSpec>::TMutex TMutex;
    options.patternMutex = new TMutex[length(readSet)];
#endif

    map<unsigned,pair< string,unsigned> > gnoToFileMap; //map to retrieve genome filename and sequence number within that file
    int error = mapReads(matches, genomeFileNames, genomeNames, gnoToFileMap, readSet, stats, options);
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

#ifdef RAZERS_PARALLEL
    delete[] options.patternMutex;
#endif

    //////////////////////////////////////////////////////////////////////////////
    // Step 3: Remove duplicates and output matches
    if (!options.spec.DONT_DUMP_RESULTS)
        dumpMatches(matches, genomeNames, genomeFileNames, gnoToFileMap, readSet, stats, readNames, readFileNames[0], errorPrbFileName, options);

    return 0;
}

void setUpArgumentParser(ArgumentParser & parser, RazerSOptions<> const & options, ParamChooserOptions const & pm_options)
{
    setAppName(parser, "razers");
    setShortDescription(parser, "Fast Read Mapping with Sensitivity Control");
    setCategory(parser, "Read Mapping");
    setVersion(parser, SEQAN_APP_VERSION " [" SEQAN_REVISION "]");
    setDate(parser, SEQAN_DATE);

    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE));
    setValidValues(parser, 0, seqan::SeqFileIn::getFileExtensions());
    setHelpText(parser, 0, "A reference genome file.");
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "READS", true));
    setValidValues(parser, 1, seqan::SeqFileIn::getFileExtensions());
    setHelpText(parser, 1, "Either one (single-end) or two (paired-end) read files.");

    addUsageLine(parser, "[\\fIOPTIONS\\fP] <\\fIGENOME FILE\\fP> <\\fIREADS FILE\\fP>");
#ifdef RAZERS_MATEPAIRS
    addUsageLine(parser, "[\\fIOPTIONS\\fP] <\\fIGENOME FILE\\fP> <\\fIMP-READS FILE1\\fP> <\\fIMP-READS FILE2\\fP>");
#endif

    addDescription(parser, "RazerS is a versatile full-sensitive read mapper based on a k-mer counting filter. "
                           "It supports single and paired-end mapping, and optimally parametrizes "
                           "the filter based on a user-defined minimal sensitivity. "
                           "See \\fIhttp://www.seqan.de/projects/razers\\fP for more information.");

    addDescription(parser, "Input to RazerS is a reference genome file and either one file with single-end reads "
                           "or two files containing left or right mates of paired-end reads. Use - to read single-end "
                           "reads from stdin.");

    addDescription(parser, "(c) Copyright 2009 by David Weese.");

    addSection(parser, "Main Options");
    addOption(parser, ArgParseOption("f", "forward", "Map reads only to forward strands."));
    addOption(parser, ArgParseOption("r", "reverse", "Map reads only to reverse strands."));
    addOption(parser, ArgParseOption("i", "percent-identity", "Percent identity threshold.", ArgParseOption::DOUBLE));
    setMinValue(parser, "percent-identity", "50");
    setMaxValue(parser, "percent-identity", "100");
    setDefaultValue(parser, "percent-identity", 100 - (100.0 * options.errorRate));
#ifndef NO_PARAM_CHOOSER
    addOption(parser, ArgParseOption("rr", "recognition-rate", "Percent recognition rate.", ArgParseOption::DOUBLE));
    setMinValue(parser, "recognition-rate", "80");
    setMaxValue(parser, "recognition-rate", "100");
    setDefaultValue(parser, "recognition-rate", 100 - (100.0 * pm_options.optionLossRate));
    addOption(parser, ArgParseOption("pd", "param-dir", "Read user-computed parameter files in the directory <\\fIDIR\\fP>.", ArgParseOption::STRING, "DIR"));
#endif
    addOption(parser, ArgParseOption("id", "indels", "Allow indels. Default: mismatches only."));
#ifdef RAZERS_MATEPAIRS
    addOption(parser, ArgParseOption("ll", "library-length", "Paired-end library length.", ArgParseOption::INTEGER));
    setMinValue(parser, "library-length", "1");
    setDefaultValue(parser, "library-length", options.libraryLength);
    addOption(parser, ArgParseOption("le", "library-error", "Paired-end library length tolerance.", ArgParseOption::INTEGER));
    setMinValue(parser, "library-error", "0");
    setDefaultValue(parser, "library-error", options.libraryError);
#endif
    addOption(parser, ArgParseOption("m", "max-hits", "Output only <\\fINUM\\fP> of the best hits.", ArgParseOption::INTEGER));
    setMinValue(parser, "max-hits", "1");
    setDefaultValue(parser, "max-hits", options.maxHits);
    addOption(parser, ArgParseOption("", "unique", "Output only unique best matches (-m 1 -dr 0 -pa)."));
    addOption(parser, ArgParseOption("tr", "trim-reads", "Trim reads to given length. Default: off.", ArgParseOption::INTEGER));
    setMinValue(parser, "trim-reads", "14");
    addOption(parser, ArgParseOption("o", "output", "Change output filename (use - to dump to stdout in razers format). Default: <\\fIREADS FILE\\fP>.razers.", ArgParseOption::OUTPUT_FILE));
    setValidValues(parser, "output", ".razers .eland .fa .fasta .gff");
    addOption(parser, ArgParseOption("v", "verbose", "Verbose mode."));
    addOption(parser, ArgParseOption("vv", "vverbose", "Very verbose mode."));

    addSection(parser, "Output Format Options");
    addOption(parser, ArgParseOption("a", "alignment", "Dump the alignment for each match (only \\fIrazer\\fP or \\fIfasta\\fP format)."));
    addOption(parser, ArgParseOption("pa", "purge-ambiguous", "Purge reads with more than <\\fImax-hits\\fP> best matches."));
    addOption(parser, ArgParseOption("dr", "distance-range", "Only consider matches with at most NUM more errors compared to the best. Default: output all.", ArgParseOption::INTEGER));
    addOption(parser, ArgParseOption("gn", "genome-naming", "Select how genomes are named (see Naming section below).", ArgParseOption::INTEGER));
    setMinValue(parser, "genome-naming", "0");
    setMaxValue(parser, "genome-naming", "1");
    setDefaultValue(parser, "genome-naming", options.genomeNaming);
    addOption(parser, ArgParseOption("rn", "read-naming", "Select how reads are named (see Naming section below).", ArgParseOption::INTEGER));
    setDefaultValue(parser, "read-naming", options.readNaming);
    setMinValue(parser, "read-naming", "0");
    setMaxValue(parser, "read-naming", "2");
    addOption(parser, ArgParseOption("so", "sort-order", "Select how matches are sorted (see Sorting section below).", ArgParseOption::INTEGER));
    setDefaultValue(parser, "sort-order", options.sortOrder);
    setMinValue(parser, "sort-order", "0");
    setMaxValue(parser, "sort-order", "1");
    addOption(parser, ArgParseOption("pf", "position-format", "Select begin/end position numbering (see Coordinate section below).", ArgParseOption::INTEGER));
    setMinValue(parser, "position-format", "0");
    setMaxValue(parser, "position-format", "1");
    setDefaultValue(parser, "position-format", options.sortOrder);

    addSection(parser, "Filtration Options");
    addOption(parser, ArgParseOption("s", "shape", "Manually set k-mer shape.", ArgParseOption::STRING, "BITSTRING"));
    setDefaultValue(parser, "shape", options.shape);
    addOption(parser, ArgParseOption("t", "threshold", "Manually set minimum k-mer count threshold.", ArgParseOption::INTEGER));
    setMinValue(parser, "threshold", "1");
    addOption(parser, ArgParseOption("oc", "overabundance-cut", "Set k-mer overabundance cut ratio.", ArgParseOption::INTEGER));
    setMinValue(parser, "overabundance-cut", "0");
    setMaxValue(parser, "overabundance-cut", "1");
    addOption(parser, ArgParseOption("rl", "repeat-length", "Skip simple-repeats of length <\\fINUM\\fP>.", ArgParseOption::INTEGER));
    setMinValue(parser, "repeat-length", "1");
    setDefaultValue(parser, "repeat-length", options.repeatLength);
    addOption(parser, ArgParseOption("tl", "taboo-length", "Set taboo length.", ArgParseOption::INTEGER));
    setMinValue(parser, "taboo-length", "1");
    setDefaultValue(parser, "taboo-length", options.tabooLength);

    addOption(parser, ArgParseOption("lm", "low-memory", "Decrease memory usage at the expense of runtime."));
#ifdef RAZERS_DIRECT_MAQ_MAPPING
    addSection(parser, "Mapping Quality Options:");
    addOption(parser, ArgParseOption("mq", "mapping-quality", "Switch on mapping quality mode."));
    addOption(parser, ArgParseOption("nbi", "no-below-id", "Do not report matches with seed identity < percent id."));
    addOption(parser, ArgParseOption("qsl", "mq-seed-length", "Set MAQ seed length." , ArgParseOption::INTEGER));
    setMinValue(parser, "mq-seed-length", "24");
    setDefaultValue(parser, "mq-seed-length", options.artSeedLength));
    addOption(parser, ArgParseOption("smq", "seed-mism-quality", "Set maximal mismatch-quality sum in the seed.", ArgParseOption::INTEGER));
    setMinValue(parser, "seed-mism-quality", "0");
    setDefaultValue(parser, "seed-mism-quality", options.maxMismatchQualSum));
    addOption(parser, ArgParseOption("tmq", "total-mism-quality", "Set total maximal mismatch-quality.", ArgParseOption::INTEGER));
    setMinValue(parser, "total-mism-quality", "0");
    setDefaultValue(parser, "total-mism-quality", options.absMaxQualSumErrors));
#endif

    addSection(parser, "Verification Options");
    addOption(parser, ArgParseOption("mN", "match-N", "N matches all other characters. Default: N matches nothing."));
    addOption(parser, ArgParseOption("ed", "error-distr", "Write error distribution to \\fIFILE\\fP.", ArgParseOption::STRING, "FILE"));

#ifdef RAZERS_SPLICED
    addOption(parser, ArgParseOption("sm", "spliced-mapping", "Set minimal match length for spliced mapping.", ArgParseOption::INTEGER));
    setMinValue(parser, "spliced-mapping", "6"));
    setDefaultValue(parser, "spliced-mapping", options.minMatchLen));
    addOption(parser, ArgParseOption("maxD", "max-distance", "Set maximal distance of prefix-suffix match.", ArgParseOption::INTEGER));
    setMinValue(parser, "max-distance", "0"));
    setDefaultValue(parser, "max-distance", options.maxDistance));
#endif

    addOption(parser, ArgParseOption("mcl", "min-clipped-len", "Set minimal read length for read clipping.", ArgParseOption::INTEGER));
    setMinValue(parser, "min-clipped-len", "0");
    setDefaultValue(parser, "min-clipped-len", options.minClippedLen);
    addOption(parser, ArgParseOption("qih", "quality-in-header", "Quality string in fasta header."));

//#ifdef RAZERS_SPLICED
//    addOption(parser, CommandLineOption("sm", "spliced-mapping",   "min. match length for prefix/suffix alignment strategy", OptionType::Int | OptionType::Label, options.minMatchLen));
//    addOption(parser, CommandLineOption("maxD", "max-distance",    "max distance of pref/suff match", OptionType::Int | OptionType::Label, options.maxDistance));
//#endif
//
//    addOption(parser, CommandLineOption("mcl", "min-clipped-len",  "min. read length for read clipping", OptionType::Int | OptionType::Label, options.minClippedLen));
//    addOption(parser, CommandLineOption("qih", "quality-in-header","quality string in fasta header", OptionType::Boolean));
//

    addTextSection(parser, "Formats, Naming, Sorting, and Coordinate Schemes");

    addText(parser, "RazerS supports various output formats. The output format is detected automatically from the file name suffix.");
    addListItem(parser, ".razers", "Razer format");
    addListItem(parser, ".fa, .fasta", "Enhanced Fasta format");
    addListItem(parser, ".eland", "Eland format");
    addListItem(parser, ".gff", "GFF format");

    addText(parser, "");
    addText(parser, "By default, reads and contigs are referred by their Fasta ids given in the input files. "
                    "With the \\fB-gn\\fP and \\fB-rn\\fP options this behaviour can be changed:");
    addListItem(parser, "0", "Use Fasta id.");
    addListItem(parser, "1", "Enumerate beginning with 1.");
    addListItem(parser, "2", "Use the read sequence (only for short reads!).");

    addText(parser, "");
    addText(parser, "The way matches are sorted in the output file can be changed with the \\fB-so\\fP option for the following formats: "
                    "\\fBrazer\\fP, \\fBfasta\\fP, \\fBsam\\fP, and \\fBamos\\fP. Primary and secondary sort keys are:");
    addListItem(parser, "0", "1. read number, 2. genome position");
    addListItem(parser, "1", "1. genome position, 2. read number");

    addText(parser, "");
    addText(parser, "The coordinate space used for begin and end positions can be changed with the "
                    "\\fB-pf\\fP option for the \\fBrazer\\fP and \\fBfasta\\fP formats:");
    addListItem(parser, "0", "Gap space. Gaps between characters are counted from 0.");
    addListItem(parser, "1", "Position space. Characters are counted from 1.");

    addTextSection(parser, "Examples");
    addListItem(parser,
                "\\fBrazers\\fP \\fBexample/genome.fa\\fP \\fBexample/reads.fa\\fP \\fB-id\\fP \\fB-a\\fP \\fB-mN\\fP \\fB-v\\fP",
                "Map single-end reads with 4% error rate, indels, and output the alignments. Ns are considered to match everything.");
    addListItem(parser,
                "\\fBrazers\\fP \\fBexample/genome.fa\\fP \\fBexample/reads.fa\\fP \\fBexample/reads2.fa\\fP \\fB-id\\fP \\fB-mN\\fP",
                "Map paired-end reads with up to 4% errors, indels, and output concordantly mapped pairs within default library size. Ns are considered to match everything.");
}

ArgumentParser::ParseResult
extractOptions(
    StringSet<CharString> & genomeFileNames, StringSet<CharString> & readFileNames,
    RazerSOptions<> & options, ParamChooserOptions & pm_options, CharString &errorPrbFileName,
    ArgumentParser const & parser)
{
    //////////////////////////////////////////////////////////////////////////////
    // Extract options

    bool stop = false;
    getOptionValue(options.forward, parser, "forward");
    getOptionValue(options.reverse, parser, "reverse");
    getOptionValue(options.errorRate, parser, "percent-identity");
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
    getOptionValue(options.hammingOnly, parser, "indels");
    options.hammingOnly = !options.hammingOnly;

#ifdef RAZERS_MATEPAIRS
    getOptionValue(options.libraryLength, parser, "library-length");
    getOptionValue(options.libraryError, parser, "library-error");
#endif
    getOptionValue(options.maxHits, parser, "max-hits");
    getOptionValue(options.purgeAmbiguous, parser, "purge-ambiguous");
    getOptionValue(options.distanceRange, parser, "distance-range");
    if (isSet(parser, "distance-range"))
        options.distanceRange++;
    getOptionValue(options.dumpAlignment, parser, "alignment");
    getOptionValue(options.sortOrder, parser, "sort-order");
    getOptionValue(options.genomeNaming, parser, "genome-naming");
    getOptionValue(options.readNaming, parser, "read-naming");
    getOptionValue(options.positionFormat, parser, "position-format");
    getOptionValue(options.shape, parser, "shape");
    getOptionValue(options.threshold, parser, "threshold");
    getOptionValue(options.abundanceCut, parser, "overabundance-cut");
    getOptionValue(options.repeatLength, parser, "repeat-length");

#ifdef RAZERS_DIRECT_MAQ_MAPPING
    getOptionValue(options.fastaIdQual, parser, "quality-in-header");
    getOptionValue(options.maqMapping, parser, "mapping-quality");
    getOptionValue(options.noBelowIdentity, parser, "no-below-id");
    getOptionValue(options.artSeedLength, parser, "mq-seed-length");
    getOptionValue(options.maxMismatchQualSum, parser, "seed-mism-quality");
    getOptionValue(options.absMaxQualSumErrors, parser, "total-mism-quality");
#endif
    getOptionValue(options.lowMemory, parser, "low-memory");
    getOptionValue(options.trimLength, parser, "trim-reads");
    getOptionValue(options.tabooLength, parser, "taboo-length");
    getOptionValue(options.matchN, parser, "match-N");
    getOptionValue(errorPrbFileName, parser, "error-distr");
#ifdef RAZERS_SPLICED
    getOptionValue(options.minMatchLen, parser, "spliced-mapping");
    getOptionValue(options.maxDistance, parser, "max-distance");
#endif
    getOptionValue(options.minClippedLen, parser, "min-clipped-len");
    getOptionValue(options.fastaIdQual, parser, "quality-in-header");
    if (isSet(parser, "verbose"))
        options._debugLevel = max(options._debugLevel, 1);
    if (isSet(parser, "vverbose"))
        options._debugLevel = max(options._debugLevel, 3);
    if (isSet(parser, "unique"))
    {
        options.maxHits = 1;
        options.distanceRange = 1;
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
    else
        options.outputFormat = 0;   // default is ".razers"

    if (isSet(parser, "shape"))
    {
        unsigned ones = 0;
        unsigned zeros = 0;
        for(unsigned i = 0; i < length(options.shape); ++i)
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
    }
    if (getArgumentValueCount(parser, 1) == 1)
        options.libraryLength = -1;        // only 1 readset -> disable mate-pair mapping
    if ((getArgumentValueCount(parser, 1) > maxReadFiles) && (stop = true))
        cerr << "More than " << maxReadFiles << " read files specified." << endl;
    if ((getArgumentValueCount(parser, 1) == 0) && (stop = true))
        cerr << "No read files specified." << endl;
    if ((options.minClippedLen < 0) && (stop = true))
        cerr << "Min. clipped read length must be a value greater 0" << endl;

    options.errorRate = (100.0 - options.errorRate) / 100.0;
    pm_options.optionLossRate = (ParamChooserOptions::TFloat)(100.0 - pm_options.optionLossRate) / 100.0;

    return (stop) ? ArgumentParser::PARSE_ERROR : ArgumentParser::PARSE_OK;
}

//////////////////////////////////////////////////////////////////////////////
// Command line parsing and parameter choosing
int main(int argc, const char *argv[])
{
    RazerSOptions<>            options;
    ParamChooserOptions        pm_options;

    StringSet<CharString>    genomeFileNames;
    StringSet<CharString>    readFileNames;
    CharString                errorPrbFileName;

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
    res = extractOptions(genomeFileNames, readFileNames, options, pm_options, errorPrbFileName, argParser);
    if (res != ArgumentParser::PARSE_OK)
    {
        cerr << "Exiting ..." << endl;
        return RAZERS_INVALID_OPTIONS;
    }

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
    if (readLength == 0) {
        cerr << "Failed to read the first read sequence." << endl;
        cerr << "Exiting ..." << endl;
        return RAZERS_READS_FAILED;
    }

    if (options.trimLength > readLength)
        options.trimLength = readLength;

#ifndef NO_PARAM_CHOOSER
    if (!(isSet(argParser, "shape") || isSet(argParser, "threshold")))
    {
        if (options.lowMemory) pm_options.maxWeight = 13;
        pm_options.verbose = (options._debugLevel >= 1);
        pm_options.optionErrorRate = (ParamChooserOptions::TFloat)options.errorRate;
        if (options.hammingOnly)
        {
            pm_options.optionProbINSERT = (ParamChooserOptions::TFloat)0.0;
            pm_options.optionProbDELETE = (ParamChooserOptions::TFloat)0.0;
        }
        else
        {
            pm_options.optionProbINSERT = (ParamChooserOptions::TFloat)0.01;    //this number is basically without meaning, any value > 0 will lead to
            pm_options.optionProbDELETE = (ParamChooserOptions::TFloat)0.01;    //edit distance parameter choosing
        }

        if (options.trimLength > 0) readLength = options.trimLength;
        if (readLength > 0)
        {
/*            if(options.maqMapping && readLength != options.artSeedLength)
                pm_options.totalN = options.artSeedLength;
            else*/
                pm_options.totalN = readLength;
#ifdef RAZERS_SPLICED
            if(options.minMatchLen>0)
                pm_options.totalN = options.minMatchLen;
#endif
            if (options._debugLevel >= 1)
                cerr << "___PARAMETER_CHOOSING__" << endl;
            if (!chooseParams(options,pm_options))
            {
                if (pm_options.verbose)
                    cerr << "Couldn't find preprocessed parameter files. Please configure manually (options --shape and --threshold)." << endl;
                cerr << "Using default configurations (shape = " << options.shape << " and q-gram lemma)." << endl;
            }
            if (options._debugLevel >= 1) cerr << endl;
        } else
        {
            cerr << "Failed to load reads" << endl;
            cerr << "Exiting ..." << endl;
            return RAZERS_READS_FAILED;
        }
    }
#endif

#ifdef RAZERS_PARALLEL
    tbb::task_scheduler_init scheduler;
#endif

    int result = mapReads(genomeFileNames, readFileNames, errorPrbFileName, options);
    if (result != 0)
        cerr << "Exiting ..." << endl;
    return result;
}
