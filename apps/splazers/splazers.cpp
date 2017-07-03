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

#define SEQAN_PROFILE					// enable time measuring
//#define SEQAN_DEBUG_SWIFT				// test SWIFT correctness and print bucket parameters
//#define RAZERS_DEBUG					// print verification regions
#define RAZERS_PRUNE_QGRAM_INDEX		// ignore highly abundant q-grams
#define RAZERS_CONCATREADS				// use <ConcatDirect> StringSet to store reads!!!
//#define RAZERS_MEMOPT					// optimize memory consumption
#define RAZERS_MASK_READS				// remove matches with max-hits optimal hits on-the-fly
//#define NO_PARAM_CHOOSER				// disable loss-rate parameter choosing
//#define RAZERS_PARALLEL				// parallelize using Intel's Threading Building Blocks
#define RAZERS_MATEPAIRS				// enable paired-end matching
//#define RAZERS_DIRECT_MAQ_MAPPING
//#define SEQAN_USE_SSE2_WORDS			// use SSE2 128-bit integers for MyersBitVector
//#define RAZERS_OPENADDRESSING
#define RAZERS_SPLICED
//#define TRY_SCORES

#ifdef RAZERS_SPLICED
  #define RAZERS_MATEPAIRS				
#endif

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

#ifdef RAZERS_PARALLEL
#include "razers_parallel.h"
#endif

#ifdef RAZERS_MATEPAIRS
#include "razers_matepairs.h"
#endif

#include "razers_spliced.h" 



#include <iostream>
#include <sstream>

using namespace std;
using namespace seqan;


template <typename TOptions>
int getGenomeFileNameList(CharString filename, StringSet<CharString> & genomeFileNames, TOptions const & options)
{
	ifstream file;
	file.open(toCString(filename),ios_base::in | ios_base::binary);
	if(!file.is_open())
		return RAZERS_GENOME_FAILED;

    DirectionIterator<std::fstream, Input>::Type reader(file);
    clear(genomeFileNames);
    // if file does not start with a fasta header --> list of multiple reference genome files
	if (!atEnd(reader) && value(reader) != '>' && value(reader) != '@')
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
	else		//if file starts with a fasta header --> regular one-genome-file input
		appendValue(genomeFileNames, filename, Exact());
	file.close();
	return 0;
}

//////////////////////////////////////////////////////////////////////////////
// Main read mapper function
template <typename TSpec>
int mapReads(
	StringSet<CharString> & genomeFileNames,
	StringSet<CharString> & readFileNames,	// NULL terminated list of one/two read files (single/mate-pairs)
	CharString & errorPrbFileName,
	RazerSOptions<TSpec> &options)
{
	TReadSet				readSet;
	TReadRegions			readRegions;
	StringSet<CharString>	genomeNames;	// genome names, taken from the Fasta file
	StringSet<CharString>	readNames;		// read names, taken from the Fasta file
	TMatches				matches;		// resulting forward/reverse matches
	String<String<unsigned short> > 	stats;		// needed for mapping quality calculation 

	// dump configuration in verbose mode
	if (options._debugLevel >= 1) 
	{
		CharString bitmap;
		Shape<Dna, GenericShape> shape;
		stringToShape(shape, options.shape);
		shapeToString(bitmap, shape);
		CharString bitmapR;
		Shape<Dna, GenericShape> shapeR;
		if(options.minMatchLen > 0)
		{
			stringToShape(shapeR, options.shapeR);
			shapeToString(bitmapR, shapeR);
		}		

		cerr << "___SETTINGS____________" << endl;
		cerr << "Genome file:                     \t" << genomeFileNames[0] << endl;
		if (length(readFileNames)<2)
			cerr << "Read file:                       \t" << std::flush << readFileNames[0] << endl;
		else
		{
			cerr << "Read files:                      \t" << std::flush << readFileNames[0] << endl;
			for (unsigned i = 1; i < length(readFileNames); ++i)
				cerr << "                                 \t" << readFileNames[i] << endl;
		}
		cerr << "Compute forward matches:         \t";
		if (options.forward)	cerr << "YES" << endl;
		else				cerr << "NO" << endl;
		cerr << "Compute reverse matches:         \t";
		if (options.reverse)		cerr << "YES" << endl;
		else				cerr << "NO" << endl;
		cerr << "Error rate:                      \t" << options.errorRate << endl;
		cerr << "Minimal threshold:               \t" << options.threshold << endl;
		cerr << "Shape:                           \t" << bitmap << endl;
		if(options.minMatchLen > 0)
		{
			cerr << "Suffix Minimal threshold:        \t" << options.thresholdR << endl;
			cerr << "Suffix Shape:                    \t" << bitmapR << endl;
		}
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
		if (!loadReads(readSet, readNames, toCString(readFileNames[0]), toCString(readFileNames[1]), options)) {
		//if (!loadReads(readSet, readQualities, readNames, readFileNames[0], readFileNames[1], options)) {
			cerr << "Failed to load reads" << endl;
			return RAZERS_READS_FAILED;
		}
	}
	else
#endif
	{
		if (options.anchored) {
			if (!loadReadsSam(readSet, readNames, readRegions, toCString(readFileNames[0]), options)) {
			//if (!loadReads(readSet, readQualities, readNames, readFileNames[0], readFileNames[1], options)) {
				cerr << "Failed to load reads" << endl;
				return RAZERS_READS_FAILED;
			}
		}
		else if (!loadReads(readSet, readNames, toCString(readFileNames[0]), options)) {
		//if (!loadReads(readSet, readQualities, readNames, readFileNames[0], readFileNames[1], options)) {
			cerr << "Failed to load reads" << endl;
			return RAZERS_READS_FAILED;
		}
	} 
	if(empty(readSet))
	{
		cerr << "Failed to load reads. File empty? " << endl;
		cerr << "Note that anchored split read mapping requires a SAM file + option -an needs to switched on.\nUnanchored split read mapping requires a Fasta/Fastq read file. " << endl;
		return RAZERS_READS_FAILED;
	}
    if(empty(readSet[0]))
    {
        // is this just the first read is a low quality read?
        unsigned i = 1;
        while(i < length(readSet))
        { 
            if(!empty(readSet[i]))
                break;
            ++i;
        }
        if(i == length(readSet)) // all read sequences are empty
        {
            cerr << "Failed to load reads. File empty? All reads low quality? " << endl;
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
	int error = mapReads(matches, genomeFileNames, genomeNames, gnoToFileMap, readSet, readRegions, stats, options);
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

	// anchored reads were all mapped onto forward strand, undo reverse complementing
	if(!empty(readRegions) && options.anchored && options.outputFormat != 4) //SAM output
	{
		for(unsigned i = 0; i < length(readSet); ++i)
			if(readRegions[i].i2.i1 > 1) // should have been > 0 before...
				reverseComplement(readSet[i]);
	}

#ifdef RAZERS_PARALLEL
    delete[] options.patternMutex;
#endif

	//////////////////////////////////////////////////////////////////////////////
	// Step 3: Remove duplicates and output matches
	if (!options.spec.DONT_DUMP_RESULTS)
		dumpMatches(matches, genomeNames, genomeFileNames, gnoToFileMap, readSet, readRegions, stats, readNames, readFileNames[0], errorPrbFileName, options);

	return 0;
}	

//////////////////////////////////////////////////////////////////////////////
// Command line parsing and parameter choosing
int main(int argc, const char *argv[]) 
{
	RazerSOptions<>			options;
	ParamChooserOptions		pm_options;

	StringSet<CharString>	genomeFileNames;
	StringSet<CharString>	readFileNames;
	CharString				errorPrbFileName;

	// Change defaults
	options.forward = false;
	options.reverse = false;
	
    ArgumentParser parser("splazers");
    addUsageLine(parser, "[\\fIOPTIONS\\fP] <\\fIGENOME FILE\\fP> <\\fIREADS FILE\\fP>");
#ifdef RAZERS_MATEPAIRS
	addUsageLine(parser, "[\\fIOPTIONS\\fP] <\\fIGENOME FILE\\fP> <\\fIREADS FILE 1\\fP> <\\fIREADS FILE 2\\fP>");
#endif
    setShortDescription(parser, "Split-map read sequences");
    addDescription(parser,
            "SplazerS uses a prefix-suffix mapping strategy to split-map read sequences."
            "If a SAM file of mapped reads is given as input, all unmapped but anchored"
            "reads are split-mapped onto anchoring target regions (specify option -an),"
            "if a Fasta/q file of reads is given, reads are split-mapped onto the whole"
            "reference sequence.");
            
    addDescription(parser, "(c) Copyright 2010 by Anne-Katrin Emde.");
    setCategory(parser, "Read Mapping");
    setVersion(parser, SEQAN_APP_VERSION " [" SEQAN_REVISION "]");
    setDate(parser, SEQAN_DATE);

    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE));
    setValidValues(parser, 0, seqan::SeqFileIn::getFileExtensions());
    setHelpText(parser, 0, "A reference genome file.");
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "READS", true));
    std::vector<std::string> exts = seqan::SeqFileIn::getFileExtensions();
    exts.push_back(".sam");
    setValidValues(parser, 1, exts);
    setHelpText(parser, 1, "Either one (single-end) or two (paired-end) read files.");

	addSection(parser, "Main Options:");
    addOption(parser, ArgParseOption("o", "output", "Change output filename. Default: <\\fIREADS FILE\\fP>.result.", ArgParseOption::OUTPUT_FILE));
	addOption(parser, ArgParseOption("f",  "forward",           "only compute forward matches"));
	addOption(parser, ArgParseOption("r",  "reverse",           "only compute reverse complement matches"));
    addOption(parser, ArgParseOption("i", "percent-identity", "Percent identity threshold.", ArgParseOption::DOUBLE));
    setMinValue(parser, "percent-identity", "50");
    setMaxValue(parser, "percent-identity", "100");
    setDefaultValue(parser, "percent-identity", 100 - (100.0 * options.errorRate));
#ifndef NO_PARAM_CHOOSER
    addOption(parser, ArgParseOption("rr", "recognition-rate",  "set the percent recognition rate", ArgParseOption::DOUBLE));
    setMinValue(parser, "recognition-rate", "80");
    setMaxValue(parser, "recognition-rate", "100");
    setDefaultValue(parser, "recognition-rate", 100 - (100.0 * pm_options.optionLossRate));
    addOption(parser, ArgParseOption("pd", "param-dir", "Read user-computed parameter files in the directory <\\fIDIR\\fP>.", ArgParseOption::STRING, "DIR"));
#endif
    addOption(parser, ArgParseOption("id", "indels", "Allow indels. Default: mismatches only."));

    addOption(parser, ArgParseOption("ll", "library-length", "Paired-end library length.", ArgParseOption::INTEGER));
    setMinValue(parser, "library-length", "1");
    setDefaultValue(parser, "library-length", options.libraryLength);
    addOption(parser, ArgParseOption("le", "library-error", "Paired-end library length tolerance.", ArgParseOption::INTEGER));
    setMinValue(parser, "library-error", "0");
    setDefaultValue(parser, "library-error", options.libraryError);

    addOption(parser, ArgParseOption("m", "max-hits", "Output only <\\fINUM\\fP> of the best hits.", ArgParseOption::INTEGER));
    setMinValue(parser, "max-hits", "1");
    setDefaultValue(parser, "max-hits", options.maxHits);

    addOption(parser, ArgParseOption("", "unique", "Output only unique best matches (-m 1 -dr 0 -pa)."));
    addOption(parser, ArgParseOption("tr", "trim-reads", "Trim reads to given length. Default: off.", ArgParseOption::INTEGER));
    setMinValue(parser, "trim-reads", "14");
	addOption(parser, ArgParseOption("mcl", "min-clipped-len",  "min. read length for read clipping", ArgParseOption::INTEGER));
    setDefaultValue(parser, "min-clipped-len", options.minClippedLen);
    setMinValue(parser, "min-clipped-len", "1");
	addOption(parser, ArgParseOption("qih", "quality-in-header","quality string in fasta header"));

	addOption(parser, ArgParseOption("ou", "outputUnmapped",    "output filename for unmapped reads", ArgParseOption::OUTPUT_FILE));
    
	addOption(parser, ArgParseOption("v",  "verbose",           "verbose mode"));
	addOption(parser, ArgParseOption("vv", "vverbose",          "very verbose mode"));
	addSection(parser, "Output Format Options:");
	addOption(parser, ArgParseOption("a",  "alignment",         "dump the alignment for each match"));
    addOption(parser, ArgParseOption("pa", "purge-ambiguous",   "purge reads with more than max-hits best matches"));
	addOption(parser, ArgParseOption("dr", "distance-range",    "only consider matches with at most NUM more errors compared to the best (default output all)", ArgParseOption::INTEGER));
	addOption(parser, ArgParseOption("of", "output-format", "Set output format. 0 = RazerS, 1 = Enhanced Fasta, 2 = Eland, 3 = GFF, 4 = SAM.", ArgParseOption::INTEGER));
    setMinValue(parser, "output-format", "0");
    setMaxValue(parser, "output-format", "4");
    addOption(parser, ArgParseOption("gn", "genome-naming", "Select how genomes are named. 0 = use Fasta id, 1 = enumerate beginning with 1.", ArgParseOption::INTEGER));
    setMinValue(parser, "genome-naming", "0");
    setMaxValue(parser, "genome-naming", "1");
    setDefaultValue(parser, "genome-naming", options.genomeNaming);
  
    addOption(parser, ArgParseOption("rn", "read-naming", "Select how reads are named. 0 = use Fasta id, 1 = enumerate beginning with 1.", ArgParseOption::INTEGER));
    setDefaultValue(parser, "read-naming", options.readNaming);
    setMinValue(parser, "read-naming", "0");
    setMaxValue(parser, "read-naming", "1");

    addOption(parser, ArgParseOption("so", "sort-order", "Select how matches are sorted. 0 = read number, 1 = genome position.", ArgParseOption::INTEGER));
    setDefaultValue(parser, "sort-order", options.sortOrder);
    setMinValue(parser, "sort-order", "0");
    setMaxValue(parser, "sort-order", "1");
    
    addOption(parser, ArgParseOption("pf", "position-format", "Select begin/end position numbering (see Coordinate section below). 0 = gap space, 1 = position space.", ArgParseOption::INTEGER));
    setMinValue(parser, "position-format", "0");
    setMaxValue(parser, "position-format", "1");
    setDefaultValue(parser, "position-format", options.positionFormat);


	addSection(parser, "Split Mapping Options:");
	addOption(parser, ArgParseOption("sm", "split-mapping",   "min. match length for prefix/suffix mapping (to disable split mapping, set to 0)", ArgParseOption::INTEGER));
    setDefaultValue(parser,"split-mapping",options.minMatchLen);
    
	addOption(parser, ArgParseOption("maxG", "max-gap",    "max. length of middle gap", ArgParseOption::INTEGER));
    setDefaultValue(parser, "max-gap", options.maxGap);
	addOption(parser, ArgParseOption("minG", "min-gap",    "min. length of middle gap (for edit distance mapping about 10% of read length is recommended)", ArgParseOption::INTEGER));
    setDefaultValue(parser, "min-gap", options.minGap);
    addOption(parser, ArgParseOption("ep", "errors-prefix",    "max. number of errors in prefix match", ArgParseOption::INTEGER));
    setDefaultValue(parser, "errors-prefix", options.maxPrefixErrors);
    
	addOption(parser, ArgParseOption("es", "errors-suffix",    "max. number of errors in suffix match", ArgParseOption::INTEGER));
    setDefaultValue(parser, "errors-suffix", options.maxSuffixErrors);
    
	addOption(parser, ArgParseOption("gl", "genome-len",    "genome length in Mb, for computation of expected number of random matches", ArgParseOption::INTEGER));
    setDefaultValue(parser, "genome-len", options.specifiedGenomeLen);
    setMaxValue(parser, "genome-len", "10000");
    
	addOption(parser, ArgParseOption("an", "anchored",           "anchored split mapping, only unmapped reads with mapped mates will be considered, requires the reads to be given in SAM format"));
	addOption(parser, ArgParseOption("pc",  "penalty-c",    "percent of read length, used as penalty for split-gap", ArgParseOption::INTEGER));
    setDefaultValue(parser, "penalty-c", options.penaltyC);

	addSection(parser, "Filtration Options:");
    addOption(parser, ArgParseOption("s", "shape", "Manually set k-mer shape.", ArgParseOption::STRING, "BITSTRING"));
    setDefaultValue(parser, "shape", options.shape);
    hideOption(parser,"shape");
    addOption(parser, ArgParseOption("t", "threshold", "Manually set minimum k-mer count threshold.", ArgParseOption::INTEGER));
    setMinValue(parser, "threshold", "1");    
    hideOption(parser,"threshold");

    addOption(parser, ArgParseOption("oc", "overabundance-cut", "Set k-mer overabundance cut ratio.", ArgParseOption::INTEGER));
    setMinValue(parser, "overabundance-cut", "0");
    setMaxValue(parser, "overabundance-cut", "1");
    addOption(parser, ArgParseOption("rl", "repeat-length", "Skip simple-repeats of length <\\fINUM\\fP>.", ArgParseOption::INTEGER));
    setMinValue(parser, "repeat-length", "1");
    setDefaultValue(parser, "repeat-length", options.repeatLength);
    addOption(parser, ArgParseOption("tl", "taboo-length", "Set taboo length.", ArgParseOption::INTEGER));
    setMinValue(parser, "taboo-length", "1");
    setDefaultValue(parser, "taboo-length", options.tabooLength);
    addOption(parser, ArgParseOption("lm", "low-memory",        "decrease memory usage at the expense of runtime"));

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


    // Parse command line.
    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
    {
        if (res == ArgumentParser::PARSE_ERROR)
            cerr << "Exiting ..." << endl;
        return (res == ArgumentParser::PARSE_ERROR) ? RAZERS_INVALID_OPTIONS : 0;
    }

    bool stop = false;

	//////////////////////////////////////////////////////////////////////////////
	// Extract options
	getOptionValue(options.forward, parser, "forward");
	getOptionValue(options.reverse, parser, "reverse");
	getOptionValue(options.errorRate, parser, "percent-identity");
#ifndef NO_PARAM_CHOOSER
	getOptionValue(pm_options.optionLossRate, parser, "recognition-rate");
	getOptionValue(pm_options.paramFolder, parser, "param-dir");
#endif
	getOptionValue(options.hammingOnly, parser, "indels");
	options.hammingOnly = !options.hammingOnly;
	getOptionValue(options.libraryLength, parser, "library-length");
	getOptionValue(options.libraryError, parser, "library-error");
	getOptionValue(options.maxHits, parser, "max-hits");
	getOptionValue(options.purgeAmbiguous, parser, "purge-ambiguous");
	getOptionValue(options.distanceRange, parser, "distance-range");
	if (isSet(parser, "distance-range")) options.distanceRange++;
	getOptionValue(options.dumpAlignment, parser, "alignment");
	getOptionValue(options.output, parser, "output");
	getOptionValue(options.outputFormat, parser, "output-format");
    getOptionValue(options.sortOrder, parser, "sort-order");
    getOptionValue(options.genomeNaming, parser, "genome-naming");
    getOptionValue(options.readNaming, parser, "read-naming");
    getOptionValue(options.positionFormat, parser, "position-format");
    getOptionValue(options.outputUnmapped, parser, "outputUnmapped");
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
	
 	getOptionValue(options.minMatchLen, parser, "split-mapping");
 	getOptionValue(options.maxGap, parser, "max-gap");
 	getOptionValue(options.minGap, parser, "min-gap");
 	getOptionValue(options.maxPrefixErrors, parser, "errors-prefix");
 	getOptionValue(options.maxSuffixErrors, parser, "errors-suffix");
 	getOptionValue(options.specifiedGenomeLen, parser, "genome-len");
	getOptionValue(options.anchored, parser, "anchored");
	getOptionValue(options.penaltyC, parser, "penalty-c");
	
	getOptionValue(options.minClippedLen, parser, "min-clipped-len");
	getOptionValue(options.fastaIdQual, parser, "quality-in-header");
	
    if (isSet(parser, "help") || isSet(parser, "version")) return 0;	// print help or version and exit
	if (isSet(parser, "verbose")) options._debugLevel = max(options._debugLevel, 1);
	if (isSet(parser, "vverbose")) options._debugLevel = max(options._debugLevel, 3);
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

    if(length(readFileNames) == 2)
        options.minMatchLen = 0;  //switch off split mapping if mate pairs are given

    //////////////////////////////////////////////////////////////////////////////
    // Check options
    if ((options.outputFormat == 4 && options.minMatchLen == 0) && (stop = true))
        cerr << "Invalid output format option. Note that SAM output is only available for split mapping." << endl;
    if ((options.outputFormat > 4 && options.outputFormat != 33 ) && (stop = true))
        cerr << "Invalid output format option." << endl;
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
		if ((ones < 7 || ones > maxOnes) && !options.anchored && !stop)
			cerr << "Warning: Shape should contain at least 7 and at most " << maxOnes << " '1's" << endl;
		else if ((ones < 5 || ones > maxOnes) && options.anchored && !stop)
			cerr << "Warning: Shape should contain at least 5 and at most " << maxOnes << " '1's" << endl;
		options.shapeR = options.shape;
		options.thresholdR = options.threshold;

	}
    if (length(readFileNames) == 1 && !options.anchored)
        options.libraryLength = -1;		// only 1 readset -> disable mate-pair mapping
    if ((getArgumentValueCount(parser, 1) > maxReadFiles) && (stop = true))
        cerr << "More than " << maxReadFiles << " read files specified." << endl;
    if ((getArgumentValueCount(parser, 1) == 0) && (stop = true))
        cerr << "No read files specified." << endl;
    if ((options.minClippedLen < 0) && (stop = true))
        cerr << "Min. clipped read length must be a value greater 0" << endl;

    options.errorRate = (100.0 - options.errorRate) / 100.0;
    pm_options.optionLossRate = (100.0 - pm_options.optionLossRate) / 100.0;
	if (stop)
	{
		cerr << "Exiting ..." << endl;
		return RAZERS_INVALID_OPTIONS;
	}

	if (options.anchored) 
		options.reverse = false; // we know which strand to look on, because the mate is anchored and orientation is known!
	if ((options.minMatchLen != 0 && options.minMatchLen < 6) && (stop = true))
		cerr << "Minimal match length in spliced mapping must be a value greater than 5" << endl;
	if ((options.maxGap <= 0) && (stop = true))
		cerr << "Max. gap length must be a value greater 0" << endl;
	if ((options.minGap < 0) && (stop = true))
		cerr << "Min. gap length must be a value greater 0" << endl;
	if(options.maxPrefixErrors == -1)
		options.maxPrefixErrors = (int)(options.minMatchLen * options.errorRate);
//	if(options.maxSuffixErrors == -1)
//		options.maxSuffixErrors = (int)(options.minMatchLen * options.errorRate);

	//////////////////////////////////////////////////////////////////////////////
	// get read length
	int readLength = options.minMatchLen;
	if(options.minMatchLen == 0)
	{
		readLength = estimateReadLength(toCString(readFileNames[0]));
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
	}
			
		
#ifndef NO_PARAM_CHOOSER
	if (!(isSet(parser, "shape") || isSet(parser, "threshold")))
	{
		if (options.lowMemory) pm_options.maxWeight = 13;
		pm_options.verbose = (options._debugLevel >= 1);
		pm_options.optionErrorRate = options.errorRate;
		if (options.hammingOnly)
		{
			pm_options.optionProbINSERT = (ParamChooserOptions::TFloat)0.0;
			pm_options.optionProbDELETE = (ParamChooserOptions::TFloat)0.0;
		}
		else
		{
			pm_options.optionProbINSERT = (ParamChooserOptions::TFloat)0.01;	//this number is basically without meaning, any value > 0 will lead to
			pm_options.optionProbDELETE = (ParamChooserOptions::TFloat)0.01;	//edit distance parameter choosing
		}

		if (empty(pm_options.paramFolder)) 
		{
			string razersFolder = argv[0];
			size_t lastPos = razersFolder.find_last_of('/') + 1;
			if (lastPos == razersFolder.npos + 1) lastPos = razersFolder.find_last_of('\\') + 1;
			if (lastPos == razersFolder.npos + 1) lastPos = 0;
			razersFolder.erase(lastPos); 
			pm_options.paramFolderPath = razersFolder;
		}
		if (options.trimLength > 0) readLength = options.trimLength;
		if (readLength > 0)
		{
/*			if(options.maqMapping && readLength != options.artSeedLength)
				pm_options.totalN = options.artSeedLength;
			else*/
				pm_options.totalN = readLength;
			if(options.minMatchLen>0)
 			{
				pm_options.totalN = options.minMatchLen;
				pm_options.optionErrorRate = (double)options.maxSuffixErrors/options.minMatchLen + 0.001;
				if (!chooseParams(options,pm_options))
				{
					if (pm_options.verbose) 
						cerr << "Couldn't find preprocessed parameter files. Please configure manually (options --shape and --threshold)." << endl;
					cerr << "Using default configurations (shape = " << options.shape << " and q-gram lemma)." << endl;
				}
				options.shapeR = options.shape;
				options.thresholdR = options.threshold;
				if (options._debugLevel >= 1) cerr << endl;
				pm_options.totalN = options.minMatchLen;
				pm_options.optionErrorRate = (double)options.maxPrefixErrors/options.minMatchLen + 0.001;
			}
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
