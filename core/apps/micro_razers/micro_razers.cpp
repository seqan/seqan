 /*==========================================================================
             MicroRazerS - Rapid Alignment of Small RNA Reads
                   http://www.seqan.de/projects/microRazers.html

 ============================================================================
  Copyright (C) 2009 by Anne-Katrin Emde

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

#define SEQAN_PROFILE			// enable time measuring
#define RAZERS_PRUNE_QGRAM_INDEX
#define RAZERS_CONCATREADS		// use <ConcatDirect> StringSet to store reads
#define RAZERS_MEMOPT			// optimize memory consumption
//#define RAZERS_MASK_READS		// remove matches with max-hits optimal hits on-the-fly
#define RAZERS_MICRO_RNA
#define RAZERS_EXTENDED_MATCH

#include "seqan/platform.h"
#ifdef PLATFORM_WINDOWS
    #define SEQAN_DEFAULT_TMPDIR "C:\\TEMP\\"
#else
    #define SEQAN_DEFAULT_TMPDIR "./"
#endif

#include "../razers/razers.h"
#include "outputFormat.h"
#include "../razers/paramChooser.h"
#include <seqan/arg_parse.h>



#include <iostream>
#include <sstream>

using namespace std;
using namespace seqan;


template<typename TOptions>
int getGenomeFileNameList(CharString const & filename, StringSet<CharString> & genomeFileNames, TOptions const & options)
{
    std::ifstream file;
    file.open(toCString(filename), std::ios_base::in | std::ios_base::binary);
    if(!file.is_open())
        return RAZERS_GENOME_FAILED;

    seqan::RecordReader<std::ifstream, seqan::SinglePass<> > reader(file);

    CharString nameStr;
    if (value(reader) != '>' && value(reader) != '@')
    {
        // If file does not start with a fasta header --> list of multiple reference genome files.
        if(options._debugLevel >=1)
            std::cout << std::endl << "Reading multiple genome files:" << std::endl;
        /*      //locations of genome files are relative to list file's location
         ::std::string tempGenomeFile(filename);
         size_t lastPos = tempGenomeFile.find_last_of('/') + 1;
         if (lastPos == tempGenomeFile.npos) lastPos = tempGenomeFile.find_last_of('\\') + 1;
         if (lastPos == tempGenomeFile.npos) lastPos = 0;
         ::std::string filePrefix = tempGenomeFile.substr(0,lastPos);*/
        unsigned i = 0;
        for (; !atEnd(reader); ++i)
        {
            clear(nameStr);
            int res = skipWhitespaces(reader);
            if (res == EOF_BEFORE_SUCCESS)
                break;  // Done, no more file name.
            if (res != 0)
                return res;  // Error reading.
            res = readGraphs(nameStr, reader);
            if (res != 0 && res != EOF_BEFORE_SUCCESS)
                return res;
            appendValue(genomeFileNames,nameStr,Generous());
            if(options._debugLevel >= 2)
                std::cout << "Genome file #" << (i + 1) << ": " << genomeFileNames[length(genomeFileNames) - 1] << std::endl;
        }
        if(options._debugLevel >=1)
            std::cout << i << " genome files total." << std::endl;
    }
    else
    {
        // If file starts with a fasta header --> regular one-genome-file input.
        appendValue(genomeFileNames,filename,Generous());
    }

    return 0;
}


//////////////////////////////////////////////////////////////////////////////
// Main read mapper function
template <typename TSpec>
int mapMicroRNAReads(
    StringSet<CharString> & genomeFileNames,
    StringSet<CharString> & readFileNames,	// NULL terminated list of one/two read files (single/mate-pairs)
    CharString              errorPrbFileName,
    RazerSOptions<TSpec> &options)
{
    MultiFasta			genomeSet;
    TReadSet			readSet;
    StringSet<CharString>		genomeNames;		// genome names, taken from the Fasta file
    StringSet<CharString>		readNames;		// read names, taken from the Fasta file
    TMatches			matches;		// resulting forward/reverse matches
    String<String<unsigned short> > stats;



    // dump configuration in verbose mode
    if (options._debugLevel >= 1)
    {
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
        ::std::cerr << "Compute forward matches:         \t";
        if (options.forward)	::std::cerr << "YES" << ::std::endl;
        else				::std::cerr << "NO" << ::std::endl;
        ::std::cerr << "Compute reverse matches:         \t";
        if (options.reverse)		::std::cerr << "YES" << ::std::endl;
        else				::std::cerr << "NO" << ::std::endl;
        ::std::cerr << "Seed length:                     \t" << options.rnaSeedLength << ::std::endl;
        if(options.exactSeed)
            ::std::cerr << "Error in seed:                   \tNO" << ::std::endl;
        else
            ::std::cerr << "Error in seed:                   \tYES" << ::std::endl;
        ::std::cerr << "Shape:                           \t" << options.shape << ::std::endl;
        ::std::cerr << "Minimal threshold:               \t" << options.threshold << ::std::endl;
        ::std::cerr << ::std::endl;
    }

    // circumvent numerical obstacles
    options.errorRate += 0.0000001;



    //////////////////////////////////////////////////////////////////////////////
    // Step 1: Load fasta files and determine genome file type
    SEQAN_PROTIMESTART(load_time);

    if (!loadReads(readSet, readNames, toCString(readFileNames[0]), options)) {
    //if (!loadReads(readSet, readQualities, readNames, readFileNames[0], readFileNames[1], options)) {
        ::std::cerr << "Failed to load reads" << ::std::endl;
        return RAZERS_READS_FAILED;
    }

    if (options._debugLevel >= 1) ::std::cerr << lengthSum(readSet) << " bps of " << length(readSet) << " reads loaded." << ::std::endl;
    options.timeLoadFiles = SEQAN_PROTIMEDIFF(load_time);

    if (options._debugLevel >= 1)
        ::std::cerr << "Loading reads took               \t" << options.timeLoadFiles << " seconds" << ::std::endl;

    StringSet<CharString> genomeFileNameList;
    int result = getGenomeFileNameList(genomeFileNames[0], genomeFileNameList, options);
    if (result == RAZERS_GENOME_FAILED)
    {
        ::std::cerr << "Failed to open genome file " << genomeFileNames[0] << ::std::endl;
        return result;
    }


    //////////////////////////////////////////////////////////////////////////////
    // Step 2: Find matches using SWIFT
    ::std::map<unsigned,::std::pair< ::std::string,unsigned> > gnoToFileMap; //map to retrieve genome filename and sequence number within that file
    int error = mapReads(matches, genomeFileNameList, genomeNames, gnoToFileMap, readSet, stats, options);
    if (error != 0)
    {
        switch (error)
        {
            case RAZERS_GENOME_FAILED:
                ::std::cerr << "Failed to load genomes" << ::std::endl;
                break;

            case RAZERS_INVALID_SHAPE:
                ::std::cerr <<	"Invalid Shape" << endl << endl;
                break;
        }
        return error;
    }

    //////////////////////////////////////////////////////////////////////////////
    // Step 3: Remove duplicates and output matches
    if (!options.spec.DONT_DUMP_RESULTS)
        dumpMatches(matches, genomeNames, genomeFileNameList, gnoToFileMap, readSet, stats, readNames, readFileNames[0], errorPrbFileName, options);


    return 0;
}



int main(int argc, const char *argv[])
{
    //////////////////////////////////////////////////////////////////////////////
    // Parse command line

    RazerSOptions<>		options;
    ParamChooserOptions	pm_options;

    bool				paramChoosing = true; //false;
    bool				paramChoosingError = false;
    StringSet<CharString>	genomeFileNames;
    StringSet<CharString>	readFileNames;
    CharString          errorPrbFileName;

    options.forward = false;
    options.reverse = false;
    options.hammingOnly = true;
    options.microRNA = true;
    pm_options.optionLossRate = 0.0;
    bool seedError = false;
    
    
    ArgumentParser parser("micro_razers");
    addUsageLine(parser, "[\\fIOPTIONS\\fP] <\\fIGENOME FILE\\fP> <\\fIREADS FILE\\fP>");
    addDescription(parser,
            "MicroRazerS uses a prefix-based mapping strategy to map "
            "small RNA reads possibly containing 3' adapter sequence.");
        
    addDescription(parser, "(c) Copyright 2009 by Anne-Katrin Emde.");
    setCategory(parser, "Read Mapping");
    setVersion(parser, "1.0.1");
    setDate(parser, "Jul 2009" );

    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE));
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE, "READS", true));

    //////////////////////////////////////////////////////////////////////////////
    // Define options

    addSection(parser, "Main Options:");
    addOption(parser, ArgParseOption("o", "output", "Change output filename. Default: <\\fIREADS FILE\\fP>.result.", ArgParseOption::OUTPUTFILE));
#ifndef NO_PARAM_CHOOSER
    addOption(parser, ArgParseOption("rr", "recognition-rate",  "set the percent recognition rate", ArgParseOption::DOUBLE));
    setMinValue(parser, "recognition-rate", "80");
    setMaxValue(parser, "recognition-rate", "100");
    setDefaultValue(parser, "recognition-rate", 100 - (100.0 * pm_options.optionLossRate));
#endif
    addOption(parser, ArgParseOption("sL", "seed-length",       "seed length", ArgParseOption::INTEGER));
    setMinValue(parser, "seed-length", "10");
    setDefaultValue(parser, "seed-length", options.rnaSeedLength);
    addOption(parser, ArgParseOption("sE", "seed-error", "allow for one error in the seed"));

    addOption(parser, ArgParseOption("f", "forward", "map reads only to forward strands."));
    addOption(parser, ArgParseOption("r", "reverse", "map reads only to reverse strands."));

    addOption(parser, ArgParseOption("mN", "match-N",           "\'N\' matches with all other characters"));
    addOption(parser, ArgParseOption("m",  "max-hits",          "output only NUM of the best hits", ArgParseOption::INTEGER));
    setMinValue(parser, "max-hits", "1");
    setDefaultValue(parser, "max-hits", options.maxHits);
    addOption(parser, ArgParseOption("s", "shape", "Manually set k-mer shape.", ArgParseOption::STRING, "BITSTRING"));
    setDefaultValue(parser, "shape", options.shape);
    hideOption(parser,"shape");
    addOption(parser, ArgParseOption("t", "threshold", "Manually set minimum k-mer count threshold.", ArgParseOption::INTEGER));
    setMinValue(parser, "threshold", "1");    
    hideOption(parser,"threshold");
    addOption(parser, ArgParseOption("pa", "purge-ambiguous",   "purge reads with more than max-hits best matches"));
    addOption(parser, ArgParseOption("lm", "low-memory",        "decrease memory usage at the expense of runtime"));
    addOption(parser, ArgParseOption("v",  "verbose",           "verbose mode"));
    addOption(parser, ArgParseOption("vv", "vverbose",          "very verbose mode"));

    addSection(parser, "Output Format Options:");
    addOption(parser, ArgParseOption("of", "output-format", "Set output format. 0 = MicroRazerS format, 1 = SAM.", ArgParseOption::INTEGER));
    setMinValue(parser, "output-format", "0");
    setMaxValue(parser, "output-format", "1");

    addOption(parser, ArgParseOption("a", "alignment",          "dump the alignment for each match"));
    
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
    getOptionValue(seedError, parser, "seed-error");
    getOptionValue(options.rnaSeedLength, parser, "seed-length");
    getOptionValue(options.forward, parser, "forward");
    getOptionValue(options.reverse, parser, "reverse");
    getOptionValue(options.maxHits, parser, "max-hits");
    getOptionValue(options.purgeAmbiguous, parser, "purge-ambiguous");
    getOptionValue(options.output, parser, "output");
    getOptionValue(options.outputFormat, parser, "output-format");   // add help line
    getOptionValue(options.dumpAlignment, parser, "alignment");
    getOptionValue(options.sortOrder, parser, "sort-order");
    getOptionValue(options.genomeNaming, parser, "genome-naming");
    getOptionValue(options.readNaming, parser, "read-naming");
    getOptionValue(options.positionFormat, parser, "position-format");

    getOptionValue(options.shape, parser, "shape");				// undocumented for own use
    getOptionValue(options.threshold, parser, "threshold");		// undocumented
    getOptionValue(options.lowMemory, parser, "low-memory");
    getOptionValue(options.matchN, parser, "match-N");

    if (isSet(parser, "help") || isSet(parser, "version")) return 0;	// print help or version and exit
    if (isSet(parser, "verbose")) options._debugLevel = max(options._debugLevel, 1);
    if (isSet(parser, "vverbose")) options._debugLevel = max(options._debugLevel, 3);
    if (!options.forward && !options.reverse)  // enable both per default
    {
        options.forward = true;
        options.reverse = true;
    }
//    appendValue(genomeFileNames, getArgumentValue(parser, 0));
//    appendValue(readFileNames, getArgumentValue(parser, 1), Generous());
   
    resize(genomeFileNames, length(genomeFileNames) + 1);
    getArgumentValue(back(genomeFileNames), parser, 0);
    resize(readFileNames, 1, Exact());
    getArgumentValue(readFileNames[0], parser, 1, 0);
    
    //////////////////////////////////////////////////////////////////////////////
    // Check options
    if (seedError) options.exactSeed = false;
    if ((pm_options.optionLossRate < 0.0 || pm_options.optionLossRate > 0.2) && (stop = true))
        cerr << "Recognition rate must be a value between 80 and 100" << endl;
    if ((options.maxHits < 1) && (stop = true))
        cerr << "Maximum hits threshold must be greater than 0" << endl;
    if ((options.outputFormat > 1) && (stop = true))
        cerr << "Invalid output format option." << endl;
    if ((options.sortOrder > 1) && (stop = true))
        cerr << "Invalid sort order options." << endl;
    if ((options.genomeNaming > 1) && (stop = true))
        cerr << "Invalid genome naming options." << endl;
    if ((options.readNaming > 2) && (stop = true))
        cerr << "Invalid read naming options." << endl;
    if ((options.positionFormat > 1) && (stop = true))
        cerr << "Invalid position format options." << endl;
    if ((options.threshold < 1) && (stop = true))
        cerr << "Threshold must be a value greater than 0" << endl;
    if ((options.rnaSeedLength < 10) && (stop = true))
        cerr << "Minimum seed length is 10" << endl << endl;
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
        if ((ones < 6 || ones > maxOnes) && !stop)
            cerr << "Warning: Shape should contain at least 6 and at most " << maxOnes << " '1's" << endl;
    }
    
    if ((getArgumentValueCount(parser, 1) > 1) && (stop = true))
        cerr << "More than one read file specified." << endl;
    if ((getArgumentValueCount(parser, 1) == 0) && (stop = true))
        cerr << "No read files specified." << endl;

    // get read length
    int readLength = estimateReadLength(toCString(readFileNames[0]));
    if (readLength == RAZERS_READS_FAILED)
    {
        ::std::cerr << "Failed to open reads file " << readFileNames[0] << ::std::endl;
        return 0;
    }
    if (readLength == 0) {
        ::std::cerr << "Failed to read the first read sequence.";
        return 0;
    }

    if(options.outputFormat == 1) options.outputFormat = 4; // of 4 is sam in all RazerS tools
    options.distanceRange = 1; // 1 actually means 0... 0 would disable distanceRange
    if(options.exactSeed) options.errorRate = 0.0001;
    else options.errorRate = (double)1.2/(options.rnaSeedLength);

    if (options.trimLength > readLength)
        options.trimLength = readLength;

    if (paramChoosing) // && !options.exactSeed)
    {
        if(options.lowMemory) pm_options.maxWeight = 13;
        pm_options.optionErrorRate = options.errorRate;
        pm_options.verbose = (options._debugLevel >= 1);
        if (options.hammingOnly)   // at the moment we do hamming only
        {
            pm_options.optionProbINSERT = 0.0;
            pm_options.optionProbDELETE = 0.0;
        }
        else
        {
            pm_options.optionProbINSERT = 0.01;	//this number is basically without meaning, any value > 0 will lead to
            pm_options.optionProbDELETE = 0.01;	//edit distance parameter choosing
        }
        string paramFolder = argv[0];
        size_t lastPos = paramFolder.find_last_of('/') + 1;
        if (lastPos == paramFolder.npos + 1) lastPos = paramFolder.find_last_of('\\') + 1;
        if (lastPos == paramFolder.npos + 1) lastPos = 0;
        paramFolder.erase(lastPos);
        pm_options.paramFolderPath = paramFolder;
        pm_options.fprefix[0] = "uniform";
        pm_options.prefixCount = true;
        pm_options.totalN = options.rnaSeedLength;
        if (readLength > 0)
        {
            if (options._debugLevel >= 1)
                cerr << "___PARAMETER_CHOOSING__" << endl;
            if (!chooseParams(options,pm_options))
            {
                paramChoosingError = true;
                if (pm_options.verbose)
                    cerr << "Couldn't find preprocessed parameter files." << endl;
                if (options._debugLevel >= 1)
                    cerr << "Possibly using suboptimal configurations." << endl;
            }
            cerr << endl;
        } else
        {
            cerr << "Failed to load reads" << endl;
            return 0;
        }
    }

    if(paramChoosingError)
    {
        if(options.exactSeed)
        {
            //shapes are designed such that span = seedLength-1 --> t = 2 for 0 errors
            if(options.rnaSeedLength==10)
                options.shape = "1111111111";
            if(options.rnaSeedLength==11)
                options.shape = "11111111111";
            if(options.rnaSeedLength==12)
                options.shape = "111111111111";
            if(options.rnaSeedLength==13)
                options.shape = "1111111111111";
            if(options.rnaSeedLength==14)
                options.shape = "11111111111111";
            if(options.rnaSeedLength==15)
                options.shape = "11111111111111";
            if(options.rnaSeedLength==16)
                options.shape = "101111111111111";
            if(options.rnaSeedLength==17)
                options.shape = "1011111111110111";
            if(options.rnaSeedLength==18)
                options.shape = "10111111011101111";
            if(options.rnaSeedLength==19)
                options.shape = "101111110111011101";
            if(options.rnaSeedLength==20)
                options.shape = "1011111010111011101";
            if(options.rnaSeedLength>20)
                options.shape = "10111110101110111001";
            options.threshold = 2;
            if(options.rnaSeedLength< 15) options.threshold = 1;
        }
        else
        {
            //best shape found by bestShape.R with t=1
            if(options.rnaSeedLength==10)
                options.shape = "11011011";
            if(options.rnaSeedLength==11)
                options.shape = "110110101";
            if(options.rnaSeedLength==12)
                options.shape = "1101101101";
            if(options.rnaSeedLength==13)
                options.shape = "11011011011";
            if(options.rnaSeedLength==14)
                options.shape = "11101110111";
            if(options.rnaSeedLength==15)
                options.shape = "1101101101101"; //fuer n=15: 13  1  2  4  5  7  8 10 11
            if(options.rnaSeedLength==16)
                options.shape = "11011011011011"; //14  1  2  4  5  7  8 10 11 13
            if(options.rnaSeedLength==17)
                options.shape = "11101110111011"; //14  1  2  3  5  6  7  9 10 11 13
            if(options.rnaSeedLength==18)
                options.shape = "111011101110111";//15  1  2  3  5  6  7  9 10 11 13 14
            if(options.rnaSeedLength==19)
                options.shape = "11011011011011011";//17  1  2  4  5  7  8 10 11 13 14 16
            if(options.rnaSeedLength==20)
                options.shape = "11101110111011101";//17  1  2  3  5  6  7  9 10 11 13 14 15
            if(options.rnaSeedLength>20)
                options.shape = "111011101110111011";//18  1  2  3  5  6  7  9 10 11 13 14 15 17
            options.threshold = 1;

        }
    }


    int result = mapMicroRNAReads(genomeFileNames, readFileNames, errorPrbFileName, options);

    if (result != 0)
        cerr << "Exiting ..." << endl;
    return result;
}
