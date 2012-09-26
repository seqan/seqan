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



#include <iostream>
#include <sstream>

using namespace std;
using namespace seqan;


template<typename TSpec>
int getGenomeFileNameList(char const * filename, StringSet<CharString> & genomeFileNames, RazerSOptions<TSpec> &options)
{
//IOREV could possibly be moved to fasta module
	::std::ifstream file;
	file.open(filename,::std::ios_base::in | ::std::ios_base::binary);
	if(!file.is_open())
		return RAZERS_GENOME_FAILED;

	char c = _streamGet(file);
	if (c != '>' && c != '@')	//if file does not start with a fasta header --> list of multiple reference genome files
	{
		if(options._debugLevel >=1)
			::std::cout << ::std::endl << "Reading multiple genome files:" <<::std::endl;
		
		unsigned i = 1;
		while(!_streamEOF(file))
		{ 
			_parseSkipWhitespace(file, c);
			appendValue(genomeFileNames,_parseReadFilepath(file,c));
			if(options._debugLevel >=2)
				::std::cout <<"Genome file #"<< i <<": " << genomeFileNames[length(genomeFileNames)-1] << ::std::endl;
			++i;
			_parseSkipWhitespace(file, c);
		}
		if(options._debugLevel >=1)
			::std::cout << i-1 << " genome files total." <<::std::endl;
	}
	else		//if file starts with a fasta header --> regular one-genome-file input
		appendValue(genomeFileNames,filename);
	file.close();
	return 0;

}


//////////////////////////////////////////////////////////////////////////////
// Main read mapper function
template <typename TSpec>
int mapMicroRNAReads(
	const char *genomeFileName,
	const char *readFileNames[],	// NULL terminated list of one/two read files (single/mate-pairs)
	const char *errorPrbFileName,
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
		::std::cerr << "___SETTINGS____________" << ::std::endl;
		::std::cerr << "Genome file:                     \t" << genomeFileName << ::std::endl;
		if (empty(readFileNames[1]))
			::std::cerr << "Read file:                       \t" << readFileNames[0] << ::std::endl;
		else
		{
			::std::cerr << "Read files:                      \t" << readFileNames[0] << ::std::endl;
			for (const char **i = readFileNames + 1; !empty(*i); ++i)
				::std::cerr << "                                 \t" << *i << ::std::endl;
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

	if (!loadReads(readSet, readNames, readFileNames[0], options)) {
	//if (!loadReads(readSet, readQualities, readNames, readFileNames[0], readFileNames[1], options)) {
		::std::cerr << "Failed to load reads" << ::std::endl;
		return RAZERS_READS_FAILED;
	}
 
	if (options._debugLevel >= 1) ::std::cerr << lengthSum(readSet) << " bps of " << length(readSet) << " reads loaded." << ::std::endl;
	options.timeLoadFiles = SEQAN_PROTIMEDIFF(load_time);

	if (options._debugLevel >= 1)
		::std::cerr << "Loading reads took               \t" << options.timeLoadFiles << " seconds" << ::std::endl;

	StringSet<CharString> genomeFileNameList;
	int result = getGenomeFileNameList(genomeFileName, genomeFileNameList, options);
	if (result == RAZERS_GENOME_FAILED)
	{
		::std::cerr << "Failed to open genome file " << genomeFileName << ::std::endl;
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


//////////////////////////////////////////////////////////////////////////////
// Print usage
void printVersion() 
{
	string rev = "$Revision: 4546 $";
	rev.replace(rev.begin(), rev.end(), '$', ' ');
	cerr << "MicroRazerS version 0.1 20090710 (prerelease)" << rev << endl;
}

template <typename TSpec>
void printHelp(int, const char *[], RazerSOptions<TSpec> &options, ParamChooserOptions &, bool longHelp = false) 
{
	if (options.printVersion)
	{
		printVersion();
		return;
	}

	cerr << "********************************************************" << endl;
	cerr << "*** MicroRazerS - Rapid Alignment of Small RNA Reads ***" << endl;
	cerr << "***    written by Anne-Katrin Emde (c) Apr 2009    � ***" << endl;
	cerr << "********************************************************" << endl << endl;
	cerr << "Usage: micro_razers [OPTION]... <GENOME FILE> <READS FILE>" << endl;
	if (longHelp) {
		cerr << endl << "Main Options:" << endl;
		cerr << "  -o,  --output FILE           \t" << "change output filename (default <READS FILE>.result)" << endl;
		cerr << "  -sL, --seed-length NUM       \t" << "seed length (default "<<options.rnaSeedLength<<")" << endl;
		cerr << "  -sE, --seed-error            \t" << "allow for one error in the seed (default off)" << endl;
		cerr << "  -rr, --recognition-rate      \t" << "set lower bound of sensitivity level for one-error matches (default 100)" << endl;
		cerr << "  -f,  --forward               \t" << "only compute forward matches" << endl;
		cerr << "  -r,  --reverse               \t" << "only compute reverse complement matches" << endl;
		cerr << "  -mN, --match-N               \t" << "\'N\' matches with all other characters" << endl;
		cerr << "  -m,  --max-hits NUM          \t" << "output only NUM of the best hits (default " << options.maxHits << ')' << endl;
		cerr << "  -pa, --purge-ambiguous       \t" << "purge reads with more than max-hits best matches" << endl;
	//	cerr << "  -tr, --trim-reads NUM        \t" << "trim reads to length NUM (default off)" << endl;
		cerr << "  -lm, --low-memory            \t" << "may decrease memory usage at the expense of runtime" << endl;
		cerr << "  -v,  --verbose               \t" << "verbose mode" << endl;
		cerr << "  -vv, --vverbose              \t" << "very verbose mode" << endl;
		cerr << "  -V,  --version               \t" << "print version number" << endl;
		cerr << "  -h,  --help                  \t" << "print this help" << endl;
		cerr << endl << "Output Format Options:" << endl;
                cerr << "  -a,  --alignment             \t" << "dump the alignment for each match" << endl;
		cerr << "  -gn, --genome-naming NUM     \t" << "select how genomes are named" << endl;
		cerr << "                               \t" << "0 = use Fasta id (default)" << endl;
		cerr << "                               \t" << "1 = enumerate beginning with 1" << endl;
		cerr << "  -rn, --read-naming NUM       \t" << "select how reads are named" << endl;
		cerr << "                               \t" << "0 = use Fasta id (default)" << endl;
		cerr << "                               \t" << "1 = enumerate beginning with 1" << endl;
		cerr << "                               \t" << "2 = use the read sequence (only for short reads!)" << endl;
		cerr << "  -so, --sort-order NUM        \t" << "select how matches are sorted" << endl;
		cerr << "                               \t" << "0 = 1. read number, 2. genome position (default)" << endl;
		cerr << "                               \t" << "1 = 1. genome position, 2. read number" << endl;
		cerr << "  -pf, --position-format       \t" << "0 = gap space (default)" << endl;
		cerr << "                               \t" << "1 = position space" << endl;
	} else {
		cerr << "Try 'micro_razers --help' for more information." << endl;
	}
}


int main(int argc, const char *argv[]) 
{
	//////////////////////////////////////////////////////////////////////////////
	// Parse command line

	RazerSOptions<>		options;
	ParamChooserOptions	pm_options;

	bool				paramChoosing = true; //false;
	bool				paramChoosingError = false;
	unsigned			fnameCount = 0;
	string				errorPrbFileName;
	const unsigned			maxFiles = 2;
	const char			*fname[maxFiles + 1] = { NULL, };

	options.forward = false;
	options.reverse = false;
	options.hammingOnly = true;
	options.microRNA = true;
	pm_options.optionLossRate = 0.0;

	for(int arg = 1; arg < argc; ++arg) {
		if (argv[arg][0] == '-') {
			// parse options

			if (strcmp(argv[arg], "-f") == 0 || strcmp(argv[arg], "--forward") == 0) {
				options.forward = true;
				continue;
			}
			if (strcmp(argv[arg], "-r") == 0 || strcmp(argv[arg], "--reverse") == 0) {
				options.reverse = true;
				continue;
			}
                        if (strcmp(argv[arg], "-sE") == 0 || strcmp(argv[arg], "--seed-error") == 0) {
                                options.exactSeed = false;
                                continue;
                        }
                        if (strcmp(argv[arg], "-sL") == 0 || strcmp(argv[arg], "--seed-length") == 0) {
                                if (arg + 1 < argc) {
                                        ++arg;
                                        istringstream istr(argv[arg]);
                                        istr >> options.rnaSeedLength;
                                        if (!istr.fail())
                                        {
                                                if (options.rnaSeedLength < 10)
                                                        cerr << "Minimum seed length is 10" << endl << endl;
                                                else
                                                        continue;
                                        }
                                }
                                printHelp(argc, argv, options, pm_options);
                                return 0;
                        }

			if (strcmp(argv[arg], "-rr") == 0 || strcmp(argv[arg], "--recognition-rate") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> pm_options.optionLossRate;
					paramChoosing = true;
					if (!istr.fail())
					{
						if (pm_options.optionLossRate < 80 || pm_options.optionLossRate > 100)
							cerr << "Recognition rate must be a value between 80 and 100" << endl << endl;
						else
						{
							pm_options.optionLossRate = 100.0-pm_options.optionLossRate;
							pm_options.optionLossRate /= 100.0;
							continue;
						}
					}
				}
				printHelp(argc, argv, options, pm_options);
				return 0;
			}
			if (strcmp(argv[arg], "-pa") == 0 || strcmp(argv[arg], "--purge-ambiguous") == 0) {
				options.purgeAmbiguous = true;
				continue;
			}
			if (strcmp(argv[arg], "-m") == 0 || strcmp(argv[arg], "--max-hits") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> options.maxHits;
					if (!istr.fail()) 
					{
						if (options.maxHits < 1)
							cerr << "Maximum hits threshold must be greater than 0" << endl << endl;
						else
							continue;
					}
				}
				printHelp(argc, argv, options, pm_options);
				return 0;
			}
			if (strcmp(argv[arg], "-a") == 0 || strcmp(argv[arg], "--alignment") == 0) {
				options.dumpAlignment = true;
				continue;
			}
			if (strcmp(argv[arg], "-o") == 0 || strcmp(argv[arg], "--output") == 0) {
				if (arg + 1 == argc) {
					printHelp(argc, argv, options, pm_options);
					return 0;
				}
				++arg;
				options.output = argv[arg];
				continue;
			}
			if (strcmp(argv[arg], "-of") == 0 || strcmp(argv[arg], "--output-format") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> options.outputFormat;
					if (!istr.fail())
					{
						if (options.outputFormat > 1)
							cerr << "Invalid output format options." << endl << endl;
						else
							continue;
					}
				}
				printHelp(argc, argv, options, pm_options);
				return 0;
			}
			if (strcmp(argv[arg], "-so") == 0 || strcmp(argv[arg], "--sort-order") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> options.sortOrder;
					if (!istr.fail())
					{
						if (options.sortOrder > 1)
							cerr << "Invalid sort order options." << endl << endl;
						else
							continue;
					}
				}
				printHelp(argc, argv, options, pm_options);
				return 0;
			}
			if (strcmp(argv[arg], "-gn") == 0 || strcmp(argv[arg], "--genome-naming") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> options.genomeNaming;
					if (!istr.fail())
					{
						if (options.genomeNaming > 1)
							cerr << "Invalid genome naming options." << endl << endl;
						else
							continue;
					}
				}
				printHelp(argc, argv, options, pm_options);
				return 0;
			}
			if (strcmp(argv[arg], "-rn") == 0 || strcmp(argv[arg], "--read-naming") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> options.readNaming;
					if (!istr.fail())
					{
						if (options.readNaming > 2)
							cerr << "Invalid read naming options." << endl << endl;
						else
							continue;
					}
				}
				printHelp(argc, argv, options, pm_options);
				return 0;
			}
			if (strcmp(argv[arg], "-pf") == 0 || strcmp(argv[arg], "--position-format") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> options.positionFormat;
					if (!istr.fail())
					{
						if (options.positionFormat > 1)
							cerr << "Invalid position format options." << endl << endl;
						else
							continue;
					}
				}
				printHelp(argc, argv, options, pm_options);
				return 0;
			}
/*			if (strcmp(argv[arg], "-oc") == 0 || strcmp(argv[arg], "--overabundance-cut") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> options.abundanceCut;
					if (!istr.fail())
					{
						if (options.abundanceCut <= 0 || options.abundanceCut > 1)
							cerr << "Overabundance cut ratio must be a value >0 and <=1. Set to 1 to disable." << endl << endl;
						else
							continue;
					}
				}
				printHelp(argc, argv, options, pm_options);
				return 0;
			}
			if (strcmp(argv[arg], "-tr") == 0 || strcmp(argv[arg], "--trim-reads") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> options.trimLength;
					if (!istr.fail()) 
					{
						if (options.trimLength < 14)
							cerr << "Minimum read length is 14" << endl << endl;
						else
							continue;
					}
				}
				printHelp(argc, argv, options, pm_options);
				return 0;
			}*/
			if (strcmp(argv[arg], "-lm") == 0 || strcmp(argv[arg], "--low-memory") == 0) {
				options.lowMemory = true;
				continue;
			}

			if (strcmp(argv[arg], "-h") == 0 || strcmp(argv[arg], "--help") == 0) {
				// print help
				printHelp(argc, argv, options, pm_options, true);
				return 0;
			}
			if (strcmp(argv[arg], "-v") == 0 || strcmp(argv[arg], "--verbose") == 0) {
				options._debugLevel = max(options._debugLevel, 1);
				continue;
			}
			if (strcmp(argv[arg], "-vv") == 0 || strcmp(argv[arg], "--vverbose") == 0) {
				options._debugLevel = 3;
				continue;
			}
			if (strcmp(argv[arg], "-V") == 0 || strcmp(argv[arg], "--version") == 0) {
				options.printVersion = true;
				continue;
			}
			if (strcmp(argv[arg], "-mN") == 0 || strcmp(argv[arg], "--match-N") == 0) {
				options.matchN = true;
				continue;
			}
			if (strcmp(argv[arg], "-a") == 0 || strcmp(argv[arg], "--alignment") == 0) {
				options.dumpAlignment = true;
				continue;
			}
			cerr << "Unknown option: " << argv[arg] << endl << endl;
			printHelp(argc, argv, options, pm_options);
			return 0;
		} else {
			// parse file name
			if (fnameCount == maxFiles) {
				cerr << "More than " << maxFiles << " input files specified." << endl << endl;
				printHelp(argc, argv, options, pm_options);
				return 0;
			}
			fname[fnameCount++] = argv[arg];
		}
	}
	if (fnameCount < 2) {
		if (argc > 1 && !options.printVersion)
			cerr << "Less than 2 input files specified." << endl << endl;
		printHelp(argc, argv, options, pm_options);
		return 0;
	}
	if (!options.forward && !options.reverse) { // eistr >> options.shape;nable both per default
		options.forward = true;
		options.reverse = true;
	}

	if (fnameCount == 2)
		options.libraryLength = -1;	// only 1 readset -> disable mate-pair mapping
	
	if (options.printVersion)
		printVersion();
		
	// get read length
	int readLength = estimateReadLength(fname[1]);
	if (readLength == RAZERS_READS_FAILED)
	{
		::std::cerr << "Failed to open reads file " << fname[1] << ::std::endl;
		return 0;
	}
	if (readLength == 0) {
		::std::cerr << "Failed to read the first read sequence.";
		return 0;
	}

	if(options.outputFormat == 1) options.outputFormat = 4; // of 4 is sam in all RazerS tools
	options.distanceRange = 1; //? 1 actually means 0...
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


	int result = mapMicroRNAReads(fname[0], fname + 1, errorPrbFileName.c_str(), options);
	
	if (result == RAZERS_INVALID_SHAPE) 
	{
		printHelp(argc, argv, options, pm_options);
		return 0;
	}
	return result;
}
