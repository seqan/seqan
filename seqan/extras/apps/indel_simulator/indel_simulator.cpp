#include <iostream>
#include <sstream>
#include <fstream>
#include <limits>



#include "indel_simulator.h"

using namespace std;
using namespace seqan;





//////////////////////////////////////////////////////////////////////////////
// Print usage
template<typename TOptions>
void printHelp(int, const char *[],TOptions & options, bool longHelp = false) 
{
	cerr << "**********************" << endl;
	cerr << "*** Indel simulator ***" << endl;
	cerr << "**********************" << endl << endl;
	cerr << "Usage: indelSimulator [OPTIONS]... <SOURCE SEQUENCE FILE>" << endl;
	cerr << "\n";
	if (longHelp) {
		cerr << "  -n,  --num-indels NUM         \t" << "total number of indels to simulate (" << options.numIndels << ")" << endl;
		cerr << "  -r,  --ranges FILE            \t" << "simulation ranges" << endl;
		cerr << "  -o,  --output FILE            \t" << "change output filename (default: <SOURCE SEQUENCE FILE>.indeled)" << endl;
		cerr << "  -oi,  --output-info FILE      \t" << "change output filename for simulation info (default: <output>.info)" << endl;
		cerr << "  -v,  --verbose                \t" << "verbose mode" << endl;
//		cerr << "  -vv, --very-verbose              \t" << "very verbose mode" << endl;
		cerr << "  -h,  --help                   \t" << "print this help" << endl;
	} else {
		cerr << "Try 'indelSimulator --help' for more information." << endl;
	}
}


//////////////////////////////////////////////////////////////////////////////
// Main part
int main(int argc, const char *argv[])
{
	srand(time(NULL));

	//unsigned fnameCount = 0;
	//const char *fname[1] = {""};
	StringSet<CharString> genomeFileNames;
	CharString rangeFile;
	IndelSimOptions<> options;


	CommandLineParser parser;

	//////////////////////////////////////////////////////////////////////////////
	// Define options
	addTitleLine(parser, "****************************************************'*");
	addTitleLine(parser, "*****    Simdel - Insertion/deletion Simulator   *****");
	addTitleLine(parser, "*****   (c) Copyright 2010 by Anne-Katrin Emde   *****");
	addTitleLine(parser, "******************************************************");
	addUsageLine(parser, "[OPTION]... <GENOME FILE> ");
	addOption(parser, CommandLineOption("n",  "num-indels",        "number of indels to simulate", OptionType::Int | OptionType::Label, options.numIndels));
	addOption(parser, addArgumentText(CommandLineOption("o",  "output",            "output filename for manipulated sequence ", OptionType::String), "FILE"));
	addOption(parser, addArgumentText(CommandLineOption("oi", "output-info",       "output filename for implanted indel information", OptionType::String), "FILE"));
	addOption(parser, addArgumentText(CommandLineOption("i",  "input-info",        "input filename that indels should be sampled from", OptionType::String), "FILE"));
	addOption(parser, addArgumentText(CommandLineOption("r",  "ranges",            "file containing ranges of indelsizes to simulate from (uniformly)", OptionType::String), "FILE"));
	addOption(parser, CommandLineOption("v",  "verbose",           "verbose mode", OptionType::Boolean));
	addOption(parser, CommandLineOption("vv", "vverbose",          "very verbose mode", OptionType::Boolean));

	bool stop = !parse(parser, argc, argv, cerr);
	if(stop) {printHelp(argc, argv, options); cout << "Exiting...\n" ;}

	getOptionValueLong(parser, "num-indels", options.numIndels);
	getOptionValueLong(parser, "output", options.output);
	getOptionValueLong(parser, "input-info", options.inputInfo);
	getOptionValueLong(parser, "output-info", options.outputInfo);
	getOptionValueLong(parser, "ranges", rangeFile);
	
	if (isSetLong(parser, "help") || isSetLong(parser, "version")) return 0;	// print help or version and exit
	if (isSetLong(parser, "verbose")) options._debugLevel = max(options._debugLevel, 1);
	if (isSetLong(parser, "vverbose")) options._debugLevel = max(options._debugLevel, 3);
	appendValue(genomeFileNames, getArgumentValue(parser, 0));


	if(argumentCount(parser) > 1) { std::cerr << "Too many arguments. Exiting... \n"; exit(0); }
	fstream file;
	clear(options.ranges);
	int absMaxValue = 0;
	file.open(toCString(rangeFile),ios_base::in | ios_base::binary);
	char c = _streamGet(file);
	while (!_streamEOF(file))
	{
		parse_skipWhitespace(file,c);
		int rangeBegin = static_cast<int>(parse_readDouble(file,c));
		if(abs(rangeBegin) > absMaxValue) absMaxValue = abs(rangeBegin);
		parse_skipWhitespace(file,c);
		int rangeEnd = static_cast<int>(parse_readDouble(file,c));
		if(abs(rangeEnd) > absMaxValue) absMaxValue = abs(rangeEnd);
		appendValue(options.ranges,Pair<int,int>(rangeBegin,rangeEnd));
	}
	options.minDistance = 2 * absMaxValue;
	options.noNsInRange = _min(30,options.minDistance);
	file.close();

	//////////////////////////////////////////////////////////////////////////////
	// Check options
	if (options.numIndels < 1 )
	{
		std::cerr << "NumIndels must be a value > 0 100" << std::endl;
		exit(0);
	}

/*
	// Command line parsing
	for(int arg = 1; arg < argc; ++arg) {
		if (argv[arg][0] == '-') {
			// parse option

			if (strcmp(argv[arg], "-n") == 0 || strcmp(argv[arg], "--num-indels") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> options.numIndels;
					if (!istr.fail())
					{
						if (options.numIndels < 1)
							cerr << "Num Indels must be a positive integer value" << endl << endl;
						else
							continue;
					}
				}
				printHelp(argc, argv, options);
				return 0;
			}
			if (strcmp(argv[arg], "-r") == 0 || strcmp(argv[arg], "--ranges") == 0) {
				if (arg + 1 < argc) {
					++arg;
					fstream file;
					clear(options.ranges);
					int absMaxValue = 0;
					file.open(argv[arg],ios_base::in | ios_base::binary);
					char c = _streamGet(file);
					while (!_streamEOF(file))
					{
						parse_skipWhitespace(file,c);
						int rangeBegin = parse_readDouble(file,c);
						if(abs(rangeBegin) > absMaxValue) absMaxValue = abs(rangeBegin);
						parse_skipWhitespace(file,c);
						int rangeEnd = parse_readDouble(file,c);
						if(abs(rangeEnd) > absMaxValue) absMaxValue = abs(rangeEnd);
						appendValue(options.ranges,Pair<int,int>(rangeBegin,rangeEnd));
					}
					options.minDistance = 2 * absMaxValue;
					options.noNsInRange = _min(30,options.minDistance);
					file.close();
					continue;
				}
				printHelp(argc, argv, options);
				return 0;
			}

			if (strcmp(argv[arg], "-o") == 0 || strcmp(argv[arg], "--output") == 0) {
				if (arg + 1 == argc) {
					printHelp(argc, argv, options);
					return 0;
				}
				++arg;
				options.output = argv[arg];
				continue;
			}
			if (strcmp(argv[arg], "-v") == 0 || strcmp(argv[arg], "--verbose") == 0) {
				options._debugLevel = max(options._debugLevel, 1);
				continue;
			}
			if (strcmp(argv[arg], "-vv") == 0 || strcmp(argv[arg], "--very-verbose") == 0) {
				options._debugLevel = max(options._debugLevel, 2);
				continue;
			}
			if (strcmp(argv[arg], "-oi") == 0 || strcmp(argv[arg], "--output-info") == 0) {
				if (arg + 1 == argc) {
					printHelp(argc, argv, options);
					return 0;
				}
				++arg;
				options.outputInfo = argv[arg];
				continue;
			}
			if (strcmp(argv[arg], "-h") == 0 || strcmp(argv[arg], "--help") == 0) {
				// print help
				printHelp(argc, argv, options, true);
				return 0;
			}
		}
		else {
			// parse file name
			if (fnameCount == 1) {
				printHelp(argc, argv, options);
				return 0;
			}
			fname[fnameCount++] = argv[arg];
		}
	}
	if (fnameCount < 1) {
		printHelp(argc, argv, options);
		return 0;
	}
*/	
	::std::map<CharString,unsigned> gIdStringToIdNumMap;
	StringSet<CharString> genomeIDs;
	StringSet<Dna5String> genomes;
	
	//if (!loadFasta(genomes, genomeIDs, fname[0])) 
	if (!loadGenomes(toCString(genomeFileNames[0]),genomes,genomeIDs,gIdStringToIdNumMap,options)) 
	{
		cerr << "Reference file " << genomeFileNames[0] << " can't be loaded." << endl;
//		cerr << "Reference file " << fname[0] << " can't be loaded." << endl;
		return 0;
	}
	if(options._debugLevel > 0 )
	{
		cout << "\n"<<lengthSum(genomes) << " bps of " << length(genomes) << " source sequence loaded." << endl;
		::std::cout << "Number of range buckets: " << length(options.ranges) << endl;
	}
	if(options._debugLevel > 1 )
	{
		for(unsigned i = 0; i < length(options.ranges); ++i)
			::std::cout << "["<<options.ranges[i].i1 << "," <<options.ranges[i].i2 <<")"   << endl;
	}
//____________________________________________________________________________
	
	// Main Part - Simulation
	
	StringSet<Dna5String>		simGenomes;
	StringSet<CharString>		simIDs;
	::std::map<int,IndelInfo>	indelInfos;


	int result = simulateIndels(genomes,genomeIDs,gIdStringToIdNumMap,simGenomes,simIDs,indelInfos,options);
	if(result > 0)
	{
		cerr << "Something went wrong.. Exiting..\n";
		return 1;
	}
//____________________________________________________________________________
	
	// output
	//saveFasta(simGenomes, simIDs, options, fname[0]); 
	//saveIndelInfos(indelInfos, genomeIDs, options, fname[0]); 
	saveFasta(simGenomes, simIDs, options, toCString(genomeFileNames[0])); 
	saveIndelInfos(indelInfos, genomeIDs, options, toCString(genomeFileNames[0])); 

	return 0;
}
