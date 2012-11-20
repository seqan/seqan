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
		cerr << "  -ns, --num-snps NUM           \t" << "total number of SNPs to simulate (" << options.numIndels << ")" << endl;
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


	ArgumentParser parser("indelSimulator");

    addUsageLine(parser, "[\\fIOPTIONS\\fP] <\\fIGENOME FILE\\fP>");
    addDescription(parser,
            "IndelSimulator simulates indels and optionally SNPs into the sequence(s) specified. "
            "SNPs are simulated randomly, indels can be simulated either randomly (uniformly within specified size ranges) or from a "
            "database of known indels (given as GFF file).");
        
    addDescription(parser, "(c) Copyright 2010 by Anne-Katrin Emde.");
    setVersion(parser, "1.0");
    setDate(parser, "May 2011" );

//    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE));
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE, "GENOME", true));
   
    addSection(parser, "Main Options:");

	addOption(parser, ArgParseOption("ns",  "num-snps",        "Number of SNPs to simulate", ArgParseOption::INTEGER));
    setMinValue(parser, "num-snps", "0");   
    setDefaultValue(parser, "num-snps", options.numSnps);

	addOption(parser, ArgParseOption("n",  "num-indels",        "Number of indels to simulate", ArgParseOption::INTEGER));
    setMinValue(parser, "num-indels", "0");   
    setDefaultValue(parser, "num-indels", options.numIndels);
 
    addOption(parser, ArgParseOption("i", "input-indel", "GFF file that indels should be sampled from", ArgParseOption::INPUTFILE));

    addOption(parser, ArgParseOption("r", "ranges", "File containing ranges of indel sizes to simulate (specifying [rangeBegin,rangeEnd) intervals, see example file ranges.txt)", ArgParseOption::INPUTFILE));

	addOption(parser, ArgParseOption("d",  "diploid",           "Simulate two haplotypes"));

    addSection(parser, "Output Options:");

    addOption(parser, ArgParseOption("o", "output", "Output filename for manipulated sequence", ArgParseOption::OUTPUTFILE));
    addOption(parser, ArgParseOption("oi", "output-indel", "Output filename for simulated indels", ArgParseOption::OUTPUTFILE));
    addOption(parser, ArgParseOption("os", "output-snp", "Output filename for simulated SNPs", ArgParseOption::OUTPUTFILE));

	addOption(parser, ArgParseOption("v",  "verbose",           "verbose mode"));
	addOption(parser, ArgParseOption("vv", "vverbose",          "very verbose mode"));

    // Parse command line.
    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
    {
        cerr << "Exiting ..." << endl;
    }


	getOptionValue(options.numIndels, parser, "num-indels");
	getOptionValue(options.numSnps, parser, "num-snps");
	getOptionValue(options.output, parser, "output");
	getOptionValue(options.inputIndel, parser, "input-indel");
	getOptionValue(options.outputIndel, parser, "output-indel");
	getOptionValue(options.outputSnp, parser, "output-snp");
	getOptionValue(rangeFile, parser, "ranges");
	getOptionValue(options.diploid, parser, "diploid");
	
	if (isSet(parser, "help") || isSet(parser, "version")) return 0;	// print help or version and exit
	if (isSet(parser, "verbose")) options._debugLevel = max(options._debugLevel, 1);
	if (isSet(parser, "vverbose")) options._debugLevel = max(options._debugLevel, 3);
	if(getArgumentValueCount(parser, 0) > 1) { std::cerr << "Too many arguments. Exiting... \n"; exit(0); }
    if(getArgumentValueCount(parser, 0) < 1) { std::cerr << "Not enough arguments. Exiting... \n"; exit(0); }
    resize(genomeFileNames, length(genomeFileNames) + 1);
    getArgumentValue(back(genomeFileNames), parser, 0, 0);
   
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
	options.minDistance = 2 * absMaxValue;                // minimum distance between indel varians, should be an option, and not hard-coded..
	options.noNsInRange = _min(30,options.minDistance);   // same here, distance of simulated indel to masked/unkonwn "N"-regions in genome
	file.close();

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
	
	StringSet<Dna5String>		simGenomes; // if diploid --> twice as many as original genomes
	StringSet<CharString>		simIDs;
	::std::map<int,IndelInfo>	indelInfos;
	::std::map<int,SnpInfo>	    snpInfos;


	int result = simulateIndels(genomes,genomeIDs,gIdStringToIdNumMap,simGenomes,simIDs,indelInfos,snpInfos,options);
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
	saveSnpInfos(snpInfos, genomeIDs, options, toCString(genomeFileNames[0])); 

	return 0;
}
