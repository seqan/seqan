#include <fstream>
#include <iostream>
#include <sstream>

#define SEQAN_PROFILE
#ifndef RELEASE	
//#define SEQAN_DEBUG			
//#define SEQAN_TEST	
#endif

#include <string>
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/map.h>
#include <seqan/refinement.h>
#include <seqan/store.h>
//#include <seqan/misc/misc_cmdparser.h>

#include "base.h"
#include "read_gff.h"
#include "read_gtf.h"
#include "create_gff.h"
#include "fusion.h"
#include "overlap_module.h"

#include <seqan/arg_parse.h>

using namespace seqan;
using namespace std;

struct InsegtOptions
{
	CharString nameSAM;
	CharString nameGFF;
	CharString outputPath;
	unsigned nTuple;
	unsigned offsetInterval;
	unsigned thresholdGaps;
	unsigned thresholdCount;
	double thresholdRPKM;
	bool maxTuple;
	bool exact_nTuple;
	bool unknownO;
	bool fusion;
	bool gtf;

    InsegtOptions() :
        outputPath(""),
        nTuple(2),
        offsetInterval(5),
        thresholdGaps(5),
        thresholdCount(1),
        thresholdRPKM(0.0),
        maxTuple(false),
        exact_nTuple(false),
        unknownO(false),
        fusion(false),
        gtf(false)
        {}
};
 

ArgumentParser::ParseResult
parseCommandLine(InsegtOptions & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    ArgumentParser parser("insegt");
    // Set short description, version, and date.
    setShortDescription(parser, "INtersecting SEcond Generation sequencing daTa with annotation");
    setVersion(parser, "1.0");
    string date = "$Date: 2012-09-11 11:21:13 +0200 (Mo, 11. Sep 2012) $";
    setDate(parser, date.substr(7, _min((int)date.size() - 8, 10)));
    // Define usage line and long description.
    addUsageLine(parser,
                 "[\\fIOPTIONS\\fP] <\\fIALIGMENTS-FILE\\fP> <\\fIANNOTATIONS-FILE\\fP> ");
    addDescription(parser,
                   "INSEGT is a tool to analyze alignments of RNA-Seq reads "
                   "(single-end or paired-end) by using gene-annotations.");
    // We require two arguments.
    addDescription(parser, "Input to INSEGT is a SAM file containing the alignments and"
                          " a file containing the annotations of the reference genome, either in GFF or GTF format (the latter has additionally to be specified with option 'z').");
    addArgument(parser, ArgParseArgument(
        ArgParseArgument::INPUTFILE, "IN"));
    addArgument(parser, ArgParseArgument(
        ArgParseArgument::INPUTFILE, "IN"));

    // Define Options -- Section Modification Options

    addSection(parser, "Options: ");
    addOption(parser, ArgParseOption("p", "output-path", "Path for output-files.", ArgParseArgument::STRING, "TEXT")); 
    setDefaultValue(parser, "output-path", "");
    addOption(parser, ArgParseOption("n", "ntuple", "ntuple", ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "ntuple", "2");
    addOption(parser, ArgParseOption("o", "offset-interval", "Offset to short alignment-intervals for search.", ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "offset-interval", "5");
    addOption(parser, ArgParseOption("t", "threshold-gaps", "Threshold for allowed gaps in alignment (not introns).", ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "threshold-gaps", "5");
    addOption(parser, ArgParseOption("c", "threshold-count", "Threshold for min. count of tuple for output.", ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "threshold-count", "1");
    addOption(parser, ArgParseOption("r", "threshold-rpkm", "Threshold for min. RPKM of tuple for output.", ArgParseArgument::DOUBLE, "DOUBLE"));
    setDefaultValue(parser, "threshold-rpkm", "0.0");
    addOption(parser, ArgParseOption("m", "max-tuple", "Create only maxTuple (which are spanned by the whole read)."));
    addOption(parser, ArgParseOption("e", "exact-ntuple", "Create only Tuple of exact length n. By default all tuple up to the given length are computed (if -m is set, -e will be ignored)."));
    addOption(parser, ArgParseOption("u", "unknown-orientation", "Orientation of reads is unknown."));
    addOption(parser, ArgParseOption("f", "fusion-genes", "Check for fusion genes and create separate output for matepair tuple."));
    addOption(parser, ArgParseOption("z", "gtf-format", "GTF format as input for annotations (instead of GFF format)."));

    // Add Examples Section.
    addTextSection(parser, "Examples");
    addListItem(parser,
                "\\fBinsegt\\fP  \\fBexample/alignments.sam\\fP \\fBexample/annotations.gff\\fP",
                "Run INSEGT on example files with default parameters.");
    addListItem(parser,
                "\\fBinsegt\\fP \\fB-m\\fP \\fBexample/alignments.sam\\fP \\fBexample/annotations.gff\\fP",
                "Run INSEGT on example files and only compute maxTuple.");
    addListItem(parser,
                "\\fBinsegt\\fP \\fB-c\\fP \\fB2\\fP \\fBexample/alignments.sam\\fP \\fBexample/annotations.gff\\fP",
                "Run INSEGT on example files and only output tuple with a min. count of 2.");


    // Parse command line.
    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    // Only extract  options if the program will continue after parseCommandLine()
    if (res != ArgumentParser::PARSE_OK)
        return res;

    // Extract option values.
    getArgumentValue(options.nameSAM, parser, 0);
    getArgumentValue(options.nameGFF, parser, 1);

    getOptionValue(options.outputPath, parser, "output-path");
    getOptionValue(options.nTuple, parser, "ntuple");
    getOptionValue(options.offsetInterval, parser, "offset-interval");
    getOptionValue(options.thresholdGaps, parser, "threshold-gaps");
    getOptionValue(options.thresholdCount, parser, "threshold-count");
    getOptionValue(options.thresholdRPKM, parser, "threshold-rpkm");

    options.maxTuple = isSet(parser, "max-tuple");
    options.exact_nTuple = isSet(parser, "exact-ntuple");
    options.unknownO = isSet(parser, "unknown-orientation");
    options.fusion = isSet(parser, "fusion-genes");
    options.gtf = isSet(parser, "gtf-format");

    // If were selected then this is an error.
    if (options.maxTuple && options.exact_nTuple)
    {
        std::cerr << "ERROR: You cannot specify both max-tuple and exact-ntuple!\n";
        return ArgumentParser::PARSE_ERROR;
    }
    return ArgumentParser::PARSE_OK;
}

///////////////////////////////////////////////////////////////////////////////
////// main
///////////////////////////////////////////////////////////////////////////////

int main( int argc, const char *argv[] ) 
{
    // Parse the command line.
    InsegtOptions options;
    ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
    // If parsing was not successful then exit with code 1 if there were errors.
    // Otherwise, exit with code 0 (e.g. help was printed).
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

	if (options.maxTuple) 
	{
		options.nTuple = 0;		// sign for maxTuple
		options.exact_nTuple = 0;	// not both possible: maxTuple is prefered over exact_nTuple and n
	}

	ngsOverlapper(options.nameSAM, options.nameGFF, options.outputPath, options.nTuple, options.exact_nTuple, options.thresholdGaps, options.offsetInterval, options.thresholdCount, options.thresholdRPKM, options.unknownO, options.fusion, options.gtf);
	return 0;
}

