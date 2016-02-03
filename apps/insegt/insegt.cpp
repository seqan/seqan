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
#include <seqan/graph_align.h>
#include <seqan/store.h>

#include "base.h"
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
	CharString readOutputFileName;
	CharString annoOutputFileName;
	CharString tupleOutputFileName;
	CharString tupleFusionOutputFileName;
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
    setCategory(parser, "Utilities");
    setVersion(parser, SEQAN_APP_VERSION " [" SEQAN_REVISION "]");
    setDate(parser, SEQAN_DATE);
    // Define usage line and long description.
    addUsageLine(parser,
                 "[\\fIOPTIONS\\fP] <\\fIALIGMENTS-FILE\\fP> <\\fIANNOTATIONS-FILE\\fP> ");
    addDescription(parser,
                   "INSEGT is a tool to analyze alignments of RNA-Seq reads "
                   "(single-end or paired-end) by using gene-annotations.");
    // We require two arguments.
    addDescription(parser, "Input to INSEGT is a SAM file containing the alignments and"
                          " a file containing the annotations of the reference genome, either in GFF or GTF format.");
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE));
    setValidValues(parser, 0, "sam");
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE));
    setValidValues(parser, 1, "gff gtf");
    // Define Options -- Section Modification Options

    addSection(parser, "Options: ");
    addOption(parser, ArgParseOption("ro", "read-output", "Output filename for read-output, which contains the mapped annotations followed by their parent annotation.", ArgParseArgument::OUTPUT_FILE)); 
    //setDefaultValue(parser, "read-output", "readOutput.gff");
    setValidValues(parser, "read-output", "gff");
    addOption(parser, ArgParseOption("ao", "anno-output", "Output filename for anno-output, which contains the annotations similar to the GFF input and additionally the counts of the mapped reads and the normalized expression levels in RPKM.", ArgParseArgument::OUTPUT_FILE)); 
    //setDefaultValue(parser, "anno-output", "annoOutput.gff");
    setValidValues(parser, "anno-output", "gff");
    addOption(parser, ArgParseOption("to", "tuple-output", "Output filename for tuple-output, which contains exon tuples connected by reads or matepairs.", ArgParseArgument::OUTPUT_FILE)); 
    //setDefaultValue(parser, "tuple-output", "tupleOutput.gff");
    setValidValues(parser, "tuple-output", "gff");
    // Check for gene fusions: currently disabled for KNIME
    addOption(parser, ArgParseOption("fo", "fusion-output", "Output filename for fusion-output, which contains exon tuple of gene fusions (Advanced option, currently no output port for KNIME).", ArgParseArgument::STRING)); 
    //setDefaultValue(parser, "fusion-output", "tupleFusionOutput.gff");
    setValidValues(parser, "fusion-output", "gff");
    //hideOption(parser, "fo");

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
    //addOption(parser, ArgParseOption("f", "fusion-genes", "Check for fusion genes and create separate output for matepair tuple."));

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
    if (endsWith(options.nameGFF, ".gff"))
    {
        options.gtf = false;
    }
    else if (endsWith(options.nameGFF, ".gtf"))
    {
        options.gtf = true;
    }
    else 
    {
        std::cerr << "ERROR: Input format of annotation file must be either GFF or GTF!\n";
        return ArgumentParser::PARSE_ERROR;
    }

    getOptionValue(options.readOutputFileName, parser, "read-output");
    getOptionValue(options.annoOutputFileName, parser, "anno-output");
    getOptionValue(options.tupleOutputFileName, parser, "tuple-output");
    getOptionValue(options.tupleFusionOutputFileName, parser, "fusion-output");
    getOptionValue(options.nTuple, parser, "ntuple");
    getOptionValue(options.offsetInterval, parser, "offset-interval");
    getOptionValue(options.thresholdGaps, parser, "threshold-gaps");
    getOptionValue(options.thresholdCount, parser, "threshold-count");
    getOptionValue(options.thresholdRPKM, parser, "threshold-rpkm");

    options.maxTuple = isSet(parser, "max-tuple");
    options.exact_nTuple = isSet(parser, "exact-ntuple");
    options.unknownO = isSet(parser, "unknown-orientation");

    options.fusion = isSet(parser, "fusion-output");
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
        options.nTuple = 0;         // sign for maxTuple
        options.exact_nTuple = 0;   // not both possible: maxTuple is prefered over exact_nTuple and n
    }

    ngsOverlapper(options);
    return 0;
}

