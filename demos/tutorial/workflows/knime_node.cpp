#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

#include <seqan/arg_parse.h>

using namespace seqan;

// ==========================================================================
// Classes
// ==========================================================================

// --------------------------------------------------------------------------
// Class KnimeNodeOptions
// --------------------------------------------------------------------------

// This struct stores the options from the command line.
//
// You might want to rename this to reflect the name of your app.

struct KnimeNodeOptions
{
    // Verbosity level.  0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    int verbosity;

    // The arguments of the program are stored here.
    CharString inputFile;
    CharString outputFile;

    KnimeNodeOptions() :
        verbosity(1)
    {}
};

// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function parseCommandLine()
// --------------------------------------------------------------------------

ArgumentParser::ParseResult
parseCommandLine(KnimeNodeOptions & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    ArgumentParser parser("knime_node");
    // Set short description, version, and date.
    setShortDescription(parser, "This is a very simple KNIME node providing an input and output port.");
    setVersion(parser, "0.1");
    setDate(parser, "Sep 2013");

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \"\\fITEXT\\fP\"");
    addDescription(parser, "This is a very simple KNIME node providing an input and output port. The code should be modified such that the node does something useful");

    // We require one argument.
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "IN"));
    setValidValues(parser, 0, SeqFileIn::getFileExtensions());

    addOption(parser, ArgParseOption("o", "outputFile", "Name of the multi-FASTA output.", ArgParseOption::OUTPUT_FILE, "OUT"));
    setValidValues(parser, "outputFile", SeqFileOut::getFileExtensions());
    setDefaultValue(parser, "outputFile", "result.fastq");

    // The verbosity option should be used to help debugging
    addOption(parser, ArgParseOption("q", "quiet", "Set verbosity to a minimum."));
    addOption(parser, ArgParseOption("v", "verbose", "Enable verbose output."));
    addOption(parser, ArgParseOption("vv", "very-verbose", "Enable very verbose output."));

    // Add Examples Section.
    addTextSection(parser, "Examples");
    addListItem(parser, "\\fBknime_node\\fP \\fB-v\\fP \\fItext\\fP",
                "Call with \\fITEXT\\fP set to \"text\" with verbose output.");

    // Parse command line.
    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != ArgumentParser::PARSE_OK)
        return res;

    // Extract option values.
    if (isSet(parser, "quiet"))
        options.verbosity = 0;
    if (isSet(parser, "verbose"))
        options.verbosity = 2;
    if (isSet(parser, "very-verbose"))
        options.verbosity = 3;

    getArgumentValue(options.inputFile, parser, 0);

    // Get output file name from command line if set.  Otherwise, autogenerate from input file name.
    if (isSet(parser, "outputFile"))
    {
        getOptionValue(options.outputFile, parser, "outputFile");
    }
    else
    {
        options.outputFile = options.inputFile;
        append(options.outputFile, ".fastq");
    }

    return ArgumentParser::PARSE_OK;
}

// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

// Program entry point.

int main(int argc, char const ** argv)
{
    // Parse the command line.
    ArgumentParser parser;
    KnimeNodeOptions options;
    ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    // If there was an error parsing or built-in argument parser functionality
    // was triggered then we exit the program.  The return code is 1 if there
    // were errors and 0 if there were none.
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    std::cout << "EXAMPLE PROGRAM\n"
              << "===============\n\n";

    // Print the command line arguments back to the user.
    if (options.verbosity > 0)
    {
        std::cout << "__OPTIONS____________________________________________________________________\n"
                  << '\n'
                  << "VERBOSITY\t" << options.verbosity << '\n'
                  << "INPUT_FILE\t" << options.inputFile << "\n\n"
                  << "OUTPUT_FILE\t" << options.outputFile << "\n\n";
    }

    // Reading the input
    SeqFileIn seqIn(toCString(options.inputFile));
    StringSet<CharString> ids;
    StringSet<Dna5String> seqs;
    StringSet<CharString> quals;


    if (atEnd(seqIn))
    {
        std::cout << "ERROR: File does not contain any sequences!\n";
        return 1;
    }
    readRecords(ids, seqs, quals, seqIn);

    // DO SOMETHING HERE
    // *
    // *
    // *

    SeqFileOut seqOut(toCString(options.outputFile));
    writeRecords(seqOut, ids, seqs, quals);

    return 0;
}
