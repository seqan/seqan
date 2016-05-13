#include <iostream>

#include <seqan/arg_parse.h>

struct ModifyStringOptions
{
    unsigned period;

    ModifyStringOptions() :
    period(1)
    {}
};

seqan::ArgumentParser::ParseResult
parseCommandLine(ModifyStringOptions & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    //seqan::ArgumentParser parser("FeatureFinderCentroided");
    seqan::ArgumentParser parser("FeatureFinderCentroided");
    setVersion(parser, "1.2.3");
    setUrl(parser, "http://open-ms.sourceforge.net/downloads/");

    // Define Options
    addOption(parser, seqan::ArgParseOption(
        "i", "period", "Period to use for the index.",
        seqan::ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "period", "1");

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    // Extract option values.    
    getOptionValue(options.period, parser, "period");

    return seqan::ArgumentParser::PARSE_OK;
}

int main(int argc, char const ** argv)
{
    // Parse the command line.
    ModifyStringOptions options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    // If parsing was not successful then exit with code 1 if there were errors.
    // Otherwise, exit with code 0 (e.g. help was printed).
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;


    return 0;
}
