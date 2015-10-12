#include <iostream>

#include <seqan/arg_parse.h>

struct ModifyStringOptions
{
    unsigned period;
    bool toUppercase;
    bool toLowercase;
    seqan::CharString text;

    ModifyStringOptions() :
    period(1), toUppercase(false), toLowercase(false)
    {}
};

seqan::ArgumentParser::ParseResult
parseCommandLine(ModifyStringOptions & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("modify_string");

    // We require one argument.
    addArgument(parser, seqan::ArgParseArgument(
        seqan::ArgParseArgument::STRING, "TEXT"));

    // Define Options
    addOption(parser, seqan::ArgParseOption(
        "i", "period", "Period to use for the index.",
        seqan::ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "period", "1");
    addOption(parser, seqan::ArgParseOption(
        "U", "uppercase", "Select to-uppercase as operation."));
    addOption(parser, seqan::ArgParseOption(
        "L", "lowercase", "Select to-lowercase as operation."));

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    // Extract option values.    
    getOptionValue(options.period, parser, "period");
    options.toUppercase = isSet(parser, "uppercase");
    options.toLowercase = isSet(parser, "lowercase");
    getArgumentValue(options.text, parser, 0);

    // If both to-uppercase and to-lowercase were selected then this is an error.
    if (options.toUppercase && options.toLowercase)
    {
        std::cerr << "ERROR: You cannot specify both to-uppercase and to-lowercase!\n";
        return seqan::ArgumentParser::PARSE_ERROR;
    }

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

    std::cout << "period   \t" << options.period << '\n'
              << "uppercase\t" << options.toUppercase << '\n'
              << "lowercase\t" << options.toLowercase << '\n'
              << "text     \t" << options.text << '\n';

    return 0;
}

