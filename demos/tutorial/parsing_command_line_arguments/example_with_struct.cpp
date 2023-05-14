//![full]
#include <iostream>

#include <seqan/arg_parse.h>

//![struct]
struct ModifyStringOptions
{
    unsigned period;
    bool toUppercase;
    bool toLowercase;
    seqan2::CharString text;

    ModifyStringOptions() :
    period(1), toUppercase(false), toLowercase(false)
    {}
};
//![struct]

int main(int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan2::ArgumentParser parser("modify_string");

    addArgument(parser, seqan2::ArgParseArgument(
        seqan2::ArgParseArgument::STRING, "TEXT"));

    addOption(parser, seqan2::ArgParseOption(
        "i", "period", "Period to use for the index.",
        seqan2::ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "period", "1");
    addOption(parser, seqan2::ArgParseOption(
        "U", "uppercase", "Select to-uppercase as operation."));
    addOption(parser, seqan2::ArgParseOption(
        "L", "lowercase", "Select to-lowercase as operation."));

    // Parse command line.
    seqan2::ArgumentParser::ParseResult res = seqan2::parse(parser, argc, argv);

    // If parsing was not successful then exit with code 1 if there were errors.
    // Otherwise, exit with code 0 (e.g. help was printed).
    if (res != seqan2::ArgumentParser::PARSE_OK)
        return res == seqan2::ArgumentParser::PARSE_ERROR;

    // Extract option values and print them.
    ModifyStringOptions options;
    getOptionValue(options.period, parser, "period");
    options.toUppercase = isSet(parser, "uppercase");
    options.toLowercase = isSet(parser, "lowercase");
    getArgumentValue(options.text, parser, 0);

    std::cout << "period   \t" << options.period << '\n'
              << "uppercase\t" << options.toUppercase << '\n'
              << "lowercase\t" << options.toLowercase << '\n'
              << "text     \t" << options.text << '\n';

    return 0;
}
//![full]
