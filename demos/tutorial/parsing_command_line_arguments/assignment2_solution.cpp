#include <iostream>

#include <seqan/arg_parse.h>

int main(int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan2::ArgumentParser parser("modify_string");

    addArgument(parser, seqan2::ArgParseArgument(
    seqan2::ArgParseArgument::STRING, "TEXT"));

    addOption(parser, seqan2::ArgParseOption(
        "i", "period", "Period to use for the index.",
        seqan2::ArgParseArgument::INTEGER, "INT"));
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
    unsigned period = 0;
    getOptionValue(period, parser, "period");
    bool toUppercase = isSet(parser, "uppercase");
    bool toLowercase = isSet(parser, "lowercase");
    seqan2::CharString text;
    getArgumentValue(text, parser, 0);

    std::cout << "period   \t" << period << '\n'
              << "uppercase\t" << toUppercase << '\n'
              << "lowercase\t" << toLowercase << '\n'
              << "text     \t" << text << '\n';

    return 0;
}
