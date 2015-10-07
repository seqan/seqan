#include <iostream>

#include <seqan/arg_parse.h>

int main(int argc, char const ** argv)
{
    //![object]
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("modify_string");
    //![object]

    //![argument]
    addArgument(parser, seqan::ArgParseArgument(
        seqan::ArgParseArgument::STRING, "TEXT"));
    //![argument]

    //![option]
    addOption(parser, seqan::ArgParseOption(
        "i", "period", "Period to use for the index.",
        seqan::ArgParseArgument::INTEGER, "INT"));
    addOption(parser, seqan::ArgParseOption(
        "U", "uppercase", "Select to-uppercase as operation."));
    //![option]

    //![parse]
    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);
    //![parse]

    // If parsing was not successful then exit with code 1 if there were errors.
    // Otherwise, exit with code 0 (e.g. help was printed).
    //![check]
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;
    //![check]

    //![print]
    // Extract option values and print them.
    unsigned period = 0;
    getOptionValue(period, parser, "period");
    bool toUppercase = isSet(parser, "uppercase");
    seqan::CharString text;
    getArgumentValue(text, parser, 0);

    std::cout << "period   \t" << period << '\n'
              << "uppercase\t" << toUppercase << '\n'
              << "text     \t" << text << '\n';
    //![print]

    return 0;
}

