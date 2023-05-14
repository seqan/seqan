#include <iostream>
#include <seqan/arg_parse.h>

int main()
{
    seqan2::ArgumentParser parser("base");

//![setMinMax]
    seqan2::addOption(parser, seqan2::ArgParseOption(
        "i", "integer-value", "An integer option",
        seqan2::ArgParseArgument::INTEGER, "INT"));

    seqan2::setMinValue(parser, "i", "10");
    seqan2::setMaxValue(parser, "integer-value", "20");
//![setMinMax]

//![setRequired]
    seqan2::addOption(parser, seqan2::ArgParseOption(
       "ir", "required-integer", "An required integer option",
       seqan2::ArgParseArgument::INTEGER, "INT"));

   	setRequired(parser, "ir");
//![setRequired]

//![setValidValues]
    seqan2::addOption(parser, seqan2::ArgParseOption(
        "", "distance-model", "Distance model, either HAMMING or EDIT.",
        seqan2::ArgParseArgument::STRING, "STR"));

    seqan2::setValidValues(parser, "distance-model", "HAMMING EDIT");
//![setValidValues]

//![addFileOption]
    seqan2::addOption(parser, seqan2::ArgParseOption(
        "I", "input-file", "Path to the input file",
        seqan2::ArgParseArgument::INPUT_FILE, "IN"));
    seqan2::addOption(parser, seqan2::ArgParseOption(
        "O", "output-file", "Path to the output file",
        seqan2::ArgParseArgument::OUTPUT_FILE, "OUT"));
//![addFileOption]

//![addFileExtension]
    seqan2::setValidValues(parser, "input-file", "txt");
    seqan2::setValidValues(parser, "output-file", "txt");
//![addFileExtension]

//![readFile]
    seqan2::CharString inputFileName, outputFileName;
    seqan2::getOptionValue(inputFileName, parser, "input-file");
    seqan2::getOptionValue(outputFileName, parser, "output-file");
//![readFile]

//![tupleOption]
    seqan2::addOption(parser, seqan2::ArgParseOption(
        "r", "range", "The range to modify.",
        seqan2::ArgParseArgument::INTEGER, "BEGIN END",
        false, 2));
//![tupleOption]

//![getTupleValue]
    unsigned rangeBegin = 0, rangeEnd = 0;
    seqan2::getOptionValue(rangeBegin, parser, "range", 0);
    seqan2::getOptionValue(rangeEnd, parser, "range", 1);
//![getTupleValue]

//![setVersion]
    seqan2::setShortDescription(parser, "String Modifier");
    seqan2::setVersion(parser, "1.0");
    seqan2::setDate(parser, "July 2012");
//![setVersion]

//![addUsageLine]
    seqan2::addUsageLine(parser,
        "[\\fIOPTIONS\\fP] \"\\fITEXT\\fP\"");
    seqan2::addDescription(parser,
        "This program allows simple character modifications to "
        "each i-th character.");
//![addUsageLine]

//![addSection]
    seqan2::addSection(parser, "Modification Options");
//![addSection]

//![addListItem]
    seqan2::addTextSection(parser, "Examples");

    seqan2::addListItem(parser,
        "\\fBmodify_string\\fP \\fB-U\\fP \\fIveryverylongword\\fP",
        "Print upper case version of \"veryverylongword\"");
    seqan2::addListItem(parser,
        "\\fBmodify_string\\fP \\fB-L\\fP \\fB-i\\fP \\fI3\\fP \\fIveryverylongword\\fP",
        "Print \"veryverylongword\" with every third character "
        "converted to upper case.");
//![addListItem]

    return 0;
}
