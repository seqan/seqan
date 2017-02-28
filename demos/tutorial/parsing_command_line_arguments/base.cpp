#include <iostream>
#include <seqan/arg_parse.h>

int main()
{
    seqan::ArgumentParser parser("base");

//![setMinMax]
    seqan::addOption(parser, seqan::ArgParseOption(
        "i", "integer-value", "An integer option",
        seqan::ArgParseArgument::INTEGER, "INT"));
 
    seqan::setMinValue(parser, "i", "10");
    seqan::setMaxValue(parser, "integer-value", "20");
//![setMinMax]

//![setRequired]
    seqan::addOption(parser, seqan::ArgParseOption(
       "ir", "required-integer", "An required integer option",
       seqan::ArgParseArgument::INTEGER, "INT"));

   	setRequired(parser, "ir");
//![setRequired]

//![setValidValues]
    seqan::addOption(parser, seqan::ArgParseOption(
        "", "distance-model", "Distance model, either HAMMING or EDIT.",
        seqan::ArgParseArgument::STRING, "STR"));

    seqan::setValidValues(parser, "distance-model", "HAMMING EDIT");
//![setValidValues]

//![addFileOption]
    seqan::addOption(parser, seqan::ArgParseOption(
        "I", "input-file", "Path to the input file",
        seqan::ArgParseArgument::INPUT_FILE, "IN"));
    seqan::addOption(parser, seqan::ArgParseOption(
        "O", "output-file", "Path to the output file",
        seqan::ArgParseArgument::OUTPUT_FILE, "OUT"));
//![addFileOption]

//![addFileExtension]
    seqan::setValidValues(parser, "input-file", "txt");
    seqan::setValidValues(parser, "output-file", "txt");
//![addFileExtension]

//![readFile]
    seqan::CharString inputFileName, outputFileName;
    seqan::getOptionValue(inputFileName, parser, "input-file");
    seqan::getOptionValue(outputFileName, parser, "output-file");
//![readFile]

//![tupleOption]
    seqan::addOption(parser, seqan::ArgParseOption(
        "r", "range", "The range to modify.",
        seqan::ArgParseArgument::INTEGER, "BEGIN END",
        false, 2));
//![tupleOption]

//![getTupleValue]
    unsigned rangeBegin = 0, rangeEnd = 0;
    seqan::getOptionValue(rangeBegin, parser, "range", 0);
    seqan::getOptionValue(rangeEnd, parser, "range", 1);
//![getTupleValue]

//![setVersion]
    seqan::setShortDescription(parser, "String Modifier");
    seqan::setVersion(parser, "1.0");
    seqan::setDate(parser, "July 2012");
//![setVersion]

//![addUsageLine]
    seqan::addUsageLine(parser,
        "[\\fIOPTIONS\\fP] \"\\fITEXT\\fP\"");
    seqan::addDescription(parser,
        "This program allows simple character modifications to "
        "each i-th character.");
//![addUsageLine]

//![addSection]
    seqan::addSection(parser, "Modification Options");
//![addSection]

//![addListItem]
    seqan::addTextSection(parser, "Examples");

    seqan::addListItem(parser,
        "\\fBmodify_string\\fP \\fB-U\\fP \\fIveryverylongword\\fP",
        "Print upper case version of \"veryverylongword\"");
    seqan::addListItem(parser,
        "\\fBmodify_string\\fP \\fB-L\\fP \\fB-i\\fP \\fI3\\fP \\fIveryverylongword\\fP",
        "Print \"veryverylongword\" with every third character "
        "converted to upper case.");
//![addListItem]

    return 0;
}
