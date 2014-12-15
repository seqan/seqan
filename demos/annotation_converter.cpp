///A tutorial about importing/exporting annotations
#include <fstream>
#include <iostream>
#include <seqan/store.h>
#include <seqan/arg_parse.h>

using namespace seqan;

int main(int argc, const char * argv[])
{
    ArgumentParser parser;

    //////////////////////////////////////////////////////////////////////////////
    // Define options
    addDescription(parser, "****************************");
    addDescription(parser, "*** Annotation Converter ***");
    addDescription(parser, "****************************");
    addDescription(parser, "");
    addOption(parser, ArgParseOption("rg", "read-gff", "Read annotation in Gff/Gtf format.", ArgParseOption::STRING));
    addOption(parser, ArgParseOption("ru", "read-ucsc", "Read annotation in Ucsc format (knownGenes.txt knownIsoforms.txt).", ArgParseOption::STRING, "TXT TXT", true, 2));
    addDescription(parser, "");
    addOption(parser, ArgParseOption("wg", "write-gff", "Write annotation in Gff/Gtf format.", ArgParseOption::STRING));
    addOption(parser, ArgParseOption("wu", "write-ucsc", "Write annotation in Ucsc format (knownGene.txt knownIsoforms.txt).", ArgParseOption::STRING, "TXT TXT", true, 2));

    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    if (res != ArgumentParser::PARSE_OK)
        return res;

    FragmentStore<> store;
    CharString fileName;

    // IMPORT
    if (isSet(parser, "read-gff"))
    {
        getOptionValue(fileName, parser, "read-gff");
        GffFileIn file(toCString(fileName));
        readRecords(store, file);
    }
    if (isSet(parser, "read-ucsc"))
    {
        getOptionValue(fileName, parser, "read-ucsc", 0);
        UcscFileIn file1(toCString(fileName));
        getOptionValue(fileName, parser, "read-ucsc", 1);
        UcscFileIn file2(toCString(fileName));
        readRecords(store, file1);
        readRecords(store, file2);
    }

    // EXPORT
    if (isSet(parser, "write-gff"))
    {
        getOptionValue(fileName, parser, "write-gff");
        GffFileOut file(toCString(fileName));
        writeRecords(file, store);
    }
    if (isSet(parser, "write-ucsc"))
    {
        getOptionValue(fileName, parser, "write-ucsc", 0);
        UcscFileOut file1(toCString(fileName));
        getOptionValue(fileName, parser, "write-ucsc", 1);
        UcscFileOut file2(toCString(fileName));
        writeRecords(file1, store);
        writeRecords(file2, store);
    }

    return 0;
}
