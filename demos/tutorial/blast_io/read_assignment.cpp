#include <iostream>
#include <seqan/basic.h>
#ifndef STDLIB_VS
#include <seqan/blast.h>

using namespace seqan;

int main(int argc, char **argv)
{
    if (argc != 2)
    {
      std::cerr << "USAGE: FILE_IN\n";
      return 0;
    }

    typedef Gaps<String<AminoAcid>, ArrayGaps> TGaps;
    typedef BlastMatch<TGaps, TGaps> TBlastMatch;
    typedef BlastRecord<TBlastMatch> TBlastRecord;
    typedef BlastIOContext<> TContext;

    BlastTabularFileIn<TContext> file(argv[1]);

    readHeader(file);

    TBlastRecord record;

    while (onRecord(file))
    {
        // read the record
        readRecord(record, file);

        // print some diagnostics
        std::cout << "Record of query sequence \"" << record.qId << "\"\n"
                  << "==========================================\n"
                  << "Number of HSPs: " << length(record.matches) << "\n";
        if  (!empty(record.matches))
            std::cout << "BitScore of best HSP: " << front(record.matches).bitScore << "\n";

        // print column composition
        std::cout << "Columns: ";
        for (auto field : context(file).fields)
            std::cout << BlastMatchField<>::optionLabels[(int)field] << " ";
        std::cout << "\n\n";

        // if there is anything unexpected, tell the user about it
        if (!empty(context(file).conformancyErrors))
        {
            std::cout << "There were non critical errors when reading the record:\n";
            write(std::cout, context(file).conformancyErrors);
            std::cout << "\n\n";
        }

        if (!empty(context(file).otherLines))
        {
            std::cout << "There were unidentified lines in the comments:\n";
            write(std::cout, context(file).otherLines);
            std::cout << "\n\n";
        }

        std::cout << "\n\n";
    }

    readFooter(file);

    std::cout << "File Format: tabular"
              << (context(file).tabularSpec == BlastTabularSpec::COMMENTS ? " with coment lines" : "")
              << '\n'
              << "Generation: "
              << (context(file).legacyFormat ? " legacy" : " BLAST+")
              << '\n'
              << "Program and version: "
              << context(file).versionString
              << '\n'
              << "Database: "
              << context(file).dbName
              << "\n\n";

    return 0;
}
#else
int main()
{
    std::cerr << "USAGE: FILE_IN\n";
    return 0;
}
#endif
