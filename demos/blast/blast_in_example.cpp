#include <iostream>
#include <seqan/blast.h>

using namespace seqan;

int main()
{
    typedef Align<String<AminoAcid>, ArrayGaps> TAlign;
    typedef BlastMatch<CharString, CharString, uint32_t, TAlign> TBlastMatch;
    typedef BlastRecord<CharString, CharString, uint32_t, TAlign> TBlastRecord;
    typedef BlastIOContext<Blosum62> TContext;

    std::string inPath = std::string(SEQAN_PATH_TO_ROOT()) + "/tests/blast/plus_header_defaults.blast";

    BlastTabularIn<TContext> in(toCString(inPath));

    TBlastRecord record;

    while (!atEnd(in))
    {
        // read the record
        readRecord(record, in);

        // print some diagnostics
        std::cout << "Record of query sequence \"" << record.qId << "\"\n"
                  << "==========================================\n"
                  << "Number of HSPs: " << length(record.matches) << "\n";
        if  (length(record.matches) > 0)
            std::cout << "E-Value of best HSP: " << front(record.matches).eValue << "\n";

        // if there is anything unexpected, tell the user about it
        if (length(context(in).conformancyErrors) > 0)
        {
            std::cout << "There were non critical errors when reading the record:\n";
            write(std::cout, context(in).conformancyErrors);
            std::cout << "\n";
        }

        std::cout << "\n";
    }

    return 0;
}