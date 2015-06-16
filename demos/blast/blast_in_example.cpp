#include <iostream>
#include <seqan/basic.h>
#ifdef SEQAN_CXX11_COMPLETE
#include <seqan/blast.h>

using namespace seqan;

int main()
{
    typedef Align<String<AminoAcid>, ArrayGaps> TAlign;
    typedef BlastMatch<TAlign> TBlastMatch;
    typedef BlastRecord<TBlastMatch> TBlastRecord;
    typedef BlastIOContext<> TContext;

    std::string inPath = std::string(SEQAN_PATH_TO_ROOT()) + "/tests/blast/plus_comments_defaults.m9";

    BlastTabularFileIn<TContext> in(toCString(inPath));

    readHeader(in);

    TBlastRecord record;

    while (onRecord(in))
    {
        // read the record
        readRecord(record, in);

        // print some diagnostics
        std::cout << "Record of query sequence \"" << record.qId << "\"\n"
                  << "==========================================\n"
                  << "Number of HSPs: " << length(record.matches) << "\n";
        if  (!empty(length(record.matches)))
            std::cout << "E-Value of best HSP: " << front(record.matches).eValue << "\n";

        // if there is anything unexpected, tell the user about it
        if (!empty(context(in).conformancyErrors))
        {
            std::cout << "There were non critical errors when reading the record:\n";
            write(std::cout, context(in).conformancyErrors);
            std::cout << "\n";
        }

        std::cout << "\n";
    }

    readFooter(in);

    return 0;
}
#else
int main()
{
    std::cerr << "Demo not run, because you don't have full C++11 support.\n";
    return 0;
}
#endif
