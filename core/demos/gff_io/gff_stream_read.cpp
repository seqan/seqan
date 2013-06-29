#include <seqan/basic.h>
#include <seqan/gff_io.h>

using namespace seqan;

// USAGE: gff_stream_read IN.gff
//
// Read GFF file and print the positions (reference, position) to stdout.

int main(int argc, char ** argv)
{
    if (argc != 2)
    {
        std::cerr << "USAGE: " << argv[1] << " IN.gff\n";
        return 1;
    }

    GffStream gffIn(argv[1]);
    if (!isGood(gffIn))
    {
        std::cerr << "ERROR: Could not open " << argv[1] << " for reading!\n";
        return 1;
    }

    GffRecord record;
    while (!atEnd(gffIn))
    {
        if (readRecord(record, gffIn) != 0)
        {
            std::cerr << "ERROR: Problem reading from " << argv[1] << "\n";
            return 1;
        }

        // Note that we print the position 1-based since we use text output
        // whereas it is 0-based in the GffRecord.
        std::cout << gffIn.sequenceNames[record.rID]
                  << "\t" << (record.beginPos + 1) << "\n";
    }

    return 0;
}
