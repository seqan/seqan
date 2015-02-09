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
        std::cerr << "USAGE: " << argv[0] << " IN.gff\n";
        return 1;
    }

    try
    {
        GffFileIn gffIn(argv[1]);
        GffRecord record;
        while (!atEnd(gffIn))
        {
            readRecord(record, gffIn);

            // Note that we print the position 1-based since we use text output
            // whereas it is 0-based in the GffRecord.
            std::cout << record.ref << "\t" << (record.beginPos + 1) << "\n";
        }
    }
    catch (Exception & e)
    {
        std::cerr << e.what() << std::endl;
        return 1;
    }

    return 0;
}
