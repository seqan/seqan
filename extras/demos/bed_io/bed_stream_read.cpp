#include <seqan/basic.h>
#include <seqan/bed_io.h>

using namespace seqan;

// USAGE: bed_stream_read IN.bed
//
// Read BED file and print the positions (reference, position) to stdout.
//
// We only read BED 3 here but the example can be easily adjusted to more complex examples.

int main(int argc, char ** argv)
{
    if (argc != 2)
    {
        std::cerr << "USAGE: " << argv[1] << " IN.bed\n";
        return 1;
    }

    BedStream bedIn(argv[1]);
    if (!isGood(bedIn))
    {
        std::cerr << "ERROR: Could not open " << argv[1] << " for reading!\n";
        return 1;
    }

    BedRecord<Bed3> record;
    while (!atEnd(bedIn))
    {
        if (readRecord(record, bedIn) != 0)
        {
            std::cerr << "ERROR: Problem reading from " << argv[1] << "\n";
            return 1;
        }

        std::cout << record.ref << "\t" << record.beginPos << "\n";
    }
    
    return 0;
}
