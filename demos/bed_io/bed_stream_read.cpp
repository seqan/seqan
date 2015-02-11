#include <seqan/basic.h>
#include <seqan/bed_io.h>

using namespace seqan;

// USAGE: bed_stream_read IN.bed
//
// Read BED file and print the positions (reference, position) to stdout.
//
// We only read BED 3 here but the example can be easily adjusted to more complex examples.

int main(int argc, char const * argv[])
{
    if (argc != 2)
    {
        std::cerr << "USAGE: " << argv[0] << " IN.bed\n";
        return 1;
    }

    BedFileIn bedIn;
    if (!open(bedIn, argv[1]))
    {
        std::cerr << "ERROR: Could not open " << argv[1] << " for reading!\n";
        return 1;
    }

    BedRecord<Bed3> record;
    while (!atEnd(bedIn))
    {
        readRecord(record, bedIn);
        std::cout << record.ref << "\t" << record.beginPos << "\n";
    }

    return 0;
}
