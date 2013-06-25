#include <seqan/basic.h>
#include <seqan/bed_io.h>

int main()
{
    // Open input stream
    seqan::BedStream bedIn("example.bed");
    if (!isGood(bedIn))
    {
        std::cerr << "ERROR: Could not open example.bed\n";
        return 1;
    }
    // Open output stream, filename "-" means stdout.
    seqan::BedStream bedOut("-", seqan::BedStream::WRITE);

    // Read the file record by record.
    seqan::BedRecord<seqan::Bed3> record;
    while (!atEnd(bedIn))
    {
        if (readRecord(record, bedIn) != 0)
        {
            std::cerr << "ERROR: Problem reading from example.bed\n";
            return 1;
        }

        // If record is on a sequence that is not known to bedOut yet then we
        // have to make it known there.
        if (record.rID >= (int)length(bedOut.sequenceNames))
            addSequenceName(bedOut, record.ref);

        if (writeRecord(bedOut, record) != 0)
        {
            std::cerr << "ERROR: Problem writing to stdout.\n";
            return 1;
        }
    }
    
    return 0;
}
