#include <seqan/basic.h>
#include <seqan/bed_io.h>

int main()
{
    // Open input stream.
    seqan::BedStream bedIn("example.bed");
    // Open output stream, filename "-" means stdout.
    seqan::BedStream bedOut("-", seqan::BedStream::WRITE);

    // Read the file record by record.
    seqan::BedRecord<seqan::Bed3> record;
    while (!atEnd(bedIn))
    {
        readRecord(record, bedIn);

        // If record is on a sequence that is not known to bedOut yet then we
        // have to make it known there.
        if (record.rID >= (int)length(bedOut.sequenceNames))
            addSequenceName(bedOut, record.ref);

        writeRecord(bedOut, record);
    }
    
    return 0;
}
