#include <seqan/basic.h>
#include <seqan/gff_io.h>

int main()
{
    // Open input stream.
    seqan::GffStream gffIn("example.gff");
    // Open output stream, filename "-" means stdout.
    seqan::GffStream gffOut("-", seqan::GffStream::WRITE);

    // Read the file record by record.
    seqan::GffRecord record;
    while (!atEnd(gffIn))
    {
        readRecord(record, gffIn);

        // If record is on a sequence that is not known to gffOut yet then we
        // have to make it known there.
        if (record.rID >= (int)length(gffOut.sequenceNames))
            addSequenceName(gffOut, record.ref);

        writeRecord(gffOut, record);
    }
    
    return 0;
}
