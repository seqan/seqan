#include <seqan/basic.h>
#include <seqan/gff_io.h>

int main()
{
    // Open input stream
    seqan::GffStream gffIn("example.gff");
    if (!isGood(gffIn))
    {
        std::cerr << "ERROR: Could not open example.gff\n";
        return 1;
    }
    // Open output stream, filename "-" means stdout.
    seqan::GffStream gffOut("-", seqan::GffStream::WRITE);

    // Read the file record by record.
    seqan::GffRecord record;
    while (!atEnd(gffIn))
    {
        if (readRecord(record, gffIn) != 0)
        {
            std::cerr << "ERROR: Problem reading from example.gff\n";
            return 1;
        }

        // If record is on a sequence that is not known to gffOut yet then we
        // have to make it known there.
        if (record.rID >= (int)length(gffOut.sequenceNames))
            addSequenceName(gffOut, record.ref);

        if (writeRecord(gffOut, record) != 0)
        {
            std::cerr << "ERROR: Problem writing to stdout.\n";
            return 1;
        }
    }
    
    return 0;
}
