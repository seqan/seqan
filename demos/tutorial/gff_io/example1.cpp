#include <seqan/gff_io.h>

int main()
{
    // Open input file.
    seqan::GffFileIn gffIn("example.gff");

    // Attach to standard output.
    seqan::GffFileOut gffOut(std::cout, seqan::Gff());

    // Copy the file record by record.
    seqan::GffRecord record;
    while (!atEnd(gffIn))
    {
        readRecord(record, gffIn);
        writeRecord(gffOut, record);
    }
    
    return 0;
}
