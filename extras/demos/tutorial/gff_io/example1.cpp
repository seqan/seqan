#include <seqan/basic.h>
#include <seqan/gff_io.h>

int main()
{
    // Open input stream.
    seqan::GffFileIn gffIn("example.gff");
    // Open output stream. If target is a ostream we must specify the format.
    seqan::GffFileOut gffOut(std::cout, seqan::Gff());

    // Read the file record by record.
    seqan::GffRecord record;
    while (!atEnd(gffIn))
    {
        readRecord(record, gffIn);
        writeRecord(gffOut, record);
    }
    
    return 0;
}
