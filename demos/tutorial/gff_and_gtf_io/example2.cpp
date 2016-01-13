#include <seqan/basic.h>
#include <seqan/gff_io.h>

using namespace seqan;

int main()
{
    // Get path to example file.
    CharString file = getAbsolutePath("demos/tutorial/gff_and_gtf_io/example.gff");

    // Open input file.
    GffFileIn gffIn(toCString(file));

    // Open output stream. If target is a ostream we must specify the format.
    GffFileOut gffOut(std::cout, Gtf());

    // Read the file record by record.
    GffRecord record;
    while (!atEnd(gffIn))
    {
        readRecord(record, gffIn);
        writeRecord(gffOut, record);
    }

    return 0;
}
