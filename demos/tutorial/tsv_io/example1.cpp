#include <seqan/tsv_io.h>

using namespace seqan;

int main()
{
    CharString tsvFileName = getAbsolutePath("demos/tutorial/tsv_io/example.tsv");

    // Open input file.
    TsvFileIn tsvFileIn(toCString(tsvFileName));

    // Open output file.
    TsvFileOut tsvFileOut(std::cout, Tsv());
    
    // Copy header.
    TsvHeader header;
    readHeader(header, tsvFileIn);
    writeHeader(tsvFileOut, header);

    // Copy records.
    TsvRecord record;
    while (!atEnd(tsvFileIn))
    {
        readRecord(record, tsvFileIn);
        writeRecord(tsvFileOut, record);
    }

    return 0;
}
