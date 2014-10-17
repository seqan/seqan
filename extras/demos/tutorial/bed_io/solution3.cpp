#include <seqan/basic.h>
#include <seqan/bed_io.h>

#include <sstream>

int main()
{
    seqan::BedFileOut out(std::cout, seqan::Bed());

    // Write out the records.
    seqan::BedRecord<seqan::Bed6> record;

    record.ref = "chr7";
    record.beginPos = 127471195;
    record.endPos = 127472363;
    record.name = "Pos1";
    record.score = "0";
    record.strand = '+';
    writeRecord(out, record);

    record.ref = "chr7";
    record.beginPos = 127472362;
    record.endPos = 127473530;
    record.name = "Pos2";
    record.score = "0";
    record.strand = '+';
    writeRecord(out, record);
    
    return 0;
}
