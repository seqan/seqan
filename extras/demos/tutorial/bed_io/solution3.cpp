#include <seqan/basic.h>
#include <seqan/bed_io.h>

#include <sstream>

int main()
{
    seqan::BedStream out("-", seqan::BedStream::WRITE);

    // Add sequence names.
    addSequenceName(out, "chr7");

    // Write out the records.
    seqan::BedRecord<seqan::Bed6> record;

    record.rID = 0;
    record.beginPos = 127471195;
    record.endPos = 127472363;
    record.name = "Pos1";
    record.score = "0";
    record.strand = '+';
    if (writeRecord(out, record) != 0)
        std::cerr << "ERROR: Problem writing output file.";

    record.rID = 0;
    record.beginPos = 127472362;
    record.endPos = 127473530;
    record.name = "Pos2";
    record.score = "0";
    record.strand = '+';
    if (writeRecord(out, record) != 0)
        std::cerr << "ERROR: Problem writing output file.";
    
    return 0;
}
