#include <seqan/basic.h>
#include <seqan/gff_io.h>

#include <sstream>

int main()
{
    seqan::GffStream out("-", seqan::GffStream::WRITE);

    // Add sequence names.
    addSequenceName(out, "ctg123");

    // Write out the records.
    seqan::GffRecord record;

    record.rID = 0;
    record.source = "";
    record.type = "gene";
    record.beginPos = 999;
    record.endPos = 9000;
    record.strand = '+';
    record.score = seqan::GffRecord::INVALID_SCORE();
    appendValue(record.tagNames, "ID");
    appendValue(record.tagValues, "gene0001");
    appendValue(record.tagNames, "Name");
    appendValue(record.tagValues, "EDEN");
    if (writeRecord(out, record) != 0)
        std::cerr << "ERROR: Problem writing output file.";

    clear(record.tagNames);
    clear(record.tagValues);

    record.rID = 0;
    record.source = "";
    record.type = "TF_binding_site";
    record.beginPos = 999;
    record.endPos = 1012;
    record.strand = '+';
    record.score = seqan::GffRecord::INVALID_SCORE();
    appendValue(record.tagNames, "Parent");
    appendValue(record.tagValues, "gene0001");
    if (writeRecord(out, record) != 0)
        std::cerr << "ERROR: Problem writing output file.";
    
    return 0;
}
