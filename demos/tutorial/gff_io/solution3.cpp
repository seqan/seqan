#include <seqan/basic.h>
#include <seqan/gff_io.h>

#include <sstream>

int main()
{
    seqan::GffFileOut out(std::cout, seqan::Gff());

    // Write out the records.
    seqan::GffRecord record;

    record.ref = "ctg123";
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
    writeRecord(out, record);

    clear(record.tagNames);
    clear(record.tagValues);

    record.ref = "ctg123";
    record.source = "";
    record.type = "TF_binding_site";
    record.beginPos = 999;
    record.endPos = 1012;
    record.strand = '+';
    record.score = seqan::GffRecord::INVALID_SCORE();
    appendValue(record.tagNames, "Parent");
    appendValue(record.tagValues, "gene0001");
    writeRecord(out, record);
    
    return 0;
}
