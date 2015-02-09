#include <seqan/gff_io.h>

using namespace seqan;

int main()
{
    GffFileOut out(std::cout, Gff());

    GffRecord record;

    // Fill and write out the first record.
    record.ref = "ctg123";
    record.source = "";
    record.type = "gene";
    record.beginPos = 999;
    record.endPos = 9000;
    record.strand = '+';
    record.score = GffRecord::INVALID_SCORE();
    appendValue(record.tagNames, "ID");
    appendValue(record.tagValues, "gene0001");
    appendValue(record.tagNames, "Name");
    appendValue(record.tagValues, "EDEN");
    writeRecord(out, record);

    // Clear the record.
    clear(record.tagNames);
    clear(record.tagValues);

    // Fill and write out the second record.
    record.ref = "ctg123";
    record.source = "";
    record.type = "TF_binding_site";
    record.beginPos = 999;
    record.endPos = 1012;
    record.strand = '+';
    record.score = GffRecord::INVALID_SCORE();
    appendValue(record.tagNames, "Parent");
    appendValue(record.tagValues, "gene0001");
    writeRecord(out, record);

    return 0;
}
