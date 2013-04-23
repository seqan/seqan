#include <iostream>
#include <sstream>
#include <seqan/bam_io.h>

int main()
{
    seqan::Dna5String ref = "CCCGATGAGCACACGATCACACGATGACA";
    std::stringstream ss;

    // --------------------------------------------------------
    // Build header.
    // --------------------------------------------------------
    seqan::BamHeader header;

    // Fill sequenceInfos.
    resize(header.sequenceInfos, 1);
    header.sequenceInfos[0].i1 = "REF";
    header.sequenceInfos[0].i2 = length(ref);

    // Fill header records.
    resize(header.records, 2);
    // @HD header.
    header.records[0].type = seqan::BAM_HEADER_FIRST;
    resize(header.records[0].tags, 1);
    // @HD header, tag/value: VN:1.4.
    header.records[0].tags[0].i1 = "VN";
    header.records[0].tags[0].i2 = "1.4";
    // @SQ header.
    header.records[1].type = seqan::BAM_HEADER_REFERENCE;
    resize(header.records[1].tags, 2);
    // @SQ header, tag/value: SN:REF
    header.records[1].tags[0].i1 = "SN";
    header.records[1].tags[0].i2 = "REF";
    // @SQ header, tag/value: LN:30
    header.records[1].tags[1].i1 = "LN";
    ss << length(ref);
    header.records[1].tags[1].i2 = ss.str();

    // --------------------------------------------------------
    // Write out records.
    // --------------------------------------------------------
    seqan::BamStream bamStream("-", seqan::BamStream::WRITE);
    bamStream.header = header;    

    seqan::BamAlignmentRecord record;

    for (unsigned i = 0; i + 12 - 1 < length(ref); ++i)
    {
        clear(record);
        // Set members that are the same for all records.
        record.rID = 0;
        record.flag = 0;
        resize(record.cigar, 1);
        record.cigar[0].operation = '=';
        record.cigar[0].count = 12;

        ss.str("");
        ss.clear();
        // The query name is REF_${START}_${END}.
        ss << "REF_" << i << "_" << (i + 12);
        record.qName = ss.str();
        // Set position.
        record.beginPos = i;
        // Set sequence.
        record.seq = infix(ref, i, i + 12);

        writeRecord(bamStream, record);
    }

    return 0;
}
