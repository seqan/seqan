#include <iostream>
#include <seqan/sequence.h>
#include <seqan/bam_io.h>

int main()
{
    seqan::Dna5String ref = "CCCGATGAGCACACGATCACACGATGACA";

    // --------------------------------------------------------
    // Build header.
    // --------------------------------------------------------
    seqan::BamFileOut bamFileOut(std::cout, seqan::Sam());

    // Fill sequenceInfos.
    assignValueById(sequenceLengths(context(bamFileOut)),
                    nameToId(nameStoreCache(context(bamFileOut)), "REF"),
                    length(ref));

    // Fill header records.
    seqan::BamHeader header;
    resize(header, 1);
    // @HD header.
    header[0].type = seqan::BAM_HEADER_FIRST;
    resize(header[0].tags, 1);
    // @HD header, tag/value: VN:1.4.
    header[0].tags[0].i1 = "VN";
    header[0].tags[0].i2 = "1.4";

    writeRecord(bamFileOut, header);

    // --------------------------------------------------------
    // Write out records.
    // --------------------------------------------------------

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

        // The query name is REF_${START}_${END}.
        record.qName = "REF_";
        appendNumber(record.qName, i);
        appendValue(record.qName, '_');
        appendNumber(record.qName, i + 12);
        // Set position.
        record.beginPos = i;
        // Set sequence.
        record.seq = infix(ref, i, i + 12);

        // Write "NH" tag.
        seqan::BamTagsDict tagsDict(record.tags);
        setTagValue(tagsDict, "NH", 1);

        writeRecord(bamFileOut, record);
    }

    return 0;
}
