#include <seqan/basic.h>
#include <seqan/bam_io.h>

//![bamRecord]
//![BamTagsDict]
using namespace seqan;

//![BamTagsDict]
class myBamAlignmentRecord
{
public:
    CharString qName;               // QNAME
    uint16_t flag;                  // FLAG
    int32_t rID;                    // REF
    int32_t beginPos;               // POS
    uint8_t mapQ;                   // MAPQ mapping quality, 255 for */invalid
    uint16_t bin;                   // bin for indexing
    String<CigarElement<> > cigar;  // CIGAR string
    int32_t rNextId;                // RNEXT (0-based)
    int32_t pNext;                  // PNEXT (0-based)
    int32_t tLen;                   // TLEN
    CharString seq;                 // SEQ, as in SAM/BAM file.
    CharString qual;                // Quality string as in SAM (Phred).
    CharString tags;                // Tags, raw as in BAM.

    // Constants for marking pos, reference id and length members invalid (== 0/*).
    static int32_t const INVALID_POS = -1;
    static int32_t const INVALID_REFID = -1;
    static int32_t const INVALID_LEN = 0;
};
//![bamRecord]

int main()
{
//![BamTagsDict]
    BamAlignmentRecord record;
    BamTagsDict tagsDict(record.tags);
//![BamTagsDict]

//![addTag]
    setTagValue(tagsDict, "NM", 2);
    // => tags: "NM:i:2"
    setTagValue(tagsDict, "NH", 1);
    // => tags: "NM:i:2 NH:i:1"
    setTagValue(tagsDict, "NM", 3);
    // => tags: "NM:i:3 NH:i:1"
//![addTag]

//![getIndex]
    unsigned tagIdx = 0;
    if (!findTagKey(tagIdx, tagsDict, "NH"))
        std::cerr << "ERROR: Unknown key!\n";
//![getIndex]

//![extractValue]
    int tagValInt = 0;
    if (!extractTagValue(tagValInt, tagsDict, tagIdx))
        std::cerr << "ERROR: There was an error extracting NH from tags!\n";
//![extractValue]

//![cast]
    short tagValShort = 0;
    extractTagValue(tagValShort, tagsDict, tagIdx);
//![cast]

    return 0;
}
