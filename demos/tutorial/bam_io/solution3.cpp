#include <seqan/bam_io.h>

using namespace seqan;

int main()
{
    // Open input file.
    seqan::BamFileIn bamFileIn;
    if (!open(bamFileIn, "example.sam"))
    {
        std::cerr << "ERROR: Could not open example.sam!" << std::endl;
        return 1;
    }

    unsigned numXXtags = 0;

    try
    {
        // Read header.
        seqan::BamHeader header;
        readRecord(header, bamFileIn);

        // Rear records.
        seqan::BamAlignmentRecord record;
        while (!atEnd(bamFileIn))
        {
            readRecord(record, bamFileIn);
            seqan::BamTagsDict tagsDict(record.tags);

            unsigned tagIdx = 0;
            if (findTagKey(tagIdx, tagsDict, "XX"))
                numXXtags += 1;
        }
    }
    catch (seqan::Exception const & e)
    {
        std::cout << "ERROR: " << e.what() << std::endl;
        return 1;
    }

    std::cout << "Number of records with the XX tag: " << numXXtags << "\n";

    return 0;
}
