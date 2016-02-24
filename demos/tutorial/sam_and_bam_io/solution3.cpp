#include <seqan/bam_io.h>

using namespace seqan;

int main()
{
    CharString bamFileName = getAbsolutePath("demos/tutorial/sam_and_bam_io/example.sam");

    // Open input file.
    BamFileIn bamFileIn;
    if (!open(bamFileIn, toCString(bamFileName)))
    {
        std::cerr << "ERROR: Could not open " << bamFileName << std::endl;
        return 1;
    }

    unsigned numXXtags = 0;

    try
    {
        // Read header.
        BamHeader header;
        readHeader(header, bamFileIn);

        // Rear records.
        BamAlignmentRecord record;
        while (!atEnd(bamFileIn))
        {
            readRecord(record, bamFileIn);
            BamTagsDict tagsDict(record.tags);

            unsigned tagIdx = 0;
            if (findTagKey(tagIdx, tagsDict, "XX"))
                numXXtags += 1;
        }
    }
    catch (Exception const & e)
    {
        std::cout << "ERROR: " << e.what() << std::endl;
        return 1;
    }

    std::cout << "Number of records with the XX tag: " << numXXtags << "\n";

    return 0;
}
