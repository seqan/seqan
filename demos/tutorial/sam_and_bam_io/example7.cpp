#include <seqan/sequence.h>
#include <seqan/bam_io.h>

using namespace seqan;

int main(int argc, char const * argv[])
{
    CharString bamFileName = getAbsolutePath("demos/tutorial/sam_and_bam_io/example.bam");
    CharString baiFileName = getAbsolutePath("demos/tutorial/sam_and_bam_io/example.bam.bai");
    CharString rName = "ref";

    // Open BamFileIn for reading.
    BamFileIn inFile;
    if (!open(inFile, toCString(bamFileName)))
    {
        std::cerr << "ERROR: Could not open " << bamFileName << " for reading.\n";
        return 1;
    }

    // Read BAI index.
    BamIndex<Bai> baiIndex;
    if (!open(baiIndex, toCString(baiFileName)))
    {
        std::cerr << "ERROR: Could not read BAI index file " << baiFileName << "\n";
        return 1;
    }

    // Read header.
    BamHeader header;
    readHeader(header, inFile);

    // Translate from reference name to rID.
    int rID = 0;
    if (!getIdByName(rID, contigNamesCache(context(inFile)), rName))
    {
        std::cerr << "ERROR: Reference sequence named " << rName << " not known.\n";
        return 1;
    }

    // Translate BEGIN and END arguments to number, 1-based to 0-based.
    int beginPos = 9, endPos = 30;

    // 1-based to 0-based.
    beginPos -= 1;
    endPos -= 1;

    // Translate number of elements to print to number.
    int num = 3;

    // Jump the BGZF stream to this position.
    bool hasAlignments = false;
    if (!jumpToRegion(inFile, hasAlignments, rID, beginPos, endPos, baiIndex))
    {
        std::cerr << "ERROR: Could not jump to " << beginPos << ":" << endPos << "\n";
        return 1;
    }
    if (!hasAlignments)
        return 0;  // No alignments here.

    // Seek linearly to the selected position.
    BamAlignmentRecord record;
    int numPrinted = 0;
    BamFileOut out(inFile, std::cout, Sam());

    while (!atEnd(inFile) && numPrinted < num)
    {
        readRecord(record, inFile);

        // If we are on the next reference or at the end already then we stop.
        if (record.rID == -1 || record.rID > rID || record.beginPos >= endPos)
            break;
        // If we are left of the selected position then we skip this record.
        if (record.beginPos < beginPos)
            continue;

        // Otherwise, we print it to the user.
        numPrinted++;
        writeRecord(out, record);
    }

    return 0;
}
