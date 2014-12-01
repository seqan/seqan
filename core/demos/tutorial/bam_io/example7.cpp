#include <iostream>
#include <fstream>

#include <seqan/sequence.h>
#include <seqan/bam_io.h>

int main(int argc, char const * argv[])
{
    if (argc != 7)
    {
        std::cerr << "USAGE: " << argv[0] << " IN.bam IN.bam.bai REF BEGIN END COUNT\n";
        return 1;
    }

    // Open BamFileIn for reading.
    seqan::BamFileIn inFile;
    if (!open(inFile, argv[1]))
    {
        std::cerr << "ERROR: Could not open " << argv[1] << " for reading.\n";
        return 1;
    }

    // Read BAI index.
    seqan::BamIndex<seqan::Bai> baiIndex;
    if (!open(baiIndex, argv[2]))
    {
        std::cerr << "ERROR: Could not read BAI index file " << argv[2] << "\n";
        return 1;
    }

    // Read header.
    seqan::BamHeader header;
    readRecord(header, inFile);

    // Translate from reference name to rID.
    int rID = 0;
    if (!getIdByName(rID, nameStoreCache(context(inFile)), argv[3]))
    {
        std::cerr << "ERROR: Reference sequence named " << argv[3] << " not known.\n";
        return 1;
    }

    // Translate BEGIN and END arguments to number, 1-based to 0-based.
    int beginPos = 0, endPos = 0;
    if (!seqan::lexicalCast(beginPos, argv[4]) || beginPos <= 0)
    {
        std::cerr << "ERROR: Begin position " << argv[4] << " is invalid.\n";
        return 1;
    }
    beginPos -= 1;  // 1-based to 0-based.
    if (!seqan::lexicalCast(endPos, argv[5]) || endPos <= 0)
    {
        std::cerr << "ERROR: End position " << argv[5] << " is invalid.\n";
        return 1;
    }
    endPos -= 1;  // 1-based to 0-based.

    // Translate number of elements to print to number.
    int num = 0;
    if (!seqan::lexicalCast(num, argv[6]))
    {
        std::cerr << "ERROR: Count " << argv[6] << " is invalid.\n";
        return 1;
    }

    // Jump the BGZF stream to this position.
    bool hasAlignments = false;
    if (!jumpToRegion(inFile, hasAlignments, rID, beginPos, endPos, baiIndex))
    {
        std::cerr << "ERROR: Could not jump to " << argv[3] << ":" << argv[4] << "\n";
        return 1;
    }
    if (!hasAlignments)
        return 0;  // No alignments here.

    // Seek linearly to the selected position.
    seqan::BamAlignmentRecord record;
    int numPrinted = 0;
    seqan::BamFileOut out(inFile, std::cout, seqan::Sam());

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
