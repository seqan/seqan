#include <iostream>
#include <fstream>

#include <seqan/sequence.h>
#include <seqan/bam_io.h>

int main(int argc, char const * argv[])
{
    if (argc != 3)
    {
        std::cerr << "USAGE: " << argv[0] << " IN.[sam|bam] OUT.[sam|bam]\n";
        return 1;
    }

    // Open BamFileIn for reading.
    seqan::BamFileIn inFile;
    if (!open(inFile, argv[1]))
    {
        std::cerr << "ERROR: Could not open " << argv[1] << " for reading.\n";
        return 1;
    }

    // Open BamFileOut for writing. Give inFile to share its BamIoContext
    seqan::BamFileOut outFile(inFile);
    if (!open(outFile, argv[2]))
    {
        std::cerr << "ERROR: Could not open " << argv[2] << " for writing.\n";
        return 1;
    }

    // Read header.
    seqan::BamHeader header;
    readRecord(header, inFile);
    writeRecord(outFile, header);

    // Copy over the alignment records.
    seqan::BamAlignmentRecord record;
    while (!atEnd(inFile))
    {
        readRecord(record, inFile);
        writeRecord(outFile, record);
    }

    return 0;
}
