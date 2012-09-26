#include <iostream>
#include <fstream>

#include <seqan/bam_io.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>

int main(int argc, char const ** argv)
{
    if (argc != 6)
    {
        std::cerr << "USAGE: " << argv[0] << " IN.bam IN.bam.bai REF POS COUNT\n";
        return 1;
    }

    // Open BGZF Stream for reading.
    seqan::Stream<seqan::Bgzf> inStream;
    if (!open(inStream, argv[1], "r"))
    {
        std::cerr << "ERROR: Could not open " << argv[1] << " for reading.\n";
        return 1;
    }

    // Read BAI index.
    seqan::BamIndex<seqan::Bai> baiIndex;
    if (read(baiIndex, argv[2]) != 0)
    {
        std::cerr << "ERROR: Could not read BAI index file " << argv[2] << "\n";
        return 1;
    }

    // Setup name store, cache, and BAM I/O context.
    typedef seqan::StringSet<seqan::CharString> TNameStore;
    typedef seqan::NameStoreCache<TNameStore>   TNameStoreCache;
    typedef seqan::BamIOContext<TNameStore>     TBamIOContext;
    TNameStore      nameStore;
    TNameStoreCache nameStoreCache(nameStore);
    TBamIOContext   context(nameStore, nameStoreCache);

    // Read header.
    seqan::BamHeader header;
    if (readRecord(header, context, inStream, seqan::Bam()) != 0)
    {
        std::cerr << "ERROR: Could not read header from BAM file " << argv[1] << "\n";
        return 1;
    }

    // Translate from reference name to rId.
    int rId = 0;
    if (!getIdByName(nameStore, argv[3], rId, nameStoreCache))
    {
        std::cerr << "ERROR: Reference sequence named " << argv[3] << " not known.\n";
        return 1;
    }

    // Translate POS argument to number, 1-based to 0-based.
    int pos = 0;
    if (!seqan::lexicalCast2(pos, argv[4]) || pos <= 0)
    {
        std::cerr << "ERROR: Position " << argv[4] << " is invalid.\n";
        return 1;
    }
    pos -= 1;  // 1-based to 0-based.

    // Translate number of elements to print to number.
    int num = 0;
    if (!seqan::lexicalCast2(num, argv[5]))
    {
        std::cerr << "ERROR: Count " << argv[5] << " is invalid.\n";
        return 1;
    }

    // Jump the BGZF stream to this position.
    bool hasAlignments = false;
    if (!jumpToPos(inStream, hasAlignments, context, rId, pos, baiIndex))
    {
        std::cerr << "ERROR: Could not jump to " << argv[3] << ":" << argv[4] << "\n";
        return 1;
    }
    if (!hasAlignments)
        return 0;  // No alignments here.

    // Seek linearly to the selected position.
    seqan::BamAlignmentRecord record;
    int numPrinted = 0;
    while (!atEnd(inStream) && numPrinted < num)
    {
        if (readRecord(record, context, inStream, seqan::Bam()) != 0)
        {
            std::cerr << "ERROR: Could not read record from BAM file.\n";
            return 1;
        }

        // If we are on the next reference or at the end already then we stop.
        if (record.rId == -1 || record.rId > rId)
            break;
        // If we are left of the selected position then we skip this record.
        if (record.pos < pos)
            continue;

        // Otherwise, we print it to the user.
        numPrinted += 1;
        if (write2(std::cout, record, context, seqan::Sam()) != 0)
        {
            std::cerr << "ERROR: Could not write record to stdout.\n";
            return 1;
        }
    }

    return 0;
}
