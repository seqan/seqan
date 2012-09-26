#include <iostream>
#include <fstream>

#include <seqan/bam_io.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>

int main(int argc, char const ** argv)
{
    if (argc != 2)
    {
        std::cerr << "USAGE: " << argv[0] << " IN.bam\n";
        return 1;
    }

    // Open BGZF Stream for reading.
    seqan::Stream<seqan::Bgzf> inStream;
    if (!open(inStream, argv[1], "r"))
    {
        std::cerr << "ERROR: Could not open " << argv[1] << " for reading.\n";
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

    // Write out header again.
    if (write2(std::cout, header, context, seqan::Sam()) != 0)
    {
        std::cerr << "ERROR: Could not write header to stdout.\n";
        return 1;
    }

    // Copy over the alignment records.
    seqan::BamAlignmentRecord record;
    while (!atEnd(inStream))
    {
        if (readRecord(record, context, inStream, seqan::Bam()) != 0)
        {
            std::cerr << "ERROR: Could not read record from BAM File " << argv[1] << "\n";
            return 1;
        }

        if (write2(std::cout, record, context, seqan::Sam()) != 0)
        {
            std::cerr << "ERROR: Could not write record to stdout.\n";
            return 1;
        }
    }

    return 0;
}
