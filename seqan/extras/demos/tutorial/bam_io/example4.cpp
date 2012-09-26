#include <iostream>
#include <fstream>

#include <seqan/bam_io.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>

int main(int argc, char const ** argv)
{
    if (argc != 3)
    {
        std::cerr << "USAGE: " << argv[0] << " IN.sam OUT.sam\n";
        return 1;
    }

    // Open std::fstream for reading.
    std::fstream inStream(argv[1], std::ios::binary | std::ios::in);
    if (!inStream.good())
    {
        std::cerr << "ERROR: Could not open " << argv[1] << " for reading.\n";
        return 1;
    }

    // Open std::fstream for writing.
    std::fstream outStream(argv[2], std::ios::binary | std::ios::out);
    if (!outStream.good())
    {
        std::cerr << "ERROR: Could not open " << argv[2] << " for writing.\n";
        return 1;
    }

    // Setup RecordReader.
    seqan::RecordReader<std::fstream, seqan::SinglePass<> > reader(inStream);

    // Setup name store, cache, and BAM I/O context.
    typedef seqan::StringSet<seqan::CharString> TNameStore;
    typedef seqan::NameStoreCache<TNameStore>   TNameStoreCache;
    typedef seqan::BamIOContext<TNameStore>     TBamIOContext;
    TNameStore      nameStore;
    TNameStoreCache nameStoreCache(nameStore);
    TBamIOContext   context(nameStore, nameStoreCache);

    // Read header.
    seqan::BamHeader header;
    if (readRecord(header, context, reader, seqan::Sam()) != 0)
    {
        std::cerr << "ERROR: Could not read header from SAM file " << argv[1] << "\n";
        return 1;
    }

    // Write out header again.
    if (write2(outStream, header, context, seqan::Sam()) != 0)
    {
        std::cerr << "ERROR: Could not write header to SAM file " << argv[2] << "\n";
        return 1;
    }

    return 0;
}
