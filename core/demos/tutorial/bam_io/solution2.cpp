#include <iostream>
#include <fstream>

#include <seqan/sequence.h>
#include <seqan/bam_io.h>

int main(int argc, char const ** argv)
{
    if (argc != 2)
    {
        std::cerr << "USAGE: " << argv[0] << " IN.bam\n";
        return 1;
    }

    // Open BGZF Stream for reading.
    typedef seqan::VirtualStream<char, seqan::Input> TInStream;
    TInStream inStream;
    if (!open(inStream, argv[1]))
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
    seqan::DirectionIterator<TInStream, seqan::Input>::Type reader = directionIterator(inStream, seqan::Input());
    readRecord(header, context, reader, seqan::Bam());

    // Write out header again.
    write(std::cout, header, context, seqan::Sam());

    return 0;
}
