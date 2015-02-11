#include <iostream>
#include <fstream>

#include <seqan/sequence.h>
#include <seqan/bam_io.h>

using namespace seqan;

int main(int argc, char const ** argv)
{
    if (argc != 2)
    {
        std::cerr << "USAGE: " << argv[0] << " IN.bam\n";
        return 1;
    }

    // Open BGZF Stream for reading.
    typedef VirtualStream<char, Input> TInStream;
    TInStream inStream;
    if (!open(inStream, argv[1]))
    {
        std::cerr << "ERROR: Could not open " << argv[1] << " for reading.\n";
        return 1;
    }

    // Setup name store, cache, and BAM I/O context.
    typedef StringSet<CharString> TNameStore;
    typedef NameStoreCache<TNameStore>   TNameStoreCache;
    typedef BamIOContext<TNameStore>     TBamIOContext;
    TNameStore      contigNames;
    TNameStoreCache contigNamesCache(contigNames);
    TBamIOContext   context(contigNames, contigNamesCache);

    // Read header.
    BamHeader header;
    DirectionIterator<TInStream, Input>::Type reader = directionIterator(inStream, Input());
    readHeader(header, context, reader, Bam());

    // Write out header again.
    write(std::cout, header, context, Sam());

    return 0;
}
