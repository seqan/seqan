#include <iostream>
#include <fstream>

#include <seqan/bam_io.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>

template <typename TOutStream, typename TInStreamOrReader, typename TTag>
int copyHeader(TOutStream & outStream,
               TInStreamOrReader & inStreamOrReader,
               seqan::BamIOContext<seqan::StringSet<seqan::CharString> > & context,
               TTag const & tag)
{

    // Read header.
    seqan::BamHeader header;
    if (readRecord(header, context, inStreamOrReader, tag) != 0)
    {
        std::cerr << "ERROR: Could not read header\n";
        return 1;
    }

    // Write out header again.
    if (write2(outStream, header, context, tag) != 0)
    {
        std::cerr << "ERROR: Could not write header.\n";
        return 1;
    }

    return 0;
}

int main(int argc, char const ** argv)
{
    if (argc != 3)
    {
        std::cerr << "USAGE: " << argv[0] << " IN.sam OUT.sam\n";
        return 1;
    }

    // Streams for SAM.
    std::fstream inStreamSam, outStreamSam;
    // Streams for BAM.
    seqan::Stream<seqan::Bgzf> inStreamBam, outStreamBam;

    if (seqan::endsWith(seqan::CharString(argv[1]), ".sam"))
    {
        inStreamSam.open(argv[1], std::ios::binary | std::ios::in);
        if (!inStreamSam.good())
        {
            std::cerr << "ERROR: Could not open " << argv[1] << " for reading.\n";
            return 1;
        }
        outStreamSam.open(argv[2], std::ios::binary | std::ios::out);
        if (!outStreamSam.good())
        {
            std::cerr << "ERROR: Could not open " << argv[2] << " for writing.\n";
            return 1;
        }
    }
    else
    {
        // Open BGZF Stream for reading.
        if (!open(inStreamBam, argv[1], "r"))
        {
            std::cerr << "ERROR: Could not open " << argv[1] << " for reading.\n";
            return 1;
        }

        // Open BGZF Stream for writing.
        if (!open(outStreamBam, argv[2], "w"))
        {
            std::cerr << "ERROR: Could not open " << argv[2] << " for writing.\n";
            return 1;
        }
    }

    // Setup name store, cache, and BAM I/O context.
    typedef seqan::StringSet<seqan::CharString> TNameStore;
    typedef seqan::NameStoreCache<TNameStore>   TNameStoreCache;
    typedef seqan::BamIOContext<TNameStore>     TBamIOContext;
    TNameStore      nameStore;
    TNameStoreCache nameStoreCache(nameStore);
    TBamIOContext   context(nameStore, nameStoreCache);

    if (endsWith(seqan::CharString(argv[1]), ".sam"))
    {
        // Stream must be open before constructing reader, thus we define it here.
        seqan::RecordReader<std::fstream, seqan::SinglePass<> > reader(inStreamSam);
        return copyHeader(outStreamSam, reader, context, seqan::Sam());
    }
    else
    {
        return copyHeader(outStreamBam, inStreamBam, context, seqan::Bam());
    }
}
