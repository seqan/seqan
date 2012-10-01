// Reading gzip-compressed FASTQ into a Dna5Q String and print it.

#include <iostream>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/seq_io.h>

int main(int argc, char const ** argv)
{
#ifdef SEQAN_HAS_ZLIB
    using namespace seqan;

    // Check arguments.
    
    if (argc != 2)
    {
        std::cerr << "USAGE: stream_read_fastq_gz IN.fq\n";
        return 1;
    }

    // Open GZ file.

    Stream<GZFile> gzStream;
    if (!open(gzStream, argv[1], "r"))
    {
        std::cerr << "ERROR: Could not open file " << argv[1] << '\n';
        
        return 1;
    }

    // Read FASTQ file from gzip-compressed file.

    RecordReader<Stream<GZFile>, SinglePass<> > reader(gzStream);
    CharString id;
    String<Dna5Q> seq;

    std::cerr << "ID\tSEQ\tQUAL\n";
    
    while (!atEnd(reader))
    {
        if (readRecord(id, seq, reader, Fastq()) != 0)
        {
            std::cerr << "Problem with your FASTQ file." << std::endl;
            return 1;
        }

        CharString qual;
        assignQualities(qual, seq);
        std::cerr << id << '\t' << seq << '\t' << qual << '\n';
    }

#endif  // #ifdef SEQAN_HAS_ZLIB
    return 0;
}
