#include <iostream>
#include <fstream>

#include <seqan/sequence.h>
#include <seqan/stream.h>

int main(int argc, char const ** argv)
{
    if (argc != 3)
    {
        std::cerr << "USAGE: " << argv[0] << " IN.bam OUT.bin\n";
        return 1;
    }

    // Open BGZF file for reading.
    seqan::Stream<seqan::Bgzf> inStream;
    if (!open(inStream, argv[1], "r"))
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

    // Copy over data.
    seqan::CharString buffer;
    resize(buffer, 1000);
    while (!seqan::atEnd(inStream) && seqan::streamError(inStream) == 0)
    {
        int num = seqan::streamReadBlock(&buffer[0], inStream, length(buffer));
        seqan::streamWriteBlock(outStream, &buffer[0], num);
    }
    
    return 0;
}
