#include <iostream>
#include <fstream>

#include <seqan/sequence.h>
#include <seqan/stream.h>

using namespace seqan;

int main(int argc, char const ** argv)
{
    if (argc != 3)
    {
        std::cerr << "USAGE: " << argv[0] << " IN.bam OUT.bin\n";
        return 1;
    }

    // Open BGZF file for reading.
    typedef VirtualStream<char, Input> TInStream;
    TInStream inStream;
    if (!open(inStream, argv[1]))
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
    DirectionIterator<TInStream, Input>::Type reader = directionIterator(inStream, Input());
    while (!atEnd(reader))
        read(outStream, reader, 1000);

    return 0;
}
