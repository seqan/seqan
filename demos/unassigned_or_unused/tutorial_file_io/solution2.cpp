#include <iostream>
#include <fstream>

#include <seqan/stream.h>

using namespace seqan;

// Copy from stream in to the stream out.

template <typename TOutStream, typename TInStream>
void copyStream(TOutStream & vout, TInStream & vin)
{
    // Create iterators to read and write.
    typedef typename DirectionIterator<TInStream, Input>::Type TReader;
    typedef typename DirectionIterator<TOutStream, Output>::Type TWriter;

    TReader reader = directionIterator(vin, Input());
    TWriter writer = directionIterator(vout, Output());

    CharString buffer;
    reserve(buffer, 1000);

    while (!atEnd(reader))
    {
        clear(buffer);
        read(buffer, reader, capacity(buffer));
        write(writer, buffer);
    }
}

// The main function parses the command line, opens the files
// and then calls either copyStream.

int main(int argc, char const ** argv)
{
    VirtualStream<char, Input> vin;
    VirtualStream<char, Output> vout;

    if (argc != 3)
    {
        CharString inExts = concat(vin.getFileExtensions(), "|", true);
        CharString outExts = concat(vout.getFileExtensions(), "|", true);
        std::cerr << "USAGE: " << argv[0] << " input[" << inExts << "] output[" << outExts << "]\n";
        return 1;
    }

    if (!open(vin, argv[1]))
    {
        std::cerr << "ERROR: Could not open file " << argv[1] << "\n";
        return 1;
    }

    if (!open(vout, argv[2]))
    {
        std::cerr << "ERROR: Could not open file " << argv[2] << "\n";
        return 1;
    }

    copyStream(vout, vin);
    return 0;
}
