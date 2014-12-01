#include <iostream>
#include <fstream>

#include <seqan/stream.h>

// Copy from stream in to the stream out.

template <typename TOutStream, typename TInStream>
void copyStream(TOutStream & vout, TInStream & vin)
{
    // Create iterators to read and write.
    typedef typename seqan::DirectionIterator<TInStream, seqan::Input>::Type TReader;
    typedef typename seqan::DirectionIterator<TOutStream, seqan::Output>::Type TWriter;

    TReader reader = directionIterator(vin, seqan::Input());
    TWriter writer = directionIterator(vout, seqan::Output());

    seqan::CharString buffer;
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
    seqan::VirtualStream<char, seqan::Input> vin;
    seqan::VirtualStream<char, seqan::Output> vout;

    if (argc != 3)
    {
        seqan::CharString inExts = seqan::concat(vin.getFileExtensions(), "|", true);
        seqan::CharString outExts = seqan::concat(vout.getFileExtensions(), "|", true);
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
