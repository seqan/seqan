#include <cstring>
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
    TReader reader = directionIterator(vin, Input());

    CharString buffer;
    reserve(buffer, 1000);

    while (!atEnd(reader))
    {
        clear(buffer);
        read(buffer, reader, capacity(buffer));
        write(vout, buffer);
    }
}

// The main function parses the command line, opens the files
// and then calls either copyStream. Read from stdin if the 1st
// argument is - or omitted and write to stdout if the 2nd argument
// is omitted.

int main(int argc, char const ** argv)
{
    VirtualStream<char, Input> vin;
    VirtualStream<char, Output> vout;

    if (argc < 1 || argc > 3)
    {
        CharString inExts = concat(vin.getFileExtensions(), "|", true);
        CharString outExts = concat(vout.getFileExtensions(), "|", true);
        std::cerr << "USAGE: " << argv[0] << " input[" << inExts << "] [output[" << outExts << "]]\n";
        return 1;
    }

    bool success;
    if (argc >= 2 && strcmp(argv[1], "-") != 0)
        success = open(vin, argv[1]);
    else
        success = open(vin, std::cin);

    if (!success)
    {
        std::cerr << "ERROR: Could not open input file " << argv[1] << "\n";
        return 1;
    }

    if (argc == 3)
        success = open(vout, argv[2]);
    else
        success = open(vout, std::cout, Nothing());  // disable compression on stdout

    if (!success)
    {
        std::cerr << "ERROR: Could not open output file " << ((argc == 3) ? argv[2] : "") << "\n";
        return 1;
    }

    copyStream(vout, vin);
    return 0;
}
