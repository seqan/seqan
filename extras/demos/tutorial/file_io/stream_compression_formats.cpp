// FRAGMENT(header)
#include <fstream>
#include <seqan/sequence.h>
#include <seqan/stream.h>

// FRAGMENT(main)
int main(int argc, char const ** argv)
{
    typedef seqan::VirtualStream<char, seqan::Input> TVStream;

    if (argc != 2)
    {
        seqan::CharString exts = seqan::concat(TVStream::getFileExtensions(), "|", true);
        std::cerr << "USAGE: " << argv[0] << " input[" << exts << "]" << std::endl;
        return 1;
    }

    TVStream vin;

    if (!open(vin, argv[1]))
    {
        std::cerr << "ERROR: Could not open " << argv[1] << std::endl;
        return 1;
    }

    // Create iterators to read and write.
    typedef seqan::DirectionIterator<TVStream, seqan::Input>::Type TReader;
    TReader reader = directionIterator(vin, seqan::Input());

    seqan::CharString buffer;
    reserve(buffer, 1000);

    while (!atEnd(reader))
    {
        clear(buffer);
        read(buffer, reader, capacity(buffer));
        write(std::cout, buffer);
    }

    return 0;
}
