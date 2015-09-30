#include <fstream>
#include <seqan/sequence.h>
#include <seqan/stream.h>

using namespace seqan;

int main()
{
    std::ifstream in("in.txt", std::ios::binary | std::ios::in);
    std::ofstream out("out.txt", std::ios::binary | std::ios::out);

    // Create iterators to read and write.
    typedef DirectionIterator<std::ifstream, Input>::Type TReader;
    typedef DirectionIterator<std::ofstream, Output>::Type TWriter;

    TReader reader = directionIterator(in, Input());
    TWriter writer = directionIterator(out, Output());

    CharString buffer;
    reserve(buffer, 1000);

    while (!atEnd(reader))
    {
        clear(buffer);
        read(buffer, reader, capacity(buffer));
        write(writer, buffer);
    }

    return 0;
}
