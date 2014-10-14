#include <fstream>
#include <seqan/sequence.h>
#include <seqan/stream.h>

int main()
{
    std::ifstream in("in.txt", std::ios::binary | std::ios::in);
    std::ofstream out("out.txt", std::ios::binary | std::ios::out);

    // Create iterators to read and write.
    typedef seqan::DirectionIterator<std::ifstream, seqan::Input>::Type TReader;
    typedef seqan::DirectionIterator<std::ofstream, seqan::Output>::Type TWriter;

    TReader reader = directionIterator(in, seqan::Input());
    TWriter writer = directionIterator(out, seqan::Output());

    seqan::CharString buffer;

    while (!atEnd(reader))
    {
        read(buffer, reader, 1000);
        write(writer, buffer);
    }

    return 0;
}
