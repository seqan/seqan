#include <iostream>
#include <fstream>

#include <seqan/sequence.h>
#include <seqan/stream.h>

int main(int argc, char const ** argv)
{
    if (argc != 2)
        return 1;
    std::fstream stream(argv[1], std::ios::binary | std::ios::in);
    if (!stream.good())
        return 1;
    
    typedef seqan::DirectionIterator<std::fstream, seqan::Input>::Type TReader;
    TReader reader = directionIterator(stream, seqan::Input());

    seqan::StringSet<seqan::CharString> result;
    
    while (!atEnd(reader))
    {
        resize(result, length(result) + 1);
        readLine(back(result), reader);
    }

    return 0;
}
