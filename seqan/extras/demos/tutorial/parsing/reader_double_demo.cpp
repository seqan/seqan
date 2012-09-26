#include <iostream>
#include <fstream>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>

using namespace seqan;

int main(int argc, char const ** argv)
{
    if (argc != 2)
        return 1;
    std::fstream stream(argv[1], std::ios::binary | std::ios::in);
    if (!stream.good())
        return 1;
    
    RecordReader<std::fstream, DoublePass<> > reader(stream);
    String<CharString> result;
    CharString buffer;
    
    while (!atEnd(reader))
    {
        startFirstPass(reader);
        clear(buffer);
        int res = readLine(buffer, reader);
        if (res != 0)
            return 1;

        startSecondPass(reader);
        resize(result, length(result) + 1);
        resize(back(result), length(buffer), Exact());
        res = readLine(back(result), reader);
        if (res != 0)
            return 1;
    }

    return 0;
}
