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
    
    seqan::RecordReader<std::fstream, seqan::SinglePass<> > reader(stream);
    seqan::StringSet<seqan::CharString> result;
    
    while (!atEnd(reader))
    {
        resize(result, length(result) + 1);
        int res = readLine(back(result), reader);
        if (res != 0)
            return 1;
    }

    return 0;
}
