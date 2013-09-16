#include <fstream>
#include <seqan/sequence.h>
#include <seqan/stream.h>

int main()
{
    std::fstream in("in.txt", std::ios::binary | std::ios::in);
    std::fstream out("out.txt", std::ios::binary | std::ios::out);

    seqan::CharString buffer;
    resize(buffer, 1000);

    while (!seqan::streamEof(in) && seqan::streamError(in) == 0)
    {
        int num = seqan::streamReadBlock(&buffer[0], in, length(buffer));
        seqan::streamWriteBlock(out, &buffer[0], num);
    }

    return 0;
}
