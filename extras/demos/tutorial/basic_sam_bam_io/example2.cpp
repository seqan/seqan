#include <iostream>
#include <seqan/bam_io.h>

int main()
{
    seqan::BamStream bamInStream("example.bam");

    for (unsigned i = 0; i < length(bamInStream.header.sequenceInfos); ++i)
        std::cout << bamInStream.header.sequenceInfos[i].i1 << '\t'
                  << bamInStream.header.sequenceInfos[i].i2 << '\n';

    return 0;
}
