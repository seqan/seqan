#include <iostream>
#include <seqan/bam_io.h>

int main()
{
    seqan::BamFileIn bamFileIn("example.bam");
    seqan::BamHeader header;
    readRecord(header, bamFileIn);

    for (unsigned i = 0; i < length(nameStore(context(bamFileIn))); ++i)
        std::cout << nameStore(context(bamFileIn))[i] << '\t'
                  << sequenceLengths(context(bamFileIn))[i] << '\n';

    return 0;
}
