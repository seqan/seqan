#include <seqan/stream.h>

using namespace seqan;

int main()
{
//![example]
    String<Dna> dnaSeq;
    // Sets the capacity of dnaSeq to 5.
    resize(dnaSeq, 4, Exact());
    // Only "TATA" is assigned to dnaSeq, since dnaSeq is limited to 4.
    assign(dnaSeq, "TATAGGGG", Limit());
    std::cout << dnaSeq << std::endl;
    // Use the default expansion strategy.
    append(dnaSeq, "GCGCGC");
    std::cout << dnaSeq << std::endl;
//![example]
    return 0;
}
