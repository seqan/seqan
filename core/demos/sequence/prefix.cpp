#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>  // for I/O


using namespace seqan;

int main()
{
    String<Dna> dnaSeq = "TATACGCGAAAA";
    Infix<String<Dna> >::Type pre = prefix(dnaSeq, 4);
    std::cout << "Prefix: " << pre << std::endl;

    return 0;
}
