#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>  // for I/O


using namespace seqan;

int main()
{
    String<Dna> dnaSeq = "TATACGCGAAAA";
    Infix<String<Dna> >::Type inf = infix(dnaSeq, 4, 8);
    std::cout << "Infix: " << inf << std::endl;

    return 0;
}
