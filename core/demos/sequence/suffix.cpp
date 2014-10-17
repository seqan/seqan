#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>  // for I/O


using namespace seqan;

int main()
{
    String<Dna> dnaSeq = "TATACGCGAAAA";
    Infix<String<Dna> >::Type suf = suffix(dnaSeq, 8);
    std::cout << "Suffix: " << suf << std::endl;

    return 0;
}
