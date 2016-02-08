#include <seqan/stream.h>
#include <seqan/sequence.h>
#include <seqan/file.h>

using namespace seqan;

int main()
{
//![main]
    String<Dna> dnaSeq = "AGTTGGCATG";
    Prefix<String<Dna> >::Type pre = prefix(dnaSeq, 4);
    std::cout << "Prefix: " << pre << std::endl;

    Infix<String<Dna> >::Type inf = infix(dnaSeq, 4, 7);
    std::cout << "Infix: " << inf << std::endl;

    Suffix<String<Dna> >::Type suf = suffix(dnaSeq, 4);
    std::cout << "Suffix: " << suf << std::endl;
//![main]
    return 0;
}
