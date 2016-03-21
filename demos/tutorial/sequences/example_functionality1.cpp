#include <seqan/sequence.h>
#include <seqan/seq_io.h>

using namespace seqan;

int main()
{
//![main]
    String<Dna> dnaSeq = "TATA";
    dnaSeq += "CGCG";
    std::cout << dnaSeq << std::endl;
//![main]
}
