#include <seqan/sequence.h>
#include <seqan/seq_io.h>

using namespace seqan;

int main()
{
//![main]
    String<Dna> dna_source = "acgtgcat";
    String<char> char_target;
    assign(char_target, dna_source);
    std::cout << char_target << std::endl;
//![main]
}
