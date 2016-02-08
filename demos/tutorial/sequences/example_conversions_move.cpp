#include <seqan/sequence.h>
#include <seqan/seq_io.h>

using namespace seqan;

int main()
{
//![main]
    String<char> char_source = "acgtgcat";
    String<Dna> dna_target;

    // The in-place move conversion.
    move(dna_target, char_source);
    std::cout << dna_target << std::endl;
//![main]
}
