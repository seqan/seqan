#include <seqan/stream.h>
#include <seqan/sequence.h>
#include <seqan/file.h>

using namespace seqan;

int main()
{
    String<Dna5> nucleotides = "AGTCGTGNNANCT";
    String<Dna5> selected;
    // Append all elements of nucleotides, apart of Gs,
    // to the list selected.
    for (unsigned i = 0; i < length(nucleotides); ++i)
    {
        appendValue(selected, nucleotides[i]);
    }
    std::cout << "Selected nucleotides: " << selected << std::endl;
    return 0;
}
