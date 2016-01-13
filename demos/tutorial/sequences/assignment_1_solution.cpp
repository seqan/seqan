// The comment lines containing ![fragment-line] are there for the
// documentation system.  You can ignore them when reading this file.
//![full]
//![top]
#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <seqan/stream.h>
#include <seqan/file.h>
#include <seqan/modifier.h>

using namespace seqan;

Dna getRevCompl(Dna const & nucleotide)
{
    if (nucleotide == 'A')
        return 'T';
    if (nucleotide == 'T')
        return 'A';
    if (nucleotide == 'C')
        return 'G';
    return 'C';
}

int main()
{
    DnaString genome = "TATATACGCGCGAGTCGT";
    DnaString revComplGenome;
//![top]
    //1.
    resize(revComplGenome, length(genome));
    //2.
    for (unsigned i = 0; i < length(genome); ++i)
        revComplGenome[length(genome) - 1 - i] = getRevCompl(genome[i]);
    //3.
    std::cout << genome << std::endl;
    std::cout << revComplGenome << std::endl;
//![bottom]

    // And to check if your output is correct,
    // use the given SeqAn function reverseComplement(),
    // which modifies the sequence in-place:
    reverseComplement(genome);
    std::cout << genome << std::endl;
    return 0;
}
//![bottom]
//![full]
