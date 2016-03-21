// The comment lines containing ![fragment-line] are there for the
// documentation system.  You can ignore them when reading this file.
//![full]
//![top]
#include <iostream>
#include <seqan/sequence.h>
#include <seqan/stream.h>

using namespace seqan;

// Function to print simple alignment between two sequences with the same length
// .. for two sequences of different types
template <typename TText1, typename TText2>
void printAlign(TText1 const & genomeFragment, TText2 const & read)
{
    std::cout <<  "Alignment " << std::endl;
    std::cout << "  genome : ";
    std::cout << genomeFragment << std::endl;
    std::cout << "  read   : ";
    std::cout << read << std::endl;
}

int main()
{
    // We have given a genome sequence
    Dna5String genome = "ATGGTTTCAACGTAATGCTGAACATGTCGCGT";
    // A read sequence
    Dna5String read = "TGGTNTCA";
    // And the begin position of a given alignment between the read and the genome
    unsigned beginPosition = 1;
//![top]

    // Create Infix of type Dna5String and get the corresponding infix sequence of genome
    Infix<Dna5String>::Type genomeFragment = infix(genome, beginPosition, beginPosition + length(read));

//![bottom]
    // Call of our function to print the simple alignment
    printAlign(genomeFragment, read);
    return 0;
}
//![bottom]
//![full]
