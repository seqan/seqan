#include <iostream>
#include <seqan/sequence.h>
#include <seqan/stream.h>

using namespace seqan;
// Function to print simple alignment between two sequences with the same length
template <typename TText1, typename TText2>
void printAlign(TText1 const & genomeFragment, TText2 const & read)
{
    std::cout <<  "Alignment " << std::endl;
    std::cout << "  genome : " << genomeFragment << std::endl;
    std::cout << "  read   : " << read << std::endl;
}

int main(int, char const **)
{
    // Build reads and genomes
    DnaString chr1 = "TATAATATTGCTATCGCGATATCGCTAGCTAGCTACGGATTATGCGC"
                     "TCTGCGATATATCGCGCTAGATGTGCAGCTCGATCGAATGCACGTGT"
                     "GTGCGATCGATTAGCGTCGATCATCGATCTATATTAGCGCGCGGTAT"
                     "CGGACGATCATATTAGCGGTCTAGCATTTAG";

    // Build List containing all reads
    typedef String<DnaString> DnaList;
    DnaList readList;
    resize(readList, 4);
    readList[0] = "TTGCTATCGCGATATCGCTAGCTAGCTACGGATTATGCGCTCTGCGATATATCGCGCT";
    readList[1] = "TCGATTAGCGTCGATCATCGATCTATATTAGCGCGCGGTATCGGACGATCATATTAGCGGTCTAGCATT";
    readList[2] = "AGCCTGCGTACGTTGCAGTGCGTGCGTAGACTGTTGCAAGCCGGGGGTTCATGTGCGCTGAAGCACACATGCACA";
    readList[3] = "CGTGCACTGCTGACGTCGTGGTTGTCACATCGTCGTGCGTGCGTACTGCTGCTGACA";

    // Append a second chromosome sequence fragment to chr1
    DnaString chr2 = "AGCCTGCGTACGTTGCAGTGCGTGCGTAGACTGTTGCAAGCCGGGGGT"
                     "TCATGTGCGCTGAAGCACACATGCACACGTCTCTGTGTTCCGACGTGT"
                     "GTCACGTGCACTGCTGACGTCGTGGTTGTCACATCGTCGTGCGTGCGT"
                     "ACTGCTGCTGACACATGCTGCTG";
    append(chr1, chr2);

    // Print readlist
    std::cout << " \n Read list: " << std::endl;
    for(unsigned i = 0; i < length(readList); ++i)
        std::cout << readList[i] << std::endl;

    // Assume we have mapped the 4 reads to chr1 (and chr2) and now have the mapping start positions (no gaps).
    // Store the start position in a String: 7, 100, 172, 272
    String<unsigned> alignPosList;
    resize(alignPosList, 4);
    alignPosList[0] = 7;
    alignPosList[1] = 100;
    alignPosList[2] = 172;
    alignPosList[3] = 272;

    // Print alignments using Segment
    std::cout << " \n Print alignment using Segment: " << std::endl;
    for(unsigned i = 0; i < length(readList); ++i)
    {
        // Begin and end position of a given alignment between the read and the genome
        unsigned beginPosition = alignPosList[i];
        unsigned endPosition = beginPosition + length(readList[i]);
        // Build infix
        Infix<DnaString>::Type genomeFragment = infix(chr1, beginPosition, endPosition);
        // Call of our function to print the simple alignment
        printAlign(genomeFragment, readList[i]);
    }

    // Iterators :)
    // Print alignments using Iterators: Do the same as above, but use Iterators to iterate over the read list.
    // First, use Standard Iterators: Build two iterators it and itEnd to traverse readList.

    std::cout << " \n Print alignment using Standard Iterators: " << std::endl;

    return 0;
}
