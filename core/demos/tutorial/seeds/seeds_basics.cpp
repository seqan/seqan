// FRAGMENT(definition)
#include <seqan/seeds.h>

using namespace seqan;

int main() {
    CharString seq0 = "drosophila";
    CharString seq1 = "phylosophy";

    typedef Seed<int, SimpleSeed> TSeed;
    TSeed seed(2,4,5);

// FRAGMENT(left-right-positions)
    std::cout << "Seed from position " << leftPosition(seed, 0);
    std::cout << " to " << rightPosition(seed, 0) << " in seq0=" << seq0 << ": ";
    std::cout << infix(seq0, leftPosition(seed, 0), rightPosition(seed, 0)+1) << std::endl;
    std::cout << "Seed from position " << leftPosition(seed, 1);
    std::cout << " to " << rightPosition(seed, 1) << " in seq1=" << seq1 << ": ";
    std::cout << infix(seq1, leftPosition(seed, 1), rightPosition(seed, 1)+1) << std::endl;

    return 0;
}
