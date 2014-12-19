#include <seqan/sequence.h>
#include <seqan/index.h>

using namespace seqan;

int main()
{
    DnaString genome = "ACGACGTGCAACGTACGACTAGCATCGGATCAGCAT";

    Shape<Dna, OneGappedShape> myShape;
    stringToShape(myShape, "1101");

    // compute hash of a search pattern
    unsigned hashedPattern = hash(myShape, "ACGA");
    std::cout << "The hash is: " << hashedPattern << std::endl;

    // compute all overlapping hashes and compare with hash of pattern
    for (unsigned i = 0; i < length(genome) - length(myShape) + 1; ++i)
        if (hash(myShape, begin(genome) + i) == hashedPattern)
            std::cout << "Hit at position: " << i << std::endl;

    return 0;
}
