#include <seqan/find.h>
#include <seqan/index.h>

using namespace seqan;

int main()
{
    DnaString genome = "ACGACGTGCAACGTACGACTAGCATCGGATCAGCAT";

    Shape<Dna, OneGappedShape> myShape;
    stringToShape(myShape, "1101");

    unsigned correctHash = hash(myShape, "ACGA");  // could also be ACAA, ACCA, ACTA
    std::cout << "The hash is: " << correctHash << std::endl;

    for (unsigned i = 0; i < length(genome) - length(myShape) + 1; ++i)
        if (hash(myShape, &genome[i]) == correctHash)
            std::cout << "Hit at position: " << i <<std::endl;
}
