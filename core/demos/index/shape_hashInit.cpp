#include <seqan/sequence.h>
#include <seqan/index.h>

using namespace seqan;

int main()
{
    DnaString genome = "AAAACACAGTTTGA";
    Shape<Dna, UngappedShape<3> > myShape;

    // loop with hash() and hashNext() starts at position 1
    std::cout << hash(myShape, begin(genome)) << '\t';
    for (unsigned i = 1; i < length(genome) - length(myShape) + 1; ++i)
        std::cout << hashNext(myShape, begin(genome) + i) << '\t';
    std::cout << std::endl;

    // loop with hashInit() and hashNext() starts at position 0
    hashInit(myShape, begin(genome));
    for (unsigned i = 0; i < length(genome) - length(myShape) + 1; ++i)
        std::cout << hashNext(myShape, begin(genome) + i) << '\t';
    std::cout << std::endl;
}
