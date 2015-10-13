#include <seqan/sequence.h>
#include <seqan/index.h>

using namespace seqan;

int main()
{
    DnaString text = "AAAACACAGTTTGA";
    Shape<Dna, UngappedShape<3> > myShape;

    // output the hash value of all overlapping 3-grams
    std::cout << hash(myShape, begin(text));
    for (unsigned i = 1; i < length(text) - length(myShape) + 1; ++i)
        std::cout << '\t' << hashNext(myShape, begin(text) + i);

    return 0;
}
