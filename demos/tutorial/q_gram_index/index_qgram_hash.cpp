//![includes]
#include <seqan/sequence.h>
#include <seqan/index.h>

using namespace seqan;

int main()
{
//![includes]
//![hash_loop1]
    DnaString text = "AAAACACAGTTTGA";
    Shape<Dna, UngappedShape<3> > myShape;

    std::cout << hash(myShape, begin(text)) << '\t';
    for (unsigned i = 1; i < length(text) - length(myShape) + 1; ++i)
        std::cout << hashNext(myShape, begin(text) + i) << '\t';
//![hash_loop1]

//![hash_loop2]
    hashInit(myShape, begin(text));
    for (unsigned i = 0; i < length(text) - length(myShape) + 1; ++i)
        std::cout << hashNext(myShape, begin(text) + i) << '\t';
//![hash_loop2]

//![end]
    return 0;
}
//![end]
