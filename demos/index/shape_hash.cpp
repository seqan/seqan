#include <seqan/sequence.h>
#include <seqan/index.h>

using namespace seqan;

int main()
{
    DnaString text = "GATTACA";
    // output all hash values as hexadecimal numbers
    std::cout << std::hex;

    // 4-gram with shape 1111 at position 0 is GATT
    // its hash value is 0b10001111 = 0x8f
    Shape<Dna, UngappedShape<4> > shape1;
    std::cout << "0x" << hash(shape1, begin(text)) << std::endl;

    // 4-gram with shape 110101 at position 0 is GATC
    // its hash value is 0b10001101 = 0x8d
    Shape<Dna, GenericShape> shape2;
    stringToShape(shape2, "110101");
    std::cout << "0x" << hash(shape2, begin(text)) << std::endl;

    // 4-gram with shape 11011 at position 0 is GATA
    // the hash value is 0b10001100 = 0x8c
    Shape<Dna, OneGappedShape> shape3;
    stringToShape(shape2, "11011");
    std::cout << "0x" << hash(shape2, begin(text)) << std::endl;

    return 0;
}
