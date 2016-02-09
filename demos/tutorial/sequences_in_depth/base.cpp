#include <seqan/sequence.h>

using namespace seqan;

int main()
{
//![default_type]
    String<Dna>           dnaSeq1; // The default string implementation: Alloc
    String<Dna, Alloc<> > dnaSeq2; // The same as above
//![default_type]

//![type_examples]
    // String with maximum length 100.
    String<char, Array<100> > myArrayString;
    // String that takes only 2 bits per nucleotide.
    String<Dna, Packed<> > myPackedString;
//![type_examples]
    return 0;
}
