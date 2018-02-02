#include <seqan/sequence.h>
#include <seqan/file.h>

using namespace seqan;

int main()
{
//![string_example]
    String<char>  myText;     // A string of characters.
    String<int>   myNumbers;  // A string of integers.
//![string_example]

//![simple_string_construction]
    myText = "Hello SeqAn!";
//![simple_string_construction]

//![string_of_strings_example]
    String<String<char> >   myStringList;   // A string of character strings.
//![string_of_strings_example]

//![special_types_example]
    String<Dna>         myGenome;   // A string of nucleotides.
    String<AminoAcid>   myProtein;  // A string of amino acids.
//![special_types_example]

//![shortcuts_example]
    // Instead of String<Dna> dnaSeq we can also write:
    DnaString dnaSeq = "TATA";
//![shortcuts_example]

//![specification_example]
    String<Dna>              myGenome1;   // A default string of nucleotides.
    String<Dna, Alloc<> >    myGenome2;    // The same type as above.
//![specification_example]

//![external_string_spec]
    // Most of the string is stored on the disk.
    String<Dna, External<> > myLargeGenome;
//![external_string_spec]

//![initialization_example]
    String<Dna> myDnaGenome = "TATACGCG";
//![initialization_example]

    Dna5String genome = "ATGGTTTCAACGTAATGCTGAACATGTCGCGT";
    Dna5String read = "TGGTNTCA";
    unsigned beginPosition = 1;
//![assignment5_code_to_change]
    // We have to create a copy of the corresponding fragment of the genome, where the read aligns to
    // Change this piece of code using an infix of the genome
    Dna5String genomeFragment;
    for (unsigned i = 0; i < length(read); ++i)
    {
        appendValue(genomeFragment, genome[beginPosition + i]);
    }    
//![assignment5_code_to_change]
} 
