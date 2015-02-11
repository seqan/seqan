#include <iostream>
#include <seqan/sequence.h>
#include <seqan/stream.h>

using namespace seqan;

int main()
{
    String<Dna5String> nucleotidesList;
    Dna5String str1 = "ATATANGCGT";
    Dna5String str2 = "AAGCATGANT";
    Dna5String str3 = "TGAAANTGAC";
    resize(nucleotidesList, 3);
    nucleotidesList[0] = str1;
    nucleotidesList[1] = str2;
    nucleotidesList[2] = str3;

    String<Dna5String> lesser;
    String<Dna5String> greater;
    Dna5String ref = "GATGCATGAT";

    // For each Dna5String of the String:
    for (unsigned i = 0; i < length(nucleotidesList); ++i)
    {
        // Compare the Dna5String with the given reference string
        // The result of the comparison is stored in comp
        Lexical<> comp(nucleotidesList[i], ref);
        // The function isLess checks only the stored result
        // without comparing the sequences again
        if (isLess(comp))
            appendValue(lesser, nucleotidesList[i]);
        else if (isGreater(comp))
            appendValue(greater, nucleotidesList[i]);
    }
    // Print the results
    std::cout << "Lesser sequences: " << std::endl;
    for (unsigned i = 0; i < length(lesser); ++i)
    {
        std::cout << lesser[i] << ", ";
    }
    std::cout << std::endl;
    std::cout << "Greater sequences: " << std::endl;
    for (unsigned i = 0; i < length(greater); ++i)
    {
        std::cout << greater[i] << ", ";
    }
}
