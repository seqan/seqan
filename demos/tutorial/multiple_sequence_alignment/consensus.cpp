//![align]
#include <iostream>
#include <seqan/align.h>
#include <seqan/graph_msa.h>

using namespace seqan;

int main()
{
    // some variangs of sonic hedgehog exon 1
    char const * strings[4] =
    {
        // gi|2440284|dbj|AB007129.1| Oryzias latipes
        "GCGGGTCACTGAGGGCTGGGATGAGGACGGCCACCACTTCGAGGAGTCCCTTCACTACGAGGGCAGGGCC"
        "GTGGACATCACCACGTCAGACAGGGACAAGAGCAAGTACGGCACCCTGTCCAGACTGGCGGTGGAAGCTG"
        "GGTTCGACTGGGTCTACTATGAGTCCAAAGCGCACATCCACTGCTCTGTGAAAGCAGAAAGCTCAGTCGC"
        "TGCAAAGTCGGGCGGTTGCTTCCCAGGATCCTCCACGGTCACCCTGGAAAATGGCACCCAGAGGCCCGTC"
        "AAAGATCTCCAACCCGGGGACAGAGTACTGGCCGCGGATTACGACGGAAACCCGGTTTATACCGACTTCA"
        "TCATGTTCAA",
        // gi|1731488|gb|U51350.1|DDU51350 Devario devario
        "CTACGGCAGAAGAAGACATCCGAAAAAGCTGACACCTCTCGCCTACAAGCAGTTCATACCTAATGTCGCG"
        "GAGAAGACCTTAGGGGCCAGCGGCAGATACGAGGGCAAGATAACGCGCAATTCGGAGAGATTTAAAGAAC"
        "TTACTCCAAATTACAATCCCGACATTATCTTTAAGGATGAGGAGAACACG",
        // gi|1731504|gb|U51352.1|PTU51352 Puntius tetrazona
        "CTACGGCAGAAGAAGACATCCCAAGAAGCTGACACCTCTCGCCTACAAGCAGTTTATACCTAATGTCGCG"
        "GAGAAGACCTTAGGGGCCAGCGGCAGATACGAGGGCAAGATCACGCGCAATTCGGAGAGATTTAAAGAAC"
        "TTACTCCAAATTACAATCCCGACATTATCTTTAAGGATGAGGAGAACACT",
        // gi|54399708|gb|AY642858.1| Bos taurus
        "TGCTGCTGCTGGCGAGATGTCTGCTGGTGCTGCTTGTCTCCTCGCTGTTGATGTGCTCGGGGCTGGCGTG"
        "CGGACCCGGCAGGGGATTTGGCAAGAGGCGGAACCCCAAAAAGCTGACCCCTTTAGCCTACAAGCAGTTT"
        "ATCCCCAACGTGGCGGAGAAGACCCTAGGGGCCAGTGGAAGATATGAGGGGAAGATCACCAGAAACTCAG"
        "AGCGATTTAAGGAACTCACCCCCAATTACAACCC"
    };

    Align<DnaString> align;
    resize(rows(align), 4);
    for (int i = 0; i < 4; ++i)
        assignSource(row(align, i), strings[i]);

    globalMsaAlignment(align, SimpleScore(5, -3, -1, -3));
    std::cout << align << "\n";
//![align]

//![profile-computation]
    // create the profile string
    String<ProfileChar<Dna> > profile;
    resize(profile, length(row(align, 0)));
    for (unsigned rowNo = 0; rowNo < 4u; ++rowNo)
        for (unsigned i = 0; i < length(row(align, rowNo)); ++i)
            profile[i].count[ordValue(getValue(row(align, rowNo), i))] += 1;
//![profile-computation]

//![consensus-calling]
    // call consensus from this string
    DnaString consensus;
    for (unsigned i = 0; i < length(profile); ++i)
    {
        int idx = getMaxIndex(profile[i]);
        if (idx < 4)  // is not gap
            appendValue(consensus, Dna(getMaxIndex(profile[i])));
    }

    std::cout << "consensus sequence is\n"
              << consensus << "\n";

    return 0;
}
//![consensus-calling]
