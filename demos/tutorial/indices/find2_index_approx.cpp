//![Complete]
#include <seqan/index.h>

using namespace seqan;

int main()
{
    DnaString genome("GAGAGGCCACTCGCAGGATTAAGTCAATAAGTTAATGGCGTGGGGTTATGGTATGTAGTACGACGCCCACAGTGACCTCATCGGTGCATTTCCTCATCGTAGGCGGAACGGTAGACACAAGGCATGATGTCAAATCGCGACTCCAATCCCAAGGTCGCAAGCCTATATAGGAACCCGCTTATGCCCTCTAATCCCGGACAGACCCCAAATATGGCATAGCTGGTTGGGGGTACCTACTAGGCACAGCCGGAAGCA");
    Index<DnaString, BidirectionalIndex<FMIndex<> > > index(genome);

    //![Delegate]
    auto delegate = [](auto & iter)
    {
        for (auto occ : getOccurrences(iter))
            std::cout << occ << std::endl;
    };
    //![Delegate]

    DnaString pattern("GGGGTTAT");
    std::cout << "Hits with up to 2 errors:" << std::endl;
    //![SinglePattern]
    find<0, 2>(delegate, index, pattern, HammingDistance());
    //![SinglePattern]

    StringSet<DnaString> patterns;
    appendValue(patterns, "GGGGTTAT");
    appendValue(patterns, "CTAGCTAA");
    std::cout << "Hits with 2-3 errors:" << std::endl;
    //![MultiplePatterns]
    find<2, 3>(delegate, index, patterns, HammingDistance(), Serial());
    //![MultiplePatterns]

    return 0;
}
//![Complete]
