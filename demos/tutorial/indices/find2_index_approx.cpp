//![Complete]
#include <set>
#include <mutex>

#include <seqan/index.h>

using namespace seqan;

int main()
{

    DnaString genome(
        "GAGAGGCCACTCGCAGGATTAAGTCAATAAGTTAATGGCGTGGGGTTATGGTATGGGGGTTCTCGCCCACAGTGACCTCATCGGT"
        "GCATTTCCTCATCGTAGGCGGAACGGTAGACACAAGGCATGATGTCAAATCGCGACTCCAATCCCAAGGTCGCAAGCCTATATAG"
        "GAACCCGCTTATGCCCTCTAATCCCGGACAGACCCCAAATATGGCATAGCTGGTTGGGGGTACCTACTAGGCACAGCCGGAAGCA");
    Index<DnaString, BidirectionalIndex<FMIndex<> > > index(genome);

    //![Delegate]
    auto delegate = [](auto & iter, DnaString const & needle, uint8_t errors)
    {
        for (auto occ : getOccurrences(iter))
            std::cout << occ << std::endl;
    };
    //![Delegate]

    DnaString pattern("GGGGTTAT");
    std::cout << "Hits with up to 2 errors (HammingDistance):" << std::endl;
    //![SinglePattern]
    find<0, 2>(delegate, index, pattern, HammingDistance());
    //![SinglePattern]

    StringSet<DnaString> patterns;
    appendValue(patterns, "GGGGTTAT");
    appendValue(patterns, "CTAGCTAA");
    std::cout << "Hits with 1-2 errors (HammingDistance):" << std::endl;
    //![MultiplePatterns]
    find<1, 2>(delegate, index, patterns, HammingDistance(), Serial());
    //![MultiplePatterns]

    //![ParallelMode]
    std::mutex mtx;
    std::set<Pair<DnaString, unsigned> > hits;
    auto delegateParallel = [&hits, &mtx](auto & iter, DnaString const & needle, uint8_t errors)
    {
        std::lock_guard<std::mutex> lck(mtx); // critical section below this line
        for (auto occ : getOccurrences(iter))
            hits.insert(Pair<DnaString, unsigned>(needle, occ));
    };
    find<2, 3>(delegateParallel, index, patterns, HammingDistance(), Parallel());
    std::cout << "Hits with 2-3 errors (HammingDistance, no duplicates):" << std::endl;
    for (auto hit : hits)
        std::cout << hit << std::endl;
    //![ParallelMode]

    return 0;
}
//![Complete]
