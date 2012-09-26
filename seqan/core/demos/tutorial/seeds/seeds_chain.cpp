// FRAGMENT(definition)
#include <seqan/seeds.h>

using namespace seqan;

// FRAGMENT(main)
int main() {
    DnaString a = "atattgatccagttactccagg";
    DnaString b = "cgattgctccgtacccatg";

    typedef Seed<int, SimpleSeed> TSeed;
	typedef Value<TSeed>::Type TPosition;
    TPosition qGramsA[] = {2, 3, 7, 8,  13, 15, 16, 17};
    TPosition qGramsB[] = {2, 3, 7, 14, 11,  6,  7, 14};

// FRAGMENT(seedSet)
    typedef int TScore;
    Score<TScore, Simple> scoreMatrix(1, -1, -1);
    TScore scoreMin = 3;
    TPosition limit = 1;
    typedef SeedSet<TPosition, SimpleSeed, DefaultScore> TSeedSet;
    TSeedSet seedSet(limit, scoreMin, scoreMatrix);

// FRAGMENT(combining-seeds)
    for (unsigned i = 0; i < 8; ++i) {
        if (!addSeed(seedSet, qGramsA[i], qGramsB[i], 3, 0, Merge())) {
            addSeed(seedSet, qGramsA[i], qGramsB[i], 3, Single());
        }
    }

    std::cout << length(seedSet) << " seeds in seedSet: " << std::endl;
    Iterator<TSeedSet>::Type it = begin(seedSet);
    while(it != end(seedSet)) {
        std::cout << " <" << leftPosition(*it, 0) << "," << rightPosition(*it, 0);
        std::cout << "> , <" << leftPosition(*it, 1) << "," << rightPosition(*it, 1) << ">" << std::endl;
        ++it;
    }

// FRAGMENT(chaining)
    String<TSeed> chain;
    globalChaining(seedSet, chain);

    std::cout << std::endl << "Chain: " << std::endl;
    for(Iterator<String<TSeed> >::Type itChain = begin(chain); itChain != end(chain); ++itChain) {
        std::cout << " <" << leftPosition(*itChain, 0) << "," << rightPosition(*itChain, 0);
        std::cout << "> , <" << leftPosition(*itChain, 1) << "," << rightPosition(*itChain, 1) << ">" << std::endl;
    }

// FRAGMENT(banded-chain-alignment)
    Align<DnaString> align;
    resize(rows(align), 2);
    assignSource(row(align, 0), a);
    assignSource(row(align, 1), b);

    int band = 2;
    bandedChainAlignment(chain, band, align, scoreMatrix);

    std::cout << std::endl << "Banded chain alignment: " << std::endl;
    std::cout << align << std::endl;

    return 0;
}
