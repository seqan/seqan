#include <seqan/stream.h>
#include <seqan/score.h>
#include <seqan/seeds.h>
#include <seqan/sequence.h>

using namespace seqan;

int main()
{
    typedef Seed<Simple>    TSeed;
    typedef SeedSet<TSeed> TSeedSet;

    Dna5String seqH;
    Dna5String seqV;
    Score<int, Simple> scoringScheme(1, -1, -1);

    String<TSeed> seeds;
    appendValue(seeds, TSeed(0, 0, 2));
    appendValue(seeds, TSeed(3, 5, 2));
    appendValue(seeds, TSeed(4, 2, 3));
    appendValue(seeds, TSeed(9, 9, 2));

    TSeedSet seedSet;
    for (unsigned i = 0; i < length(seeds); ++i)
    {
        if (!addSeed(seedSet, seeds[i], 2, 2, scoringScheme,
                     seqH, seqV, Chaos()))
            addSeed(seedSet, seeds[i], Single());
    }

    std::cout << "Resulting seeds.\n";
    typedef Iterator<TSeedSet>::Type TIter;
    for (TIter it = begin(seedSet, Standard());
         it != end(seedSet, Standard()); ++it)
        std::cout << "(" << beginPositionH(*it) << ", " << endPositionH(*it)
                  << ", " << beginPositionV(*it) << ", " << endPositionV(*it)
                  << ", " << lowerDiagonal(*it) << ", " << upperDiagonal(*it)
                  << ")\n";

    return 0;
}
