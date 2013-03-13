// FRAGMENT(header)
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/seeds.h>

int main()
{
    // FRAGMENT(example)
    typedef seqan::Seed<seqan::Simple>    TSeed;
    typedef seqan::SeedSet<seqan::Simple> TSeedSet;
    
    TSeedSet seedSet;
    addSeed(seedSet, TSeed(0, 0, 2), seqan::Single());
    addSeed(seedSet, TSeed(3, 5, 2), seqan::Single());
    addSeed(seedSet, TSeed(4, 2, 3), seqan::Single());
    addSeed(seedSet, TSeed(9, 9, 2), seqan::Single());

    std::cout << "Resulting seeds.\n";
    typedef seqan::Iterator<TSeedSet>::Type TIter;
    for (TIter it = begin(seedSet, seqan::Standard()); it != end(seedSet, seqan::Standard()); ++it)
        std::cout << "(" << beginPositionH(*it) << ", " << endPositionH(*it)
                  << ", " << beginPositionV(*it) << ", " << endPositionV(*it)
                  << ", " << lowerDiagonal(*it) << ", " << upperDiagonal(*it)
                  << ")\n";

    // FRAGMENT(footer)
    return 0;
}
