//![header]
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/seeds.h>

using namespace seqan;

int main()
{
//![header]
//![example]
    typedef Seed<Simple>    TSeed;
    typedef SeedSet<TSeed> TSeedSet;

    TSeedSet seedSet;
    addSeed(seedSet, TSeed(0, 0, 2), Single());
    addSeed(seedSet, TSeed(3, 5, 2), Single());
    addSeed(seedSet, TSeed(4, 2, 3), Single());
    addSeed(seedSet, TSeed(9, 9, 2), Single());

    std::cout << "Resulting seeds.\n";
    typedef Iterator<TSeedSet>::Type TIter;
    for (TIter it = begin(seedSet, Standard()); it != end(seedSet, Standard()); ++it)
        std::cout << "(" << beginPositionH(*it) << ", " << endPositionH(*it)
                  << ", " << beginPositionV(*it) << ", " << endPositionV(*it)
                  << ", " << lowerDiagonal(*it) << ", " << upperDiagonal(*it)
                  << ")\n";
//![example]

//![footer]
    return 0;
}
//![footer]
