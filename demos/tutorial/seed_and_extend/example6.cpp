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

    String<TSeed> result;
    chainSeedsGlobally(result, seedSet, SparseChaining());
//![example]

//![footer]
    return 0;
}
//![footer]
