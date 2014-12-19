#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/seeds.h>

using namespace seqan;

int main()
{
    typedef Seed<Simple>    TSeed;
    typedef SeedSet<TSeed> TSeedSet;

    TSeedSet seedSet;
    addSeed(seedSet, TSeed(1, 1, 3), Single());
    addSeed(seedSet, TSeed(6, 9, 2), Single());
    addSeed(seedSet, TSeed(10, 13, 3), Single());
    addSeed(seedSet, TSeed(20, 22, 5), Single());

    String<TSeed> result;
    chainSeedsGlobally(result, seedSet, SparseChaining());

    return 0;
}
