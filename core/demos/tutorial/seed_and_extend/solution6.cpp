#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/seeds.h>

int main()
{
    typedef seqan::Seed<seqan::Simple>    TSeed;
    typedef seqan::SeedSet<seqan::Simple> TSeedSet;
    
    TSeedSet seedSet;
    addSeed(seedSet, TSeed(1, 1, 3), seqan::Single());
    addSeed(seedSet, TSeed(6, 9, 2), seqan::Single());
    addSeed(seedSet, TSeed(10, 13, 3), seqan::Single());
    addSeed(seedSet, TSeed(20, 22, 5), seqan::Single());

    seqan::String<TSeed> result;
    chainSeedsGlobally(result, seedSet, seqan::SparseChaining());

    return 0;
}
