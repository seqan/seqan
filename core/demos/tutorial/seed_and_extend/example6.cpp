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

    seqan::String<TSeed> result;
    chainSeedsGlobally(result, seedSet, seqan::SparseChaining());

    // FRAGMENT(footer)
    return 0;
}
