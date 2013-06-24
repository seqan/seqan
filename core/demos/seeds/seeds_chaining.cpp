//An example for global seed chaining.
#include <seqan/seeds.h>

using namespace seqan;

int main()
{
    SeedSet<Simple, Unordered> seedSet;
    addSeed(seedSet, Seed<Simple>(0, 93, 281, 342), Single());
    addSeed(seedSet, Seed<Simple>(3, 237, 127, 364), Single());
    addSeed(seedSet, Seed<Simple>(3, 284, 86, 368), Single());
    addSeed(seedSet, Seed<Simple>(5, 146, 239, 374), Single());

    String<Seed<Simple> > chain;
    chainSeedsGlobally(chain, seedSet, SparseChaining());

    return 0;
}
