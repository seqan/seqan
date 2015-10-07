///A tutorial about seed extension.
#include <iostream>
#include <seqan/seeds.h>

using namespace seqan;

int main()
{
///Example 1: three algorithms for seed extension.
    String<char> a = "SEEDabXcdXefXXX";
    String<char> b = "SEEDabYcdefYYYY";

    Seed<Simple> seed1(0, 0, 4);          //left=0; length=4
    extendSeed(seed1, a, b, EXTEND_BOTH, MatchExtend());
    std::cout << endPositionH(seed1) << std::endl;  //output: 6
    std::cout << endPositionV(seed1) << std::endl;  //output: 6

    Seed<Simple> seed2(0, 0, 4);          //left=0; length=4
    Score<> scoring(1, -1, -1);
    extendSeed(seed2, a, b, EXTEND_BOTH, scoring, 2, UnGappedXDrop());
    std::cout << endPositionH(seed2) << std::endl;  //output: 9
    std::cout << endPositionV(seed2) << std::endl;  //output: 9

    Seed<Simple> seed3(0, 0, 4);          //left=0; length=4
    extendSeed(seed3, a, b, EXTEND_BOTH, scoring, 2, GappedXDrop());
    std::cout << endPositionH(seed3) << std::endl;  //output: 14
    std::cout << endPositionV(seed3) << std::endl;  //output: 13

///Example 2: global chaining.
    SeedSet<Seed<Simple>, Unordered> seedSet;
    addSeed(seedSet, Seed<Simple>(0, 93, 281, 342), Single());
    addSeed(seedSet, Seed<Simple>(3, 237, 127, 364), Single());
    addSeed(seedSet, Seed<Simple>(3, 284, 86, 368), Single());
    addSeed(seedSet, Seed<Simple>(5, 146, 239, 374), Single());

    String<Seed<Simple> > chain;
    chainSeedsGlobally(chain, seedSet, SparseChaining());

    return 0;
}
