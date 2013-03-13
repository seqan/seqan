// FRAGMENT(header)
#include <seqan/file.h>
#include <seqan/score.h>
#include <seqan/seeds.h>
#include <seqan/sequence.h>

int main()
{
    // FRAGMENT(example)
    typedef seqan::Seed<seqan::Simple>    TSeed;
    typedef seqan::SeedSet<seqan::Simple> TSeedSet;
    
    seqan::Dna5String seqH;
    seqan::Dna5String seqV;
    seqan::Score<int, seqan::Simple> scoringScheme(1, -1, -1);
    
    seqan::String<TSeed> seeds;
    appendValue(seeds, TSeed(0, 0, 2));
    appendValue(seeds, TSeed(3, 5, 2));
    appendValue(seeds, TSeed(4, 2, 3));
    appendValue(seeds, TSeed(9, 9, 2));
    
    TSeedSet seedSet;
    for (unsigned i = 0; i < length(seeds); ++i)
    {
        if (!addSeed(seedSet, seeds[i], 2, 2, scoringScheme,
                     seqH, seqV, seqan::SimpleChain()))
            addSeed(seedSet, seeds[i], seqan::Single());
    }

    std::cout << "Resulting seeds.\n";
    typedef seqan::Iterator<TSeedSet>::Type TIter;
    for (TIter it = begin(seedSet, seqan::Standard());
         it != end(seedSet, seqan::Standard()); ++it)
        std::cout << "(" << beginPositionH(*it) << ", " << endPositionH(*it)
                  << ", " << beginPositionV(*it) << ", " << endPositionV(*it)
                  << ", " << lowerDiagonal(*it) << ", " << upperDiagonal(*it)
                  << ")\n";

    // FRAGMENT(footer)
    return 0;
}
