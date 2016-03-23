#include <seqan/sequence.h>
#include <seqan/align.h>
#include <seqan/score.h>
#include <seqan/seeds.h>

using namespace seqan;

int main()
{
    typedef Seed<Simple>    TSeed;

    Dna5String sequenceH = "CGAATCCATCCCACACA";
    Dna5String sequenceV = "GGCGATNNNCATGGCACA";
    Score<int, Simple> scoringSchemeAnchor(0, -1, -1);
    Score<int, Simple> scoringSchemeGap(2, -1, -1, -2);

    String<TSeed> seedChain;
    appendValue(seedChain, TSeed(0, 2, 5, 6));
    appendValue(seedChain, TSeed(6, 9, 9, 12));
    appendValue(seedChain, TSeed(11, 14, 17, 16));

    Align<Dna5String, ArrayGaps> alignment;
    resize(rows(alignment), 2);
    assignSource(row(alignment, 0), sequenceH);
    assignSource(row(alignment, 1), sequenceV);
    AlignConfig<true, false, false, true> alignConfig;

    int result = bandedChainAlignment(alignment, seedChain, scoringSchemeAnchor, scoringSchemeGap, alignConfig, 2);

    std::cout << "Score: " << result << std::endl;
    std::cout << alignment << std::endl;

    return 0;
}
