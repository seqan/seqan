#include <seqan/sequence.h>
#include <seqan/align.h>
#include <seqan/score.h>
#include <seqan/seeds.h>

int main()
{
    typedef seqan::Seed<seqan::Simple>    TSeed;

    seqan::Dna5String sequenceH = "CGAATCCATCCCACACA";
    seqan::Dna5String sequenceV = "GGCGATNNNCATGGCACA";
    seqan::Score<int, seqan::Simple> scoringSchemeAnchor(0, -1, -1);
    seqan::Score<int, seqan::Simple> scoringSchemeGap(2, -1, -1, -2);

    seqan::String<TSeed> seedChain;
    seqan::appendValue(seedChain, TSeed(0, 2, 5, 6));
    seqan::appendValue(seedChain, TSeed(6, 9, 9, 12));
    seqan::appendValue(seedChain, TSeed(11, 14, 17, 16));

    seqan::Align<seqan::Dna5String, seqan::ArrayGaps> alignment;
    seqan::resize(seqan::rows(alignment), 2);
    seqan::assignSource(seqan::row(alignment, 0), sequenceH);
    seqan::assignSource(seqan::row(alignment, 1), sequenceV);
    seqan::AlignConfig<true, false, false, true> alignConfig;

    int result = seqan::bandedChainAlignment(alignment, seedChain, scoringSchemeAnchor, scoringSchemeGap, alignConfig, 2);

    std::cout << "Score: " << result << std::endl;
    std::cout << alignment << std::endl;

    return 0;
}
