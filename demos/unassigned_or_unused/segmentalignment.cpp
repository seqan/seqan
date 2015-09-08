#include <iostream>
#include <seqan/graph_msa.h>
#include <seqan/index.h>
#include <seqan/seq_io.h>

using namespace seqan;

int main(int argc, const char * argv[])
{
    typedef String<Dna> TSequence;
    typedef StringSet<TSequence> TSequenceSet;

    if (argc < 3)
    {
        std::cout << "\nAt least two Fasta sequence files need to be specified.\n\n";
        std::cout << "Usage: ./segmentalignment <seq1.fa> <seq2.fa> ...\n\n";
        return 1;
    }

    CharString id;
    TSequence seq;
    TSequenceSet seqs;
    for (int i = 1; i < argc; ++i)
    {
        SeqFileIn seqFileIn(argv[i]);
        readRecord(id, seq, seqFileIn);
        appendValue(seqs, seq);
    }


    typedef Fragment<> TMatch;
    String<TMatch> matches;

    typedef Index<TSequenceSet> TIndex;
    TIndex index(seqs);

    Iterator<TIndex, Mums>::Type mumIt(index, 5);
    String<SAValue<TIndex>::Type> occs;

    while (!atEnd(mumIt))
    {
        occs = getOccurrences(mumIt);
        for (unsigned i = 0; i < length(occs); ++i)
        {
            for (unsigned j = i + 1; j < length(occs); ++j)
            {
                TMatch m(getValueI1(occs[i]), getValueI2(occs[i]), getValueI1(occs[j]), getValueI2(occs[j]), repLength(mumIt));
                appendValue(matches, m);
            }
        }
        ++mumIt;
    }

    typedef StringSet<TSequence, Dependent<> > TDepSequenceSet;
    typedef Graph<Alignment<TDepSequenceSet> > TAlignmentGraph;
    TAlignmentGraph g(seqs);
    matchRefinement(matches, seqs, g);
    clear(matches);

    tripletLibraryExtension(g);

    typedef String<double> TDistanceMatrix;
    TDistanceMatrix distanceMatrix;
    getDistanceMatrix(g, distanceMatrix, KmerDistance());

    typedef Graph<Tree<double> > TGuideTree;
    TGuideTree  guideTree;
    upgmaTree(distanceMatrix, guideTree);

    TAlignmentGraph gOut(seqs);
    progressiveAlignment(g, guideTree, gOut);
    std::cout << gOut;

    return 0;
}
