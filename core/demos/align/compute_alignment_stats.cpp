#include <iostream>

#include <seqan/align.h>
#include <seqan/sequence.h>

int main()
{
    // Create an alignment between subject and query.
    seqan::Peptide subject =
            "MGLSDGEWQLVLNVWGKVEADIPGHGQEVLIRLFKGHPETLEKFDKFKHLKSEDEMKASE"
            "DLKKHGATVLTALGGILKKKGHHEAEIKPLAQSHATKHKIPVKYLEFISECIIQVLQSKH"
            "PGDFGADAQGAMNKALELFRKDMASNYK";
    seqan::Peptide query =
            "MSLTKTERTIIVSMWAKISTQADTIGTETLERLFLSHPQTKTYFPHFDLHPGSA"
            "QLRAHGSKVVAAVGDAVKSIDDIGGALSKLSELHAYILRVDPVNFKLLSHCLLVTLAARF"
            "PADFTAEAHAAWDKFLSVTEKYR";

    seqan::Align<seqan::Peptide> align;
    resize(rows(align), 2);
    setSource(row(align, 0), subject);
    setSource(row(align, 1), query);

    seqan::Blosum62 scoringScheme(-1, -12);
    globalAlignment(align, scoringScheme);

    // Compute the statistics of the alignment.
    seqan::AlignmentStats stats;
    int scoreVal = computeAlignmentStats(stats, align, scoringScheme);
    SEQAN_ASSERT_EQ(scoreVal, stats.alignmentScore);
    std::cout << align
              << "gap opens:           " << stats.numGapOpens << "\n"
              << "gap extensions:      " << stats.numGapExtensions << "\n"
              << "num matches:         " << stats.numMatches << "\n"
              << "num mismatches:      " << stats.numMismatches << "\n"
              << "num positive scores: " << stats.numPositiveScores << "\n"
              << "num negative scores: " << stats.numNegativeScores << "\n\n\n";

    // Clip alignment rows and compute score of this view.
    setClippedEndPosition(row(align, 0), 100);
    setClippedEndPosition(row(align, 1), 100);
    setClippedBeginPosition(row(align, 0), 5);
    setClippedBeginPosition(row(align, 1), 5);

    scoreVal = computeAlignmentStats(stats, align, scoringScheme);
    SEQAN_ASSERT_EQ(scoreVal, stats.alignmentScore);
    std::cout << "Clipping alignment to (5, 100)\n"
              << align
              << "gap opens:           " << stats.numGapOpens << "\n"
              << "gap extensions:      " << stats.numGapExtensions << "\n"
              << "num matches:         " << stats.numMatches << "\n"
              << "num mismatches:      " << stats.numMismatches << "\n"
              << "num positive scores: " << stats.numPositiveScores << "\n"
              << "num negative scores: " << stats.numNegativeScores << "\n";

    return 0;
}
