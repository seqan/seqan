#include <iostream>

#include <seqan/align.h>
#include <seqan/sequence.h>

int main()
{
    using namespace seqan;

    // Create an alignment between subject and query.
    seqan::Peptide subject =
            "MGLSDGEWQLVLNVWGKVEADIPGHGQEVLIRLFKGHPETLEKFDKFKHLKSEDEMKASE"
            "DLKKHGATVLTALGGILKKKGHHEAEIKPLAQSHATKHKIPVKYLEFISECIIQVLQSKH"
            "PGDFGADAQGAMNKALELFRKDMASNYK";
    seqan::Peptide query =
            "MSLTKTERTIIVSMWAKISTQADTIGTETLERLFLSHPQTKTYFPHFDLHPGSA"
            "QLRAHGSKVVAAVGDAVKSIDDIGGALSKLSELHAYILRVDPVNFKLLSHCLLVTLAARF"
            "PADFTAEAHAAWDKFLSVTEKYR";

    typedef seqan::Align<seqan::Peptide> TAlign;
    typedef typename seqan::Row<TAlign>::Type TRow;      
    
    seqan::Align<seqan::Peptide> align;
    resize(rows(align), 2);
    setSource(row(align, 0), subject);
    setSource(row(align, 1), query);

    seqan::Blosum62 scoringScheme(-1, -12);
    globalAlignment(align, scoringScheme);

    TRow queryRow = row(align, 0);
    TRow refRow = row(align, 1);

    // Compute the statistics of the alignment.
    std::cout << "Unclipped alignment (Default Protein Similarity)\n"
              << row(align, 0) << "\n"
              << alignmentSimilarityString(align, scoringScheme) << "\n"
              << row(align, 1) << "\n\n";

    std::cout << "Unclipped alignment (Custom Similarity)\n"
              << row(align, 0) << "\n"
              << alignmentSimilarityString(queryRow, refRow, scoringScheme, ':', '.', '*', '*') << "\n"
              << row(align, 1) << "\n\n";

    std::cout << "Unclipped alignment (Default Identity)\n"
              << row(align, 0) << "\n"
              << alignmentIdentityString(align) << "\n"
              << row(align, 1) << "\n\n";

    std::cout << "Unclipped alignment (Custom Identity)\n"
              << row(align, 0) << "\n"
              << alignmentIdentityString(queryRow, refRow, ':', '.', ' ') << "\n"
              << row(align, 1) << "\n\n";

    // Clip alignment rows and compute score of this view.
    setClippedEndPosition(row(align, 0), 100);
    setClippedEndPosition(row(align, 1), 100);
    setClippedBeginPosition(row(align, 0), 5);
    setClippedBeginPosition(row(align, 1), 5);

    TRow clipQueryRow = row(align, 0);
    TRow clipRefRow = row(align, 1);

    std::cout << "Clipped alignment (Default Protein Similarity)\n"
              << row(align, 0) << "\n"
              << alignmentSimilarityString(align, scoringScheme) << "\n"
              << row(align, 1) << "\n\n";

    std::cout << "Clipped alignment (Custom Similarity)\n"
              << row(align, 0) << "\n"
              << alignmentSimilarityString(clipQueryRow, clipRefRow, scoringScheme, ':', '.', '*', '*') << "\n"
              << row(align, 1) << "\n\n";

    std::cout << "Clipped alignment (Default Identity)\n"
              << row(align, 0) << "\n"
              << alignmentIdentityString(align) << "\n"
              << row(align, 1) << "\n\n";

    std::cout << "Clipped alignment (Custom Identity)\n"
              << row(align, 0) << "\n"
              << alignmentIdentityString(clipQueryRow, clipRefRow, ':', '.', ' ') << "\n"
              << row(align, 1) << "\n";

    return 0;
}
