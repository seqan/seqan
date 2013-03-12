// FRAGMENT(header)
#include <seqan/align.h>
#include <seqan/file.h>
#include <seqan/score.h>
#include <seqan/seeds.h>
#include <seqan/sequence.h>

int main()
{
    // FRAGMENT(example)
    // The horizontal and vertical sequence (database and query).
    seqan::CharString seqH = "The quick BROWN fox jumped again!";
    seqan::CharString seqV =     "thick BROWN boxes of brownies!";
                                     //  ^^^
    // Create seed and print the seeed sequence.
    seqan::Seed<seqan::Simple> seed(11, 7, 14, 10);

    // Perform match extension.
    seqan::Score<int, seqan::Simple> scoringScheme(1, -1, -1);
    extendSeed(seed, seqH, seqV, seqan::EXTEND_BOTH, scoringScheme, 3,
               seqan::UnGappedXDrop());

    // Perform a banded alignment.
    seqan::Align<seqan::CharString> align;
    resize(rows(align), 2);
    assignSource(row(align, 0), infix(seqH, beginPositionH(seed),
                                      endPositionH(seed)));
    assignSource(row(align, 1), infix(seqV, beginPositionV(seed),
                                      endPositionV(seed)));

    // TODO(holtgrew): Use seed diagonals as bands.
    globalAlignment(align, scoringScheme);
    std::cerr << "Resulting alignment\n" << align << "\n";

    // FRAGMENT(footer)
    return 0;
}
