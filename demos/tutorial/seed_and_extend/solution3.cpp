#include <seqan/align.h>
#include <seqan/stream.h>
#include <seqan/score.h>
#include <seqan/seeds.h>
#include <seqan/sequence.h>

using namespace seqan2;

int main()
{
    // The horizontal and vertical sequence (subject and query sequences).
    CharString seqH = "The quick BROWN fox jumped again!";
    CharString seqV =     "thick BROWN boxes of brownies!";
    // Create the seed sequence.
    Seed<Simple> seed(11, 7, 14, 10);

    // Perform match extension.
    Score<int, Simple> scoringScheme(1, -1, -2, -2);
    extendSeed(seed, seqH, seqV, EXTEND_BOTH, scoringScheme, 3,
               GappedXDrop());

    // Perform a banded alignment.
    Align<CharString> align;
    resize(rows(align), 2);
    assignSource(row(align, 0), infix(seqH, beginPositionH(seed),
                                      endPositionH(seed)));
    assignSource(row(align, 1), infix(seqV, beginPositionV(seed),
                                      endPositionV(seed)));

    globalAlignment(align, scoringScheme);
    std::cout << "Resulting alignment\n" << align << "\n";

    return 0;
}
