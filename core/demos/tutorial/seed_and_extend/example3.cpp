// FRAGMENT(header)
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/score.h>
#include <seqan/seeds.h>

int main()
{
    // FRAGMENT(example)
    // The horizontal and vertical sequence (database and query).
    seqan::CharString seqH = "The quick BROWN fox jumped again!";
    seqan::CharString seqV =     "thick BROWN boxes of brownies!";
                                     //  ^^^
    // Create seed and print the seeed sequence.
    seqan::Seed<seqan::Simple> seed(11, 7, 14, 10);
    std::cout << "original\n"
              << "seedH: " << infix(seqH, beginPositionH(seed),
                                    endPositionH(seed)) << "\n"
              << "seedV: " << infix(seqV, beginPositionV(seed),
                                    endPositionV(seed)) << "\n";

    // Perform match extension.
    seqan::Score<int, seqan::Simple> scoringScheme(1, -1, -1);
    extendSeed(seed, seqH, seqV, seqan::EXTEND_BOTH, scoringScheme, 3,
               seqan::UnGappedXDrop());
    // Print the resulting seed.
    std::cout << "result\n"
              << "seedH: " << infix(seqH, beginPositionH(seed),
                                    endPositionH(seed)) << "\n"
              << "seedV: " << infix(seqV, beginPositionV(seed),
                                    endPositionV(seed)) << "\n";

    // FRAGMENT(footer)
    return 0;
}
