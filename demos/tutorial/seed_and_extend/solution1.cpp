#include <seqan/stream.h>
#include <seqan/seeds.h>
#include <seqan/sequence.h>

using namespace seqan;

int main()
{
    // Default-construct seed.
    Seed<Simple> seed1;

    // Construct seed with begin and end position in both sequences.
    Seed<Simple> seed2(3, 10, 7, 14);
    setUpperDiagonal(seed2, -7);
    setLowerDiagonal(seed2, -9);

    // Update seed1.
    setBeginPositionH(seed1, 2 * beginPositionH(seed2));
    setEndPositionH(seed1, 2 * endPositionH(seed2));
    setBeginPositionV(seed1, 2 * beginPositionV(seed2));
    setEndPositionV(seed1, 2 * endPositionV(seed2));
    setLowerDiagonal(seed1, 2 * lowerDiagonal(seed2));
    setUpperDiagonal(seed1, 2 * upperDiagonal(seed2));

    // Print resulting seed1.
    std::cout << "seed1\n"
              << "beginPositionH == " << beginPositionH(seed1) << "\n"
              << "endPositionH == " << endPositionH(seed1) << "\n"
              << "beginPositionV == " << beginPositionV(seed1) << "\n"
              << "endPositionV == " << endPositionV(seed1) << "\n"
              << "lowerDiagonal == " << lowerDiagonal(seed1) << "\n"
              << "upperDiagonal == " << upperDiagonal(seed1) << "\n\n";

    return 0;
}
