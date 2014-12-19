//![header]
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/seeds.h>

using namespace seqan;

int main()
{
//![header]
//![example]
    // Default-construct seed.
    Seed<Simple> seed1;
    std::cout << "seed1\n"
              << "beginPositionH == " << beginPositionH(seed1) << "\n"
              << "endPositionH == " << endPositionH(seed1) << "\n"
              << "beginPositionV == " << beginPositionV(seed1) << "\n"
              << "endPositionV == " << endPositionV(seed1) << "\n"
              << "lowerDiagonal == " << lowerDiagonal(seed1) << "\n"
              << "upperDiagonal == " << upperDiagonal(seed1) << "\n\n";

    // Construct seed with begin and end position in both sequences.
    Seed<Simple> seed2(3, 10, 7, 14);
    setUpperDiagonal(seed2, -7);
    setLowerDiagonal(seed2, -9);
    std::cout << "seed2\n"
              << "beginPositionH == " << beginPositionH(seed2) << "\n"
              << "endPositionH == " << endPositionH(seed2) << "\n"
              << "beginPositionV == " << beginPositionV(seed2) << "\n"
              << "endPositionV == " << endPositionV(seed2) << "\n"
              << "lowerDiagonal == " << lowerDiagonal(seed2) << "\n"
              << "upperDiagonal == " << upperDiagonal(seed2) << "\n\n";
//![example]

//![footer]
    return 0;
}
//![footer]
