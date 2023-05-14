#include <iostream>
#include <seqan/sequence.h>
#include <seqan/stream.h>

using namespace seqan2;

int main()
{
    Dna5String genome = "ANTGGTTNCAACNGTAANTGCTGANNNACATGTNCGCGTGTA";
    for (unsigned i = 0; i < length(genome); ++i)
    {
        if (genome[i] == 'N')
            genome[i] = 'A';
    }
    std::cout << "Modified genome: " << genome << std::endl;
    return 0;
}
