#include <iostream>
#include <seqan/sequence.h>
#include <seqan/stream.h>

using namespace seqan;

int main()
{
    Dna5String genome = "ANTGGTTNCAACNGTAANTGCTGANNNACATGTNCGCGTGTA";

    Iterator<Dna5String, Rooted>::Type it = begin(genome);

    for (; !atEnd(it); goNext(it))
    {
        if (getValue(it) == 'N')
            value(it) = 'A';
    }
    std::cout << "Modified genome: " << genome << std::endl;
    return 0;
}
