// The comment lines containing ![fragment-line] are there for the
// documentation system.  You can ignore them when reading this file.i
// This is the draft for the new iterator tutorial
//![includes]
#include <iostream>
#include <seqan/sequence.h>
#include <seqan/stream.h>

using namespace seqan;

int main()
{
//![includes]
//![metafunctions]
    Dna5String genome = "TATANNNGCGCG";
    Iterator<Dna5String>::Type it = begin(genome);
    Iterator<Dna5String>::Type itEnd = end(genome);
//![metafunctions]
//![iterators]
    while (it != itEnd)
    {
        std::cout << *it;
        ++it;
    }
    std::cout << std::endl;
//![iterators]
//![standard-iterators]
    for (goBegin(it, genome); !atEnd(it, genome); goNext(it))
    {
        std::cout << *it;
    }
    std::cout << std::endl;
//![standard-iterators]
//![rooted-iterators]
    Iterator<Dna5String, Rooted>::Type it2 = begin(genome);
    for (goBegin(it2); !atEnd(it2); goNext(it2))
    {
        if (getValue(it2) == 'A')
            std::cout << 'T';
        else if (getValue(it2) == 'T')
            std::cout << 'A';
        else if (getValue(it2) == 'G')
            std::cout << 'C';
        else if (getValue(it2) == 'C')
            std::cout << 'G';
        else
            std::cout << 'N';
    }
    std::cout << std::endl;
//![rooted-iterators]
//![iterators-reverse]
    goEnd(it2);
    while (!atBegin(it2))
    {
        goPrevious(it2);
        std::cout << getValue(it2);
    }
    std::cout << std::endl;
//![iterators-reverse]
//![assign-value]
    assignValue(begin(genome), 'N');
    std::cout << genome << std::endl;

    return 0;
}
//![assign-value]
