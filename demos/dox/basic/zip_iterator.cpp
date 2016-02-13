#include <iostream>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>

using namespace seqan;

int main()
{
    StringSet<CharString> names;
    StringSet<CharString> surnames;
    String<int>           ages;

    appendValue(names, "John");
    appendValue(names, "Johnny");
    appendValue(names, "Garfield");

    appendValue(surnames, "Doe");
    appendValue(surnames, "Donny");
    appendValue(surnames, "the Cat");

    appendValue(ages, 98);
    appendValue(ages, 20);
    appendValue(ages, 42);

    auto it = makeZipIterator(begin(names), begin(surnames), begin(ages));
    auto itEnd = makeZipIterator(end(names), end(surnames), end(ages));

    for (; it != itEnd; ++it)
        std::cout << std::get<1>(*it) << ", " << std::get<0>(*it) << ": " << std::get<2>(*it) << std::endl;

    return 0;
}
