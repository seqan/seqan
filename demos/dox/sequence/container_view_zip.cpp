#include <iostream>

#include <seqan/sequence.h>
#include <seqan/stream.h>  // for output

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

    auto pRegister = makeZipView(names, surnames, ages);

    for (auto person : pRegister)
        std::cout << std::get<1>(person) << ", " << std::get<0>(person) << ": " << std::get<2>(person) << std::endl;
    return 0;
}
