#include <iostream>
#include <seqan/sequence.h>
#include <seqan/stream.h>

#include <time.h>

using namespace seqan;

int main()
{
    unsigned num = 1000000;
    time_t start;

    String<Dna> str;
    clear(str);
    start = time(NULL);
    for (unsigned i = 0; i < num; ++i)
    {

        appendValue(str, 'A', Exact());
    }
    std::cout << "Strategy Exact() took: " << time(NULL) - start << " s\n\n";

    clear(str);
    start = time(NULL);
    for (unsigned i = 0; i < num; ++i)
    {

        appendValue(str, 'A', Generous());
    }
    std::cout << "Strategy Generous() took: " << time(NULL) - start << " s\n\n";

    return 0;
}
