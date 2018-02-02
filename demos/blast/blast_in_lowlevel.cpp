#include <iostream>
#include <seqan/basic.h>
#ifndef STDLIB_VS
#include <seqan/blast.h>

using namespace seqan;

int main()
{
    std::string inPath = getAbsolutePath("/tests/blast/plus_comments_defaults.m9");

    std::ifstream fin(toCString(inPath), std::ios_base::in | std::ios_base::binary);
    auto fit = directionIterator(fin, Input());

    typedef std::pair<std::string, std::string> THsp;
    std::vector<THsp> hsps;

    while (!atEnd(fit))
    {
        // skip any comment lines
        if (!onMatch(fit, BlastTabularLL()))
        {
            skipUntilMatch(fit, BlastTabularLL());
            if (atEnd(fit))
                break;
        }

        // resize output list
        resize(hsps, length(hsps)+1);

        // read only the first two fields into our variables
        readMatch(fit, BlastTabularLL(), back(hsps).first, back(hsps).second);
    }

    std::sort(std::begin(hsps), std::end(hsps));
    std::unique(std::begin(hsps), std::end(hsps));

    for (THsp const & hsp : hsps)
        std::cout << '(' << hsp.first << ", " << hsp.second << ")\n";

    return 0;
}
#else
int main()
{
    std::cerr << "Demo not run, because of a bug in Microsoft Visual Studio 2015.\n";
    return 0;
}
#endif
