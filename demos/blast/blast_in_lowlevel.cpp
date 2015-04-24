#include <iostream>
#include <seqan/blast.h>

using namespace seqan;

int main()
{
    std::string inPath = std::string(SEQAN_PATH_TO_ROOT()) + "/tests/blast/plus_header_defaults.blast";

    std::ifstream fin(toCString(inPath), std::ios_base::in | std::ios_base::binary);
    auto fit = directionIterator(fin, Input());

    typedef std::pair<std::string, std::string> THsp;
    std::vector<THsp> hsps;

    while (!atEnd(fit))
    {
        // skip any headers
        if (!onMatch(fit, BlastTabular()))
        {
            skipUntilMatch(fit, BlastTabular());
            if (atEnd(fit))
                break;
        }

        // resize output list
        resize(hsps, length(hsps)+1);

        // read only the first two fields into our variables
        readMatch0(fit, BlastTabular(), back(hsps).first, back(hsps).second);
    }

    std::sort(std::begin(hsps), std::end(hsps));
    std::unique(std::begin(hsps), std::end(hsps));

    for (THsp const & hsp : hsps)
        std::cout << '(' << hsp.first << ", " << hsp.second << ")\n";

    return 0;
}