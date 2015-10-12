//![main]
#include <iostream>
#include <seqan/align.h>

using namespace seqan;

int main()
{
    typedef String<char> TSequence;

    TSequence text =    "MISSISSIPPIANDMISSOURI";
    TSequence pattern = "SISSI";

    String<int> locations;
    for (unsigned i = 0; i < length(text) - length(pattern); ++i)
    {
        // Compute the MyersBitVector in current window of text.
        TSequence tmp = infix(text, i, i + length(pattern));

        // Report hits with at most 2 errors.
        if (globalAlignmentScore(tmp, pattern, MyersBitVector()) >= -2)
        {
            appendValue(locations, i);
        }
    }
    return 0;
}
//![main]
