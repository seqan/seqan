//![main]
#include <iostream>
#include <seqan/align.h>

using namespace seqan;

int main()
{
    typedef String<char> TSequence;
    typedef Gaps<TSequence, ArrayGaps> TGaps;
    typedef Iterator<String<int> >::Type TIterator;

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

    TGaps gapsText;
    TGaps gapsPattern;
    assignSource(gapsPattern, pattern);
    std::cout << "Text: " << text << "\tPattern: " << pattern << std::endl;
    for (TIterator it = begin(locations); it != end(locations); ++it)
    {
        // Clear previously computed gaps.
        clearGaps(gapsText);
        clearGaps(gapsPattern);

        // Only recompute the area within the current window over the text.
        TSequence textInfix = infix(text, *it, *it + length(pattern));
        assignSource(gapsText, textInfix);

        // Use semi-global alignment since we do not want to track leading/trailing gaps in the pattern.
        // Restirct search space using a band allowing at most 2 errors in vertical/horizontal direction.
        int score = globalAlignment(gapsText, gapsPattern, Score<int>(0, -1, -1), AlignConfig<true, false, false, true>(), -2, 2);
        std::cout << "Hit at position " << *it << "\ttotal edits: " << abs(score) << std::endl;
    }

    return 0;
}
//![main]
