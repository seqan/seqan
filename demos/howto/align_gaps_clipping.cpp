// Demo program for clipping with Gaps objects.

#include <iostream>

#include <seqan/sequence.h>
#include <seqan/align.h>

using namespace seqan;

int main()
{
    // Create sequence variable and gaps basd on sequence.
    CharString seq("ABCDEFGHIJ");
    Gaps<CharString> gaps(seq);

    // Insert gaps, the positions are in (clipped) view space.
    insertGaps(gaps, 2, 2);
    insertGap(gaps, 6);
    insertGap(gaps, 10);

    // Print to stdout.
    std::cout << "gaps\t" << gaps << "\n"
              << "seq \t" << seq << "\n\n";

    // Print the begin and end positions in sequence space and the clipped
    // begin and end positions in gap space.  We have no clipping, so no
    // surprises here.
    std::cout << "beginPosition(gaps)        == " << beginPosition(gaps) << "\n"
              << "endPosition(gaps)          == " << endPosition(gaps) << "\n"
              << "clippedBeginPosition(gaps) == " << clippedBeginPosition(gaps) << "\n"
              << "clippedEndPosition(gaps)   == " << clippedEndPosition(gaps) << "\n\n";

    // Now, clip the alignment and again print the gaps, sequence and begin/end
    // positions.  Note that the clipping positions are relative to the unclipped
    // view.
    setClippedBeginPosition(gaps, 3);
    setClippedEndPosition(gaps, 10);

    std::cout << "gaps\t" << gaps << "\n"
              << "seq \t" << infix(seq, beginPosition(gaps), endPosition(gaps)) << "\n\n";

    std::cout << "beginPosition(gaps)        == " << beginPosition(gaps) << "\n"
              << "endPosition(gaps)          == " << endPosition(gaps) << "\n"
              << "clippedBeginPosition(gaps) == " << clippedBeginPosition(gaps) << "\n"
              << "clippedEndPosition(gaps)   == " << clippedEndPosition(gaps) << "\n\n";

    // We can translate between the (clipped) gapped position (aka view) and
    // the unclipped ungapped positions (aka) source using toSourcePosition()
    // and toViewPosition().  Note that because of projection to the right of
    // gaps, these operations are not symmetric.
    std::cout << "4 view position => " << toSourcePosition(gaps, 4) << " source position\n"
              << "2 source position => " << toViewPosition(gaps, 2) << " view position\n\n";

    // Translating between clipped gapped and unclipped gapped position can
    // be done by adding/subtracting clippedBeginPosition(gaps).
    std::cout << "3 clipped gapped => " << 3 + clippedBeginPosition(gaps) << " unclipped gapped\n"
              << "6 unclipped gapped => " << 5 - clippedBeginPosition(gaps) << " clipped gapped\n\n";

    // Translating between clipped ungapped and unclipped ungapped position can
    // be done by adding/subtracing beginPosition(gaps).  Since there are no
    // gaps, this operation is symmetric.
    std::cout << "3 clipped ungapped => " << 3 + beginPosition(gaps) << " unclipped ungapped\n"
              << "5 unclipped ungapped => " << 5 - beginPosition(gaps) << " clipped ungapped\n\n";

    // Translating between gapped clipped position and ungapped clipped
    // position and between gapped unclipped and ungapped unclipped positions
    // has to be done using the translations above.
    std::cout << "3 clipped gapped => " << toSourcePosition(gaps, 3) - beginPosition(gaps) << " clipped ungapped\n"
              << "4 unclipped ungapped => " << toViewPosition(gaps, 4) + clippedBeginPosition(gaps) << " unclipped gapped\n";

    return 0;
}
