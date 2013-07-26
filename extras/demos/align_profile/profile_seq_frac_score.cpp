#include <iostream>
#include <seqan/align_profile.h>

int main()
{
    typedef seqan::ProfileChar<seqan::Dna, int> TDnaProfile;
    typedef seqan::String<TDnaProfile> TProfileString;

    TProfileString profile = "CGAT";
    seqan::DnaString seq = "CGGAAT";

    seqan::Gaps<TProfileString> gapsH(profile);
    seqan::Gaps<seqan::DnaString> gapsV(seq);

    seqan::Score<int, seqan::ProfileSeqFracScore> sScheme(profile);

    int val = globalAlignment(gapsH, gapsV, sScheme, seqan::NeedlemanWunsch());
    std::cout << "score value = " << val << "\n";

    std::cout << "gaps in profile/sequence\n"
              << "pos\tG\tS\n";
    for (unsigned i = 0; i < length(gapsH); ++i)
        std::cerr << i << "\t" << isGap(gapsH, i) << "\t" << isGap(gapsV, i) << "\n";

    return 0;
}
