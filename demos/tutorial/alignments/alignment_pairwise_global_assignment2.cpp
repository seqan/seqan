//![main]
#include <iostream>
#include <seqan/align.h>

using namespace seqan;

int main()
{
    typedef String<char>                TSequence;  // sequence type
    typedef StringSet<TSequence, Dependent<> >   TDepStringSet;
    typedef Graph<Alignment<TDepStringSet> >    TAlignGraph;
//![main]

//![init]
    StringSet<TSequence> seq;
    appendValue(seq, "blablablu");
    appendValue(seq, "abab");

    TAlignGraph alignG(seq);
//![init]

//![alignment]
    AlignConfig<true, false, false, true> ac;
    int score = globalAlignment(alignG, Score<int>(1, -1, -1, -1), ac, Gotoh());
    std::cout << "Score = " << score << std::endl;
    std::cout << alignG;

    return 0;
}
//![alignment]
