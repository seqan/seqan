///A tutorial about accessing the lcp table of an extended suffix array
#include <iostream>
#include <seqan/index.h>

using namespace seqan;

int main()
{
///We begin with a @Class.StringSet@ that stores multiple strings.
    StringSet<String<char> > mySet;
    resize(mySet, 3);
    mySet[0] = "SeqAn is a library for sequence analysis.";
    mySet[1] = "The String class is the fundamental sequence type in SeqAn.";
    mySet[2] = "Subsequences can be handled with SeqAn's Segment class.";

///Then we create an @Class.Index@ of this @Class.StringSet@.
    typedef Index<StringSet<String<char> > > TMyIndex;
    TMyIndex myIndex(mySet);

///Now we require that the index has a suffix array and a lcp table using the function @Function.indexRequire@.

    indexRequire(myIndex, EsaLcp());
    indexRequire(myIndex, EsaSA());

///We iterate over the lcp table and output the position, the value of the lcp table and the corresponding suffix using the functions @Function.saAt@ and @Function.IndexEsa#lcpAt@.

    for (Size<TMyIndex>::Type i = 0; i < length(myIndex); ++i)
    {
        SAValue<TMyIndex>::Type p = saAt(i, myIndex);
        std::cout << i << " " << lcpAt(i, myIndex) << " " << p << " " << suffix(mySet[getSeqNo(p, stringSetLimits(mySet))], getSeqOffset(p, stringSetLimits(mySet))) << std::endl;
    }

    return 0;
}
