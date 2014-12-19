//![includes]
#include <iostream>
#include <seqan/index.h>

using namespace seqan;

//![includes]
//![initialization]
int main()
{
    typedef Index<DnaString, IndexQGram<UngappedShape<3> > > TIndex;
    TIndex index("CATGATTACATA");

//![initialization]
//![output]
    hash(indexShape(index), "CAT");
    for (unsigned i = 0; i < length(getOccurrences(index, indexShape(index))); ++i)
        std::cout << getOccurrences(index, indexShape(index))[i] << std::endl;

    return 0;
}
//![output]
