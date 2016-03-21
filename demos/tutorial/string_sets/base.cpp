#include <seqan/sequence.h>

using namespace seqan;

int main()
{
//![example_sets]
    StringSet<DnaString>               ownerSet;
    StringSet<DnaString, Owner<> >     ownerSet2;      // same as above
    StringSet<DnaString, Dependent<> > dependentSet;
//![example_sets]

//![appendValue]
    StringSet<DnaString> stringSet;
    DnaString str0 = "TATA";
    DnaString str1 = "CGCG";
    appendValue(stringSet, str0);
    appendValue(stringSet, str1);
//![appendValue]
    return 0;
}
