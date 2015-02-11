#include <seqan/sequence.h>
#include <seqan/index.h>

using namespace seqan;

int main()
{
    // One possible solution to the first sub assignment
    String<Dna> text = "ACGTTTGACAGCT";
    Index<String<Dna>, IndexEsa<> > index(text);

    // One possible solution to the second sub assignment
    StringSet<String<Dna> > stringSet;
    appendValue(stringSet, "ACGTCATCAT");
    appendValue(stringSet, "ACTTTG");
    appendValue(stringSet, "CACCCCCCTATTT");

    Index<StringSet<String<Dna> >, IndexEsa<> > indexSet(stringSet);

    return 0;
}
