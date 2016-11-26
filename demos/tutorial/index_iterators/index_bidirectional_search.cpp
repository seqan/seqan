#include <iostream>
#include <seqan/index.h>

using namespace seqan;

int main()
{
//![Create]
String<Dna> text = "ACTTTGACAGCT";
typedef FastFMIndexConfig<void, uint32_t> TFastConfig;
typedef Index<String<Dna>, BidirectionalIndex<FMIndex<void, TFastConfig> > > TIndex;
TIndex index(text);
Iter<TIndex, VSTree<TopDown<ParentLinks<> > > > iter(index);
//![Create]

//![Search]
goDown(iter, DnaString("TTTC"), Fwd()); // search CTTT in the prefix trie
goDown(iter, Dna('G'), Rev()); // extend to CTTTG
goUp(iter);

std::cout << representative(iter, Fwd()) << std::endl;
std::cout << representative(iter, Rev()) << std::endl;
//![Search]

//![output]
// if we get here the pattern was found
// output match positions
for (unsigned i = 0; i < length(getOccurrences(iter, Fwd())); ++i)
    std::cout << getOccurrences(iter, Fwd())[i] << std::endl;

for (unsigned i = 0; i < length(getOccurrences(iter, Rev())); ++i)
    std::cout << getOccurrences(iter, Rev())[i] << std::endl;
//![output]

return 0;
}
