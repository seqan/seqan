#include <seqan/sequence.h>
#include <seqan/index.h>

using namespace std;
using namespace seqan;

//![FMIndexConfigs]
typedef FMIndexConfig<void, uint32_t>     TConfig;
typedef FastFMIndexConfig<void, uint32_t> TFastConfig;

Index<String<Dna>, FMIndex<void, TFastConfig> > index(genome);
//![FMIndexConfigs]

//![FMIndexConfigs2]
typedef FastFMIndexConfig<void, uint32_t, 2, 1> TMyFastConfig;
//![FMIndexConfigs2]


int main()
{
//![BidirectionalIndex]
String<Dna> text = "ACTTTGACAGCT";
typedef Index<String<Dna>, BidirectionalIndex<FMIndex<void, TMyFastConfig> > > TIndex;
TIndex index(text);
//![BidirectionalIndex]

return 0;
}
