#include <seqan/sequence.h>
#include <seqan/index.h>

using namespace seqan;

int main()
{
//![FMIndexConfigs]
typedef FMIndexConfig<void, uint32_t>     TConfig;
typedef FastFMIndexConfig<void, uint32_t> TFastConfig;

String<Dna> genome = "ACTTTGACAGCT";
Index<String<Dna>, FMIndex<void, TConfig> > index(genome);
//![FMIndexConfigs]

Index<String<Dna>, FMIndex<void, TFastConfig> > index2(genome);

//![FMIndexConfigs2]
typedef FastFMIndexConfig<void, uint32_t, 2, 1> TMyFastConfig;
//![FMIndexConfigs2]

//![BidirectionalIndex]
String<Dna> text = "ACTTTGACAGCT";
typedef Index<String<Dna>, BidirectionalIndex<FMIndex<void, TMyFastConfig> > > TIndex;
TIndex bidirectionalIndex(text);
//![BidirectionalIndex]

return 0;
}
