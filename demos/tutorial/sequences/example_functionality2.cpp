#include <seqan/sequence.h>
#include <seqan/seq_io.h>

using namespace seqan;

int main()
{
//![main]
    String<Dna5String> readList;
    resize(readList, 2);
    readList[0] = "GGTTTCGACG";
    readList[1] = "AAGATGTCGC";
    appendValue(readList, "TATGCATGAT");
//![main]

//![print]
    std::cout << length(readList) << std::endl;
    std::cout << length(readList[0]) << std::endl;
//![print]

//![clear]
    clear(readList);
//![clear]
}
