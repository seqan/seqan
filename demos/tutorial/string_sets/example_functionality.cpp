//![appendValue]
#include <seqan/sequence.h>
#include <seqan/stream.h>

using namespace seqan;

int main()
{
    StringSet<DnaString> stringSet;
    DnaString str0 = "TATA";
    DnaString str1 = "CGCG";
    appendValue(stringSet, str0);
    appendValue(stringSet, str1);
//![appendValue]
//![retrieve_id]
    // (1) Access by position
    std::cout << "Owner: " << '\n';
    std::cout << "Position 0: " << value(stringSet, 0) << '\n';

    // Get the corresponding ids
    unsigned id0 = positionToId(stringSet, 0);
    unsigned id1 = positionToId(stringSet, 1);

    // (2) Access by id
    std::cout << "Id 0:  " << valueById(stringSet, id0) << '\n';
    std::cout << "Id 1:  " << valueById(stringSet, id1) << '\n';

    return 0;
}
//![retrieve_id]
