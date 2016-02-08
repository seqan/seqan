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

    // Get the corresponding ids
    unsigned id0 = positionToId(stringSet, 0);
    unsigned id1 = positionToId(stringSet, 1);
    std::cout << "//![main]" << '\n';
//![main]
    // Let's create a string set of type dependent to represent strings,
    // which are stored in the StringSet of type Owner
    StringSet<DnaString, Dependent<Tight> > depSet;
    // We assign the first two strings of the owner string set to the dependent StringSet,
    // but in a reverse order
    assignValueById(depSet, stringSet, id1);
    assignValueById(depSet, stringSet, id0);

    std::cout << "Dependent: " << '\n';
    // (1) Access by position
    std::cout << "Pos 0: " << value(depSet, 0) << '\n';
    // (2) Access by id
    std::cout << "Id 0:  " << valueById(depSet, id0) << '\n';
//![main]
    std::cout << "//![main]" << '\n';

    std::cout << "//![difference]" << '\n';
//![difference]
    std::cout << "Position 0: Id " << positionToId(depSet, 0) << '\n';
    std::cout << "Position 1: Id " << positionToId(depSet, 1) << '\n';
//![difference]
    std::cout << "//![difference]" << '\n';

    return 0;
}
