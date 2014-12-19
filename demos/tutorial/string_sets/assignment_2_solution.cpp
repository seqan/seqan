#include <iostream>
#include <seqan/sequence.h>
#include <seqan/stream.h>

using namespace seqan;

// Check whether the string set contains the string with the given id,
// without comparing the actual sequences
template <typename TStringSet, typename TId>
bool isElement(TStringSet & stringSet1, TId & id)
{

    for (unsigned i = 0; i < length(stringSet1); ++i)
    {
        // Get the id of the element at position i
        if (positionToId(stringSet1, i) == id)
            return true;
    }
    return false;
}

int main()
{
    // Build strings
    DnaString str0 = "TATA";
    DnaString str1 = "CGCG";
    DnaString str2 = "TTAAGGCC";
    DnaString str3 = "ATGC";
    DnaString str4 = "AGTGTCA";
    // Build owner string set and append strings
    StringSet<DnaString> stringSetOw;
    appendValue(stringSetOw, str0);
    appendValue(stringSetOw, str1);
    appendValue(stringSetOw, str2);
    appendValue(stringSetOw, str3);
    appendValue(stringSetOw, str4);
    // Get corresponding ids for positions
    unsigned id0 = positionToId(stringSetOw, 0);
    unsigned id1 = positionToId(stringSetOw, 1);
    unsigned id2 = positionToId(stringSetOw, 2);
    unsigned id3 = positionToId(stringSetOw, 3);
    // Build dependent string set and assigns strings by id
    StringSet<DnaString, Dependent<Generous> > stringSetDep;
    assignValueById(stringSetDep, stringSetOw, id0);
    assignValueById(stringSetDep, stringSetOw, id1);
    assignValueById(stringSetDep, stringSetOw, id3);
    // Call function to check if a string is contained and output result
    std::cout << "Does the string set contain the string with the id 'id3'? " <<  isElement(stringSetDep, id3) << std::endl;
    std::cout << "Does the string set contain the string with the id 'id2'? " <<  isElement(stringSetDep, id2) << std::endl;

    return 0;
}
