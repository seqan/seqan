#include <seqan/basic.h>
#include <seqan/sequence.h>

using namespace seqan;

int main()
{
//![iterators]
    String<char> str = "ACME";
    Iterator<String<char> >::Type it1 = begin(str);                   // a standard iterator
    Iterator<String<char>, Standard>::Type it2 = begin(str);          // same as above
    Iterator<String<char>, Rooted>::Type it3 = begin(str);            // a rooted iterator
    Iterator<String<char>, Rooted>::Type it4 = begin(str, Rooted());  // same as above
//![iterators]

    ignoreUnusedVariableWarning(it1);
    ignoreUnusedVariableWarning(it2);
    ignoreUnusedVariableWarning(it3);
    ignoreUnusedVariableWarning(it4);

    return 0;
}
