#include <seqan/basic.h>
#include <seqan/sequence.h>

using namespace seqan;

//![findPatternInReferenceHint]
template <typename TString, typename TPattern>
void findPatternInReference(String<int> & hits,
                            TString const & reference,
                            TPattern const & pattern)
{
    // [A] Check whether pattern fits into the sequence.

    // [B] Iterate over all positions at which the pattern might occur.

    // [C] Evaluate all positions of the pattern until you find a mismatch or you have found a hit.

    // [D] Report begin position at which pattern matches the sequence.
}
//![findPatternInReferenceHint]

//![main]
int main()
{
    return 0;
}
//![main]
