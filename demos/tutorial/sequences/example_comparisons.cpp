#include <seqan/sequence.h>

using namespace seqan;

int main()
{
//![main]
    String<char> a = "beta";
    String<char> b = "alpha";

    std::cout << (a != b) << std::endl;
    std::cout << (a < b) << std::endl;
    std::cout << (a > b) << std::endl;
//![main]

//![first]
    if (a < b)      { /* code for case "a < b"  */ }
    else if (a > b) { /* code for case "a > b"  */ }
    else            { /* code for case "a == b" */ }
//![first]

//![second]
    // Compare a and b and store the result in comp
    Lexical<> comp(a, b);

    if (isLess(comp))         { /* code for case "a < b"  */ }
    else if (isGreater(comp)) { /* code for case "a > b"  */ }
    else                      { /* code for case "a == b" */ }
//![second]
}
