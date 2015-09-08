///A tutorial about the use of the reverse modifier.
#include <iostream>
#include <seqan/stream.h>
#include <seqan/modifier.h>

using namespace seqan;


int main()
{
///A reverse modifier applied to a string.
    String<char> myString = "A man, a plan, a canal-Panama";
    ModifiedString<String<char>, ModReverse> myModifier(myString);

    std::cout << myString << std::endl;
    std::cout << myModifier << std::endl;
    replace(myString, 9, 9, "master ");
    std::cout << myString << std::endl;
    std::cout << myModifier << std::endl;
    return 0;
}
