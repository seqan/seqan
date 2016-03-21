//![main]
#include <iostream>
#include <seqan/stream.h>
#include <seqan/modifier.h>

using namespace seqan;

int main()
{
    String<char> myString = "A man, a plan, a canal-Panama";
//![main]
//![modifier]
    ModifiedString<String<char>, ModReverse> myModifier(myString);
//![modifier]

    std::cout << "//![output1]\n";
//![output1]
    std::cout << myString << std::endl;
    std::cout << myModifier << std::endl;
//![output1]
    std::cout << "//![output1]\n";
    std::cout << "//![output2]\n";
//![output2]
    replace(myString, 9, 9, "master ");
    std::cout << myString << std::endl;
    std::cout << myModifier << std::endl;
    
//![output2]
    std::cout << "//![output2]\n";
//![output2]
    return 0;
}
//![output2]
