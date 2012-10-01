#include <iostream>

#include <seqan/basic.h>
#include <seqan/file.h>

int main(int argc, char const ** argv)
{
    if (argc != 2)
        return 1;
    
    // Memory mapped string, automatically backed by temporary file.
    seqan::String<char, seqan::MMap<> > str1;
    str1 = "This is the first mapped string!";
    std::cout << str1 << std::endl;

    // Open file as memory mapped string.
    seqan::String<char, seqan::MMap<> > str2;
    if (!open(str2, argv[1], seqan::OPEN_RDONLY))
    {
        std::cerr << "Could not open file " << argv[1] << std::endl;
        return 1;
    }
    std::cout << str2 << std::endl;

    return 0;
}
