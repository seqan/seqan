///A tutorial about the use of a user-defined modifier.
#include <iostream>
#include <seqan/stream.h>
#include <seqan/modifier.h>

using namespace seqan;

///A user-defined modifier that transforms all characters to upper case.
struct MyFunctor :
    public std::unary_function<char, char>
{
    inline char operator()(char x) const
    {
        if (('a' <= x) && (x <= 'z')) return x + ('A' - 'a');

        return x;
    }

};


int main()
{
///The modifier is applied to a string.
    String<char> myString = "A man, a plan, a canal-Panama";
    ModifiedString<String<char>, ModView<MyFunctor> > myModifier(myString);

    std::cout << myString << std::endl;
    std::cout << myModifier << std::endl;
    replace(myString, 9, 9, "master ");
    std::cout << myString << std::endl;
    std::cout << myModifier << std::endl;

    return 0;
}
