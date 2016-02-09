///A tutorial about the use of a user-defined modifier.
#include <iostream>
#include <seqan/stream.h>
#include <seqan/modifier.h>

using namespace seqan;

//![functor]
struct MyFunctor :
    public std::unary_function<char, char>
{
    inline char operator()(char x) const
    {
        if (('a' <= x) && (x <= 'z'))
            return x + ('A' - 'a');

        return x;
    }

};
//![functor]

int main()
{
///The modifier is applied to a string.
//![mod_str]
    String<char> myString = "A man, a plan, a canal-Panama";
    ModifiedString<String<char>, ModView<MyFunctor> > myModifier(myString);
//![mod_str]
//![output]
    std::cout << myString << std::endl;
    std::cout << myModifier << std::endl;
//![output]

//![predefined]
    ModifiedString< String<char>, ModView<FunctorUpcase<char> > > myPredefinedModifier(myString);
//![predefined]
    return 0;
}
