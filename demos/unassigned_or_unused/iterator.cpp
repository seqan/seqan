///A tutorial about the use of iterators.
#include <iostream>
#include <seqan/sequence.h>
#include <seqan/stream.h>

using namespace seqan;

int main()
{
///The metafunction @Metafunction.Iterator@ returns the iterator type for a given container type.
    String<char> str = "admn";
    Iterator<String<char> >::Type it = begin(str);
    Iterator<String<char> >::Type itEnd = end(str);
///We can use iterators to iterate over the elements of a container.
    while (it != itEnd)
    {
        std::cout << *it;
        ++it;
    }
    std::cout << std::endl;
///Rooted iterators know their container (@Concept.RootedIteratorConcept|Rooted Iterator@).
///Hence, the functions @Function.goBegin@ and @Function.atEnd@ do
///not get $str$ as an argument.
///The following loop increments each character in $str$.
    Iterator<String<char>, Rooted>::Type it2 = begin(str);
    for (goBegin(it2); !atEnd(it2); goNext(it2))
    {
        ++value(it2);
    }
///Some iterators support an iteration in reverse order.
///Note that @Function.goPrevious@ is called before the value of $it2$ is accessed.
///Remember that the end position of a container is always the position behind the last item in the container.
    goEnd(it2);
    while (!atBegin(it2))
    {
        goPrevious(it2);
        std::cout << getValue(it2);
    }
    std::cout << std::endl;
///@Function.assignValue@ can be used to change the value of an iterator.
    assignValue(begin(str), 'X');
    std::cout << str << std::endl;

    return 0;
}
