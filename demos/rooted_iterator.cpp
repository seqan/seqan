///A tutorial about the use of rooted iterators.
#include <iostream>
#include <seqan/sequence.h>
#include <seqan/stream.h>

using namespace seqan;

int main()
{
    typedef String<char> TText;
    TText str = "abcdefg";
///The rooted iterator can be specified by passing the @Tag.Iterator Spec.iterator spec Rooted@ as a second argument.
    Iterator<TText, Rooted>::Type it = begin(str);
///The same iterator spec can be used as a tag for functions that return an iterator, e.g. @Function.begin@ or @Function.end@.
    it = begin(str, Rooted());
///A rooted iterator "knows" its container, so it supports the function @Function.container@.
    std::cout << container(it);          //output: "abcdefg"
    goNext(it);
///It also supports functions such as @Function.goBegin@ or @Function.position@.
    std::cout << position(it);           //output: 7
    return 0;
}
