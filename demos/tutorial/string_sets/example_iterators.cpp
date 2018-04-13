#include <seqan/sequence.h>
#include <seqan/stream.h>

using namespace seqan;

int main()
{
    StringSet<DnaString> stringSet;
    DnaString str0 = "TATA";
    DnaString str1 = "CGCG";
    appendValue(stringSet, str0);
    appendValue(stringSet, str1);

    std::cout << "//![simple_example]" << '\n';
//![simple_example]
    typedef Iterator<StringSet<DnaString> >::Type TStringSetIterator;
    for (TStringSetIterator it = begin(stringSet); it != end(stringSet); ++it)
    {
        std::cout << *it << '\n';
    }
//![simple_example]
    std::cout << "//![simple_example]" << '\n';

    std::cout << "//![concatenator]" << '\n';
//![concatenator]
    typedef Concatenator<StringSet<DnaString> >::Type TConcat;
    TConcat concatSet = concat(stringSet);

    Iterator<TConcat>::Type it = begin(concatSet);
    Iterator<TConcat>::Type itEnd = end(concatSet);
    for (; it != itEnd; goNext(it))
    {
        std::cout << *it << " ";
    }
    std::cout << '\n';
//![concatenator]
    std::cout << "//![concatenator]" << '\n';

    return 0;
}
