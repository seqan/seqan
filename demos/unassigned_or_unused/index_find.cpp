///A tutorial about find algorithms using an index.
#include <iostream>
#include <seqan/index.h>

using namespace seqan;

int main()
{
///The following code creates an @Class.Index@ of $"tobeornottobe"$.
///As there is no second template parameter given to $Index<..>$,
///the default index based on an enhanced suffix array is used.
    Index<String<char> > index_esa("tobeornottobe");
    Finder<Index<String<char> > > finder_esa(index_esa);

    std::cout << "hit at ";
    while (find(finder_esa, "be"))
        std::cout << position(finder_esa) << " ";
    std::cout << std::endl;

///Now we explicitly create a q-gram index using an ungapped 2-gram and do the same.
///Instead of this @Spec.UngappedShape.fixed-size shape@ you can use arbitrary
///@Class.Shape.shapes@ and assign them before calling @Function.find@ via @Function.indexShape@.
    typedef Index<String<char>, IndexQGram<UngappedShape<2> > > TQGramIndex;
    TQGramIndex index_2gram("tobeornottobe");
    Finder<TQGramIndex> finder_2gram(index_2gram);

    std::cout << "hit at ";
    while (find(finder_2gram, "be"))
        std::cout << position(finder_2gram) << " ";
    std::cout << std::endl;

    return 0;
}
