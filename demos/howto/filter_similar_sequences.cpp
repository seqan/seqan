//![includes]
#include <seqan/stream.h>
#include <seqan/index.h>
#include <seqan/store.h>
#include <iostream>

using namespace seqan;
//![includes]

//![load_reads]
int main(int argc, char const * argv[])
{
    FragmentStore<> fragStore;
    if (argc < 2 || !loadReads(fragStore, argv[1]))
    {
        std::cerr << "ERROR: Coud not load reads." << std::endl;
        return 0;
    }
//![load_reads]

//![filter]
    typedef FragmentStore<>::TReadSeqStore TReadSeqStore;
    typedef GetValue<TReadSeqStore>::Type TReadSeq;
    typedef Index<TReadSeqStore, IndexQGram<Shape<Dna, UngappedShape<11> >, OpenAddressing> > TIndex;
    typedef Pattern<TIndex, Swift<SwiftSemiGlobal> > TPattern;
    typedef Finder<TReadSeq, Swift<SwiftSemiGlobal> > TFinder;

    TIndex index(fragStore.readSeqStore);
    TPattern pattern(index);
    for (unsigned i = 0; i < length(fragStore.readSeqStore); ++i)
    {
        if ((i % 1000) == 0)
            std::cout << "." << std::flush;
        TFinder finder(fragStore.readSeqStore[i]);
        while (find(finder, pattern, 0.1))
        {
            if (i == position(pattern).i1)
                continue;
            // do further alignment here
/*			std::cout << "Found possible overlap of " << std::endl;
            std::cout << "\t" << infix(finder) << std::endl;
            std::cout << "\t" << seqs[position(pattern).i1] << std::endl;
*/      }
    }

    return 0;
}
//![filter]
