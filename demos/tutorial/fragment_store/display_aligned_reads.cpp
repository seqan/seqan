//![includes]
#include <iostream>
#include <seqan/store.h>
#include <seqan/misc/svg.h>

using namespace seqan;

int main()
{
    CharString fastaFileName = getAbsolutePath("demos/tutorial/fragment_store/example.fa");
    CharString samFileName = getAbsolutePath("demos/tutorial/fragment_store/example.sam");

    typedef FragmentStore<> TStore;

    TStore store;
    loadContigs(store, toCString(fastaFileName));
    BamFileIn file(toCString(samFileName));
    readRecords(store, file);

//![includes]

//![ascii]
    AlignedReadLayout layout;
    layoutAlignment(layout, store);
    printAlignment(std::cout, layout, store, 1, 0, 150, 0, 36);
//![ascii]

//![svg]
    SVGFile svg("layout.svg");
    printAlignment(svg, layout, store, 1, 0, 150, 0, 36);

    return 0;
}
//![svg]
