// FRAGMENT(includes)
#include <iostream>
#include <seqan/store.h>
#include <seqan/misc/misc_svg.h>

using namespace seqan;

int main ()
{
    FragmentStore<> store;
    loadContigs(store, "ex1.fa");
    BamFileIn file("ex1.sam");
    readRecords(store, file);

// FRAGMENT(ascii)
    AlignedReadLayout layout;
    layoutAlignment(layout, store);
    printAlignment(std::cout, layout, store, 1, 0, 150, 0, 36);

// FRAGMENT(svg)
    SVGFile svg("layout.svg");
    printAlignment(svg, layout, store, 1, 0, 150, 0, 36);

    return 0;
}
