#include <seqan/store.h>
#include <seqan/realign.h>

using namespace seqan;

int main()
{
    // Build paths to reference and to alignment.
    std::string refPath = getAbsolutePath("demos/tutorial/realignment/ref.fa");
    std::string samPath = getAbsolutePath("demos/tutorial/realignment/reads.sam");

    // Declare fragment store.
    FragmentStore<> store;

    // Load contigs and read alignment.
    loadContigs(store, refPath.c_str());
    BamFileIn bamFileIn(samPath.c_str());
    readRecords(store, bamFileIn);

    // Perform the realignment.
    reAlignment(store, /*contigID=*/ 0, /*method=*/ 1, /*bandwidth=*/ 20,
                /*includeReference=*/ true);

    // Layout alignment and print.
    AlignedReadLayout layout;
    layoutAlignment(layout, store);
    printAlignment(std::cout, layout, store, /*contigID=*/ 0, /*posBegin=*/ 1060,
                   /*posEnd=*/ 1140, /*lineBegin=*/ 0, /*lineEnd=*/ 100);

    return 0;
}
