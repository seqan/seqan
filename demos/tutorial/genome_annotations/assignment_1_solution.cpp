#include <iostream>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/store.h>

using namespace seqan;
int main()
{
    CharString fileName = getAbsolutePath("demos/tutorial/genome_annotations/assignment_annotations.gtf");
    GffFileIn file(toCString(fileName));

    FragmentStore<> store;
    readRecords(store, file);
    // Create AnnotationTree iterator
    Iterator<FragmentStore<>, AnnotationTree<> >::Type it;
    it = begin(store, AnnotationTree<>());
    unsigned count = 0;
    // Go down to the first leaf (first child of the first mRNA)
    while (goDown(it))
    {}
    std::cout << "Is leaf: " << isLeaf(it) << std::endl;

    ++count;
    // Iterate over all siblings and count
    while (goRight(it))
        ++count;
    std::cout << "No. of children of the first mRNA: " << count << std::endl;
    return 0;
}
