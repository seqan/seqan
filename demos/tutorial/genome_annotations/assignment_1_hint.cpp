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
    // Move iterator one node down
    goDown(it);

    std::cout << "Is leaf: " << isLeaf(it) << std::endl;
    return 0;
}
