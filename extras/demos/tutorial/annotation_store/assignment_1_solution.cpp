#include <fstream>
#include <iostream>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/store.h>

using namespace seqan;
int main()
{
    FragmentStore<> store;
    std::ifstream file("assignment_annotations.gtf", std::ios_base::in | std::ios_base::binary);
    read(file, store, Gtf());
    // Create AnnotationTree iterator
    Iterator<FragmentStore<>, AnnotationTree<> >::Type it;
    it = begin(store, AnnotationTree<>());
    unsigned count = 0;
    // Go down to the first leaf (first child of the first mRNA)
    while (goDown(it)) ;
    ++count;
    // Iterate over all siblings and count
    while (goRight(it))
        ++count;
    std::cout << "No. of children of the first mRNA: " << count << std::endl;
    return 0;
}
