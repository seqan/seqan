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
    // Create iterator
    Iterator<FragmentStore<>, AnnotationTree<> >::Type it;
    it = begin(store, AnnotationTree<>());
    // Iterate to the first annotation of type "exon"
    while (!atEnd(it) && getType(it) != "exon") goNext(it);
    // Output:
    std::cout << "  type: " << getType(it) << std::endl;
    std::cout << "  begin position: " << getAnnotation(it).beginPos << std::endl;
    std::cout << "  end position: " << getAnnotation(it).endPos << std::endl;
    std::cout << "  id: " << value(it) << std::endl;
    std::cout << "  parent id: " << getAnnotation(it).parentId << std::endl;
    std::cout << "  parent name: " << getParentName(it) << std::endl;
    return 0;
}
