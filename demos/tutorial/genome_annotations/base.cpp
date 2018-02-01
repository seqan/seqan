//![INCLUDE]
#include <seqan/store.h>
//![INCLUDE]
//![MAIN]
#include <iostream>
#include <seqan/sequence.h>
#include <seqan/stream.h>

using namespace seqan;

int main()
{
//![MAIN]
//![SIGNATURE]
    FragmentStore<> store;
//![SIGNATURE]
//![LOAD]
    CharString fileName = getAbsolutePath("demos/tutorial/genome_annotations/example.gtf");
    GffFileIn file(toCString(fileName));

    readRecords(store, file);
//![LOAD]
//![ITERATOR]
    Iterator<FragmentStore<>, AnnotationTree<> >::Type it;
    it = begin(store, AnnotationTree<>());
//![ITERATOR]
//![MOVE]
    // Move the iterator down to a leaf
    while (goDown(it))
    {}
    // Create a new iterator and if possible move it to the right sibling of the first iterator
    Iterator<FragmentStore<>, AnnotationTree<> >::Type it2;
    if (isLastChild(it))
        it2 = nodeRight(it);
//![MOVE]
//![DFS]
    // Move the iterator back to the beginning
    goBegin(it);
    // Iterate over the nodes in preorder DFS while the end is not reached and
    // output if the current node is a leaf
    while (!atEnd(it))
    {
        if (isLeaf(it))
            std::cout << " current node is leaf" << std::endl;
        goNext(it);
    }
//![DFS]
//![ACCESS]
    // Move the iterator to the begin of the annotation tree
    it = begin(store, AnnotationTree<>());
    // Go down to exon level
    while (goDown(it)) ;
    std::cout << "type: " <<  getType(it) << std::endl;
    std::cout << "id: " << value(it) << std::endl;
    std::cout << "begin position: " <<  getAnnotation(it).beginPos << std::endl;
//![ACCESS]
//![CREATE]
    Iterator<FragmentStore<>, AnnotationTree<> >::Type it3;
    // Create a right sibling of the current node and return an iterator to this new node
    it3 = createSibling(it);
//![CREATE]
//![OUT]
    // Open output stream
    GffFileOut fileOut("example_out.gtf");
    // Write annotations to GTF file
    writeRecords(fileOut, store);
//![OUT]
    return 0;
}
