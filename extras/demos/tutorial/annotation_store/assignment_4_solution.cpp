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
    unsigned countGenes = 0;
    unsigned countmRNAs = 0;
    unsigned countExons = 0;
    unsigned length = 0;
    // Iterate over annotation tree and count different elements and compute exon lengths
    while (!atEnd(it)){
        if (getType(it) == "gene") ++countGenes;
        else if (getType(it) == "mRNA") ++countmRNAs;
        else if (getType(it) == "exon") {
            ++countExons;
            length += abs((int)getAnnotation(it).endPos - (int)getAnnotation(it).beginPos);
        }
        goNext(it);
    }
    // Ouput some stats:
    std::cout << "Average number of mRNAs for genes: " << (float)countmRNAs/(float)countGenes << std::endl;
    std::cout << "Average number of exons for mRNAs: " << (float)countExons/(float)countmRNAs << std::endl;
    std::cout << "Average length of exons: " << (float)length/(float)countExons << std::endl;
    return 0;
}
