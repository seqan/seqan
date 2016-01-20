#include <iostream>
#include <seqan/store.h>
#include <seqan/misc/interval_tree.h>
#include <seqan/parallel.h>

using namespace seqan;


// define used types
typedef FragmentStore<> TStore;

//![solution]
//
// 1. Load annotations and alignments from files
//
bool loadFiles(TStore & store, std::string const & annotationFileName,  std::string const & alignmentFileName)
{
    BamFileIn alignmentFile;
    if (!open(alignmentFile, alignmentFileName.c_str()))
    {
        std::cerr << "Couldn't open alignment file " << alignmentFileName << std::endl;
        return false;
    }
    std::cout << "Loading read alignments ..... " << std::flush;
    readRecords(store, alignmentFile);
    std::cout << "[" << length(store.alignedReadStore) << "]" << std::endl;

    // load annotations
    GffFileIn annotationFile;
    if (!open(annotationFile, toCString(annotationFileName)))
    {
        std::cerr << "Couldn't open annotation file" << annotationFileName << std::endl;
        return false;
    }
    std::cout << "Loading genome annotation ... " << std::flush;
    readRecords(store, annotationFile);
    std::cout << "[" << length(store.annotationStore) << "]" << std::endl;

    return true;
}
//![solution]

//![main]
int main(int argc, char const * argv[])
{
    TStore store;
    std::string annotationFileName = getAbsolutePath("demos/tutorial/simple_rna_seq/example.gtf");
    std::string alignmentFileName = getAbsolutePath("demos/tutorial/simple_rna_seq/example.sam");

    if (!loadFiles(store, annotationFileName, alignmentFileName))
        return 1;

    return 0;
}
//![main]
