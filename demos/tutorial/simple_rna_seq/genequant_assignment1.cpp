#include <iostream>
#include <seqan/store.h>
#include <seqan/misc/interval_tree.h>
#include <seqan/parallel.h>

using namespace seqan2;


// define used types
typedef FragmentStore<> TStore;
//
// 2. Load annotations and alignments from files
//
bool loadFiles(TStore & store, std::string const & annotationFileName,  std::string const & alignmentFileName)
{
    (void) store;
    (void) annotationFileName;
    (void) alignmentFileName;
    // INSERT YOUR CODE HERE ...
    //

    return true;
}

//![main]
int main(int, char const **)
{
    TStore store;
    std::string annotationFileName = getAbsolutePath("demos/tutorial/simple_rna_seq/example.gtf");
    std::string alignmentFileName = getAbsolutePath("demos/tutorial/simple_rna_seq/example.sam");

    if (!loadFiles(store, annotationFileName, annotationFileName))
        return 1;

    return 0;
}
//![main]
