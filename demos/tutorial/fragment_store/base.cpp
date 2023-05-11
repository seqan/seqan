#include <seqan/basic.h>
#include <seqan/store.h>

using namespace seqan2;

int main()
{
    typedef FragmentStore<> TStore;
    TStore store;

//![iterator]
    Iterator<FragmentStore<>, AnnotationTree<> >::Type it;
    it = begin(store, AnnotationTree<>());
//![iterator]

    return 0;
}
