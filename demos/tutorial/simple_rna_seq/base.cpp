#include <seqan/store.h>

using namespace seqan;

int main()
{
//![store]
    FragmentStore<> store;
//![store]

//![typedefs]
    typedef FragmentStore<> TStore;
    typedef Value<TStore::TAnnotationStore>::Type TAnnotation;
    typedef TAnnotation::TId TId;
    typedef TAnnotation::TId TPos;
    typedef IntervalAndCargo<TPos, TId> TInterval;
//![typedefs]
//![interval]
    String<String<TInterval> > intervals;
//![interval]

//![resize]
    resize(intervals, length(store.contigStore));
//![resize]

//![tree]
    typedef IntervalTree<TPos, TId> TIntervalTree;
    String<TIntervalTree> intervalTrees;
//![tree]

//![resize_tree]
    resize(intervals, length(store.contigStore));
//![resize_tree]

//![reads]
    String<unsigned> readsPerGene;
//![reads]

//![read_alignment_type]
    typedef Value<TStore::TAlignedReadStore>::Type TAlignedRead;
//![read_alignment_type]
    TAlignedRead tmp; // sp no unused-Variable warning is triggered

//![resize_reads]
    resize(readsPerGene, length(store.annotationStore), 0);
//![resize_reads]

//![result]
    String<TId> result;
//![result]

    return 0;
}
