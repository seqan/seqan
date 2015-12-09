// Demo for the updated interval trees with the GenomicRegion class
// that implements the IntervalConcept.
//
// In this demo, we use Interval class with unsigned cargo and int position
// type.

#include <vector>

#include <seqan/basic.h>
#include <seqan/interval_tree.h>
#include <seqan/sequence.h>
#include <seqan/seq_io/genomic_region.h>
#include <seqan/stream.h>

using namespace seqan;

int main()
{
    // Using vector here, but any container will do.
    std::vector<GenomicRegion> intervals;
    // Note that GenomicRegion parses 1-based coordinates.
    intervals.push_back(GenomicRegion("chr1:10-20"));  // => [9,20)
    intervals.push_back(GenomicRegion("chr1:15-25"));  // => [14,25)

    // Build interval tree, not modifiable from here.
    IntervalTree<GenomicRegion, Static> tree(intervals);

    std::cout << "length(tree) == " << length(tree) << "\n\n";

    // Query interval tree, results are iterators into the tree, each points at
    // a TValue of the IntervalTree.  We are using vectors here, but any
    // container will do.
    typedef Iterator<IntervalTree<GenomicRegion, Static> const>::Type TIterator;

    // point query
    {
        std::vector<TIterator> result;
        findOverlappingWithPoint(tree, 15, result);
        std::cout << "query for point 15\n"
                  << "=>";
        for (unsigned i = 0; i < result.size(); ++i)
        {
            CharString out;
            result[i]->toString(out);
            std::cout << " " << out;
        }
        std::cout << "\n\n";
    }

    // interval query
    {
        std::vector<TIterator> result;
        findOverlappingWithInterval(tree, 20, 30, result);
        std::cout << "query for interval [20, 30)\n"
                  << "=>";
        for (unsigned i = 0; i < result.size(); ++i)
        {
            CharString out;
            result[i]->toString(out);
            std::cout << " " << out;
        }
        std::cout << "\n\n";
    }

    return 0;
}
