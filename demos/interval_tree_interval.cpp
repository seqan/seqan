// Demo for the updated interval trees with the Interval<> class.
//
// In this demo, we use Interval class with unsigned cargo and int position
// type.

#include <vector>

#include <seqan/basic.h>
#include <seqan/interval_tree.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>

using namespace seqan;

int main()
{
    typedef Interval<unsigned, int> TInterval;

    // Using vector here, but any container will do.
    std::vector<TInterval > intervals;
    intervals.push_back(TInterval(0, 10, 20));
    intervals.push_back(TInterval(1, 15, 25));

    // Build interval tree, not modifiable from here.
    IntervalTree<Interval<unsigned, int>, Static> tree(intervals);

    std::cout << "length(tree) == " << length(tree) << "\n\n";

    // Query interval tree, results are iterators into the tree, each points at
    // a TValue of the IntervalTree.  We are using vectors here, but any
    // container will do.
    typedef Iterator<IntervalTree<TInterval, Static> const>::Type TIterator;

    // point query
    {
        std::vector<TIterator> result;
        findOverlappingWithPoint(tree, 15, result);
        std::cout << "query for point 15\n"
                  << "=>";
        for (unsigned i = 0; i < result.size(); ++i)
            std::cout << " " << cargo(*result[i]);
        std::cout << "\n\n";
    }

    // interval query
    {
        std::vector<TIterator> result;
        findOverlappingWithInterval(tree, 20, 30, result);
        std::cout << "query for interval [20, 30)\n"
                  << "=>";
        for (unsigned i = 0; i < result.size(); ++i)
            std::cout << " " << cargo(*result[i]);
        std::cout << "\n\n";
    }

    return 0;
}
