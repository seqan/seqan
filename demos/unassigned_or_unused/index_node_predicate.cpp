///A tutorial showing how to extent an index with a node predicate.
#include <iostream>
#include <seqan/index.h>

using namespace seqan;

/// constraint parameters
struct TMyConstraints
{
    double p_min;
    unsigned int replen_max;
    bool _cachedPred;
};

/// SeqAn extensions
namespace seqan {
// custom TSpec for our customized wotd-Index
struct TMyIndex;

template <typename TText>
struct Cargo<Index<TText, IndexWotd<TMyIndex> > >
{
    typedef TMyConstraints Type;
};

// node predicate
template <typename TText, typename TSpec>
bool nodePredicate(Iter<Index<TText, IndexWotd<TMyIndex> >, TSpec> & it)
{
    TMyConstraints & cons = cargo(container(it));
    unsigned int delta = countSequences(container(it)) * repLength(it);
    unsigned int textLen = length(container(it));

    if (textLen <= delta)
        return cons._cachedPred = true;

    return cons._cachedPred =
        ((double)countOccurrences(it) / (double)(textLen - delta)) > cons.p_min;
}

// monotonic hull
template <typename TText, typename TSpec>
bool nodeHullPredicate(Iter<Index<TText, IndexWotd<TMyIndex> >, TSpec> & it)
{
    TMyConstraints const & cons = cargo(container(it));
    unsigned int textLen = length(container(it));

    if (repLength(it) > cons.replen_max)
        return false;

    unsigned int delta = countSequences(container(it)) * cons.replen_max;
    if (textLen <= delta)
        return true;

    return ((double)countOccurrences(it) /
            (double)(textLen - delta)) > cons.p_min;
}

}

int main()
{
///We begin with a @Class.String@ to store our sequence.
    String<char> myString = "How many wood would a woodchuck chuck.";

///Then we create our customized index which is a specialization
///of the deferred @Class.Index.wotd-Index@
    typedef Index<String<char>, IndexWotd<TMyIndex> > TMyIndex;
    TMyIndex myIndex(myString);

    cargo(myIndex).replen_max = 10;
    cargo(myIndex).p_min = 0.05;

///To find all strings that fulfill our constraints,
///we simply do a dfs-traversal via @Function.goBegin@ and @Function.goNext@.
    typedef Iterator<TMyIndex, TopDown<ParentLinks<> > >::Type TConstrainedIterator;
    TConstrainedIterator myConstrainedIterator(myIndex);

    goBegin(myConstrainedIterator);
    while (!atEnd(myConstrainedIterator))
    {

///@Function.countOccurrences@ returns the number of hits of the representative.
        std::cout << countOccurrences(myConstrainedIterator) << "x  ";

///The representative string can be determined with @Function.representative@
        std::cout << "\t\"" << representative(myConstrainedIterator) << '\"' << std::endl;

        goNext(myConstrainedIterator);
    }

    return 0;
}
