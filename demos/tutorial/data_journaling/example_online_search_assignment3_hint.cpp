#include <seqan/basic.h>
#include <seqan/sequence.h>

using namespace seqan;

//![findInOriginalNode]
template <typename TJournalEntriesIterator, typename TPattern>
void _findInOriginalNode(String<int> & hitTarget,
                         TJournalEntriesIterator & entriesIt,
                         TPattern const & pattern,
                         String<int> const & refHits)
{
    // [A] Check if hits exist in the reference.

    // [B] Find upper bound to current physical position in sorted refHits using std::upper_bound.

    // [C] Make sure we do not miss hits that begin at physical position of current node.

    // [D] Store all hits that are found in the region of the reference which is covered by this node.

    // [E] Store the correct virtual position and check next hit.
}
//![findInOriginalNode]

//![main]
int main()
{
    return 0;
}
//![main]
