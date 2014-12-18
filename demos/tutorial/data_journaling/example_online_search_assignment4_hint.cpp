#include <seqan/basic.h>
#include <seqan/sequence.h>

using namespace seqan;

//![searchAtBorder]
template <typename TJournalEntriesIterator, typename TJournal, typename TPattern>
void _searchAtBorder(String<int> & hitTarget,
                     TJournalEntriesIterator & entriesIt,
                     TJournal const & journal,
                     TPattern const & pattern)
{
    // [A] Determine first position of the at which pattern crosses the border of current node.

    // [B] Determine last position before pattern exits the current node or reaches the end of the sequence.

    // [C] Move step by step over search region.

    // [D] Scan pattern in current window and report possible hits.
}
//![searchAtBorder]

//![main]
int main()
{
    return 0;
}
//![main]
