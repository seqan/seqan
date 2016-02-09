#include <seqan/basic.h>
#include <seqan/sequence.h>

using namespace seqan;

//![assignment2]
template <typename TString, typename TPattern>
void findPatternInReference(String<int> & hits,
                            TString const & reference,
                            TPattern const & pattern)
{
// [A] Check whether pattern fits into the sequence.

// [B] Iterate over all positions at which the pattern might occur.

// [C] Evaluate all positions of the pattern until you find a mismatch or you have found a hit.

// [D] Report begin position at which pattern matches the sequence.
}
//![assignment2]

//![assignment3]
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
//![assignment3]

//![assignment4]
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
//![assignment4]

int main()
{
	return 0;
}
