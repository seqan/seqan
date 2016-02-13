//![include]
#include <iostream>
#include <seqan/seq_io.h>
#include <seqan/journaled_set.h>

using namespace seqan;
//![include]

//![findInOriginalNode]
template <typename TJournalEntriesIterator, typename TPattern>
void _findInOriginalNode(String<int> & hitTarget,
                         TJournalEntriesIterator & entriesIt,
                         TPattern const & pattern,
                         String<int> const & refHits)
{
    // Define an Iterator which iterates over the reference hit set.
    typedef typename Iterator<String<int> const, Standard>::Type THitIterator;

    // [A] Check if hits exist in the reference.
    if (!empty(refHits))
    {
        // [B] Find upper bound to current physical position in sorted refHits using std::upper_bound.
        THitIterator itHit = std::upper_bound(begin(refHits), end(refHits), (int)entriesIt->physicalPosition);
        // [C] Make sure we do not miss hits that begin at physical position of current node.
        if (itHit != begin(refHits) && *(itHit - 1) >= (int)entriesIt->physicalPosition)
            --itHit;
        // [D] Store all hits that are found in the region of the reference which is covered by this node.
        while ((int)*itHit < ((int)entriesIt->physicalPosition + (int)entriesIt->length - (int)length(pattern) + 1) && itHit != end(refHits))
        {
            // [E] Store the correct virtual position and check next hit.
            appendValue(hitTarget, entriesIt->virtualPosition + (*itHit - (int)entriesIt->physicalPosition));
            ++itHit;
        }
    }
}
//![findInOriginalNode]

//![findPatternInJournalString]
template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec, typename TPattern>
void findPatternInJournalString(String<int> & hitTarget,
                                String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > const & journal,
                                TPattern const & pattern,
                                String<int> const & refHits)
{
    typedef String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > const TJournal;
    typedef typename JournalType<TJournal>::Type TJournalEntries;
    typedef typename Iterator<TJournalEntries>::Type TJournalEntriesIterator;

    if (length(pattern) > length(journal))
        return;

    TJournalEntriesIterator it = begin(journal._journalEntries);
    TJournalEntriesIterator itEnd = findInJournalEntries(journal._journalEntries, length(journal) - length(pattern) + 1) + 1;

    while (it != itEnd)
    {
        if (it->segmentSource == SOURCE_ORIGINAL) // Find a possible hit in the current source vertex.
        {
            _findInOriginalNode(hitTarget, it, pattern, refHits);
        }
        if (it->segmentSource == SOURCE_PATCH) // Search for pattern within the patch node.
        { //            _findInPatchNode(hitTarget, it, journal, pattern);
        }
        // Scan the border for a possible match.
//        _searchAtBorder(hitTarget, it, journal, pattern);
        ++it;
    }
}
//![findPatternInJournalString]

//![findPatternInReference]
template <typename TString, typename TPattern>
void findPatternInReference(String<int> & hits,
                            TString const & reference,
                            TPattern const & pattern)
{
    // Check whether the pattern fits into the sequence.
    if (length(pattern) > length(reference))
        return;

    //
    for (unsigned pos = 0; pos < length(reference) - length(pattern) + 1; ++pos)
    {
        bool isHit = true;

        for (unsigned posPattern = 0; posPattern < length(pattern); ++posPattern)
        {
            if (pattern[posPattern] != reference[posPattern + pos])
            {
                isHit = false;
                break;
            }
        }
        // Report the position if found a hit.
        if (isHit)
            appendValue(hits, pos);
    }
}
//![findPatternInReference]

//![searchPattern]
template <typename TString, typename TPattern>
void searchPattern(StringSet<String<int> > & hitSet,
                   StringSet<TString, Owner<JournaledSet> > const & journalSet,
                   TPattern const & pattern)
{
    typedef StringSet<TString, Owner<JournaledSet> > TJournalSet;
    typedef typename Host<TJournalSet const>::Type THost;

    // Check for valid initial state.
    if (empty(host(journalSet)))
    {
        std::cout << "No reference set. Aborted search!" << std::endl;
        return;
    }

    // Reset the hitSet to avoid phantom hits coming from a previous search.
    clear(hitSet);
    resize(hitSet, length(journalSet) + 1);
    // Access the reference sequence.
    THost & globalRef = host(journalSet);
    // Search for pattern in the reference sequence.
    findPatternInReference(hitSet[0], globalRef, pattern);
    // Search for pattern in the journaled sequences.
    for (unsigned i = 0; i < length(journalSet); ++i)
        findPatternInJournalString(hitSet[i + 1], journalSet[i], pattern, hitSet[0]);
}
//![searchPattern]

//![laodAndJoin]
template <typename TString, typename TSpec>
inline int
loadAndJoin(StringSet<TString, Owner<JournaledSet> > & journalSet,
            SeqFileIn & databaseFile,
            JoinConfig<TSpec> const & joinConfig)
{
    typedef typename Host<TString>::Type THost;

    clear(journalSet);

    String<char> seqId;
    THost sequence;

    // No sequences in the fasta file!
    if (atEnd(databaseFile))
    {
        std::cerr << "Empty FASTA file." << std::endl;
        return -1;
    }
    // First read sequence for reference sequence.
    readRecord(seqId, sequence, databaseFile);

    // We have to create the global reference sequence otherwise we loose the information after this function terminates.
    createHost(journalSet, sequence);

    // If there are more
    while (!atEnd(databaseFile))
    {
        readRecord(seqId, sequence, databaseFile);
        appendValue(journalSet, TString(sequence));
        join(journalSet, length(journalSet) - 1, joinConfig);
    }
    return 0;
}
//![laodAndJoin]

//![main]
int main()
{
    // Definition of the used types.
    typedef String<Dna, Alloc<> > TSequence;
    typedef String<Dna, Journaled<Alloc<>, SortedArray, Alloc<> > > TJournal;
    typedef StringSet<TJournal, Owner<JournaledSet> > TJournaledSet;

    // Open the stream to the file containing the sequences.
    CharString seqDatabasePath = getAbsolutePath("demos/tutorial/journaled_set/sequences.fasta");
    SeqFileIn databaseFile(toCString(seqDatabasePath));

    // Reading each sequence and journal them.
    TJournaledSet journalSet;
    JoinConfig<GlobalAlign<JournaledCompact> > joinConfig;
    loadAndJoin(journalSet, databaseFile, joinConfig);

    // Define a pattern and start search.
    StringSet<String<int> > hitSet;
    TSequence pattern = "GTGGT";
    std::cout << "Search for: " << pattern << ":\n";
    searchPattern(hitSet, journalSet, pattern);
//![main]

//![printResultReference]
    if (empty(hitSet[0]))
    {
        std::cout << "No hit in reference " << std::endl;
    }
    else
    {
        std::cout << "Hit in reference " << " at ";
        for (unsigned j = 0; j < length(hitSet[0]); ++j)
            std::cout << hitSet[0][j] << ": " << infix(host(journalSet), hitSet[0][j], hitSet[0][j] + length(pattern)) << "\t";
    }
    std::cout << std::endl;
//![printResultReference]

//![printResultJournalSequence]
    for (unsigned i = 1; i < length(hitSet); ++i)
    {
        if (empty(hitSet[i]))
        {
            std::cout << "No hit in sequence " << i - 1 << std::endl;
        }
        else
        {
            std::cout << "Hit in sequence " << i - 1 << " at ";
            for (unsigned j = 0; j < length(hitSet[i]); ++j)
                std::cout << hitSet[i][j] << ": " << infix(value(journalSet, i - 1), hitSet[i][j], hitSet[i][j] + length(pattern)) << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << "Done!" << std::endl;
    return 0;
}
//![printResultJournalSequence]
