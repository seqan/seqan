//![include]
#include <iostream>
#include <seqan/seq_io.h>
#include <seqan/journaled_set.h>

using namespace seqan;
//![include]

//![findPatternInReference]
template <typename TString, typename TPattern>
void findPatternInReference(String<int> & hits,
                            TString const & reference,
                            TPattern const & pattern)
{
    // [A] Check whether pattern fits into the sequence.
    if (length(pattern) > length(reference))
        return;

    // [B] Iterate over all positions at which the pattern might occur.
    for (unsigned pos = 0; pos < length(reference) - length(pattern) + 1; ++pos)
    {
        bool isHit = true;
        // [C] Evaluate all positions of the pattern until you find a mismatch or you have found a hit.
        for (unsigned posPattern = 0; posPattern < length(pattern); ++posPattern)
        {
            if (pattern[posPattern] != reference[posPattern + pos])
            {
                isHit = false;
                break;
            }
        }
        // [D] Report begin position at which pattern matches the sequence.
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
    {
//        findPatternInJournalString(hitSet[i+1], journalSet[i], pattern, hitSet[0]);
    }
}
//![searchPattern]

//![loadAndJoin]
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
//![loadAndJoin]

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

//![printResult]
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

    std::cout << "Done!" << std::endl;
    return 0;
}
//![printResult]
