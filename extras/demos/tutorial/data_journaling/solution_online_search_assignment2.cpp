// FRAGMENT(include)
#include <iostream>
#include <seqan/seq_io.h>
#include <seqan/journaled_set.h>

using namespace seqan;

// FRAGMENT(findPatternInReference)
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
            if(pattern[posPattern] != reference[posPattern + pos])
            {
                isHit = false;
                break;
            }
        }
        // [D] Report begin position at which pattern matches the sequence.
        if(isHit)
            appendValue(hits, pos);
    }
}

// FRAGMENT(searchPattern)
template <typename TString, typename TPattern>
void searchPattern(StringSet<String<int> > & hitSet,
                   StringSet<TString, Owner<JournaledSet> > const & journalSet,
                   TPattern const & pattern)
{
    typedef typename Host<TString>::Type THost;

    // Check for valid initial state.
    if (empty(globalReference(journalSet)))
    {
        std::cout << "No reference set. Aborted search!" << std::endl;
        return;
    }

    // Reset the hitSet to avoid phantom hits coming from a previous search.
    clear(hitSet);
    resize(hitSet, length(journalSet) + 1);
    // Access the reference sequence.
    THost & globalRef = globalReference(journalSet);
    // Search for pattern in the reference sequence.
    findPatternInReference(hitSet[0], globalRef, pattern);

    // Search for pattern in the journaled sequences.
    for (unsigned i = 0; i < length(journalSet); ++i)
    {
//        findPatternInJournalString(hitSet[i+1], journalSet[i], pattern, hitSet[0]);
    }
}

// FRAGMENT(loadAndJoin)
template <typename TString, typename TStream, typename TSpec>
inline int
loadAndJoin(StringSet<TString, Owner<JournaledSet> > & journalSet,
            TStream & stream,
            JoinConfig<TSpec> const & joinConfig)
{
    typedef typename Host<TString>::Type THost;

    RecordReader<std::ifstream, SinglePass<> > reader(stream);

    clear(journalSet);

    String<char> tempSeqId;
    THost sequence;

    // No sequences in the fasta file!
    if (atEnd(reader))
    {
        std::cerr << "Empty FASTA file." << std::endl;
        return -1;
    }
    // First read sequence for reference sequence.
    if (readRecord(tempSeqId, sequence, reader, Fasta()) != 0)
    {
        std::cerr << "ERROR reading FASTA." << std::endl;
        return 1;
    }
    // We have to create the global reference sequence otherwise we loose the information after this function terminates.
    createGlobalReference(journalSet, sequence);

    // If there are more
    while (!atEnd(reader))
    {
        if (readRecord(tempSeqId, sequence, reader, Fasta()) != 0)
        {
            std::cerr << "ERROR reading FASTA." << std::endl;
            return 1;
        }
        appendValue(journalSet, TString(sequence));
        join(journalSet, length(journalSet) - 1, joinConfig);
    }
    return 0;
}

// FRAGMENT(main)
int main()
{
    // Definition of the used types.
    typedef String<Dna,Alloc<> > TSequence;
    typedef String<Dna,Journaled<Alloc<>,SortedArray,Alloc<> > > TJournal;
    typedef StringSet< TJournal, Owner<JournaledSet> > TJournaledSet;

    // Open the stream to the file containing the sequences.
    String<char> seqDatabasePath =  "/Users/rahn_r/Downloads/sequences.fasta";
    std::ifstream databaseFile(toCString(seqDatabasePath), std::ios_base::in);
    if(!databaseFile.good())
    {
        std::cerr << "Cannot open file <" << seqDatabasePath << ">!" << std::endl;
    }


    // Reading each sequence and journal them.
    TJournaledSet journalSet;
    JoinConfig<GlobalAlign<JournaledCompact> > joinConfig;
    loadAndJoin(journalSet, databaseFile, joinConfig);
    databaseFile.close();

    // Define a pattern and start search.
    StringSet<String<int> > hitSet;
    TSequence pattern = "GTGGT";
    std::cout << "Search for: " << pattern << ":\n";
    searchPattern(hitSet, journalSet, pattern);

    // FRAGMENT(printResult)
    if (empty(hitSet[0]))
    {
        std::cout << "No hit in reference " << std::endl;
    }
    else
    {
        std::cout << "Hit in reference " << " at ";
        for (unsigned j = 0; j < length(hitSet[0]); ++j)
            std::cout << hitSet[0][j] << ": " << infix(globalReference(journalSet), hitSet[0][j],hitSet[0][j] + length(pattern)) << "\t";
    }
    std::cout << std::endl;

    return 0;
}
