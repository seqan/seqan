//![main]
#include <iostream>
#include <seqan/seq_io.h>
#include <seqan/journaled_set.h>

using namespace seqan;

template <typename TString, typename TSpec>
inline int
loadAndJoin(StringSet<TString, Owner<JournaledSet> > & journalSet,
            SeqFileIn & databaseFile,
            JoinConfig<TSpec> const & joinConfig)
{
    typedef typename Host<TString>::Type THost;
    // [A]
    clear(journalSet);

    // Construct the temporary buffers for the read id and sequence.
    String<char> tempSeqId;
    THost sequence;

    // No sequences in the fasta file!
    if (atEnd(databaseFile))
    {
        std::cerr << "Empty FASTA file." << std::endl;
        return -1;
    }
    // First read sequence for reference sequence.
    readRecord(tempSeqId, sequence, databaseFile);
    // [B]
    createHost(journalSet, sequence);  // When using create we copy the reference instead of storing a pointer.

    // Read remaining sequences.
    while (!atEnd(databaseFile))
    {
        readRecord(tempSeqId, sequence, databaseFile);
        // [C]
        appendValue(journalSet, TString(sequence)); // First we append the sequence to the set.
        join(journalSet, length(journalSet) - 1, joinConfig); // Second we join it to the set.
    }
    return 0;
}

int main()
{
    // Definition of the used types.
    typedef String<Dna, Journaled<Alloc<>, SortedArray, Alloc<> > > TJournal;
    typedef StringSet<TJournal, Owner<JournaledSet> > TJournaledSet;

    // Open the stream to the file containing the sequences.
    CharString seqDatabasePath = getAbsolutePath("demos/tutorial/journaled_set/sequences.fasta");
    SeqFileIn databaseFile(toCString(seqDatabasePath));

    // Reading each sequence and journal them.
    TJournaledSet journalSet;
    JoinConfig<GlobalAlign<JournaledCompact> > joinConfig;
    loadAndJoin(journalSet, databaseFile, joinConfig);

    std::cout << "Done!" << std::endl;
    return 0;
}
//![main]
