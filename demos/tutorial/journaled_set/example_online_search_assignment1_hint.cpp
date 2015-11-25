//![main]
#include <iostream>
#include <seqan/seq_io.h>
#include <seqan/journaled_set.h>

using namespace seqan;

template <typename TString, typename TSpec>
inline int
loadAndJoin(StringSet<TString, Owner<JournaledSet> > & /*journalSet*/,
            SeqFileIn & databaseFile,
            JoinConfig<TSpec> const & /*joinConfig*/)
{
    typedef typename Host<TString>::Type THost;

    // [A] Ensure the Journal Set is not occupied by other sequences.

    // Construct the temporary buffers for the read id and sequence.
    String<char> tempSeqId;
    THost tempSeq;

    // No sequences in the fasta file!
    if (atEnd(databaseFile))
    {
        std::cerr << "Empty FASTA file." << std::endl;
        return -1;
    }
    // First read sequence for reference sequence.
    readRecord(tempSeqId, tempSeq, databaseFile);

    // [B] Set the reference sequence to the Journal Set

    // Read remaining sequences.
    while (!atEnd(databaseFile))
    {
        readRecord(tempSeqId, tempSeq, databaseFile);
        // [C] Append and join the current read sequence.
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
    // [D] Construct Journaled Set and call loadAndJoin
    loadAndJoin(journalSet, databaseFile, joinConfig);

    std::cout << "Done!" << std::endl;

    return 0;
}
//![main]
