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
    // Open the stream to the file containing the sequences.
    CharString seqDatabasePath = "/path/to/your/fasta/file/sequences.fasta";
    SeqFileIn databaseFile(toCString(seqDatabasePath));

    // Reading each sequence and journal them.
    // [D] Construct Journaled Set and call loadAndJoin

    return 0;
}
//![main]
