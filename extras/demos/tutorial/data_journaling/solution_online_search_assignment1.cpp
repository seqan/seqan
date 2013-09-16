// FRAGMENT(main)
#include <iostream>
#include <seqan/seq_io.h>
#include <seqan/journaled_set.h>

using namespace seqan;

template <typename TString, typename TStream, typename TSpec>
inline int
loadAndJoin(StringSet<TString, Owner<JournaledSet> > & journalSet,
            TStream & stream,
            JoinConfig<TSpec> const & joinConfig)
{
    typedef typename Host<TString>::Type THost;
    // Define the RecordReader.
    RecordReader<std::ifstream, SinglePass<> > reader(stream);
    // [A]
    clear(journalSet);

    // Construct the temporary buffers for the read id and sequence.
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
    // [B]
    createGlobalReference(journalSet, sequence);  // When using create we copy the reference instead of storing a pointer.

    // Read remaining sequences.
    while (!atEnd(reader))
    {
        if (readRecord(tempSeqId, sequence, reader, Fasta()) != 0)
        {
            std::cerr << "ERROR reading FASTA." << std::endl;
            return 1;
        }
        // [C]
        appendValue(journalSet, TString(sequence)); // First we append the sequence to the set.
        join(journalSet, length(journalSet) - 1, joinConfig); // Second we join it to the set.
    }
    return 0;
}


int main()
{
    // Definition of the used types.
    typedef String<Dna,Alloc<> > TSequence;
    typedef String<Dna,Journaled<Alloc<>,SortedArray,Alloc<> > > TJournal;
    typedef StringSet< TJournal, Owner<JournaledSet> > TJournaledSet;

    // Open the stream to the file containing the sequences.
    String<char> seqDatabasePath = "/Users/rahn_r/Downloads/sequences.fasta";
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

    return 0;
}
