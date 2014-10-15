// FRAGMENT(main)
#include <iostream>
#include <seqan/seq_io.h>
#include <seqan/journaled_set.h>

using namespace seqan;

template <typename TString, typename TSpec>
void loadAndJoin(StringSet<TString, Owner<JournaledSet> > & journalSet,
                 SeqFileIn & databaseFile,
                 JoinConfig<TSpec> const & joinConfig)
{
    typedef typename Host<TString>::Type THost;

    clear(journalSet);

    String<char> seqId;
    THost sequence;

    // No sequences in the fasta file!
    if (atEnd(databaseFile))
        throw IOError("empty FASTA file");

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
}

int main()
{
    // Definition of the used types.
    typedef String<Dna,Alloc<> > TSequence;
    typedef String<Dna,Journaled<Alloc<>,SortedArray,Alloc<> > > TJournal;
    typedef StringSet< TJournal, Owner<JournaledSet> > TJournaledSet;

    // Open the stream to the file containing the sequences.
    CharString seqDatabasePath = "/path/to/your/fasta/file/sequences.fasta";
    SeqFileIn databaseFile(toCString(seqDatabasePath));

    // Reading each sequence and journal them.
    TJournaledSet journalSet;
    JoinConfig<GlobalAlign<JournaledCompact> > joinConfig;
    try
    {
        loadAndJoin(journalSet, databaseFile, joinConfig);
    }
    catch (ParseError const & err)
    {
        std::cerr << "Problem with parsing the files: " << err.what();
        return 1;
    }
    catch (IOError const & err)
    {
        std::cerr << "Problem during I/O: " << err.what();
        return 1;
    }

    return 0;
}
