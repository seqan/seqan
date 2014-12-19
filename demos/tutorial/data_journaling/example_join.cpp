//![main]
#include <iostream>
#include <seqan/file.h>
#include <seqan/journaled_set.h>

using namespace seqan;

int main()
{
    typedef String<char, Journaled<Alloc<>, SortedArray, Alloc<> > > TJournalString;
    typedef Host<TJournalString>::Type THost;
    typedef StringSet<TJournalString, Owner<JournaledSet> > TJournaledSet;

    TJournaledSet journaledSet;

    THost reference = "DPKKPRGKMSSYAFFVQTSREEHKKKHPDASVNFSEFSKKCSERWKTMSAKEKGKFEDMAKADKARYEREMKTYIPPKGE";
    THost seq0 = "DPKKPRGKMVNSPPAFFVQTSREEHKKKHPDASVFSKKCSERWKTMSAKEKGKFEDMAKARYEREMKTTYIPKGETYIPPKGE";
    THost seq1 = "DPHHPPKPRGKMVNSPPAFFVQTSREEHKPDASVFSKKCSERRMPNHHTMSAKEKGKFEDMAKARYEREMKTTYIPKGETYIPPKGE";
    THost seq2 = "DPKKPRGKMSSYAFFVQTSREEHKKKHPKKCDEFSKKCSERWKTMSAKEKGKFEDARYEREMKTYIPPKGE";
//![main]

//![init]
    setHost(journaledSet, reference);
    appendValue(journaledSet, TJournalString(seq0));
    appendValue(journaledSet, TJournalString(seq1));
    appendValue(journaledSet, TJournalString(seq2));
//![init]

//![join]
    join(journaledSet, 0, JoinConfig<GlobalAlign<JournaledManhatten> >());  // Simply inserts the
    join(journaledSet, 1, JoinConfig<GlobalAlign<JournaledCompact> >());    // Uses default scoring scheme to compute compact journal.
    JoinConfig<GlobalAlign<JournaledCompact> > joinConfig;
    setScoringScheme(joinConfig, Score<int, BiAffine>(0, -1, -1));    // Note the mismatch score is forbidden internally when used in the context of journaling.
    join(journaledSet, 2, joinConfig);  // Compute journal using Levenshtein distance.

    std::cout << "Reference: " << host(journaledSet) << std::endl;
    for (unsigned i = 0; i < length(journaledSet); ++i)
        std::cout << "Journaled Sequence " << i << ": " << value(journaledSet, i) << std::endl;

    return 0;
}
//![join]
