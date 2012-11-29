// ==========================================================================
//                                mini_bowtie
// ==========================================================================

#include <iostream>
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/index.h>
#include <seqan/store.h>

using namespace seqan;

void search() {}

int main(int argc, char *argv[]) 
{
    // type definitions
    typedef String<Dna5> TString;
    typedef StringSet<TString> TStringSet;
    typedef Index<StringSet<TString>, FMIndex<> > TIndex;
    typedef Iterator<TIndex, TopDown<ParentLinks<> > >::Type TIter;

    // reading the command line arguments
    if (argc < 3) {
        std::cerr << "Invalid number of arguments." << std::endl
                  << "USAGE: minimapper GENOME.fasta READS.fasta OUT.sam" << std::endl;
        return 1;
    }

    // declaration and initialization of the fragment store
    FragmentStore<> fragStore;
    if (!loadContigs(fragStore, argv[1])) return 1;
    if (!loadReads(fragStore, argv[2])) return 1;

    StringSet<TString> text;
    for (unsigned i = 0; i < length(fragStore.contigStore); ++i)
        appendValue(text, fragStore.contigStore[i].seq);
        
    TIndex fmIndex(text);
    TIter it(fmIndex);
    search();
    clear(fmIndex);
    clear(it);

    reverse(text);
    reverse(fragStore.readSeqStore);

    fmIndex = TIndex(text);
    it = TIter(fmIndex);
    search();
    clear(fmIndex);
    clear(it);

    reverse(text);
    reverse(fragStore.readSeqStore);
    clear(fmIndex);
    clear(it);

    return 0;
}
