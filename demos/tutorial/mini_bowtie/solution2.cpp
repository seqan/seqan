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

template <typename TIter, typename TStringSet>
void search(TIter & it, TStringSet const & pattern)
{
    typedef typename Iterator<TStringSet const, Standard>::Type TPatternIter;

    for (TPatternIter patternIt = begin(pattern, Standard()); patternIt != end(pattern, Standard()); ++patternIt)
    {
        unsigned startApproxSearch = length(value(patternIt)) / 2;
        goDown(it, infix(value(patternIt), 0, startApproxSearch - 1));
        goRoot(it);
    }
}

int main(int argc, char *argv[]) 
{
    typedef String<Dna5> TString;
    typedef StringSet<String<Dna5> > TStringSet;
    typedef Index<StringSet<TString>, FMIndex<> > TIndex;
    typedef Iterator<TIndex, TopDown<ParentLinks<> > >::Type TIter;

    /*String<Dna> text = "ACGTACGT";
    Index<String<Dna>, FMIndex<> > index(text);

    Finder<Index<String<Dna> > > finder(index);

    find(finder, "AC");
    std::cerr << position("AC") << std::endl;
*/
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

    // combining the contigs of the reference into one string set
    StringSet<TString> text;
    for (unsigned i = 0; i < length(fragStore.contigStore); ++i)
        appendValue(text, fragStore.contigStore[i].seq);
        
    // forward search
    reverse(text);
    TIndex fmIndex(text);
    TIter it(fmIndex);
    search(it, fragStore.readSeqStore);
    clear(fmIndex);
    clear(it);

    // reversing the sequences for backward search
    reverse(text);
    reverse(fragStore.readSeqStore);

    // backward search
    fmIndex = TIndex(text);
    it = TIter(fmIndex);
    search(it, fragStore.readSeqStore);
    clear(fmIndex);
    clear(it);

    return 0;
}
