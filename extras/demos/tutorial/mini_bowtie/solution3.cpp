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
        if (goDown(it, infix(value(patternIt), startApproxSearch + 1, length(value(patternIt)))))
        {
            for (unsigned i = startApproxSearch; ; --i)
            {
                for (Dna5 c = MinValue<Dna>::VALUE; c < +ValueSize<Dna>::VALUE; ++c)
                {
                    TIter localIt = it;
                    if (goDown(localIt, c))
                    {
                        if (goDown(localIt, infix(value(patternIt), 0, i)))
                        {
                            // HIT
                        }
                    }
                }
                if (i == 0 || !goDown(it, getValue(patternIt)[i]))
                {
                    break;
                }
            }
        }
        goRoot(it);
    }
}

int main(int argc, char *argv[]) 
{
    typedef String<Dna5> TString;
    typedef StringSet<String<Dna5> > TStringSet;
    typedef Index<StringSet<TString>, FMIndex<> > TIndex;
    typedef Iterator<TIndex, TopDown<> >::Type TIter;

    // 0) Handle command line arguments.
    if (argc < 3) {
        std::cerr << "Invalid number of arguments." << std::endl
                  << "USAGE: minimapper GENOME.fasta READS.fasta OUT.sam" << std::endl;
        return 1;
    }
    // 1) Load contigs and reads.
    FragmentStore<> fragStore;
    if (!loadContigs(fragStore, argv[1])) return 1;
    if (!loadReads(fragStore, argv[2])) return 1;

    StringSet<TString> text;
    for (unsigned i = 0; i < length(fragStore.contigStore); ++i)
        appendValue(text, fragStore.contigStore[i].seq);
        
    TIndex fmIndex(text);
    TIter it(fmIndex);
    search(it, fragStore.readSeqStore);
    clear(fmIndex);
    clear(it);

    reverse(text);
    reverse(fragStore.readSeqStore);

    fmIndex = TIndex(text);
    it = TIter(fmIndex);
    search(it, fragStore.readSeqStore);
    clear(fmIndex);
    clear(it);

    return 0;
}
