/*
  Usage: seq_sample N IN.fastq OUT.fastq
         seq_sample N IN.1.fastq IN.2.fastq OUT.1.fastq OUT.2.fastq

  Samples N reads from IN.fastq and writes them to OUT.fastq without
  duplicates, alternatively from mate-paired reads.
 */

#include <fstream>
#include <iostream>
#include <set>

#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/misc/misc_random.h>
#include <seqan/sequence.h>
#include <seqan/store.h>
#include <seqan/random.h>

using namespace seqan;


template <typename TStream, typename TStringSetSpec, typename TStringSetSpec2>
void write(TStream & stream,
           StringSet<String<Dna5Q>, TStringSetSpec> const & sequences,
           StringSet<CharString, TStringSetSpec2> const & seqIds,
           Fastq const &) {
    typedef StringSet<String<Dna5Q>, TStringSetSpec> TStringSet;
    typedef typename Position<TStringSet>::Type TPosition;

    CharString qualBuffer;
    for (TPosition i = 0; i < length(sequences); ++i) {
        stream << "@" << seqIds[i] << std::endl;
        stream << sequences[i] << std::endl;
        stream << "+" << seqIds[i] << std::endl;
        resize(qualBuffer, length(sequences[i]), Exact());
        for (TPosition j = 0; j < length(sequences[i]); ++j) {
            qualBuffer[j] = getQualityValue(sequences[i][j]) + '!';
        }
        stream << qualBuffer << std::endl;
    }
}


int main(int argc, char **argv) {
    if (argc != 4 && argc != 6) {
        std::cerr << "ERROR: Wrong number of parameters!" << std::endl;
        std::cerr << "USAGE: seq_sample N IN.fastq OUT.fastq" << std::endl;
        std::cerr << "       seq_sample N IN.1.fastq IN.2.fastq OUT.1.fastq OUT.2.fastq" << std::endl;
        return 1;
    }

    const unsigned SEED = 42;
	Rng<MersenneTwister> rng(SEED);

    unsigned n = atoi(argv[1]);

    // Load reads using the fragment store.
    FragmentStore<> fragStore;
    bool hasMatePairs = (argc == 6);
    unsigned matePairFactor = hasMatePairs ? 2 : 1;
    if (!hasMatePairs) {
        if (!loadReads(fragStore, argv[2])) return 1;
    } else {
        if (!loadReads(fragStore, argv[2], argv[3])) return 1;
    }
    unsigned matePairCount = length(fragStore.readSeqStore) / matePairFactor;

    // The factor 10 is arbitrarily chosen to make sure the generation
    // is sufficiently fast.
    if (matePairCount < 2 * n) {
        std::cerr << "Number of reads/mate pairs in file/s should be >= 2*n" << std::endl;
        return 1;
    }

    // Randomly build the resulting string set.
    // Pick ids to select.
    std::set<unsigned> selectedIds;
    while (length(selectedIds) < n) {
		Pdf<Uniform<unsigned> > pdf(0, matePairCount - 1);
		unsigned x = pickRandomNumber(rng, pdf);
        selectedIds.insert(x);
    }
    // Actually build result.
    StringSet<String<Dna5Q> > resultSeqsL, resultSeqsR;
    StringSet<CharString> resultIdsL, resultIdsR;
    reserve(resultSeqsL, n);
    reserve(resultIdsL, n);
    if (hasMatePairs) {
        reserve(resultIdsR, n);
        reserve(resultSeqsR, n);
    }
    for (std::set<unsigned>::const_iterator it = selectedIds.begin(); it != selectedIds.end(); ++it) {
        if (!hasMatePairs) {
            appendValue(resultSeqsL, fragStore.readSeqStore[*it]);
            appendValue(resultIdsL, fragStore.readNameStore[*it]);
        } else {
            appendValue(resultSeqsL, fragStore.readSeqStore[2 * (*it)]);
            appendValue(resultIdsL, fragStore.readNameStore[2 * (*it)]);
            appendValue(resultSeqsR, fragStore.readSeqStore[2 * (*it) + 1]);
            appendValue(resultIdsR, fragStore.readNameStore[2 * (*it) + 1]);
        }
    }

    // Write out the result.
    if (!hasMatePairs) {
        if (CharString("-") == argv[3]) {
            write(std::cout, resultSeqsL, resultIdsL, Fastq());
        } else {
            std::fstream fstrm(argv[3], std::ios_base::out);
            if (! fstrm.is_open()) {
                std::cerr << "Could not open out file \"" << argv[3] << "\"" << std::endl;
                return 1;
            }
            write(fstrm, resultSeqsL, resultIdsL, Fastq());
        }
    } else {
        if (CharString("-") == argv[4]) {
            write(std::cout, resultSeqsL, resultIdsL, Fastq());
        } else {
            std::fstream fstrm(argv[4], std::ios_base::out);
            if (! fstrm.is_open()) {
                std::cerr << "Could not open out file \"" << argv[4] << "\"" << std::endl;
                return 1;
            }
            write(fstrm, resultSeqsL, resultIdsL, Fastq());
        }
        if (CharString("-") == argv[5]) {
            write(std::cout, resultSeqsR, resultIdsR, Fastq());
        } else {
            std::fstream fstrm(argv[5], std::ios_base::out);
            if (! fstrm.is_open()) {
                std::cerr << "Could not open out file \"" << argv[5] << "\"" << std::endl;
                return 1;
            }
            write(fstrm, resultSeqsR, resultIdsR, Fastq());
        }
    }

    return 0;
}
