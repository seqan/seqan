#ifndef BENCHMARKS_READ_MAPPERS_FILE_HELPERS_H_
#define BENCHMARKS_READ_MAPPERS_FILE_HELPERS_H_

#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/sequence.h>


template <typename TString>
inline
void loadFileIntoStringSet_assignQualities_doAssign(TString & sequence, CharString const & qualities) {
    // Convert ascii to values from 0..62 store dna and quality
    // together in Dna5Q.
    for (unsigned j = 0; j < length(qualities) && j < length(sequence); ++j)
        assignQualityValue(sequence[j], static_cast<int>(ordValue(qualities[j]) - 33));
}

template <typename TStringSpec>
inline
void loadFileIntoStringSet_assignQualities(String<Dna5Q, TStringSpec> & sequence,
                                           CharString const & qualities) {
    // Forward to function that actually performs loading.
    loadFileIntoStringSet_assignQualities_doAssign(sequence, qualities);
}


template <typename TStringSpec>
inline
void loadFileIntoStringSet_assignQualities(String<DnaQ, TStringSpec> & sequence,
                                           CharString const & qualities) {
    // Forward to function that actually performs loading.
    loadFileIntoStringSet_assignQualities_doAssign(sequence, qualities);
}


template <typename TAlphabet, typename TStringSpec>
inline
void loadFileIntoStringSet_assignQualities(String<TAlphabet, TStringSpec> &,
                                           CharString const &) {
    // No quality assignment for non-DNA5Q.
}


/**
.Function.truncateAfterFirst
..signature:truncateAfterFirst(str, c)
..param.str:String to truncate.
..param.c:Alphabet character to truncate after.
..summary:Truncate a string after the first occurence of c.
 */
template <typename TString>
void truncateAfterFirst(TString & str, typename Value<TString>::Type const & c) {
    typedef typename Size<TString>::Type TSize;
    typedef typename Iterator<TString, Standard>::Type TIterator;
    TSize firstPos = 0;
    for (TIterator it = begin(str, Standard()); it != end(str, Standard()); ++it, ++firstPos) {
        if (*it == c)
            break;
    }
    if (firstPos == length(str))
        return;
    resize(str, firstPos, Exact());
}


template <typename TSequenceStringSet, typename TNameStringSet>
void loadFileIntoStringSet(TSequenceStringSet & sequences, TNameStringSet & names, const CharString & filename) {
    typedef typename Value<TSequenceStringSet>::Type TSequence;
    
    MultiSeqFile multiSeqFile;
    if (!open(multiSeqFile.concat, toCString(filename), OPEN_RDONLY)) {
        std::cerr << "Could not open file " << filename << std::endl;
        exit(1);
    }

    AutoSeqFormat format;
    guessFormat(multiSeqFile.concat, format);
    split(multiSeqFile, format);

    unsigned seqCount = length(multiSeqFile);
    reserve(sequences, seqCount, Exact());
    reserve(names, seqCount, Exact());

    TSequence seq;
    CharString qual;
    CharString id;
    for (unsigned i = 0; i < seqCount; ++i) {
        assignSeq(seq, multiSeqFile[i], format);    // read sequence
        assignQual(qual, multiSeqFile[i], format);  // read ascii quality values
        assignSeqId(id, multiSeqFile[i], format);   // read sequence id
        truncateAfterFirst(id, ' ');

        loadFileIntoStringSet_assignQualities(seq, qual);

        // we use reserve and append, as assign is not supported
        // by StringSet<..., Owner<ConcatDirect<> > >
        appendValue(sequences, seq, Generous());
        appendValue(names, id, Generous());
    }
}

#endif  // BENCHMARKS_READ_MAPPERS_FILE_HELPERS_H_
