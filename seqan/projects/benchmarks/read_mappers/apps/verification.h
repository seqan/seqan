/*
  Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>

  Code for the verification of the wit_builder results.
 */

#ifndef WIT_BUILDER_VERIFICATION_H_
#define WIT_BUILDER_VERIFICATION_H_

#include <seqan/modifier.h>
#include <seqan/store.h>

#include "curve_smoothing.h"
#include "wit_builder_options.h"
#include "find_myers_ukkonen_reads.h"
#include "find_myers_ukkonen_ext.h"
#include "find_approx_dp_quality.h"
#include "find_hamming_simple_quality.h"

using namespace seqan;


// For each read, an error curve is a list of weighted matches.
typedef std::map<int, TWeightedMatches> TErrorCurves;


// Represents an interval in a contig.
struct ContigInterval {
    // Contig the interval is in.
    size_t contigId;

    // true iff the interval is on the forward strand.
    bool isForward;

    // First value of interval
    size_t first;

    // Last value of interval.
    size_t last;

    // Default constructor so it can be used in containers.
    ContigInterval() {}

    // Constructor for the record.
    ContigInterval(size_t _contigId, bool _isForward, size_t _first, size_t _last)
            : contigId(_contigId), isForward(_isForward), first(_first), last(_last) {}
};


// Function to use for finding reads with a finder only -- MyersUkkonenReads version.
//
// We search for the read (resp. reverse-complemented read), excluding
// the last character.  We then compute the final score by computing
// the match score between the last characters.  This enforces the
// last character of the read to actually align with the string.
//
// There also is FindMyersUkkonenReads pattern but we do it manually
// here.  This seems safer for verification purposes.
template <typename TString, typename TContigIdType>
void verifyMatchestoErrorFunctionResults_FindReads(
        TString /*const*/ & contig,
        bool const & isForward,
        TString /*const*/ & read,
        const Options &options,
        int maxError,
        TContigIdType contigId,
        String<WeightedMatch> &foundMatches,
        const MyersUkkonenReads &) {
    Finder<TString> finder(contig);
    Pattern<Segment<TString, InfixSegment>, Myers<FindInfix> > pattern(infix(read, 0, length(read) - 1));
    setScoreLimit(pattern, -maxError);
    // std::cout << read << std::endl; 
    while (find(finder, pattern)) {
        // std::cout << "endpos = " << endPosition(finder) << std::endl;
        if (endPosition(finder) >= length(contig))
            continue;  // Skip if aligning beyond the contig.
        int score = getScore(pattern);
        // Now, apply the cost for a mismatch at the end
        // of the read.
        score -= back(read) != contig[endPosition(finder)];
        // std::cout << "contig = " << contig << ", score = " << score << ", endpos = " << endPosition(finder) << std::endl;
        // Skip this hit if the score is not good enough.
        if (score < -maxError)
            continue;
        int relativeScore = (int)ceilAwayFromZero(100.0 * score / length(read));
        SEQAN_ASSERT_GEQ(relativeScore, -options.maxError);
        (void)options;  // Supress warnings in non-debug mode.
        // std::cout << "append value(expected, WeightedMatch(" << contigId << ", " << endPosition(finder) + 1 << ", " << relativeScore << ") delta = " << (back(reversedAndComplementedRead) != contig[endPosition(finder) - 1]) << ", score " << score << std::endl;
        bool ret = findBegin(finder, pattern, getScore(pattern));  // Compute begin position for smoothing.
        (void)ret;  // Supress warning in non-debug mode.
        SEQAN_ASSERT(ret);
        appendValue(foundMatches, WeightedMatch(contigId, isForward, endPosition(finder), relativeScore, beginPosition(finder)));
    }
}


// MyersUkkonen (*not* Reads) version.
template <typename TString, typename TContigIdType>
void verifyMatchestoErrorFunctionResults_FindReads(
        TString /*const*/ & contig,
        bool const & isForward,
        TString /*const*/ & read,
        const Options & options,
        int maxError,
        TContigIdType contigId,
        String<WeightedMatch> &foundMatches,
        const Myers<FindInfix> &) {
    Finder<TString> finder(contig);
    Pattern<TString, Myers<FindInfix> > pattern(read);
    setScoreLimit(pattern, -maxError);
    while (find(finder, pattern)) {
        if (endPosition(finder) >= length(contig))
            continue;  // Skip if aligning beyond the contig.
        int score = getScore(pattern);
        if (score < -maxError)
            continue;
        int relativeScore = ceilAwayFromZero(100.0 * score / length(read));
        SEQAN_ASSERT_GEQ(relativeScore, -options.maxError);
        (void)options;  // Supress warnings in non-debug mode.
        bool ret = findBegin(finder, pattern, getScore(pattern));  // Compute begin position for smoothing.
        (void)ret;  // Supress warning in non-debug mode.
        SEQAN_ASSERT(ret);
        appendValue(foundMatches, WeightedMatch(contigId, isForward, endPosition(finder) - 1, relativeScore, beginPosition(finder)));
    }
}


// Approximate DP search -- quality based version.
template <typename TString, typename TContigIdType>
void verifyMatchestoErrorFunctionResults_FindReads(
        TString /*const*/ & contig,
        bool const & isForward,
        TString /*const*/ & read,
        const Options &options,
        int maxError,
        TContigIdType contigId,
        String<WeightedMatch> & foundMatches,
        const QualityDpSearch<FindInfix> &) {
    Finder<TString> finder(contig);
    Pattern<TString, QualityDpSearch<FindInfix> > pattern(read, -maxError);
    //    setScoreLimit(pattern, -maxError);
    while (find(finder, pattern)) {
        if (endPosition(finder) >= length(contig))
            continue;  // Skip if aligning beyond the contig.
        int score = getScore(pattern);
        if (score < -maxError)
            continue;
        int relativeScore = ceilAwayFromZero(100.0 * score / length(read));
        SEQAN_ASSERT_GEQ(relativeScore, -options.maxError);
        bool ret = findBegin(finder, pattern, getScore(pattern));  // Compute begin position for smoothing.
        SEQAN_ASSERT(ret);
        appendValue(foundMatches, WeightedMatch(contigId, isForward, endPosition(finder) - 1, relativeScore, beginPosition(finder)));
    }
}


// Function to use for finding reads with a finder only -- Hamming Simple Finder version.
template <typename TString, typename TContigIdType>
void verifyMatchestoErrorFunctionResults_FindReads(
        TString /*const*/ & contig,
        bool const & isForward,
        TString /*const*/ & read,
        const Options & options,
        int maxError,
        TContigIdType contigId,
        String<WeightedMatch> &foundMatches,
        const HammingSimple &) {
    Finder<TString> finder(contig);
    Pattern<TString, HammingSimple> pattern(read);
    setScoreLimit(pattern, -maxError);
    // std::cout << read << std::endl; 
    while (find(finder, pattern)) {
        int score = getScore(pattern);
        if (score < -maxError)
            continue;
        int relativeScore = (int)ceilAwayFromZero(100.0 * score / length(read));
        SEQAN_ASSERT_GEQ(relativeScore, -options.maxError);
        (void)options;  // Supress warnings in non-debug mode.
        bool ret = findBegin(finder, pattern, getScore(pattern));  // Compute begin position for smoothing.
        (void)ret;  // Supress warning in non-debug mode.
        SEQAN_ASSERT(ret);
        appendValue(foundMatches, WeightedMatch(contigId, isForward, endPosition(finder) - 1, relativeScore, beginPosition(finder)));
    }
}


// Function to use for finding reads with a finder only -- Hamming Simple Finder With Qualities version.
template <typename TString, typename TContigIdType>
void verifyMatchestoErrorFunctionResults_FindReads(
        TString /*const*/ & contig,
        bool const & isForward,
        TString /*const*/ & read,
        const Options &options,
        int maxError,
        TContigIdType contigId,
        String<WeightedMatch> &foundMatches,
        const HammingSimpleQuality &) {
    Finder<TString> finder(contig);
    Pattern<TString, HammingSimpleQuality> pattern(read);
    setScoreLimit(pattern, -maxError);
    // std::cout << read << std::endl; 
    while (find(finder, pattern)) {
        int score = getScore(pattern);
        if (score < -maxError)
            continue;
        int relativeScore = ceilAwayFromZero(100.0 * score / length(read));
        SEQAN_ASSERT_GEQ(relativeScore, -options.maxError);
        bool ret = findBegin(finder, pattern, getScore(pattern));  // Compute begin position for smoothing.
        SEQAN_ASSERT(ret);
        appendValue(foundMatches, WeightedMatch(contigId, isForward, endPosition(finder) - 1, relativeScore, beginPosition(finder)));
    }
}


// Verify the result of matchesToErrorFunction().  Writes to stderr
// for logging and on errors.  Returns true iff there were no errors.
//
// Only use this when testing on small samples: We will perform a DP
// search for all reads on all contigs with all error values.
//
// The verification is performed by simply using the Myers<FindInfix>
// approximate search algorithm for each read over all contigs.
template <typename TFragmentStore, typename TPatternSpec>
bool verifyMatchesToErrorFunctionResults(TFragmentStore /*const*/ & fragments,
                                         TErrorCurves const & errorCurves,
                                         Options const & options,
                                         TPatternSpec const &) {
    bool valid = true;
    std::cerr << "Verifying error curves." << std::endl;

    typedef typename TFragmentStore::TAlignedReadStore                   TAlignedReadStore;
    typedef typename TFragmentStore::TContigSeq                          TContigSeq;
    typedef typename TFragmentStore::TReadSeq                            TReadSeq;
    typedef typename TFragmentStore::TContigStore                        TContigStore;
    typedef typename TFragmentStore::TReadStore                          TReadStore;
    typedef typename TFragmentStore::TReadSeqStore                       TReadSeqStore;
    typedef typename Value<TReadStore>::Type                             TRead;
    typedef typename TRead::TId                                          TReadId;
    typedef typename Value<TContigStore>::Type                           TContig;
    typedef typename TContig::TId                                        TContigId;
    typedef Gaps<TContigSeq, AnchorGaps<typename TContig::TGapAnchors> > TContigGaps;

    TContigStore /*const*/ & contigs = fragments.contigStore;
    TReadStore /*const*/ & reads = fragments.readStore;
    TReadSeqStore /*const*/ & readSeqs = fragments.readSeqStore;

    // For all contigs:
    //   For all reads:
    //     For direction in {forward, reverse}:
    //       Search $read on the $direction strand of $contig.
    //       Fill gaps in gathered point list and smooth curve.
    //       Compare error curve points found above to the error curve for the read
    //           from errorCurves.

    for (TContigId contigId = 0; contigId < length(contigs); ++contigId) {
        std::cerr << "  contig " << contigId << "/" << length(contigs) << std::endl;
        TContigSeq /*const*/ & contig = contigs[contigId].seq;
        // Get reverse-complement of the contig.
        TContigSeq rcContig(contig);
        reverseComplement(rcContig);

        for (TReadId readId = 0; readId < length(reads); ++readId) {
            std::cerr << "    read " << readId << "/" << length(reads) << std::endl;
            String<WeightedMatch> expected;
            TReadSeq /*const*/ read = readSeqs[readId];
            int maxError = (int)floor(options.maxError / 100.0 * length(read));

            // Search read in contig and reverse complemented contig.
            verifyMatchestoErrorFunctionResults_FindReads(contig, true, read, options, maxError, contigId, expected, TPatternSpec());
            verifyMatchestoErrorFunctionResults_FindReads(rcContig, false, read, options, maxError, contigId, expected, TPatternSpec());

            // Fill the gaps in the "expected" error curve and smooth it so it
            // really is what we expect.
            std::sort(begin(expected, Standard()), end(expected, Standard()));
            fillGaps(expected);
            smoothErrorCurve(expected);

            // Get a string of "result" matches, i.e. the error curve
            // points for the current read on the current contig.
            SEQAN_ASSERT(errorCurves.find(readId) != errorCurves.end());
            String<WeightedMatch> const & allResults = errorCurves.find(readId)->second;
            String<WeightedMatch> result;
            for (size_t i = 0; i < length(allResults); ++i) {
                if (allResults[i].contigId != contigId) continue;
                appendValue(result, allResults[i]);
            }

            // Compare resulting and expected weighted matches.  We do this by
            // computing the set differences.
            //
            // Now get the sequences of superflous and missing weighted matches.
            String<WeightedMatch> superflous, missing;
            resize(superflous, _max(length(result), length(expected)));
            resize(missing, _max(length(result), length(expected)));
            typedef typename Iterator<String<WeightedMatch> >::Type TIter;
            TIter end_superflous, end_missing;
            end_superflous = std::set_difference(begin(result), end(result),
                                                 begin(expected), end(expected),
                                                 begin(superflous));
            end_missing = std::set_difference(begin(expected), end(expected),
                                              begin(result), end(result),
                                              begin(missing));
            // And check.
            if (end_missing != begin(missing) || end_superflous != begin(superflous)) {
                std::cerr << "read id = " << readId << ", read name = " << fragments.readNameStore[readId] << std::endl;
                std::cerr << "Matches are:" << std::endl;
                for (TIter it = begin(result); it != end(result); ++it)
                    std::cerr << value(it) << std::endl;
                std::cerr << "Matches should be:" << std::endl;
                for (TIter it = begin(expected); it != end(expected); ++it)
                    std::cerr << value(it) << std::endl;
            }
            if (end_superflous != begin(superflous)) {
                valid = false;
                std::cerr << "VERIFICATION ERROR: Superflous weighted matches:" << std::endl;
                for (TIter it = begin(superflous); it != end_superflous; ++it) {
                    TReadSeq const & read = fragments.readSeqStore[readId];
                    TContigSeq const & contig = contigs[it->contigId].seq;
                    typedef typename Value<TContigSeq>::Type TAlphabet;
                    typedef ModifiedString<ModifiedString<TContigSeq, ModView<FunctorComplement<TAlphabet> > >, ModReverse> TModRevContigSeq;
                    TModRevContigSeq rcContig(contig);
                    std::cerr << "read: " << read << std::endl;
                    std::cerr << "contigId: " << it->contigId << (it->isForward ? " (F)" : " (R)") << std::endl;
                    std::cerr << "contig length: " << length(contig) << std::endl;
                    if (it->isForward)
                        std::cerr << "contig[it->pos - length(read):it->pos] == " << infix(contig, it->pos - length(read), it->pos) << std::endl;
                    else
                        std::cerr << "rcContig[it->pos - length(read):it->pos] == " << infix(rcContig, it->pos - length(read), it->pos) << std::endl;
                    std::cerr << value(it) << std::endl;
                }
            }
            if (end_missing != begin(missing)) {
                valid = false;
                std::cerr << "VERIFICATION ERROR: Missing weighted matches:" << std::endl;
                for (TIter it = begin(missing); it != end_missing; ++it) {
                    TReadSeq const & read = fragments.readSeqStore[readId];
                    TContigSeq const & contig = contigs[it->contigId].seq;
                    typedef typename Value<TContigSeq>::Type TAlphabet;
                    typedef ModifiedString<ModifiedString<TContigSeq, ModView<FunctorComplement<TAlphabet> > >, ModReverse> TModRevContigSeq;
                    TModRevContigSeq rcContig(contig);
                    std::cerr << "read: " << read << std::endl;
                    std::cerr << "contigId: " << it->contigId << (it->isForward ? " (F)" : " (R)") << std::endl;
                    std::cerr << "contig length: " << length(contig) << std::endl;
                    if (it->isForward)
                        std::cerr << "contig[it->pos - length(read):it->pos] == " << infix(contig, it->pos - length(read), it->pos) << std::endl;
                    else
                        std::cerr << "rcContig[it->pos - length(read):it->pos] == " << infix(rcContig, it->pos - length(read), it->pos) << std::endl;
                    std::cerr << value(it) << std::endl;
                }
            }
            if (end_superflous - begin(superflous) != 0 ||
                end_missing - begin(missing) != 0)
                std::cerr << "-----" << std::endl;
        }
    }

    std::cerr << "Done verifying error curves." << std::endl;
    return valid;
}

#endif  // WIT_BUILDER_VERIFICATION_H_
