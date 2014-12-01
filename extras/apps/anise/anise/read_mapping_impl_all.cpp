// ==========================================================================
//                                   ANISE
// ==========================================================================
// Copyright (c) 2006-2014, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

#include "read_mapping_impl.h"

#include <iostream>
#include <iterator>

#include <seqan/bam_io.h>
#include <seqan/index.h>
#include <seqan/index/find_pigeonhole.h>

namespace {  // anonymous namespace

// --------------------------------------------------------------------------
// Class SingleEndMatch
// --------------------------------------------------------------------------

// Represent a single-end match.

struct SingleEndMatch
{
    unsigned readID;
    char strand;  // '+', '-', 'N'
    __int32 contigID;
    __int32 beginPos;
    __int32 endPos;
    __int32 score;  // Hamming/Levenshtein distance if negative/score if positive.

    SingleEndMatch() :
            readID(seqan::maxValue<unsigned>()), strand('N'), contigID(-1), beginPos(0), endPos(0), score(0)
    {}

    SingleEndMatch(unsigned _readID, char _strand, __int32 _contigID, __int32 _beginPos,
                   __int32 _endPos, __int32 _score) :
            readID(_readID), strand(_strand), contigID(_contigID), beginPos(_beginPos), endPos(_endPos),
            score(_score)
    {}

    // Returns true if this match was assigned.
    bool isValid() const { return (strand != 'N'); }

    // Mark as invalid.
    void setInvalid() { strand = 'N'; }

    static bool ltReadID(SingleEndMatch const & lhs, SingleEndMatch const & rhs)
    { return (lhs.readID < rhs.readID); }

    static bool ltReadIDScore(SingleEndMatch const & lhs, SingleEndMatch const & rhs)
    { return (std::make_pair(lhs.readID, -lhs.score) < std::make_pair(rhs.readID, -rhs.score)); }

    static bool ltCoord(SingleEndMatch const & lhs, SingleEndMatch const & rhs)
    { return (std::make_pair(lhs.contigID, lhs.beginPos) < std::make_pair(rhs.contigID, rhs.beginPos)); }
};

std::ostream & operator<<(std::ostream & out, SingleEndMatch const & match)
{
    return out << "SingleEndMatch(readID=" << match.readID << ", strand=" << match.strand
               << ", contigID=" << match.contigID << ", beginPos=" << match.beginPos
               << ", endPos=" << match.endPos << ", score=" << match.score << ")";
}

// --------------------------------------------------------------------------
// Class CompactingStore
// --------------------------------------------------------------------------

// Wrapper for std::vector<> that allows to grow with hysteresis and filtration on resize.

template <typename T>
class CompactingStore
{
public:
    // Next cleanup count.
    size_t nextCleanup { 1 };
    // Hysteresis factor.
    double factor { 1.5 };
    // The stored elements.
    std::vector<T> elements;
    // The function to use for cleaning up.
    std::function<void(std::vector<T>&)> cleanupFunc;

    CompactingStore() : cleanupFunc([](std::vector<T>&) {}) {}
    CompactingStore(std::function<void(std::vector<T>&)> fun) : cleanupFunc(fun) {}

    void push_back(T const & elem) { elements.push_back(elem); checkCleanup(); }

    typename std::vector<T>::const_iterator begin() const { return elements.begin(); }
    typename std::vector<T>::const_iterator end()   const { return elements.end();   }

private:

    void checkCleanup()
    {
        if (elements.size() < nextCleanup)
            return;
        cleanupFunc(elements);
        nextCleanup *= factor;
    }
};

// --------------------------------------------------------------------------
// Function cleanupMatchStore()
// --------------------------------------------------------------------------

// Helper function used for cleaning up MatchStore::matches.

void cleanupMatchStore(std::vector<SingleEndMatch> & /*matches*/)
{
    // if (matches.empty())
    //     return;  // nothing to do

    // // std::cerr << "CLEANING UP STORE (#" << matches.size() << "\n";
    // // std::copy(matches.begin(), matches.end(), std::ostream_iterator<SingleEndMatch>(std::cerr, "\n"));

    // // Sort by (readID, edit distance), scanning over this then allows us to purge the non-optimal matches.
    // std::sort(matches.begin(), matches.end(), SingleEndMatch::ltReadIDScore);

    // auto it = matches.begin(), itNext = matches.begin() + 1;
    // while (itNext != matches.end())
    //     if (it->readID != itNext->readID)  // write out *itNext
    //     {
    //         ++it;
    //         *it = *itNext++;
    //     }
    //     else  // write out *itNext if score is equal to current (the best)
    //     {
    //         if (it->score == itNext->score)
    //         {
    //             ++it;
    //             *it = *itNext;
    //         }
    //         ++itNext;
    //     }
    // ++it;  // keep in any case
    // matches.resize(it - matches.begin());

    // // std::cerr << "DONE CLEANING UP STORE (#" << matches.size() << "\n";
    // // std::copy(matches.begin(), matches.end(), std::ostream_iterator<SingleEndMatch>(std::cerr, "\n"));
}

// --------------------------------------------------------------------------
// Class MatchStore
// --------------------------------------------------------------------------

// Manages SingleEndMatches

class MatchStore
{
public:
    struct Entry
    {
        // Number of matches for this read, regardless of error rate.
        unsigned count { 0 };

        Entry() = default;
        Entry(unsigned count) : count(count) {}
    };

    // Maximal number of matches
    unsigned maxMatches { 0 };
    // Number of reads.
    unsigned numReads { 0 };
    // Lowest match error rate for each read and the number of matches with this error count.
    std::vector<Entry> matchInfos;
    // Buffer of SingleEndMatch records.
    CompactingStore<SingleEndMatch> matches;

    size_t size() const { return matches.elements.size(); }
    std::vector<SingleEndMatch>::const_iterator begin() const { return matches.begin(); }
    std::vector<SingleEndMatch>::const_iterator end()   const { return matches.end();   }

    MatchStore(unsigned maxMatches, unsigned numReads) :
            maxMatches(maxMatches), numReads(numReads), matchInfos(numReads), matches(cleanupMatchStore)
    {}

    // Returns the match information for the given read.
    Entry matchInfo(unsigned readID) const
    {
        SEQAN_ASSERT_LT(readID, numReads);
        return matchInfos[readID];
    }

    // Try to register a new match with the given error rate for the given read.  If this match is of sufficient quality
    // then return true and add register it in the store (maybe resetting the count).  Otherwise, return false.
    bool registerMatch(unsigned readID)
    {
        SEQAN_ASSERT_LT(readID, numReads);
        if (limitReached(readID))
            return false;
        matchInfos[readID].count += 1;
        return true;
    }

    bool registerMatch(SingleEndMatch const & match)
    {
        if (DEBUG)
            std::cerr << "Registering match\t" << match << "\n"
                      << "\t" << limitReached(match.readID) << "\n";
        if (!registerMatch(match.readID))
            return false;  // match not to be registered
        matches.push_back(match);
        return true;
    }

    // Force running cleanup.
    void runCleanup() { cleanupMatchStore(matches.elements); }

    // Print store.
    void print(std::ostream & out) const
    {
        out << "MATCH STORE\n"
            << "INFO\n";
        for (unsigned i = 0; i < matchInfos.size(); ++i)
            out << i << ":\tcount=" << matchInfos[i].count << "\n";
        out << "MATCHES\n";
        for (auto const & match : matches)
        out << "\t" << match << "\n";
    }

private:

    // Returns whether the limit for the given read has been reached or not.
    bool limitReached(unsigned readID) const
    {
        SEQAN_ASSERT_LT(readID, matchInfos.size());
        return (matchInfos[readID].count >= maxMatches);
    }


    static const bool DEBUG = false;
};

}  // anonymous namespace

// --------------------------------------------------------------------------
// Class FilteringAllMapperImpl
// --------------------------------------------------------------------------

class FilteringAllMapperImpl
{
    // Typedes for the filtration.
    typedef seqan::StringSet<seqan::Dna5String>                    TReadSet;
    typedef seqan::Shape<seqan::Dna5, seqan::OneGappedShape>       TShape;
    typedef seqan::IndexQGram<TShape, seqan::OpenAddressing>       TIndexSpec;
    typedef seqan::Index<TReadSet, TIndexSpec>                     TIndex;
    typedef seqan::Pattern<TIndex, seqan::Pigeonhole<>>            TFilterPattern;
    typedef seqan::Finder<seqan::Dna5String, seqan::Pigeonhole<>>  TFilterFinder;

    // Typedefs for the verification.
    typedef typename seqan::Value<TReadSet>::Type /*const*/        TReadSeq;
    typedef typename seqan::Prefix<TReadSeq>::Type                 TReadPrefix;
    typedef seqan::ModifiedString<TReadPrefix, seqan::ModReverse>  TRevReadPrefix;

    typedef typename seqan::Infix<seqan::Dna5String>::Type         TGenomeInfix;
    typedef typename seqan::Position<TGenomeInfix>::Type           TPosition;
    typedef seqan::ModifiedString<TGenomeInfix, seqan::ModReverse> TGenomeInfixRev;
    typedef seqan::Finder<TGenomeInfix>                            TMyersFinder;
    typedef seqan::Finder<TGenomeInfixRev>                         TMyersFinderRev;
    typedef seqan::AlignTextBanded<seqan::FindInfix,
                                   seqan::NMatchesNone_,
                                   seqan::NMatchesNone_>            TAlignTextBanded;
    typedef seqan::PatternState_<TReadPrefix,
                                 seqan::Myers<TAlignTextBanded,
                                              seqan::True, void> >  TPatternState;
    typedef seqan::PatternState_<TRevReadPrefix,
                                 seqan::Myers<TAlignTextBanded,
                                              seqan::True, void> >  TRPatternState;

public:
    typedef FilteringAllMapper::Options Options;

    FilteringAllMapperImpl(Options options) :
    rightClip(0), contigLength(0), options(options), numPHHits(0), numMatches(0)
    {}

    void run(std::vector<seqan::BamAlignmentRecord> & result,
             seqan::StringSet<seqan::Dna5String> & contigs,
             seqan::StringSet<seqan::Dna5String> /*const*/ & readSeqs,
             seqan::StringSet<seqan::CharString> /*const*/ & quals,
             seqan::StringSet<seqan::CharString> const & readNames)
    {
        if (options.debug)
            std::cerr << "RUNNING ALL-MAPPING WITH ERROR RATE\t" << options.errorRate << "\n";

        if (empty(readSeqs) || empty(contigs))
            return;  // no matches possible

        // Setup match store.
        MatchStore matchStore(options.maxMatches, length(readSeqs));

        // Build pigeonhole pattern for the reads.
        TIndex index(readSeqs);
        TFilterPattern filterPattern(index);
        // Initialize shape according to error rate and already build q-gram index.
        _patternInit(filterPattern, options.errorRate);

        // Remove bad reads.
        maskBadReads(filterPattern, options.errorRate);

        // Map against both strands of each contig, this fills matches with verified matches.
        for (unsigned contigID = 0; contigID < length(contigs); ++contigID)
        {
            seqan::Dna5String & contig = contigs[contigID];

            mapAgainstContig(matchStore, contig, filterPattern, contigID, /*forward=*/true);
            reverseComplement(contig);
            mapAgainstContig(matchStore, contig, filterPattern, contigID, /*forward=*/false);
            reverseComplement(contig);
        }

        // Run cleanup on matchStore to get rid of suboptimal matches.
        matchStore.runCleanup();

        // Convert matches into BAM alignment records in result.
        buildBamRecords(result, matchStore, contigs, readSeqs, quals, readNames);

        if (options.debug)
            fprintf(stderr, "\n#reads: %u, PH hits: %u, Matches: %u\n", (unsigned)length(readSeqs),
                    numPHHits, numMatches);
    }

private:

    // Mask out read with too many Ns.
    void maskBadReads(TFilterPattern & filterPattern, double /*errorRate*/)  // error rate not required for PH
    {
        auto const & reads = indexText(host(filterPattern));
        for (unsigned readID = 0; readID < length(reads); ++readID)
        {
            unsigned count = 0;
            for (auto c : reads[readID])
                count += (c == 'N');
            if (count > options.errorRate * length(reads[readID]))
            {
                if (options.debug)
                    std::cerr << "MASKING OUT\t" << readID << "\t" << reads[readID] << "\n";
                maskPatternSequence(filterPattern, readID, false);
            }
        }
    }

    void buildBamRecords(std::vector<seqan::BamAlignmentRecord> & result,
                         MatchStore const & matchStore,
                         seqan::StringSet<seqan::Dna5String> & contigs,
                         seqan::StringSet<seqan::Dna5String> /*const*/ & readSeqs,
                         seqan::StringSet<seqan::CharString> /*const*/ & quals,
                         seqan::StringSet<seqan::CharString> const & readNames) const
    {
        result.resize(matchStore.size());

        seqan::Align<seqan::Dna5String> align;
        resize(rows(align), 2);
        seqan::Score<short, seqan::EditDistance> scoringScheme;

        for (auto itM = matchStore.begin(); itM != matchStore.end(); ++itM)
        {
            auto & record = result[itM - matchStore.begin()];
            clear(record);
            auto const & match = *itM;

            record.flag = 0;
            record._qID = match.readID;
            record.qName = readNames[match.readID];
            record.seq = readSeqs[match.readID];
            if (match.strand == '-')  // use RC sequence in case of reverse match
                reverseComplement(record.seq);
            record.qual = quals[match.readID];
            reverse(record.qual);

            record.rID = match.contigID;
            record.beginPos = match.beginPos;
            if (match.strand == '-')
                record.flag |= seqan::BAM_FLAG_RC;

            assignSource(row(align, 0), infix(contigs[match.contigID], match.beginPos, match.endPos));
            assignSource(row(align, 1), readSeqs[match.readID]);
            if (match.strand == '-')  // use RC sequence in case of reverse match
                reverseComplement(source(row(align, 1)));
            int band = abs(match.score);
            if (!(-match.score == 0 || (-match.score == 1 && length(record.seq)) == length(source(row(align, 0)))))
                globalAlignment(align, scoringScheme, -band - 1, band + 1);
            getCigarString(record.cigar, row(align, 0), row(align, 1));

            seqan::BamTagsDict tags(record.tags);
            setTagValue(tags, "NM", static_cast<__int32>(-match.score));
        }
    }

    // Map against a single contig.
    void mapAgainstContig(MatchStore & matchStore,
                          seqan::Dna5String & contig,
                          TFilterPattern & filterPattern,
                          int contigID,
                          bool forward)
    {
        SingleEndMatch match;  // match stored in verification

        // Update state.
        rightClip = 0;
        contigLength = length(contig);

        // Enumerate the pigeonhole hits.
        TFilterFinder filterFinder(contig, /*minRepeatLength=*/1000, /*periodSize=*/1);
        while (find(filterFinder, filterPattern, options.errorRate))
        {
            ++numPHHits;
            int readID = position(filterPattern).i1;
            if (options.debug)
                std::cerr << "PH hit for read=" << readID << ", contig=" << contigID << ", strand="
                          << (forward ? '+' : '-') << ", beginPos="
                          << (forward ? beginPosition(infix(filterFinder)) :
                              contigLength - beginPosition(infix(filterFinder))) << ", endPos="
                          << (forward ? endPosition(infix(filterFinder)) :
                              contigLength - endPosition(infix(filterFinder))) << "\n"
                          << "  CONTIG\t" << infix(filterFinder) << "\n"
                          << "  READ  \t" << indexText(host(filterPattern))[readID] << "\n";

            // Get id of read as shortcut and verify candidate.
            verifyMatch(match, filterFinder, filterPattern, contigID, forward);
            // If the candidate could be verified then update it in case of better score.
            if (match.isValid())
            {
                matchStore.registerMatch(match);
                if (options.debug)
                    std::cerr << "\t=> VALID " << match << "\n";
                if (false && options.debug)
                    matchStore.print(std::cerr);
            }
        }
    }

    void verifyMatch(SingleEndMatch & match,
                     TFilterFinder & filterFinder,
                     TFilterPattern & filterPattern,
                     int contigID,
                     bool isForward)
    {
        match.setInvalid();  // mark as invalid

        // Compute clipping for banded Myers if match it is left/right of the text.
        patternState.leftClip = (beginPosition(filterFinder) >= 0) ? 0 : -beginPosition(filterFinder);
        rightClip = (endPosition(filterFinder) <= contigLength) ? 0 : endPosition(filterFinder) - contigLength;

        // Verify pigeonhole hit.
        //
        // We use a much simpler version of verification in RazerS 3.

        // First, we search for the most promising end position.
        unsigned readID = position(filterPattern).i1;
        int ndlLength = sequenceLength(readID, indexText(host(filterPattern)));
        int minScore = -ndlLength * options.errorRate;

        TGenomeInfix inf = infix(filterFinder);
        TReadPrefix readPrefix(indexText(host(filterPattern))[readID], ndlLength);  // TODO(holtgrew): Should take full seq!
        TMyersFinder myersFinder(inf);

        if (false && options.debug)
            std::cerr << "\nVERIFYING\t" << readPrefix << "\t" << inf << "\n";

        int bestScore = seqan::minValue<int>();
        unsigned bestPos = 0;
        // Enumerate all start positions to find best.
        int score = 0;
        while (find(myersFinder, readPrefix, patternState, minScore))
            if ((score = getScore(patternState)) > bestScore)
            {
                bestScore = score;
                bestPos = position(hostIterator(myersFinder));
            }

        if (bestScore == seqan::minValue<int>())
            return;  // No Myers hit, look for next pigeonhole hit.

        // Second, search for the leftmost start position for the hit.
        auto infEndPos = endPosition(inf);
        auto newInfEndPos = beginPosition(inf) + bestPos + 1;
        revPatternState.leftClip = infEndPos - newInfEndPos + rightClip;
        setEndPosition(inf, newInfEndPos);
        if (endPosition(inf) > (unsigned)(ndlLength - bestScore))
            setBeginPosition(inf, endPosition(inf) - ndlLength + bestScore);
        else
            setBeginPosition(inf, 0);

        TRevReadPrefix readRev(readPrefix);
        TGenomeInfixRev infRev(inf);
        TMyersFinderRev myersFinderRev(infRev);
        if (false && options.debug)
            std::cerr << "\nVERIFYING REV\t" << readPrefix << "\t" << inf << "\n"
                      << "   => " << readRev << "\t" << infRev << "\n";
        auto beginPos = newInfEndPos;
        while (find(myersFinderRev, readRev, revPatternState, bestScore))
            beginPos = newInfEndPos - position(myersFinderRev);
        if (beginPos == newInfEndPos)
            return;  // No banded Myers hit, skip.
        if (beginPos > (decltype(beginPos))0)
            --beginPos;

        // Write out match.
        match.readID = readID;
        match.strand = isForward ? '+' : '-';
        match.contigID = contigID;
        match.beginPos = isForward ? beginPos : contigLength - newInfEndPos;
        match.endPos = isForward ? newInfEndPos : contigLength - beginPos;
        match.score = bestScore;

        SEQAN_CHECK(match.beginPos <= match.endPos, "Invalid match!");

        if (options.debug)
            std::cerr << "\nBEGIN\t" << match.beginPos << "\t" << match.endPos << "\t" << match.score << "\t"
                      << "STRAND=" << match.strand << "\n"
                      << "    GENOME=" << infix(host(infix(filterFinder)), match.beginPos, match.endPos) << "\n"
                      << "    READ=  " << readPrefix << "\n";
    }

    // Variables/state for the verification.
    TPatternState	patternState;
    TRPatternState  revPatternState;
    unsigned        rightClip;
    unsigned        contigLength;
    // The configuration.
    Options options;
    // Number of PH hits and matches.
    unsigned numPHHits, numMatches;
};

// --------------------------------------------------------------------------
// Class FilteringBestMapper
// --------------------------------------------------------------------------

FilteringAllMapper::FilteringAllMapper() : impl(new FilteringAllMapperImpl(FilteringAllMapper::Options()))
{}

FilteringAllMapper::FilteringAllMapper(FilteringAllMapper::Options options) : impl(new FilteringAllMapperImpl(options))
{}

FilteringAllMapper::~FilteringAllMapper()
{}

void FilteringAllMapper::run(std::vector<seqan::BamAlignmentRecord> & result,
                                 seqan::StringSet<seqan::Dna5String> & contigs,
                                 seqan::StringSet<seqan::Dna5String> /*const*/ & reads,
                                 seqan::StringSet<seqan::CharString> /*const*/ & quals,
                                 seqan::StringSet<seqan::CharString> const & readNames)
{
    impl->run(result, contigs, reads, quals, readNames);
}
