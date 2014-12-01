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
    bool isValid() const
    {
        return (strand != 'N');
    }

    // Mark as invalid.
    void setInvalid()
    {
        strand = 'N';
    }
};


}  // anonymous namespace

// --------------------------------------------------------------------------
// Class FilteringBestMapperImpl
// --------------------------------------------------------------------------

class FilteringBestMapperImpl
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

    FilteringBestMapperImpl(double errorRate) :
            rightClip(0), contigLength(0), errorRate(errorRate), numPHHits(0), numMatches(0)
    {}

    void run(std::vector<seqan::BamAlignmentRecord> & result,
             seqan::StringSet<seqan::Dna5String> & contigs,
             seqan::StringSet<seqan::Dna5String> /*const*/ & readSeqs,
             seqan::StringSet<seqan::CharString> /*const*/ & quals,
             seqan::StringSet<seqan::CharString> const & readNames)
    {
        if (empty(readSeqs) || empty(contigs))
            return;  // no matches possible

        // Prepare single end matches, we only need one per read for best mapping.
        std::vector<SingleEndMatch> matches(length(readSeqs));

        // Build pigeonhole pattern for the reads.
        TIndex index(readSeqs);
        TFilterPattern filterPattern(index);
        // Initialize shape according to error rate and already build suffix array.x
        _patternInit(filterPattern, 0.01 * errorRate);

        // Remove bad reads.
        maskBadReads(filterPattern, 0.01 * errorRate);
        // TODO(holtgrew): Mask against stretches of Ns.

        // Map against both strands of each contig, this fills matches with verified matches.
        for (unsigned contigID = 0; contigID < length(contigs); ++contigID)
        {
            seqan::Dna5String & contig = contigs[contigID];

            mapAgainstContig(matches, contig, filterPattern, contigID, /*forward=*/true);
            reverseComplement(contig);
            mapAgainstContig(matches, contig, filterPattern, contigID, /*forward=*/false);
            reverseComplement(contig);
        }

        // Convert matches into BAM alignment records in result.
        buildBamRecords(result, matches, contigs, readSeqs, quals, readNames);

        if (DEBUG)
            fprintf(stderr, "\n#reads: %u, PH hits: %u, Matches: %u\n", (unsigned)length(readSeqs),
                    numPHHits, numMatches);
    }

private:

    void maskBadReads(TFilterPattern & filterPattern, double errorRate)
    {
        auto const & reads = indexText(host(filterPattern));
        for (unsigned readID = 0; readID < length(reads); ++readID)
        {
            unsigned count = 0;
            for (auto c : reads[readID])
                count += (c == 'N');
            if (count > errorRate * length(reads[readID]))
            {
                if (DEBUG)
                    std::cerr << "MASKING OUT\t" << readID << "\t" << reads[readID] << "\n";
                maskPatternSequence(filterPattern, readID, false);
            }
        }
    }

    void buildBamRecords(std::vector<seqan::BamAlignmentRecord> & result,
                         std::vector<SingleEndMatch> const & matches,
                         seqan::StringSet<seqan::Dna5String> & contigs,
                         seqan::StringSet<seqan::Dna5String> /*const*/ & readSeqs,
                         seqan::StringSet<seqan::CharString> /*const*/ & quals,
                         seqan::StringSet<seqan::CharString> const & readNames) const
    {
        result.resize(length(readNames));

        seqan::Align<seqan::Dna5String> align;
        resize(rows(align), 2);
        seqan::Score<short, seqan::EditDistance> scoringScheme;

        for (unsigned readID = 0; readID < length(readNames); ++readID)
        {
            seqan::BamAlignmentRecord & record = result[readID];
            SingleEndMatch const & match = matches[readID];

            if (match.isValid() && match.strand == '-')  // use RC sequence in case of reverse match
                reverseComplement(readSeqs[readID]);

            clear(record);
            record.flag = seqan::BAM_FLAG_UNMAPPED;
            record._qID = readID;
            record.rID = match.contigID;
            record.qName = readNames[readID];
            record.seq = readSeqs[readID];
            record.qual = quals[readID];
            if (match.isValid())
            {
                record.flag = 0;
                record.rID = match.contigID;
                record.beginPos = match.beginPos;
                if (match.strand == '-')
                    record.flag |= seqan::BAM_FLAG_RC;

                assignSource(row(align, 0), infix(contigs[match.contigID], match.beginPos, match.endPos));
                assignSource(row(align, 1), readSeqs[readID]);
                int band = abs(match.score);
                if (!(-match.score == 0 || (-match.score == 1 && length(record.seq)) == length(source(row(align, 0)))))
                    globalAlignment(align, scoringScheme, -band - 1, band + 1);
                getCigarString(record.cigar, row(align, 0), row(align, 1));

                seqan::BamTagsDict tags(record.tags);
                setTagValue(tags, "NM", static_cast<__int32>(-match.score));
            }
        }
    }

    // Map against a single contig.
    void mapAgainstContig(std::vector<SingleEndMatch> & matches,
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
        while (find(filterFinder, filterPattern, 0.01 * errorRate))
        {
            ++numPHHits;
            int readID = position(filterPattern).i1;
            if (DEBUG)
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
            if (match.isValid() && (!matches[readID].isValid() || match.score > matches[readID].score))
            {
                ++numMatches;
                matches[readID] = match;
                // Mark read as ignored if perfect match was found;  this can save verifications.
                if (match.score == 0)
                    maskPatternSequence(filterPattern, readID, false);
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
        int minScore = -ndlLength * 0.01 * errorRate;

        TGenomeInfix inf = infix(filterFinder);
        TReadPrefix readPrefix(indexText(host(filterPattern))[readID], ndlLength);  // TODO(holtgrew): Should take full seq!
        TMyersFinder myersFinder(inf);

        if (DEBUG)
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
        if (DEBUG)
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

        if (DEBUG)
            std::cerr << "\nBEGIN\t" << match.beginPos << "\t" << match.endPos << "\t" << match.score << "\t"
                      << "STRAND=" << match.strand << "\n"
                      << "    GENOME=" << infix(host(infix(filterFinder)), match.beginPos, match.endPos) << "\n"
                      << "    READ=  " << readPrefix << "\n";
    }

    static const bool DEBUG = false;

    // Variables/state for the verification.
    TPatternState	patternState;
    TRPatternState  revPatternState;
    unsigned        rightClip;
    unsigned        contigLength;
    // The error rate to use for the mapping, in percent.
    double errorRate;
    // Number of PH hits and matches.
    unsigned numPHHits, numMatches;
};

// --------------------------------------------------------------------------
// Class FilteringBestMapper
// --------------------------------------------------------------------------

FilteringBestMapper::FilteringBestMapper(double errorRate) : impl(new FilteringBestMapperImpl(errorRate))
{}

FilteringBestMapper::~FilteringBestMapper()
{}

void FilteringBestMapper::run(std::vector<seqan::BamAlignmentRecord> & result,
                              seqan::StringSet<seqan::Dna5String> & contigs,
                              seqan::StringSet<seqan::Dna5String> /*const*/ & reads,
                              seqan::StringSet<seqan::CharString> /*const*/ & quals,
                              seqan::StringSet<seqan::CharString> const & readNames)
{
    impl->run(result, contigs, reads, quals, readNames);
}
