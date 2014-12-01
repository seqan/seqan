// ==========================================================================
//                      RABEMA Read Alignment Benchmark
// ==========================================================================
// Copyright (C) 2010-1012 Manuel Holtgrewe, FU Berlin
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License along
// with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// ==========================================================================
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Defines the FindMyersUkkonenReads pattern.  The aim is to provide a version
// of the Myers-Ukkonen pattern that creates alignments that make sense for
// read mapping: We do not want to find matches where there is a delete at the
// end of the pattern when finding an approximate match.
//
// We use the following trick to achieve this:  Do an approximate search for
// the needle, excluding the last characters, using Myers-Ukkonen.  The
// maximal score is decremented by one.  Then, compare the last character of
// the original needle and the character behind the last matched one in the
// pattern.  Add the score of this comparison to the score of Myers-Ukkonen
// search.  This is the final score of the alignment that is guaranteed to
// end in a non-delete for the pattern.
//
// Note that the finder currently only works with very large absolute score
// values (i.e. < -(length(needle))).
// ==========================================================================

// TODO(holtgrew): Integrate into SeqAn (with tests) or give as an example on how to extend SeqAn in your own code.

#ifndef SEQAN_APPS_FIND_FIND_MYERS_UKKONEN_READS_H_
#define SEQAN_APPS_FIND_FIND_MYERS_UKKONEN_READS_H_

#include <seqan/sequence.h>
#include <seqan/find.h>

namespace seqan {

struct MyersUkkonenReads_;
typedef Tag<MyersUkkonenReads_> MyersUkkonenReads;

template <typename TNeedle>
class Pattern<TNeedle, MyersUkkonenReads>
{
public:
    // Shortcut to the MyersUkkonen-Pattern we will use internally.
    typedef Pattern<Segment<TNeedle, InfixSegment>, Myers<FindInfix> > TMyersUkkonen_;

    // The type of the needle's characters.
    typedef typename Value<TNeedle>::Type TAlphabet_;

    // The wrapped pattern.
    TMyersUkkonen_ _wrappedPattern;

    // The last character of the pattern.
    TAlphabet_ _lastCharacter;

    // The score of the match at the current position.
    int _score;

    // The score limit.
    int _scoreLimit;

    // A holder for the needle.
    Holder<TNeedle> data_host;

    // Boolean flag "find() has not been called".
    bool _firstFind;

    // Mask that determines whether N is matched against all/none.
    String<int> _matchNMask;

    // The whole needle, wrapped has one char less.
    //Holder<TNeedle> data_host;

    Pattern() :
        _score(0),
        _scoreLimit(0),
        _firstFind(true)
    {}

    template <typename TNeedle2>
    Pattern(TNeedle2 & ndl, int score = -1) :
        _wrappedPattern(infix(ndl, 0, length(ndl) - 1), score),
        _lastCharacter(back(ndl)),
        _score(score),
        _scoreLimit(score),
        data_host(ndl),
        _firstFind(true)
    {
        typedef typename Value<TNeedle>::Type TAlphabet;
        resize(_matchNMask, ValueSize<TAlphabet>::VALUE * ValueSize<TAlphabet>::VALUE);

        for (unsigned i = 0; i < ValueSize<TAlphabet>::VALUE; ++i)
        {
            for (unsigned j = 0; j < ValueSize<TAlphabet>::VALUE; ++j)
            {
                if (i == ordValue(TAlphabet('N')) || j == ordValue(TAlphabet('N')))
                    _matchNMask[ValueSize < TAlphabet > ::VALUE * i + j] = -1;
                else if (i == j)
                    _matchNMask[ValueSize < TAlphabet > ::VALUE * i + j] = 0;
                else
                    _matchNMask[ValueSize < TAlphabet > ::VALUE * i + j] = -1;
            }
        }
    }

};


template <typename TNeedle>
void _patternMatchNOfPattern(Pattern<TNeedle, MyersUkkonenReads> & me, bool match)
{
    _patternMatchNOfPattern(me._wrappedPattern, match);

    typedef typename Value<TNeedle>::Type TAlphabet;
    for (unsigned i = 0; i < ValueSize<TAlphabet>::VALUE; ++i)
        me._matchNMask[ValueSize < TAlphabet > ::VALUE * i + ordValue(TAlphabet('N'))] = match ? 0 : -1;
}

template <typename TNeedle>
void _patternMatchNOfFinder(Pattern<TNeedle, MyersUkkonenReads> & me, bool match)
{
    _patternMatchNOfFinder(me._wrappedPattern, match);

    typedef typename Value<TNeedle>::Type TAlphabet;
    for (unsigned i = 0; i < ValueSize<TAlphabet>::VALUE; ++i)
        me._matchNMask[ValueSize < TAlphabet > ::VALUE * ordValue(TAlphabet('N')) + i] = match ? 0 : -1;
}

template <typename TNeedle>
inline int
scoreLimit(Pattern<TNeedle, MyersUkkonenReads> const & me)
{
    // The score limit for us is "one more" than for the wrapped
    // pattern, meaning "one more negative."
    return scoreLimit(me._wrappedPattern) - 1;
}

template <typename TNeedle, typename TScoreValue>
inline void
setScoreLimit(Pattern<TNeedle, MyersUkkonenReads> & me,
              TScoreValue _limit)
{
    // The score limit for the wrapped pattern is "one less", meaning
    // "one less negative."
    setScoreLimit(me._wrappedPattern, _limit + 1);
}

template <typename TNeedle>
int getScore(const Pattern<TNeedle, MyersUkkonenReads> & me)
{
    return static_cast<int>(me._score);
}

template <typename TNeedle, typename TNeedle2>
inline void setHost(Pattern<TNeedle, MyersUkkonenReads> & me,
                    const TNeedle2 & needle)
{
    setHost(me._wrappedPattern, prefix(needle, length(needle) - 1));
    me._lastCharacter = back(needle);
}

template <typename TNeedle, typename TNeedle2>
inline void setHost(Pattern<TNeedle, MyersUkkonenReads> & me, TNeedle2 & needle)
{
    setHost(me, const_cast<TNeedle2 const &>(needle));
}

template <typename TNeedle>
inline void _patternInit(Pattern<TNeedle, MyersUkkonenReads> & me)
{
    _patternInit(me._wrappedPattern);
}

template <typename TNeedle>
inline int getScoreLimit(Pattern<TNeedle, MyersUkkonenReads> & me)
{
    return me._scoreLimit;
}

template <typename TFinder, typename TNeedle>
inline bool find(TFinder & finder, Pattern<TNeedle, MyersUkkonenReads> & me)
{
    // Go back in finder by one if this is not the first find call
    // because the wrapped pattern on the infix trick.
    if (me._firstFind)
    {
        me._firstFind = false;
    }
    else
    {
        // std::cout << "before go previous " << endPosition(finder);
        goPrevious(finder);
        _setFinderEnd(finder);
        // std::cout << ", after go previous " << endPosition(finder) << std::endl;
    }

    do
    {
        // Make sure we find something with the wrapped pattern.
        if (!find(finder, me._wrappedPattern))
            return false;

        // Make sure we are not at the end of the haystack with the
        // wrapped pattern yet.
        size_t endPos = endPosition(finder);
        if (endPos == length(haystack(finder)))
            return false;

        // Compute the current score.
        typedef typename Value<TNeedle>::Type TAlphabet;
        int lastCharScore = me._matchNMask[ordValue(me._lastCharacter) + ValueSize < TAlphabet > ::VALUE * ordValue(haystack(finder)[endPos])];
        me._score = getScore(me._wrappedPattern) + lastCharScore;
    }
    while (me._score < scoreLimit(me));

    // Advance finder by one, will wind back by one before calling
    // wrapped finder in next call.
    // std::cout << "before go next " << endPosition(finder);
    goNext(finder);
    _setFinderEnd(finder);
    // std::cout << ", after go next " << endPosition(finder) << std::endl;

    return true;
}

template <typename TFinder, typename TNeedle>
inline bool find(TFinder & finder,
                 Pattern<TNeedle, MyersUkkonenReads> & me,
                 int const k)
{
    setScoreLimit(me, k);
    return find(finder, me);
}

template <typename TFinder, typename TNeedle, typename TLimit>
inline bool findBegin(TFinder & finder,
                      Pattern<TNeedle, MyersUkkonenReads> & me,
                      TLimit limit)
{
    // Compute score of last character.
    size_t endPos = endPosition(finder) - 1;
    int lastCharScore = (me._lastCharacter == haystack(finder)[endPos]) ? 0 : -1;
    // Go back one with the finder because of the "last character is compared manually" trick.
    goPrevious(finder);
    finder.data_endPos -= 1;
    finder.data_length -= 1;
    // We have to incorporate the last character's score into limit
    // used for the wrapped pattern's findBegin().
    bool res = findBegin(finder, me._wrappedPattern, limit - lastCharScore);
    goNext(finder);
    finder.data_endPos += 1;
    finder.data_length += 1;
    return res;
}

template <typename TFinder, typename TNeedle>
inline bool findBegin(TFinder & finder,
                      Pattern<TNeedle, MyersUkkonenReads> & me)
{
    return findBegin(finder, me._wrappedPattern, getScoreLimit(me));
}

// Set the end position of the pattern in the finder.
template <typename THaystack, typename TNeedle, typename TPosition>
inline bool setEndPosition(Finder<THaystack, void> & finder,
                           Pattern<TNeedle, MyersUkkonenReads> & pattern,
                           const TPosition & pos)
{
    // Compute delta, such that we start searching at pos - delta.
    TPosition delta = length(needle(pattern)) + _min(length(needle(pattern)), static_cast<size_t>(-getScoreLimit(pattern)));
    if (delta > pos)
        delta = pos;

    // Set end position in the finder to pos - delta.
    setPosition(finder, pos - delta);

    // Clear the pattern, and search until we are at pos.
    _patternInit(pattern._wrappedPattern, finder);
    bool result;
    while ((result = find(finder, pattern)) &&
           endPosition(finder) < pos)
        continue;
    return result;
}

}  // namespace seqan

#endif  // SEQAN_APPS_FIND_MYERS_UKKONEN_READS_H_
