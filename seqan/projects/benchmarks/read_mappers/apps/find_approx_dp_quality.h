/*
  Approximate string search that takes scores into considerations.

  The needle sequence must support the getQualityValue() function for
  this.  Qualities are expected to be positive integers.
*/

// TODO(holtgrew): Banding of prefix search.

#include <seqan/basic.h>
#include <seqan/find.h>
#include <seqan/modifier.h>
#include <seqan/sequence.h>

#ifndef WIT_BUILDER_APPROX_DP_QUALITY_H_
#define WIT_BUILDER_APPROX_DP_QUALITY_H_

namespace seqan {

//.Tag.QualityDpSearch
//..param.TSpec: One of FindInfix, FindPrefix.
template <typename TSpec, bool HasFindBeginSupport = true>
struct QualityDpSearch {};


// Empty find begin pattern struct with empty base class optimization trick.
struct EmptyBase_ {};
struct EmptyFindBegin_ : EmptyBase_ {
    template <typename T1, typename T2> EmptyFindBegin_(T1 const &, T2 const &) {}
};


template <typename TNeedle, typename TSpec, bool HasFindBeginSupport>
class Pattern<TNeedle, QualityDpSearch<TSpec, HasFindBeginSupport> > {
public:
    typedef int TScoreValue;

    Holder<TNeedle> data_host;
    TScoreValue data_limit;
    String<TScoreValue> data_tab;
    TScoreValue data_maxscore;  // Score of the needle matching itself

    typedef typename If<HasFindBeginSupport, Pattern<ModifiedString<TNeedle, ModReverse>, QualityDpSearch<FindPrefix, false> >, EmptyFindBegin_>::Type TFindBeginPattern;
    TFindBeginPattern _findBeginPattern;

    Pattern() {}

    template <typename TNeedle2>
	Pattern(TNeedle2 & _needle, TScoreValue _limit) : data_limit(_limit), _findBeginPattern(_needle, _limit)
	{
        SEQAN_CHECKPOINT;
		setHost(*this, _needle);
	}

	Pattern(TScoreValue _limit) : data_limit(_limit)
	{
        SEQAN_CHECKPOINT;
	}
};


template <typename TNeedle, typename TSpec, bool HasFindBeginSupport>
inline typename Host<Pattern<TNeedle, QualityDpSearch<TSpec, HasFindBeginSupport> > >::Type &
host(Pattern<TNeedle, QualityDpSearch<TSpec, HasFindBeginSupport> > & me)
{
    SEQAN_CHECKPOINT;
	return value(me.data_host);
}


template <typename TNeedle, typename TSpec, bool HasFindBeginSupport>
inline typename Host<Pattern<TNeedle, QualityDpSearch<TSpec, HasFindBeginSupport> > const>::Type &
     host(Pattern<TNeedle, QualityDpSearch<TSpec, HasFindBeginSupport> > const & me)
{
    SEQAN_CHECKPOINT;
	return value(me.data_host);
}


template <typename TNeedle, typename TNeedle2, typename TSpec, bool HasFindBeginSupport>
void
setHost(Pattern<TNeedle, QualityDpSearch<TSpec, HasFindBeginSupport> > & me,
		TNeedle2 & ndl)
{
	me.data_host = ndl;
	clear(me.data_tab);
}


template <typename TNeedle, typename TNeedle2, typename TSpec, bool HasFindBeginSupport>
void
setHost(Pattern<TNeedle, QualityDpSearch<TSpec, HasFindBeginSupport> > & me,
		TNeedle2 const & ndl)
{
	me.data_host = ndl;
	clear(me.data_tab);
}


template <typename TNeedle, typename TSpec, bool HasFindBeginSupport>
inline int
scoreLimit(Pattern<TNeedle, QualityDpSearch<TSpec, HasFindBeginSupport> > const & me)
{
    SEQAN_CHECKPOINT;
	return me.data_limit;
}


template <typename TNeedle, typename TSpec, bool HasFindBeginSupport>
inline void
setScoreLimit(Pattern<TNeedle, QualityDpSearch<TSpec, HasFindBeginSupport> > & me,
			  int _limit)
{
    SEQAN_CHECKPOINT;
	me.data_limit = _limit;
}


template <typename TNeedle, typename TSpec, bool HasFindBeginSupport>
inline int
getScore(Pattern<TNeedle, QualityDpSearch<TSpec, HasFindBeginSupport> > & me)
{
	return back(me.data_tab);
}


template <typename TNeedle, typename TSpec, bool HasFindBeginSupport>
inline void _patternInit(Pattern<TNeedle, QualityDpSearch<TSpec, HasFindBeginSupport> > & me)
{
	typedef Pattern<TNeedle, QualityDpSearch<FindInfix, HasFindBeginSupport> > TPattern;
	typedef int TScoreValue;
	typedef typename Size<TPattern>::Type TSize;
	typedef String<TScoreValue> TTab;
	typedef typename Iterator<TTab, Standard>::Type TIterator;
	typedef typename Iterator<TNeedle, Standard>::Type TNeedleIterator;

    SEQAN_ASSERT_GT_MSG(length(needle(me)), 0u, "Empty needles not allowed.");

    TNeedle const & ndl = needle(me);

    // Resize matrix column and fill with default values.
    resize(me.data_tab, length(ndl));
    me.data_tab[0] = -getQualityValue(ndl[0]);
    for (TSize i = 1; i < length(ndl); ++i)
        me.data_tab[i] = me.data_tab[i - 1] - getQualityValue(ndl[i]);
}


template <typename TFinder, typename TNeedle, typename TSpec, bool HasFindBeginSupport>
inline bool
find(TFinder & finder,
	 Pattern<TNeedle, QualityDpSearch<TSpec, HasFindBeginSupport> > & me)
{
    typedef typename Haystack<TFinder>::Type THaystack;
	typedef typename Size<THaystack>::Type TSize;
    typedef typename Value<THaystack>::Type THaystackValue;

    TSize prefixBeginPosition;  // For prefix search.

    if (empty(finder)) {
        _patternInit(me);
        _finderSetNonEmpty(finder);
        prefixBeginPosition = position(finder);
    } else {
        goNext(finder);
        prefixBeginPosition = beginPosition(finder);
    }
    // TODO(holtgrew): This is HACKY. Prefix search will ONLY work correctly for non-begin search if has find begin support is turned on. Maybe this bug is in DPSearch, too?
    if (HasFindBeginSupport)
        prefixBeginPosition = 0;

//     {
//         std::cerr << "me.data_tab == ";
//         for (unsigned i = 0; i < length(me.data_tab); ++i) {
//             std::cerr << me.data_tab[i] << " ";
//         }
//         std::cerr << std::endl;
//     }

    TSize haystack_length = length(container(hostIterator(finder)));
    TNeedle const & ndl = needle(me);
    for (; position(finder) < haystack_length; ++finder) {
        THaystackValue c = *finder;

        // Update the column...
        int northWest = me.data_tab[0];  // "Western" fields are one step back.
        // First row.
        if (IsSameType<TSpec, FindInfix>::VALUE) {
            me.data_tab[0] = (ndl[0] == *finder) ? 0 : -getQualityValue(ndl[0]);
        } else {
//             std::cerr << "prefixBeginPosition = " << prefixBeginPosition << std::endl;
            int const north = -(position(finder) + 1 - prefixBeginPosition) * getQualityValue(ndl[0]);
            int const west = me.data_tab[0];
            if (length(ndl) > 1)
                me.data_tab[0] = west - _max(1, static_cast<int>(0.5 * (getQualityValue(ndl[0]) + getQualityValue(ndl[1]))));
            else
                me.data_tab[0] = west - getQualityValue(ndl[0]);
            me.data_tab[0] = _max(me.data_tab[0], north - getQualityValue(ndl[0]));
            me.data_tab[0] = _max(me.data_tab[0], -static_cast<int>(position(finder) - prefixBeginPosition) * getQualityValue(ndl[0]) + ((ndl[0] == *finder) ? 0 : -getQualityValue(ndl[0])));
        }
        // All but the last row.
        TSize i;
        for (i = 1; i + 1 < length(ndl); ++i) {
            int const & north = me.data_tab[i - 1];
            int const & west = me.data_tab[i];

            // Compute value for the i-th field of the current column.
            int newValue = north - getQualityValue(ndl[i]);
            newValue = _max(newValue, west - _max(1, static_cast<int>(0.5 * (getQualityValue(ndl[i]) + getQualityValue(ndl[i + 1])))));
            newValue = _max(newValue, northWest - getQualityValue(ndl[i]) * (ndl[i] != *finder));

            // Save current value for next iteration, then write new value into the field.
            northWest = me.data_tab[i];
            me.data_tab[i] = newValue;
        }
        if (length(ndl) > 1) {
            // Last row.
            int const north = me.data_tab[i - 1];
            int const west = me.data_tab[i];
            back(me.data_tab) = north - getQualityValue(ndl[i]);
            back(me.data_tab) = _max(back(me.data_tab), west - _max(1, getQualityValue(ndl[i])));
            back(me.data_tab) = _max(back(me.data_tab), northWest - getQualityValue(ndl[i]) * (ndl[i] != *finder));
        }
//         {
//             std::cerr << "me.data_tab == ";
//             for (unsigned i = 0; i < length(me.data_tab); ++i) {
//                 std::cerr << me.data_tab[i] << " ";
//             }
//             std::cerr << std::endl;
//         }

        // Check whether we can report a hit.
        if (getScore(me) >= scoreLimit(me)) {
			_setFinderEnd(finder);
            return true;
        }
    }

    // Report that there are no more hits.
    return false;
}


template <typename THaystack, typename TFinderSpec, typename TNeedle, typename TSpec>
inline
bool findBegin(Finder<THaystack, TFinderSpec> & finder,
               Pattern<TNeedle, QualityDpSearch<TSpec, true> > & pattern,
               int scoreLimit) {
    typedef Pattern<TNeedle, QualityDpSearch<TSpec, true> > TPattern;
    typedef typename TPattern::TFindBeginPattern TFindBeginPattern;
    typedef ModifiedString<THaystack, ModReverse> TReverseHaystack;
    typedef Finder<TReverseHaystack, TFinderSpec> TFindBeginFinder;
    typedef typename Position<THaystack>::Type TPosition;

    TFindBeginPattern & findBeginPattern = pattern._findBeginPattern;
    setScoreLimit(findBeginPattern, scoreLimit);

    // Build beginFinder.
    TFindBeginFinder beginFinder;
    THaystack & hstck = haystack(finder);
    setContainer(host(hostIterator(beginFinder)), hstck);
    TPosition beginFinderBeginPosition = position(finder);
    TPosition beginFinderPosition;

//     std::cerr << "reverse needle = " << needle(findBeginPattern) << std::endl;
//     std::cerr << "forward haystack = " << haystack(finder) << std::endl;

    if (!finder._beginFind_called) {
        // Start finding.
        finder._beginFind_called = true;
        _setFinderLength(finder, 0);
        clear(beginFinder);
        beginFinderPosition = beginFinderBeginPosition;
    } else {
        // Resume finding.
        _finderSetNonEmpty(beginFinder);
        SEQAN_ASSERT_GT(length(finder), 0u);
        beginFinderPosition = endPosition(finder) - length(finder);
    }

    setPosition(host(hostIterator(beginFinder)), beginFinderPosition);
    _setFinderEnd(beginFinder);
    _setFinderLength(beginFinder, length(finder));

    bool beginFound = find(beginFinder, findBeginPattern);
    if (beginFound) // New begin found:  Report in finder.
        _setFinderLength(finder, beginFinderBeginPosition - position(host(hostIterator(beginFinder))) + 1);

    return beginFound;
}


template <typename THaystack, typename TFinderSpec, typename TNeedle, typename TSpec>
inline
bool findBegin(Finder<THaystack, TFinderSpec> & finder,
               Pattern<TNeedle, QualityDpSearch<TSpec, true> > & pattern) {
    SEQAN_CHECKPOINT;
    return findBegin(finder, pattern, getScoreLimit(pattern));
}


template <typename TNeedle, typename TSpec>
inline
int getFindBeginScore(Pattern<TNeedle, QualityDpSearch<TSpec, true> > & pattern) {
    SEQAN_CHECKPOINT;
    return getScore(pattern._findBeginPattern);
}


template <typename THaystack, typename TFinderSpec, typename TNeedle, typename TSpec, bool HasFindBeginSupport, typename TPosition>
inline bool setEndPosition(Finder<THaystack, TFinderSpec> & finder,
                           Pattern<TNeedle, QualityDpSearch<TSpec, HasFindBeginSupport> > & pattern,
                           const TPosition & pos) {
    // Compute delta, such that we start searching at pos - delta.
    int minQualityValue = 1000;
    for (unsigned i = 0; i < length(needle(pattern)); ++i) {
        minQualityValue = _min(minQualityValue, getQualityValue(needle(pattern)[i]));
    }
    TPosition x = ceil(-1.0 * scoreLimit(pattern) / minQualityValue);
    TPosition delta = length(needle(pattern)) + _min(length(needle(pattern)), x) + 1;
    if (delta > pos)
        delta = pos;

    // Set end position in the finder to pos - delta.
    setPosition(finder, pos - delta);
    if (pos == delta)
        clear(finder);

    // Clear the pattern, and search until we are at pos.
    _patternInit(pattern);
    bool result;
    while ((result = find(finder, pattern)) and endPosition(finder) < pos)
        continue;
    return result;
}

}  // namespace seqan

#endif  // WIT_BUILDER_APPROX_DP_QUALITY_H_
