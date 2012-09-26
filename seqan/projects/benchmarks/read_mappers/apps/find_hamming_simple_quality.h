/*
  Approximate Hamming-distance string search that takes scores into
  considerations.

  The needle sequence must support the getQualityValue() function for
  this.  Qualities are expected to be positive integers.
*/

// TODO(holtgrew): copy and paste code from HammingSimple!

#include <seqan/basic.h>
#include <seqan/find.h>
#include <seqan/modifier.h>
#include <seqan/sequence.h>

#ifndef WIT_BUILDER_HAMMING_SIMPLE_QUALITY_H_
#define WIT_BUILDER_HAMMING_SIMPLE_QUALITY_H_

namespace seqan {

struct _HammingSimpleQuality;
typedef Tag<_HammingSimpleQuality> HammingSimpleQuality;


template <typename TNeedle>
class Pattern<TNeedle, HammingSimpleQuality> {
public:
    // The holder for the needle.
    Holder<TNeedle> data_host;

    // The maximal distance.  Must be >= 0, i.e. -score.
    int maxDistance;

    // The current distance, >= 0, i.e. -current score.
    int distance;

    Pattern() : maxDistance(-1), distance(0) {}

    template <typename TNeedle2>
    Pattern(const TNeedle2 &ndl, int k = -1) {
        SEQAN_CHECKPOINT;
        setHost(*this, ndl, k);
    }
};


template <typename TNeedle, typename TNeedle2>
void setHost (Pattern<TNeedle, HammingSimpleQuality> & me, 
              const TNeedle2 & needle, int k) {
    SEQAN_CHECKPOINT;

    SEQAN_ASSERT_NOT(empty(needle));
    SEQAN_ASSERT_LEQ_MSG(k, 0, "Are you confusing distances and scores?");

    setValue(me.data_host, needle);
    me.maxDistance = -k;
}


template <typename TNeedle, typename TNeedle2>
void
setHost(Pattern<TNeedle, HammingSimpleQuality> &horsp, TNeedle2 &ndl, int k) {
    SEQAN_CHECKPOINT;
    setHost(horsp, reinterpret_cast<const TNeedle2&>(ndl), k);
}


template <typename TNeedle>
inline void _finderInit(Pattern<TNeedle, HammingSimpleQuality> & me) {
    SEQAN_CHECKPOINT;
    (void) me;  // Suppress unused variable warning.
}


template <typename TNeedle>
inline int score(const Pattern<TNeedle, HammingSimpleQuality> &me) {
    SEQAN_CHECKPOINT;
    return -me.distance;
}


template <typename TNeedle>
inline int getScore(const Pattern<TNeedle, HammingSimpleQuality> &me) {
    SEQAN_CHECKPOINT;
    return -me.distance;
}


template <typename TNeedle>
inline void setScoreLimit(Pattern<TNeedle, HammingSimpleQuality> & me, int _limit) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_LEQ(_limit, 0);
    me.maxDistance = -_limit;
}


template <typename TFinder, typename TNeedle>
inline bool find(TFinder &finder, 
                 Pattern<TNeedle, HammingSimpleQuality> &me) {
    SEQAN_CHECKPOINT;

    typedef typename Haystack<TFinder>::Type THaystack;
    typedef typename Size<THaystack>::Type TSize;

    // Shortcuts to haystack and needle.
    const THaystack &hstk = haystack(finder);
    const TNeedle &ndl = needle(me);

    // If the needle is longer than the haystack then we cannot find anything.
    if (length(hstk) < length(ndl))
        return false;

    // Initialize or advance finder, depending whether it has been
    // initialized before.
    if (empty(finder)) {
        _finderInit(me);
        _setFinderLength(finder, length(needle(me)));
        _finderSetNonEmpty(finder);
    } else {
        finder += 1;
    }

    // Check whether we are beyond the last possible match position.
    if (position(finder) > length(hstk) - length(ndl))
        return false;

    // TODO(holtgrew): Switch from indices to iterators to improve performance.

    // Perform a naive search for the needle in the haystack such that
    // the difference is <= me.maxDistance.
    TSize i;
    for (i = position(finder); i <= length(hstk) - length(ndl); ++i) {
        me.distance = 0;  // Reset mismatch count.
        for (TSize j = 0; j < length(ndl); ++j) {
            me.distance += (ndl[j] != hstk[i + j]) * getQualityValue(ndl[j]);
            if (me.distance > me.maxDistance)
                break;
        }
        if (me.distance <= me.maxDistance)
            break;
    }

    // Return false if we did not break out of the for-loop but it
    // stopped normally.
    if (i > length(hstk) - length(ndl))
        return false;

    _setFinderEnd(finder, i + length(ndl));
    setPosition(finder, beginPosition(finder));
    return true; 
}


// Set the end position of the pattern in the finder.
template <typename THaystack, typename TNeedle, typename TPosition>
inline bool setEndPosition(Finder<THaystack, void> & finder,
                           Pattern<TNeedle, HammingSimpleQuality> & pattern,
                           const TPosition & pos) {
//     std::cerr << "setEndPosition(finder, pattern, " << pos << ")" << std::endl;
    // Compute delta, such that we start searching at pos - delta.
    TPosition delta = length(needle(pattern));
    if (delta > pos)
        delta = pos;
//     std::cerr << "delta == " << delta << std::endl;

    // Set end position in the finder to pos - delta.
    finder.data_length = length(needle(pattern));
    setPosition(finder, pos - delta);
    finder.data_endPos = pos - delta;
//     std::cerr << "beginPosition(finder) == " << beginPosition(finder) << std::endl;
//     std::cerr << "endPosition(finder) == " << endPosition(finder) << std::endl;

    // Clear the pattern, and search until we are at pos.
    bool result;
    while ((result = find(finder, pattern)) && endPosition(finder) < pos) {
//         std::cerr << "Skipping over end pos " << endPosition(finder) << std::endl;
        continue;
    }
//     std::cerr << "XXX beginPosition(finder) == " << beginPosition(finder) << std::endl;
//     std::cerr << "XXX endPosition(finder) == " << endPosition(finder) << std::endl;
    return result;
}

}  // namespace seqan

#endif  // WIT_BUILDER_HAMMING_SIMPLE_QUALITY_H_
