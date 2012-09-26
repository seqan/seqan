/* Tests for find_approx_dp_quality.h. */

#include <seqan/basic.h>
#include <seqan/find.h>
#include "find_approx_dp_quality.h"

using namespace seqan;

// A simple test, search for "AGGA" in "TAGAGA" with simple,
// increasing qualities, high limit.
SEQAN_DEFINE_TEST(test_dp_quality_find_infix_simple_high_limit) {
    String<Dna5> hstck = "TAGAGA";
    String<Dna5Q> ndl = "AGGA";
    assignQualityValue(ndl[0], 5);
    assignQualityValue(ndl[1], 10);
    assignQualityValue(ndl[2], 15);
    assignQualityValue(ndl[3], 20);

    Finder<String<Dna5> > finder(hstck);
    Pattern<String<Dna5Q>, QualityDpSearch<FindInfix, false> > pattern(ndl, -1000);

    // The following has been verified on paper :)

    bool ret = find(finder, pattern);
    SEQAN_ASSERT(ret);
    SEQAN_ASSERT_EQ(1u, endPosition(finder));
    SEQAN_ASSERT_EQ(-50, getScore(pattern));

    ret = find(finder, pattern);
    SEQAN_ASSERT(ret);
    SEQAN_ASSERT_EQ(2u, endPosition(finder));
    SEQAN_ASSERT_EQ(-30, getScore(pattern));

    ret = find(finder, pattern);
    SEQAN_ASSERT(ret);
    SEQAN_ASSERT_EQ(3u, endPosition(finder));
    SEQAN_ASSERT_EQ(-30, getScore(pattern));

    ret = find(finder, pattern);
    SEQAN_ASSERT(ret);
    SEQAN_ASSERT_EQ(4u, endPosition(finder));
    SEQAN_ASSERT_EQ(-10, getScore(pattern));

    ret = find(finder, pattern);
    SEQAN_ASSERT(ret);
    SEQAN_ASSERT_EQ(5u, endPosition(finder));
    SEQAN_ASSERT_EQ(-30, getScore(pattern));

    ret = find(finder, pattern);
    SEQAN_ASSERT(ret);
    SEQAN_ASSERT_EQ(6u, endPosition(finder));
    SEQAN_ASSERT_EQ(-10, getScore(pattern));

    ret = find(finder, pattern);
    SEQAN_ASSERT_NOT(ret);
}


// A simple test, search for "AGGA" in "TAGAGA" with simple,
// increasing qualities, not so high limit.
SEQAN_DEFINE_TEST(test_dp_quality_find_infix_simple) {
    String<Dna5> hstck = "TAGAGA";
    String<Dna5Q> ndl = "AGGA";
    assignQualityValue(ndl[0], 5);
    assignQualityValue(ndl[1], 10);
    assignQualityValue(ndl[2], 15);
    assignQualityValue(ndl[3], 20);

    Finder<String<Dna5> > finder(hstck);
    Pattern<String<Dna5Q>, QualityDpSearch<FindInfix> > pattern(ndl, -15);

    // The following has been verified on paper :)

    bool ret = find(finder, pattern);
    SEQAN_ASSERT(ret);
    SEQAN_ASSERT_EQ(4u, endPosition(finder));
    SEQAN_ASSERT_EQ(-10, getScore(pattern));

    ret = find(finder, pattern);
    SEQAN_ASSERT(ret);
    SEQAN_ASSERT_EQ(6u, endPosition(finder));
    SEQAN_ASSERT_EQ(-10, getScore(pattern));

    ret = find(finder, pattern);
    SEQAN_ASSERT_NOT(ret);
}


// A simple test, search for "AAG" in "TAGAGA" with "random" scores,
// high limit.
SEQAN_DEFINE_TEST(test_dp_quality_find_infix_not_increasing_high_limit) {
    String<Dna5> hstck = "TAGAGA";
    String<Dna5Q> ndl = "AAG";
    assignQualityValue(ndl[0], 5);
    assignQualityValue(ndl[1], 10);
    assignQualityValue(ndl[2], 2);

    Finder<String<Dna5> > finder(hstck);
    Pattern<String<Dna5Q>, QualityDpSearch<FindInfix> > pattern(ndl, -1000);

    // The following has been verified on paper :)

    bool ret = find(finder, pattern);
    SEQAN_ASSERT(ret);
    SEQAN_ASSERT_EQ(1u, endPosition(finder));
    SEQAN_ASSERT_EQ(-17, getScore(pattern));

    ret = find(finder, pattern);
    SEQAN_ASSERT(ret);
    SEQAN_ASSERT_EQ(2u, endPosition(finder));
    SEQAN_ASSERT_EQ(-7, getScore(pattern));

    ret = find(finder, pattern);
    SEQAN_ASSERT(ret);
    SEQAN_ASSERT_EQ(3u, endPosition(finder));
    SEQAN_ASSERT_EQ(-5, getScore(pattern));

    ret = find(finder, pattern);
    SEQAN_ASSERT(ret);
    SEQAN_ASSERT_EQ(4u, endPosition(finder));
    SEQAN_ASSERT_EQ(-7, getScore(pattern));

    ret = find(finder, pattern);
    SEQAN_ASSERT(ret);
    SEQAN_ASSERT_EQ(5u, endPosition(finder));
    SEQAN_ASSERT_EQ(-5, getScore(pattern));

    ret = find(finder, pattern);
    SEQAN_ASSERT(ret);
    SEQAN_ASSERT_EQ(6u, endPosition(finder));
    SEQAN_ASSERT_EQ(-7, getScore(pattern));

    ret = find(finder, pattern);
    SEQAN_ASSERT_NOT(ret);
}


// A simple test, search for "AAG" in "TAGAGA" with "random" scores,
// not so high limit.
SEQAN_DEFINE_TEST(test_dp_quality_find_infix_not_increasing) {
    String<Dna5> hstck = "TAGAGA";
    String<Dna5Q> ndl = "AAG";
    assignQualityValue(ndl[0], 5);
    assignQualityValue(ndl[1], 10);
    assignQualityValue(ndl[2], 2);

    Finder<String<Dna5> > finder(hstck);
    Pattern<String<Dna5Q>, QualityDpSearch<FindInfix> > pattern(ndl, -7);

    // The following has been verified on paper :)

    bool ret = find(finder, pattern);
    SEQAN_ASSERT(ret);
    SEQAN_ASSERT_EQ(2u, endPosition(finder));
    SEQAN_ASSERT_EQ(-7, getScore(pattern));

    ret = find(finder, pattern);
    SEQAN_ASSERT(ret);
    SEQAN_ASSERT_EQ(3u, endPosition(finder));
    SEQAN_ASSERT_EQ(-5, getScore(pattern));

    ret = find(finder, pattern);
    SEQAN_ASSERT(ret);
    SEQAN_ASSERT_EQ(4u, endPosition(finder));
    SEQAN_ASSERT_EQ(-7, getScore(pattern));

    ret = find(finder, pattern);
    SEQAN_ASSERT(ret);
    SEQAN_ASSERT_EQ(5u, endPosition(finder));
    SEQAN_ASSERT_EQ(-5, getScore(pattern));

    ret = find(finder, pattern);
    SEQAN_ASSERT(ret);
    SEQAN_ASSERT_EQ(6u, endPosition(finder));
    SEQAN_ASSERT_EQ(-7, getScore(pattern));

    ret = find(finder, pattern);
    SEQAN_ASSERT_NOT(ret);
}


// A simple test, search for "AAG" in "TAGAGA" with "random" scores,
// not so high limit.  Prefix Search.
SEQAN_DEFINE_TEST(test_dp_quality_find_prefix_not_increasing) {
    String<Dna5> hstck = "TAGAGA";
    String<Dna5Q> ndl = "AAG";
    assignQualityValue(ndl[0], 5);
    assignQualityValue(ndl[1], 10);
    assignQualityValue(ndl[2], 2);

    Finder<String<Dna5> > finder(hstck);
    Pattern<String<Dna5Q>, QualityDpSearch<FindPrefix> > pattern(ndl, -1000);

    // The following has been verified on paper :)

    bool ret = find(finder, pattern);
    SEQAN_ASSERT(ret);
    SEQAN_ASSERT_EQ(1u, endPosition(finder));
    SEQAN_ASSERT_EQ(-17, getScore(pattern));

    ret = find(finder, pattern);
    SEQAN_ASSERT(ret);
    SEQAN_ASSERT_EQ(2u, endPosition(finder));
    SEQAN_ASSERT_EQ(-7, getScore(pattern));

    ret = find(finder, pattern);
    SEQAN_ASSERT(ret);
    SEQAN_ASSERT_EQ(3u, endPosition(finder));
    SEQAN_ASSERT_EQ(-5, getScore(pattern));

    ret = find(finder, pattern);
    SEQAN_ASSERT(ret);
    SEQAN_ASSERT_EQ(4u, endPosition(finder));
    SEQAN_ASSERT_EQ(-7, getScore(pattern));

    ret = find(finder, pattern);
    SEQAN_ASSERT(ret);
    SEQAN_ASSERT_EQ(5u, endPosition(finder));
    SEQAN_ASSERT_EQ(-9, getScore(pattern));

    ret = find(finder, pattern);
    SEQAN_ASSERT(ret);
    SEQAN_ASSERT_EQ(6u, endPosition(finder));
    SEQAN_ASSERT_EQ(-11, getScore(pattern));

    ret = find(finder, pattern);
    SEQAN_ASSERT_NOT(ret);
}


// A simple test, search for "AAG" in "TAGAGA" with "random" scores,
// not so high limit.  Also testing findBegin().
SEQAN_DEFINE_TEST(test_dp_quality_find_begin_not_increasing) {
    String<Dna5> hstck = "TAGAGA";
    String<Dna5Q> ndl = "AAG";
    assignQualityValue(ndl[0], 5);
    assignQualityValue(ndl[1], 10);
    assignQualityValue(ndl[2], 2);

    Finder<String<Dna5> > finder(hstck);
    Pattern<String<Dna5Q>, QualityDpSearch<FindInfix> > pattern(ndl, -7);

    // The following has been verified on paper :)

    bool ret = find(finder, pattern);
    SEQAN_ASSERT(ret);
    SEQAN_ASSERT_EQ(2u, endPosition(finder));
    SEQAN_ASSERT_EQ(-7, getScore(pattern));

    ret = findBegin(finder, pattern, getScore(pattern));
    SEQAN_ASSERT(ret);
    SEQAN_ASSERT_EQ(-7, getFindBeginScore(pattern));
    SEQAN_ASSERT_EQ(1u, beginPosition(finder));
    ret = findBegin(finder, pattern, getScore(pattern));
    SEQAN_ASSERT(ret);
    SEQAN_ASSERT_EQ(-7, getFindBeginScore(pattern));
    SEQAN_ASSERT_EQ(0u, beginPosition(finder));
    ret = findBegin(finder, pattern, getScore(pattern));
    SEQAN_ASSERT_NOT(ret);

    ret = find(finder, pattern);
    SEQAN_ASSERT(ret);
    SEQAN_ASSERT_EQ(3u, endPosition(finder));
    SEQAN_ASSERT_EQ(-5, getScore(pattern));

    ret = findBegin(finder, pattern, getScore(pattern));
    SEQAN_ASSERT(ret);
    SEQAN_ASSERT_EQ(-5, getFindBeginScore(pattern));
    SEQAN_ASSERT_EQ(1u, beginPosition(finder));
    ret = findBegin(finder, pattern, getScore(pattern));
    SEQAN_ASSERT(ret);
    SEQAN_ASSERT_EQ(-5, getFindBeginScore(pattern));
    SEQAN_ASSERT_EQ(0u, beginPosition(finder));
    ret = findBegin(finder, pattern, getScore(pattern));
    SEQAN_ASSERT_NOT(ret);

    ret = find(finder, pattern);
    SEQAN_ASSERT(ret);
    SEQAN_ASSERT_EQ(4u, endPosition(finder));
    SEQAN_ASSERT_EQ(-7, getScore(pattern));

    ret = findBegin(finder, pattern, getScore(pattern));
    SEQAN_ASSERT(ret);
    SEQAN_ASSERT_EQ(-7, getFindBeginScore(pattern));
    SEQAN_ASSERT_EQ(3u, beginPosition(finder));
    ret = findBegin(finder, pattern, getScore(pattern));
    SEQAN_ASSERT(ret);
    SEQAN_ASSERT_EQ(-7, getFindBeginScore(pattern));
    SEQAN_ASSERT_EQ(2u, beginPosition(finder));
    ret = findBegin(finder, pattern, getScore(pattern));
    SEQAN_ASSERT(ret);
    SEQAN_ASSERT_EQ(-7, getFindBeginScore(pattern));
    SEQAN_ASSERT_EQ(1u, beginPosition(finder));
    ret = findBegin(finder, pattern, getScore(pattern));
    SEQAN_ASSERT(ret);
    SEQAN_ASSERT_EQ(-7, getFindBeginScore(pattern));
    SEQAN_ASSERT_EQ(0u, beginPosition(finder));
    ret = findBegin(finder, pattern, getScore(pattern));
    SEQAN_ASSERT_NOT(ret);

    ret = find(finder, pattern);
    SEQAN_ASSERT(ret);
    SEQAN_ASSERT_EQ(5u, endPosition(finder));
    SEQAN_ASSERT_EQ(-5, getScore(pattern));

    ret = findBegin(finder, pattern, getScore(pattern));
    SEQAN_ASSERT(ret);
    SEQAN_ASSERT_EQ(-5, getFindBeginScore(pattern));
    SEQAN_ASSERT_EQ(3u, beginPosition(finder));
    ret = findBegin(finder, pattern, getScore(pattern));
    SEQAN_ASSERT(ret);
    SEQAN_ASSERT_EQ(-5, getFindBeginScore(pattern));
    SEQAN_ASSERT_EQ(2u, beginPosition(finder));
    ret = findBegin(finder, pattern, getScore(pattern));
    SEQAN_ASSERT_NOT(ret);

    ret = find(finder, pattern);
    SEQAN_ASSERT(ret);
    SEQAN_ASSERT_EQ(6u, endPosition(finder));
    SEQAN_ASSERT_EQ(-7, getScore(pattern));

    ret = findBegin(finder, pattern, getScore(pattern));
    SEQAN_ASSERT(ret);
    SEQAN_ASSERT_EQ(-7, getFindBeginScore(pattern));
    SEQAN_ASSERT_EQ(5u, beginPosition(finder));
    ret = findBegin(finder, pattern, getScore(pattern));
    SEQAN_ASSERT(ret);
    SEQAN_ASSERT_EQ(-7, getFindBeginScore(pattern));
    SEQAN_ASSERT_EQ(4u, beginPosition(finder));
    ret = findBegin(finder, pattern, getScore(pattern));
    SEQAN_ASSERT(ret);
    SEQAN_ASSERT_EQ(-7, getFindBeginScore(pattern));
    SEQAN_ASSERT_EQ(3u, beginPosition(finder));
    ret = findBegin(finder, pattern, getScore(pattern));
    SEQAN_ASSERT(ret);
    SEQAN_ASSERT_EQ(-7, getFindBeginScore(pattern));
    SEQAN_ASSERT_EQ(2u, beginPosition(finder));
    ret = findBegin(finder, pattern, getScore(pattern));
    SEQAN_ASSERT_NOT(ret);

    ret = find(finder, pattern);
    SEQAN_ASSERT_NOT(ret);
}


// Test setEndPosition when searching for "AAG" in "TAGAGA" with
// "random" scores.
SEQAN_DEFINE_TEST(test_dp_quality_set_end_position_not_increasing) {
    String<Dna5> hstck = "TAGAGA";
    String<Dna5Q> ndl = "AAG";
    assignQualityValue(ndl[0], 5);
    assignQualityValue(ndl[1], 10);
    assignQualityValue(ndl[2], 2);

    Finder<String<Dna5> > finder(hstck);
    Pattern<String<Dna5Q>, QualityDpSearch<FindPrefix> > pattern(ndl, -7);

    // The following has been verified on paper :)

    bool ret = find(finder, pattern);
    SEQAN_ASSERT(ret);
    SEQAN_ASSERT_EQ(2u, endPosition(finder));
    SEQAN_ASSERT_EQ(-7, getScore(pattern));

    ret = setEndPosition(finder, pattern, 3u);
    SEQAN_ASSERT(ret);
    SEQAN_ASSERT_EQ(3u, endPosition(finder));
    SEQAN_ASSERT_EQ(-5, getScore(pattern));

    ret = setEndPosition(finder, pattern, 2u);
    SEQAN_ASSERT(ret);
    SEQAN_ASSERT_EQ(2u, endPosition(finder));
    SEQAN_ASSERT_EQ(-7, getScore(pattern));

    ret = find(finder, pattern);
    SEQAN_ASSERT(ret);
    SEQAN_ASSERT_EQ(3u, endPosition(finder));
    SEQAN_ASSERT_EQ(-5, getScore(pattern));
}


SEQAN_BEGIN_TESTSUITE(test_dp_quality) {
    SEQAN_CALL_TEST(test_dp_quality_find_infix_simple_high_limit);
    SEQAN_CALL_TEST(test_dp_quality_find_infix_simple);
    SEQAN_CALL_TEST(test_dp_quality_find_infix_not_increasing_high_limit);
    SEQAN_CALL_TEST(test_dp_quality_find_infix_not_increasing);
    SEQAN_CALL_TEST(test_dp_quality_find_prefix_not_increasing);
    SEQAN_CALL_TEST(test_dp_quality_find_begin_not_increasing);
    SEQAN_CALL_TEST(test_dp_quality_set_end_position_not_increasing);
}
SEQAN_END_TESTSUITE
