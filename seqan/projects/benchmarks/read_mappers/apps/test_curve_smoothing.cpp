/* Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>

   Tests for curve_smoothing.h.
*/

#include <seqan/basic.h>

#include "curve_smoothing.h"


// Test with one beginPos interval, non-decreasing only.
SEQAN_DEFINE_TEST(test_curve_smoothing_one_non_decreasing) {
    // Build string of weighted matches.
    String<WeightedMatch> matches;
    appendValue(matches, WeightedMatch(0, true, 0, -6, 0));
    appendValue(matches, WeightedMatch(0, true, 1, -5, 0));
    appendValue(matches, WeightedMatch(0, true, 2, -5, 0));
    appendValue(matches, WeightedMatch(0, true, 3, -4, 0));

    // Make copy of the string's state before for later comparison.
    const String<WeightedMatch> matchesBefore(matches);

    // Call the function under test -- smoothErrorCurve().
    smoothErrorCurve(matches);

    // Test state after running the algorithm.
    SEQAN_ASSERT_EQ(length(matches), length(matchesBefore));
    SEQAN_ASSERT_EQ(matches[0], matchesBefore[0]);
    SEQAN_ASSERT_EQ(matches[1], matchesBefore[1]);
    SEQAN_ASSERT_EQ(matches[2], matchesBefore[2]);
    SEQAN_ASSERT_EQ(matches[3], matchesBefore[3]);
}


// Test with one beginPos interval, non-increasing only.
SEQAN_DEFINE_TEST(test_curve_smoothing_one_non_increasing) {
    // Build string of weighted matches.
    String<WeightedMatch> matches;
    appendValue(matches, WeightedMatch(0, true, 0, -4, 0));
    appendValue(matches, WeightedMatch(0, true, 1, -5, 0));
    appendValue(matches, WeightedMatch(0, true, 2, -5, 0));
    appendValue(matches, WeightedMatch(0, true, 3, -6, 0));

    // Make copy of the string's state before for later comparison.
    const String<WeightedMatch> matchesBefore(matches);

    // Call the function under test -- smoothErrorCurve().
    smoothErrorCurve(matches);

    // Test state after running the algorithm.
    SEQAN_ASSERT_EQ(length(matches), length(matchesBefore));
    SEQAN_ASSERT_EQ(matches[0], matchesBefore[0]);
    SEQAN_ASSERT_EQ(matches[1], matchesBefore[1]);
    SEQAN_ASSERT_EQ(matches[2], matchesBefore[2]);
    SEQAN_ASSERT_EQ(matches[3], matchesBefore[3]);
}


// Test with one beginPos interval, hill-shaped.
SEQAN_DEFINE_TEST(test_curve_smoothing_one_hill) {
    // Build string of weighted matches.
    String<WeightedMatch> matches;
    appendValue(matches, WeightedMatch(0, true, 0, -4, 0));
    appendValue(matches, WeightedMatch(0, true, 1, -5, 0));
    appendValue(matches, WeightedMatch(0, true, 2, -6, 0));
    appendValue(matches, WeightedMatch(0, true, 3, -4, 0));

    // Make copy of the string's state before for later comparison.
    const String<WeightedMatch> matchesBefore(matches);

    // Call the function under test -- smoothErrorCurve().
    smoothErrorCurve(matches);

    // Test state after running the algorithm.
    SEQAN_ASSERT_EQ(length(matches), length(matchesBefore));
    SEQAN_ASSERT_EQ(matches[0], WeightedMatch(0, true, 0, -4, 0));
    SEQAN_ASSERT_EQ(matches[1], WeightedMatch(0, true, 1, -4, 0));
    SEQAN_ASSERT_EQ(matches[2], WeightedMatch(0, true, 2, -4, 0));
    SEQAN_ASSERT_EQ(matches[3], WeightedMatch(0, true, 3, -4, 0));
}


// Test with one beginPos interval, valley-shaped.
SEQAN_DEFINE_TEST(test_curve_smoothing_one_valley) {
    // Build string of weighted matches.
    String<WeightedMatch> matches;
    appendValue(matches, WeightedMatch(0, true, 0, -6, 0));
    appendValue(matches, WeightedMatch(0, true, 1, -5, 0));
    appendValue(matches, WeightedMatch(0, true, 2, -4, 0));
    appendValue(matches, WeightedMatch(0, true, 3, -6, 0));

    // Make copy of the string's state before for later comparison.
    const String<WeightedMatch> matchesBefore(matches);

    // Call the function under test -- smoothErrorCurve().
    smoothErrorCurve(matches);

    // Test state after running the algorithm.
    SEQAN_ASSERT_EQ(length(matches), length(matchesBefore));
    SEQAN_ASSERT_EQ(matches[0], matchesBefore[0]);
    SEQAN_ASSERT_EQ(matches[1], matchesBefore[1]);
    SEQAN_ASSERT_EQ(matches[2], matchesBefore[2]);
    SEQAN_ASSERT_EQ(matches[3], matchesBefore[3]);
}


// Test with two intervals, first one-element, second longer hill shaped.
SEQAN_DEFINE_TEST(test_curve_smoothing_two_short_long) {
    // Build string of weighted matches.
    String<WeightedMatch> matches;
    appendValue(matches, WeightedMatch(0, true, 0, -4, 0));
    appendValue(matches, WeightedMatch(0, true, 1, -5, 1));
    appendValue(matches, WeightedMatch(0, true, 2, -6, 1));
    appendValue(matches, WeightedMatch(0, true, 3, -4, 1));
    appendValue(matches, WeightedMatch(0, true, 4, -4, 1));
    appendValue(matches, WeightedMatch(0, true, 5, -4, 1));

    // Make copy of the string's state before for later comparison.
    const String<WeightedMatch> matchesBefore(matches);

    // Call the function under test -- smoothErrorCurve().
    smoothErrorCurve(matches);

    // Test state after running the algorithm.
    SEQAN_ASSERT_EQ(length(matches), length(matchesBefore));
    SEQAN_ASSERT_EQ(matches[0], WeightedMatch(0, true, 0, -4, 0));
    SEQAN_ASSERT_EQ(matches[1], WeightedMatch(0, true, 1, -5, 1));
    SEQAN_ASSERT_EQ(matches[2], WeightedMatch(0, true, 2, -5, 1));
    SEQAN_ASSERT_EQ(matches[3], WeightedMatch(0, true, 3, -4, 1));
    SEQAN_ASSERT_EQ(matches[4], WeightedMatch(0, true, 4, -4, 1));
    SEQAN_ASSERT_EQ(matches[5], WeightedMatch(0, true, 5, -4, 1));
}


// Test with two intervals, first longer hill shaped, second short.
SEQAN_DEFINE_TEST(test_curve_smoothing_two_long_short) {
    // Build string of weighted matches.
    String<WeightedMatch> matches;
    appendValue(matches, WeightedMatch(0, true, 0, -4, 0));
    appendValue(matches, WeightedMatch(0, true, 1, -5, 0));
    appendValue(matches, WeightedMatch(0, true, 2, -6, 0));
    appendValue(matches, WeightedMatch(0, true, 3, -5, 0));
    appendValue(matches, WeightedMatch(0, true, 4, -5, 0));
    appendValue(matches, WeightedMatch(0, true, 5, -3, 1));

    // Make copy of the string's state before for later comparison.
    const String<WeightedMatch> matchesBefore(matches);

    // Call the function under test -- smoothErrorCurve().
    smoothErrorCurve(matches);

    // Test state after running the algorithm.
    SEQAN_ASSERT_EQ(length(matches), length(matchesBefore));
    SEQAN_ASSERT_EQ(matches[0], WeightedMatch(0, true, 0, -4, 0));
    SEQAN_ASSERT_EQ(matches[1], WeightedMatch(0, true, 1, -5, 0));
    SEQAN_ASSERT_EQ(matches[2], WeightedMatch(0, true, 2, -5, 0));
    SEQAN_ASSERT_EQ(matches[3], WeightedMatch(0, true, 3, -5, 0));
    SEQAN_ASSERT_EQ(matches[4], WeightedMatch(0, true, 4, -5, 0));
    SEQAN_ASSERT_EQ(matches[5], WeightedMatch(0, true, 5, -3, 1));
}


// Test with two short intervals.
SEQAN_DEFINE_TEST(test_curve_smoothing_two_short_short) {
    // Build string of weighted matches.
    String<WeightedMatch> matches;
    appendValue(matches, WeightedMatch(0, true, 0, -4, 0));
    appendValue(matches, WeightedMatch(0, true, 1, -5, 1));

    // Make copy of the string's state before for later comparison.
    const String<WeightedMatch> matchesBefore(matches);

    // Call the function under test -- smoothErrorCurve().
    smoothErrorCurve(matches);

    // Test state after running the algorithm.
    SEQAN_ASSERT_EQ(length(matches), length(matchesBefore));
    SEQAN_ASSERT_EQ(matches[0], matchesBefore[0]);
    SEQAN_ASSERT_EQ(matches[1], matchesBefore[1]);
}


// Test with two long intervals.
SEQAN_DEFINE_TEST(test_curve_smoothing_two_long_long) {
    // Build string of weighted matches.
    String<WeightedMatch> matches;
    appendValue(matches, WeightedMatch(0, true, 0, -4, 0));
    appendValue(matches, WeightedMatch(0, true, 1, -5, 0));
    appendValue(matches, WeightedMatch(0, true, 2, -6, 0));
    appendValue(matches, WeightedMatch(0, true, 3, -5, 0));
    appendValue(matches, WeightedMatch(0, true, 4, -4, 0));
    appendValue(matches, WeightedMatch(0, true, 5, -2, 1));
    appendValue(matches, WeightedMatch(0, true, 6, -4, 1));
    appendValue(matches, WeightedMatch(0, true, 7, -5, 1));
    appendValue(matches, WeightedMatch(0, true, 8, -3, 1));

    // Make copy of the string's state before for later comparison.
    const String<WeightedMatch> matchesBefore(matches);

    // Call the function under test -- smoothErrorCurve().
    smoothErrorCurve(matches);

    // Test state after running the algorithm.
    SEQAN_ASSERT_EQ(length(matches), length(matchesBefore));
    SEQAN_ASSERT_EQ(matches[0], WeightedMatch(0, true, 0, -4, 0));
    SEQAN_ASSERT_EQ(matches[1], WeightedMatch(0, true, 1, -4, 0));
    SEQAN_ASSERT_EQ(matches[2], WeightedMatch(0, true, 2, -4, 0));
    SEQAN_ASSERT_EQ(matches[3], WeightedMatch(0, true, 3, -4, 0));
    SEQAN_ASSERT_EQ(matches[4], WeightedMatch(0, true, 4, -4, 0));
    SEQAN_ASSERT_EQ(matches[5], WeightedMatch(0, true, 5, -2, 1));
    SEQAN_ASSERT_EQ(matches[6], WeightedMatch(0, true, 6, -3, 1));
    SEQAN_ASSERT_EQ(matches[7], WeightedMatch(0, true, 7, -3, 1));
    SEQAN_ASSERT_EQ(matches[8], WeightedMatch(0, true, 8, -3, 1));
}


// Test fillGaps without any gaps.
SEQAN_DEFINE_TEST(test_curve_smoothing_fill_gaps_no_gap) {
    // Build string of weighted match objects.
    String<WeightedMatch> matches;

    appendValue(matches, WeightedMatch(0, true, 10, -1, 0));
    appendValue(matches, WeightedMatch(0, true, 11, -1, 0));
    appendValue(matches, WeightedMatch(0, true, 12, -1, 0));
    appendValue(matches, WeightedMatch(1, true, 10, -1, 0));
    appendValue(matches, WeightedMatch(1, true, 11, -1, 0));
    appendValue(matches, WeightedMatch(1, true, 12, -1, 0));

    // Make copy of the string's state before for later comparison.
    const String<WeightedMatch> matchesBefore(matches);

    // Call the function under test -- fillGaps().
    fillGaps(matches);

    // Test state after running the algorithm.
    SEQAN_ASSERT_EQ(length(matches), length(matchesBefore));
    SEQAN_ASSERT_EQ(matches[0], matches[0]);
    SEQAN_ASSERT_EQ(matches[1], matches[1]);
    SEQAN_ASSERT_EQ(matches[2], matches[2]);
    SEQAN_ASSERT_EQ(matches[3], matches[3]);
    SEQAN_ASSERT_EQ(matches[4], matches[4]);
    SEQAN_ASSERT_EQ(matches[5], matches[5]);
    SEQAN_ASSERT_EQ(matches[6], matches[6]);
}


SEQAN_DEFINE_TEST(test_curve_smoothing_fill_gaps_one_short_gap) {
    // Build string of weighted match objects.
    String<WeightedMatch> matches;

    appendValue(matches, WeightedMatch(0, true, 10, -1, 0));
    appendValue(matches, WeightedMatch(0, true, 12, -1, 0));
    appendValue(matches, WeightedMatch(0, true, 13, -1, 0));

    // Make copy of the string's state before for later comparison.
    const String<WeightedMatch> matchesBefore(matches);

    // Call the function under test -- fillGaps().
    fillGaps(matches);

    // Test state after running the algorithm.
    SEQAN_ASSERT_EQ(length(matches), length(matchesBefore) + 1);
    SEQAN_ASSERT_EQ(matches[0], WeightedMatch(0, true, 10, -1, 0));
    SEQAN_ASSERT_EQ(matches[1], WeightedMatch(0, true, 11, -1, 0));
    SEQAN_ASSERT_EQ(matches[2], WeightedMatch(0, true, 12, -1, 0));
    SEQAN_ASSERT_EQ(matches[3], WeightedMatch(0, true, 13, -1, 0));
}


SEQAN_DEFINE_TEST(test_curve_smoothing_fill_gaps_one_long_gap) {
    // Build string of weighted match objects.
    String<WeightedMatch> matches;

    appendValue(matches, WeightedMatch(0, true, 10, -1, 0));
    appendValue(matches, WeightedMatch(0, true, 16, -1, 0));
    appendValue(matches, WeightedMatch(0, true, 17, -1, 0));

    // Make copy of the string's state before for later comparison.
    const String<WeightedMatch> matchesBefore(matches);

    // Call the function under test -- fillGaps().
    fillGaps(matches);

    // Test state after running the algorithm.
    SEQAN_ASSERT_EQ(length(matches), length(matchesBefore) + 5);
    SEQAN_ASSERT_EQ(matches[0], WeightedMatch(0, true, 10, -1, 0));
    SEQAN_ASSERT_EQ(matches[1], WeightedMatch(0, true, 11, -1, 0));
    SEQAN_ASSERT_EQ(matches[2], WeightedMatch(0, true, 12, -1, 0));
    SEQAN_ASSERT_EQ(matches[3], WeightedMatch(0, true, 13, -1, 0));
    SEQAN_ASSERT_EQ(matches[4], WeightedMatch(0, true, 14, -1, 0));
    SEQAN_ASSERT_EQ(matches[5], WeightedMatch(0, true, 15, -1, 0));
    SEQAN_ASSERT_EQ(matches[6], WeightedMatch(0, true, 16, -1, 0));
    SEQAN_ASSERT_EQ(matches[7], WeightedMatch(0, true, 17, -1, 0));
}


SEQAN_DEFINE_TEST(test_curve_smoothing_fill_gaps_long_short_gap) {
    // Build string of weighted match objects.
    String<WeightedMatch> matches;

    appendValue(matches, WeightedMatch(0, true, 8, -1, 0));
    appendValue(matches, WeightedMatch(0, true, 14, -1, 0));
    appendValue(matches, WeightedMatch(0, true, 16, -1, 0));
    appendValue(matches, WeightedMatch(0, true, 17, -1, 0));

    // Make copy of the string's state before for later comparison.
    const String<WeightedMatch> matchesBefore(matches);

    // Call the function under test -- fillGaps().
    fillGaps(matches);

    // Test state after running the algorithm.
    SEQAN_ASSERT_EQ(length(matches), length(matchesBefore) + 6);
    SEQAN_ASSERT_EQ(matches[0], WeightedMatch(0, true, 8, -1, 0));
    SEQAN_ASSERT_EQ(matches[1], WeightedMatch(0, true, 9, -1, 0));
    SEQAN_ASSERT_EQ(matches[2], WeightedMatch(0, true, 10, -1, 0));
    SEQAN_ASSERT_EQ(matches[3], WeightedMatch(0, true, 11, -1, 0));
    SEQAN_ASSERT_EQ(matches[4], WeightedMatch(0, true, 12, -1, 0));
    SEQAN_ASSERT_EQ(matches[5], WeightedMatch(0, true, 13, -1, 0));
    SEQAN_ASSERT_EQ(matches[6], WeightedMatch(0, true, 14, -1, 0));
    SEQAN_ASSERT_EQ(matches[7], WeightedMatch(0, true, 15, -1, 0));
    SEQAN_ASSERT_EQ(matches[8], WeightedMatch(0, true, 16, -1, 0));
    SEQAN_ASSERT_EQ(matches[9], WeightedMatch(0, true, 17, -1, 0));
}


SEQAN_DEFINE_TEST(test_curve_smoothing_fill_gaps_short_long_gap)  {
    // Build string of weighted match objects.
    String<WeightedMatch> matches;

    appendValue(matches, WeightedMatch(0, true, 8, -1, 0));
    appendValue(matches, WeightedMatch(0, true, 10, -1, 0));
    appendValue(matches, WeightedMatch(0, true, 16, -1, 0));
    appendValue(matches, WeightedMatch(0, true, 17, -1, 0));

    // Make copy of the string's state before for later comparison.
    const String<WeightedMatch> matchesBefore(matches);

    // Call the function under test -- fillGaps().
    fillGaps(matches);

    // Test state after running the algorithm.
    SEQAN_ASSERT_EQ(length(matches), length(matchesBefore) + 6);
    SEQAN_ASSERT_EQ(matches[0], WeightedMatch(0, true, 8, -1, 0));
    SEQAN_ASSERT_EQ(matches[1], WeightedMatch(0, true, 9, -1, 0));
    SEQAN_ASSERT_EQ(matches[2], WeightedMatch(0, true, 10, -1, 0));
    SEQAN_ASSERT_EQ(matches[3], WeightedMatch(0, true, 11, -1, 0));
    SEQAN_ASSERT_EQ(matches[4], WeightedMatch(0, true, 12, -1, 0));
    SEQAN_ASSERT_EQ(matches[5], WeightedMatch(0, true, 13, -1, 0));
    SEQAN_ASSERT_EQ(matches[6], WeightedMatch(0, true, 14, -1, 0));
    SEQAN_ASSERT_EQ(matches[7], WeightedMatch(0, true, 15, -1, 0));
    SEQAN_ASSERT_EQ(matches[8], WeightedMatch(0, true, 16, -1, 0));
    SEQAN_ASSERT_EQ(matches[9], WeightedMatch(0, true, 17, -1, 0));
}


// Test fillGaps() with real data that failed.
SEQAN_DEFINE_TEST(test_curve_smoothing_fill_gaps_real_data) {
    // Build string of weighted match objects.
    String<WeightedMatch> matches;

    appendValue(matches, WeightedMatch(6, true, 1074876, -6, 1074842));
    appendValue(matches, WeightedMatch(6, true, 1074878, -3, 1074842));
    appendValue(matches, WeightedMatch(6, true, 1074879, -6, 1074842));

    // Make copy of the string's state before for later comparison.
    const String<WeightedMatch> matchesBefore(matches);

    // Call the function under test -- fillGaps().
    fillGaps(matches);

    // Test state after running the algorithm.
    SEQAN_ASSERT_EQ(length(matches), length(matchesBefore) + 1);
    SEQAN_ASSERT_EQ(matches[0], WeightedMatch(6, true, 1074876, -6, 1074842));
    SEQAN_ASSERT_EQ(matches[1], WeightedMatch(6, true, 1074877, -6, 1074842));
    SEQAN_ASSERT_EQ(matches[2], WeightedMatch(6, true, 1074878, -3, 1074842));
    SEQAN_ASSERT_EQ(matches[3], WeightedMatch(6, true, 1074879, -6, 1074842));
}


// Test smoothErrorCurve() with another real data set that failed.
SEQAN_DEFINE_TEST(test_curve_smoothing_real_data) {
    // Build string of weighted match objects.
    String<WeightedMatch> matches;

    appendValue(matches, WeightedMatch(6, true, 8028062, -3, 8028031));
    appendValue(matches, WeightedMatch(6, true, 8028063, -2, 8028031));
    appendValue(matches, WeightedMatch(6, true, 8028064, -3, 8028031));
    appendValue(matches, WeightedMatch(6, true, 8028065, -2, 8028031));
    appendValue(matches, WeightedMatch(6, true, 8028066, -3, 8028031));

    // Make copy of the string's state before for later comparison.
    const String<WeightedMatch> matchesBefore(matches);

    // Call the function under test -- fillGaps().
    smoothErrorCurve(matches);

    // Test state after running the algorithm.
    SEQAN_ASSERT_EQ(length(matches), length(matchesBefore));
    SEQAN_ASSERT_EQ(matches[0], WeightedMatch(6, true, 8028062, -3, 8028031));
    SEQAN_ASSERT_EQ(matches[1], WeightedMatch(6, true, 8028063, -2, 8028031));
    SEQAN_ASSERT_EQ(matches[2], WeightedMatch(6, true, 8028064, -2, 8028031));
    SEQAN_ASSERT_EQ(matches[3], WeightedMatch(6, true, 8028065, -2, 8028031));
    SEQAN_ASSERT_EQ(matches[4], WeightedMatch(6, true, 8028066, -3, 8028031));
}


SEQAN_BEGIN_TESTSUITE(test_curve_smoothing) {
    // Test the error curve smoothing algorithm, increasing/decreasing
    // is related to the distance not the score.
    SEQAN_CALL_TEST(test_curve_smoothing_one_non_decreasing);
    SEQAN_CALL_TEST(test_curve_smoothing_one_non_increasing);
    SEQAN_CALL_TEST(test_curve_smoothing_one_hill);
    SEQAN_CALL_TEST(test_curve_smoothing_one_valley);

    // Test the error curve smoothing algorith with two intervals.
    SEQAN_CALL_TEST(test_curve_smoothing_two_short_long);
    SEQAN_CALL_TEST(test_curve_smoothing_two_long_short);
    SEQAN_CALL_TEST(test_curve_smoothing_two_short_short);
    SEQAN_CALL_TEST(test_curve_smoothing_two_long_long);

    // Test the gap filling algorithm with various test cases.
    SEQAN_CALL_TEST(test_curve_smoothing_fill_gaps_no_gap);
    SEQAN_CALL_TEST(test_curve_smoothing_fill_gaps_one_short_gap);
    SEQAN_CALL_TEST(test_curve_smoothing_fill_gaps_one_long_gap);
    SEQAN_CALL_TEST(test_curve_smoothing_fill_gaps_long_short_gap);
    SEQAN_CALL_TEST(test_curve_smoothing_fill_gaps_short_long_gap);

    SEQAN_CALL_TEST(test_curve_smoothing_fill_gaps_real_data);
    SEQAN_CALL_TEST(test_curve_smoothing_real_data);
}
SEQAN_END_TESTSUITE
