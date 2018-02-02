// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
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

SEQAN_DEFINE_TEST(Infix)
{
//____________________________________________________________________________
// infix

    Infix<String<char> >::Type infix_1;
    String<char> str_1 = "this is a string";
    setHost(infix_1, str_1);
    SEQAN_ASSERT_EQ(length(infix_1), 0u);
    SEQAN_ASSERT_EQ(getObjectId(infix_1), getObjectId(str_1));

    setEnd(infix_1, end(str_1));
    SEQAN_ASSERT_EQ(infix_1, str_1);
    SEQAN_ASSERT_EQ(length(infix_1), length(infix_1));

    setEndPosition(infix_1, 9);
    SEQAN_ASSERT_EQ(infix_1, infix(str_1, 0, 9));

    Infix<String<char> >::Type infix_2(infix_1);
    SEQAN_ASSERT_EQ(infix_2, infix(str_1, 0, 9));
    SEQAN_ASSERT_EQ(infix_2, infix_1);
    SEQAN_ASSERT_EQ(getObjectId(infix_1), getObjectId(infix_2));

    setBeginPosition(infix_2, 5);
    SEQAN_ASSERT_EQ(infix_2, "is a");

    setBegin(infix_2, begin(str_1));
    SEQAN_ASSERT_EQ(infix_2, "this is a");

    Infix<String<char> >::Type infix_3(str_1);
    SEQAN_ASSERT_EQ(infix_3, str_1);
    SEQAN_ASSERT_EQ(getObjectId(infix_3), getObjectId(str_1));

    Infix<String<char> >::Type infix_4(str_1, 5, 9);
    SEQAN_ASSERT_EQ(infix_4, "is a");

    SEQAN_ASSERT_EQ(capacity(infix_4), capacity(str_1) - length(str_1) + length(infix_4));

    Infix<String<char> >::Type infix_5(str_1, begin(str_1), end(str_1));
    SEQAN_ASSERT_EQ(infix_5, str_1);

    SEQAN_ASSERT(begin(infix_5) == begin(str_1));
    SEQAN_ASSERT_EQ(beginPosition(infix_5), 0u);
    SEQAN_ASSERT(end(infix_5) == end(str_1));
    SEQAN_ASSERT_EQ(endPosition(infix_5), length(str_1));

    SEQAN_ASSERT(begin(infix(str_1, 0, length(str_1))) == begin(str_1));
    SEQAN_ASSERT(end(infix(str_1, 0, length(str_1))) == end(str_1));
    SEQAN_ASSERT_EQ(length(infix(str_1, 0, length(str_1))), length(str_1));

    // str_1 = "begin middle end";
    // assign(infix(str_1, 6, 12),  "to");
    // SEQAN_ASSERT_EQ(str_1, "begin to end");

    // assign(infix(str_1, 6, 8), "the test", 14);
    // SEQAN_ASSERT_EQ(str_1, "begin the test");

//    setEnd(infix_1);
//    SEQAN_ASSERT(infix_1 == "");

//____________________________________________________________________________
// test infix iteration

    str_1 = "begin middle end";
    goBegin(infix_1);
    SEQAN_ASSERT_EQ(infix_1, "b");

    goBegin(infix_1, str_1);
    SEQAN_ASSERT_EQ(infix_1, "b");

    goPrevious(infix_1);
    SEQAN_ASSERT(atBegin(infix_1));

    goEnd(infix_1);
    SEQAN_ASSERT_EQ(infix_1, str_1);

    goEnd(infix_1, str_1);
    SEQAN_ASSERT_EQ(infix_1, str_1);

    goNext(infix_1);
    SEQAN_ASSERT(atEnd(infix_1));
//____________________________________________________________________________
// compare operators

    str_1 = "hellohello";
    Infix<String<char> >::Type infix_6(str_1, 0, 10);

//     SEQAN_ASSERT_EQ(infix_6, str_1);
//
//     infix_6 += str_1;
    SEQAN_ASSERT(isEqual(infix_6, "hellohello"));

    SEQAN_ASSERT_NEQ(infix_6, "bla");
    SEQAN_ASSERT_NOT(isNotEqual(infix_6, "hellohello"));

    SEQAN_ASSERT_NOT(infix_6 < "hello");
    SEQAN_ASSERT_NOT(isLess(infix_6, "hello"));

    SEQAN_ASSERT_NOT(infix_6 <= "hello");
    SEQAN_ASSERT_NOT(isLessOrEqual(infix_6, "hello"));

    SEQAN_ASSERT_GT(infix_6, "hello");
    SEQAN_ASSERT(isGreater(infix_6, "hello"));

    SEQAN_ASSERT_GEQ(infix_6, "hello");
    SEQAN_ASSERT(isGreaterOrEqual(infix_6, "hello"));
//____________________________________________________________________________

//     clear(infix_6);
//     SEQAN_ASSERT_EQ(infix_6, "");
//____________________________________________________________________________
}

//////////////////////////////////////////////////////////////////////////////

SEQAN_DEFINE_TEST(Suffix)
{
//____________________________________________________________________________
// suffix

    String<char> str_1 = "this is a string";

    Suffix<String<char> >::Type suffix_1;
    setHost(suffix_1, str_1);
    SEQAN_ASSERT_EQ(length(suffix_1), length(str_1));
    SEQAN_ASSERT_EQ(getObjectId(suffix_1), getObjectId(str_1));

    Suffix<String<char> >::Type suffix_2(suffix_1);
    SEQAN_ASSERT_EQ(suffix_2, suffix_1);
    SEQAN_ASSERT_EQ(getObjectId(suffix_1), getObjectId(suffix_2));

    setBeginPosition(suffix_2, 5);
    SEQAN_ASSERT_EQ(suffix_2, "is a string");

    setBegin(suffix_2, begin(str_1));
    SEQAN_ASSERT_EQ(suffix_2, "this is a string");

    Suffix<String<char> >::Type suffix_3(str_1);
    SEQAN_ASSERT_EQ(suffix_3, str_1);
    SEQAN_ASSERT_EQ(getObjectId(suffix_3), getObjectId(str_1));

    Suffix<String<char> >::Type suffix_4(str_1, 5);
    SEQAN_ASSERT_EQ(suffix_4, "is a string");

    SEQAN_ASSERT_EQ(capacity(suffix_4), capacity(str_1) - length(str_1) + length(suffix_4));

    Suffix<String<char> >::Type suffix_5(str_1, begin(str_1));
    SEQAN_ASSERT_EQ(suffix_5, str_1);

    SEQAN_ASSERT(begin(suffix_5) == begin(str_1));
    SEQAN_ASSERT_EQ(beginPosition(suffix_5), 0u);
    SEQAN_ASSERT(end(suffix_5) == end(str_1));
    SEQAN_ASSERT_EQ(endPosition(suffix_5), length(str_1));

    SEQAN_ASSERT(begin(suffix(str_1, 0)) == begin(str_1));
    SEQAN_ASSERT(end(suffix(str_1, 3)) == end(str_1));
    SEQAN_ASSERT_EQ(length(suffix(str_1, 0)), length(str_1));

    // str_1 = "begin middle end";
    // assign(suffix(str_1, 6), "to panic");
    // SEQAN_ASSERT_EQ(str_1, "begin to panic");

    // assign(suffix(str_1, 6), "the test", 9);
    // SEQAN_ASSERT_EQ(str_1, "begin the");

    // char str_2[200] = "begin middle end";
    // assign(suffix(str_2, 6), "to panic");
    // SEQAN_ASSERT(isEqual(str_2, "begin to panic"));

    // assign(suffix(str_2, 6), "the test", 9);
    // SEQAN_ASSERT(isEqual(str_2, "begin the"));

//____________________________________________________________________________
// test suffix iteration
/*
    str_1 = "begin middle end";
    goBegin(suffix_1);
    SEQAN_ASSERT(suffix_1 == str_1);

    goBegin(suffix_1, str_1);
    SEQAN_ASSERT(suffix_1 == str_1);
    SEQAN_ASSERT(atBegin(suffix_1));

    goEnd(suffix_1);
    SEQAN_ASSERT(suffix_1 == "d");

    goEnd(suffix_1, str_1);
    SEQAN_ASSERT(suffix_1 == "d");

    goNext(suffix_1);
    SEQAN_ASSERT(atEnd(suffix_1));

    goPrevious(suffix_1);
    SEQAN_ASSERT(atEnd(suffix_1));
*/

}

SEQAN_DEFINE_TEST(ticket317)
{
    // http://trac.mi.fu-berlin.de/seqan/ticket/317

    CharString text = "thisisatext";
    Infix<CharString>::Type sub1 = infix(text, 2, 8);
    Infix<Infix<CharString>::Type>::Type sub2 = infix(sub1, 1, 3);

    // Check our assumptions.
    SEQAN_ASSERT_EQ(sub1, "isisat");
    SEQAN_ASSERT_EQ(sub2, "si");

    // operator - (triggered the ticket)
    SEQAN_ASSERT_EQ(begin(sub1) - begin(text), 2);
    SEQAN_ASSERT_EQ(begin(sub2) - begin(text), 3);
    SEQAN_ASSERT_EQ(begin(sub2) - begin(sub1), 1);

    // operator == and !=
#   define TEST_IT(a, b) \
    SEQAN_ASSERT_EQ(a, b); \
    SEQAN_ASSERT_LEQ(a, b); \
    SEQAN_ASSERT_GEQ(a, b);

    TEST_IT(begin(sub1), begin(text) + 2);
    TEST_IT(begin(sub2), begin(text) + 3);
    TEST_IT(begin(sub2), begin(sub1) + 1);
#   undef TEST_IT

    SEQAN_ASSERT_NEQ(begin(sub1), begin(text));
    SEQAN_ASSERT_NEQ(begin(sub2), begin(text));
    SEQAN_ASSERT_NEQ(begin(sub2), begin(sub1));

    // operator <, >, <=, >=
#   define TEST_IT(a, b) \
    SEQAN_ASSERT_LT(a, b); \
    SEQAN_ASSERT_GT(b, a); \
    SEQAN_ASSERT_LEQ(a, b); \
    SEQAN_ASSERT_GEQ(b, a);

    TEST_IT(begin(sub1), end(text));
    TEST_IT(begin(text), begin(sub1));
    TEST_IT(begin(sub2), end(sub1));
    TEST_IT(begin(sub1), begin(sub2));
    TEST_IT(begin(sub2), end(text));
    TEST_IT(begin(text), begin(sub2));
#   undef TEST_IT
}

SEQAN_DEFINE_TEST(ticket848)
{
    // http://trac.mi.fu-berlin.de/seqan/ticket/848

    CharString text = "012345";

    Infix<CharString>::Type inf1 = infix(text, begin(text) + 1, begin(text) + 5);
    SEQAN_ASSERT_EQ(CharString("1234"), CharString(inf1));

    Infix<Infix<CharString>::Type>::Type inf2 = infix(inf1, begin(inf1) + 1, begin(inf1) + 3);
    SEQAN_ASSERT_EQ(CharString("23"), CharString(inf2));
}
