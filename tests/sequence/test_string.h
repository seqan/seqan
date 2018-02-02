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
// Test sequence default interface: Non-container objects are treated like
// containers of length 1
// ==========================================================================

struct DummyClass
{
};

SEQAN_DEFINE_TEST(Sequence_Interface)
{
//* ???Anti Default Sequences
    DummyClass const c = DummyClass();
    DummyClass d;
    DummyClass e;

    void* const nullPtr = 0;
    SEQAN_ASSERT_NEQ(getObjectId(c), getObjectId(e));
    SEQAN_ASSERT_NEQ(getObjectId(d), getObjectId(e));
    SEQAN_ASSERT_NEQ(getObjectId(c), nullPtr);
    SEQAN_ASSERT_NEQ(getObjectId(d), nullPtr);

    SEQAN_ASSERT_EQ(begin(c, Standard()), & c);        //begin
    SEQAN_ASSERT_EQ(begin(d, Standard()), & d);

    SEQAN_ASSERT_EQ(end(c, Standard()), & c + 1);    //end
    SEQAN_ASSERT_EQ(end(d, Standard()), & d + 1);

    SEQAN_ASSERT_EQ(beginPosition(c), 0u);            //beginPosition
    SEQAN_ASSERT_EQ(beginPosition(d), 0u);
    SEQAN_ASSERT_EQ(endPosition(c), (size_t)(end(c) - begin(c))); //endPosition
    SEQAN_ASSERT_EQ(endPosition(d), (size_t)(end(d) - begin(d)));

    SEQAN_ASSERT_EQ(iter(c, endPosition(c)), end(c));    //iter
    SEQAN_ASSERT_EQ(iter(d, endPosition(d)), end(d));

    SEQAN_ASSERT_EQ(& getValue(c, 10), & c);
    SEQAN_ASSERT_EQ(& getValue(d, 10), & d);

    SEQAN_ASSERT_EQ(length(c), 1u);
    SEQAN_ASSERT_EQ(length(d), 1u);

    SEQAN_ASSERT_EQ(capacity(c), 1u);
    SEQAN_ASSERT_EQ(capacity(d), 1u);

    SEQAN_ASSERT_NOT(empty(c));
    SEQAN_ASSERT_NOT(empty(d));

    Iterator<DummyClass, Standard>::Type it_1;
    it_1 = begin(d, Standard());

    goNext(it_1);
    goPrevious(it_1);
    SEQAN_ASSERT_EQ(it_1, begin(d, Standard()));
//*/

//____________________________________________________________________________
//test interfaces for assignment functions:
//semantics are tested somewhere else

    char str[100] = "hallo";

    // assign(str, str);
    // assign(str, infix(str, 0, 5));
    // assign(infix(str, 0, 5), str);
    // assign(infix(str, 0, 5), infix(str, 0, 5));
    // SEQAN_ASSERT(isEqual(str, "hallo"));

    // assign(str, str, 5);
    // assign(str, infix(str, 0, 5), 5);
    // assign(infix(str, 0, 5), str, 5);
    // assign(infix(str, 0, 5), infix(str, 0, 5), 5);

    // append(str, str);
    // append(str, infix(str, 0, 5));
    // append(infix(str, 0, 5), str);
    // append(infix(str, 0, 5), infix(str, 0, 5));

    // append(str, str, 5);
    // append(str, infix(str, 0, 5), 5);
    // append(infix(str, 0, 5), str, 5);
    // append(infix(str, 0, 5), infix(str, 0, 5), 5);

    // replace(str, 2, 3, str);
    // replace(str, 2, 3, infix(str, 0, 5));
    // replace(infix(str, 0, 5), 2, 3, str);
    // replace(infix(str, 0, 5), 2, 3, infix(str, 0, 5));

    // replace(str, 2, 3, str, 5);
    // replace(str, 2, 3, infix(str, 0, 5), 5);
    // replace(infix(str, 0, 5), 2, 3, str, 5);
    // replace(infix(str, 0, 5), 2, 3, infix(str, 0, 5), 5);

//____________________________________________________________________________

    String<char, Alloc<> > str2;

    Size<String<char, Alloc<> > >::Type cap = reserve(str2, 200);
    SEQAN_ASSERT_EQ(cap, capacity(str2));
    SEQAN_ASSERT_GEQ(cap, 200u);

    cap = reserve(str2, 400, Insist());
    SEQAN_ASSERT_EQ(cap, 400u);

    cap = reserve(str2, 1000, Limit());
    SEQAN_ASSERT_EQ(cap, capacity(str2));


    Size<String<char, Alloc<> > >::Type len = resize(str2, 100);
    SEQAN_ASSERT_LEQ(len, length(str2));
    SEQAN_ASSERT_LEQ(len, 100u);

    len = resize(str2, 150, 'C');
    SEQAN_ASSERT_LEQ(len, length(str2));
    SEQAN_ASSERT_LEQ(len, 150u);


    len = resizeSpace(str2, 100, 50, 100);
    SEQAN_ASSERT_EQ(len, 100u);
    SEQAN_ASSERT_EQ(length(str2), 200u);

    len = resizeSpace(str2, 100, 150, 200, 220);
    SEQAN_ASSERT_EQ(len, 70u);
    SEQAN_ASSERT_EQ(length(str2), 220u);

//____________________________________________________________________________
//test interfaces for appending different string types

    assign(str, "test mixed append");
    String<char> str3;
    clear(str3);
    append(str3, str);
    SEQAN_ASSERT(isEqual(str3, str));
    SEQAN_ASSERT(isEqual(str, str3));

    const char *str4 = "test const char*";
    clear(str3);
    append(str3, str4);
    SEQAN_ASSERT(isEqual(str4, str3));
    SEQAN_ASSERT(isEqual(str3, str4));

    clear(str3);
    append(str3, "test");
    SEQAN_ASSERT(isEqual("test", str3));
    SEQAN_ASSERT(isEqual(str3, "test"));
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//Test semantics of assignment functions like assign, append, ...

template <typename TExpand, typename TMe>
void Test_String_Base_Assignments(TMe & str)
{
    assign(str, "goldfishgoldfishgoldfishgoldfishgoldfish", TExpand());
    SEQAN_ASSERT(isEqual(str, "goldfishgoldfishgoldfishgoldfishgoldfish"));

    assign(str, "x", TExpand());
    SEQAN_ASSERT(isEqual(str,"x"));

    assign(str, "goldfishgoldfishgoldfishgoldfishgoldfish", 8, TExpand());
    SEQAN_ASSERT(isEqual(str,"goldfish"));

    assign(str, "x", 8, TExpand());
    SEQAN_ASSERT(isEqual(str,"x"));

    //test append
    assign(str, "goldfish");

    append(str, "goldfishgoldfishgoldfishgoldfish", TExpand());
    SEQAN_ASSERT(isEqual(str,"goldfishgoldfishgoldfishgoldfishgoldfish"));

    assign(str, "goldfish");
    append(str, "", TExpand());
    SEQAN_ASSERT(isEqual(str,"goldfish"));

    append(str, "goldfish", 4, TExpand());
    SEQAN_ASSERT(isEqual(str,"gold"));

    append(str, "goldfish", 8, TExpand());
    SEQAN_ASSERT(isEqual(str,"goldgold"));

    append(str, "goldfish", 16, TExpand());
    SEQAN_ASSERT(isEqual(str,"goldgoldgoldfish"));

    //test replace

    str = "goldfish";
    replace(str, 2, 4, "ndel chips and ", TExpand());
    SEQAN_ASSERT_EQ(str, "gondel chips and fish");

    replace(str, 7, 16, "is", TExpand());
    SEQAN_ASSERT_EQ(str, "gondel is fish");

    replace(str, 6, 10, "", TExpand());
    SEQAN_ASSERT_EQ(str, "gondelfish");

    replace(str, 2, 2, "ld in my ma", TExpand());
    SEQAN_ASSERT_EQ(str, "gold in my mandelfish");

    replace(str, (size_t) 0, 0, "here is ", TExpand());
    SEQAN_ASSERT_EQ(str, "here is gold in my mandelfish");

    replace(str, length(str), length(str), " and silver", TExpand());
    SEQAN_ASSERT_EQ(str, "here is gold in my mandelfish and silver");

    replace(str, 8, length(str), "nothing", TExpand());
    SEQAN_ASSERT_EQ(str, "here is nothing");


    assignValue(str, 1, 'a');        //assignValue
    SEQAN_ASSERT_EQ(str, "hare is nothing");

    moveValue(str, 1, 'e');            //moveValue
    SEQAN_ASSERT_EQ(str, "here is nothing");

    appendValue(str, '!');            //appendValue
    SEQAN_ASSERT_EQ(str, "here is nothing!");

    insertValue(str, 7, 't');        //insertValue
    SEQAN_ASSERT_EQ(str, "here ist nothing!");

    erase(str, 7);                    //erase
    SEQAN_ASSERT_EQ(str, "here is nothing!");
    erase(str, 1, 3);
    SEQAN_ASSERT_EQ(str, "he is nothing!");

    TMe str2 = "another string";
    move(str, str2);                //move
    SEQAN_ASSERT_EQ(str, "another string");
}

SEQAN_DEFINE_TEST(String_Base)
{
    String<char> str1("hello");
    String<char> const str2("HELLO");

    SEQAN_ASSERT_EQ(getValue(str1, 0), 'h');
    SEQAN_ASSERT_EQ(getValue(str1, 1), 'e');
    SEQAN_ASSERT_EQ(getValue(str2, 1), 'E');

    clear(str1);
    SEQAN_ASSERT_EQ(length(str1), 0u);

//____________________________________________________________________________
//note: the following tests destroy the same items several times
//note: most of the tests only check syntax not semantic.
//semantic should be tested via the calling assignment functions
/*
    str1 = "hello";
    _clearSpace(str1, 200, Exact());
    SEQAN_ASSERT(length(str1) == 200);
    SEQAN_ASSERT(capacity(str1) == 200);

    _clearSpace(str1, 200, Limit());
    _clearSpace(str1, 300, Generous());
    _clearSpace(str1, 200, Insist());

    _clearSpace(str1, 300, 100, Exact());
    SEQAN_ASSERT(length(str1) == 100);
    _clearSpace(str1, 300, 10000, Limit());
    _clearSpace(str1, 400, 100, Generous());
    _clearSpace(str1, 300, 100, Insist());

    _clearSpace(str1, 300, 10, 20, Exact());
    _clearSpace(str1, 300, 10, 20, Limit());
    _clearSpace(str1, 300, 10, 20, Generous());
    _clearSpace(str1, 100, 10, 200, Insist());

    _clearSpace(str1, 300, begin(str1), end(str1), Exact());
    _clearSpace(str1, 300, begin(str1), end(str1), Limit());
    _clearSpace(str1, 300, begin(str1), end(str1), Generous());
    _clearSpace(str1, 100, begin(str1), end(str1), Insist());

    _clearSpace(str1, 300, 10, 20, 100, Exact());
    _clearSpace(str1, 300, 10, 20, 100, Limit());
    _clearSpace(str1, 300, 10, 20, 100, Generous());
    _clearSpace(str1, 100, 10, 50, 100, Insist());

    _clearSpace(str1, 300, 10, 20, 20000, Exact());
    _clearSpace(str1, 300, 10, 20, 20000, Limit());
    _clearSpace(str1, 300, 10, 20, 20000, Generous());
    _clearSpace(str1, 100, 10, 50, 20000, Insist());
*/

/*
    _clearSpace(str1, 300, begin(str1), end(str1), 100, Exact());
    _clearSpace(str1, 300, begin(str1), end(str1), 100, Limit());
    _clearSpace(str1, 300, begin(str1), end(str1), 100, Generous());
    _clearSpace(str1, 100, begin(str1), end(str1), 100, Insist());
*/

//____________________________________________________________________________
// test assignment functions

    Test_String_Base_Assignments<Exact>(str1);
    Test_String_Base_Assignments<Generous>(str1);
    reserve(str1, 10000);
    Test_String_Base_Assignments<Insist>(str1);

    str1 += str2;
    str1 += str1;

//* ???Anti Default Sequences
    str1 += 'x';
//*/
//____________________________________________________________________________

    String<char> str3;
    resize(str3, 0);
    SEQAN_ASSERT_EQ(length(str3), 0u);

    resize(str3, 20);
    SEQAN_ASSERT_EQ(length(str3), 20u);
    resize(str3, 10);
    SEQAN_ASSERT_EQ(length(str3), 10u);
    resize(str3, 20);
    SEQAN_ASSERT_EQ(length(str3), 20u);
    resize(str3, 200);
    SEQAN_ASSERT_EQ(length(str3), 200u);

    resize(str3, 100, 'x');
    SEQAN_ASSERT_EQ(length(str3), 100u);
    resize(str3, 200, 'x');
    SEQAN_ASSERT_EQ(length(str3), 200u);
    resize(str3, 400, 'y');
    SEQAN_ASSERT_EQ(length(str3), 400u);

//____________________________________________________________________________

    str3 = "abc";
    SEQAN_ASSERT_LT(str3, "abcd");
}

//////////////////////////////////////////////////////////////////////////////
//test some basic string features

template <typename TMe>
void TestStringBasics()
{
    //test default ctor
    TMe str1;
    SEQAN_ASSERT_EQ(str1, "");
    SEQAN_ASSERT_EQ(length(str1), 0u);
    SEQAN_ASSERT_LEQ(length(str1), capacity(str1));

    //test assignment ctor (with length == 0)
    TMe str2 = str1;
    SEQAN_ASSERT_EQ(str1, str2);
    SEQAN_ASSERT_EQ(length(str1), 0u);
    SEQAN_ASSERT_LEQ(length(str1), capacity(str1));

    //test assignment with char const []
    str1 = "hamster";
    SEQAN_ASSERT_EQ(str1, "hamster");
    SEQAN_ASSERT_EQ(length(str1), 7u);
    SEQAN_ASSERT_LEQ(length(str1), capacity(str1));

    //test assignment with char *
    char * s1 = (char *) "goldfish";
    str1 = s1;
    SEQAN_ASSERT_EQ(str1, "goldfish");
    SEQAN_ASSERT_EQ(length(str1), 8u);
    SEQAN_ASSERT_LEQ(length(str1), capacity(str1));

    //test assignment ctor (with length > 0)
    TMe str3 = str1;
    SEQAN_ASSERT_EQ(str3, str1);
    SEQAN_ASSERT_EQ(length(str3), 8u);
    SEQAN_ASSERT_LEQ(length(str3), capacity(str3));

    //test independency
    {
        TMe str4 = "hamster";
        str3 = str4;
        str4 = "...";
    }
    SEQAN_ASSERT_EQ(str3, "hamster");
    SEQAN_ASSERT_EQ(length(str3), 7u);
    SEQAN_ASSERT_LEQ(length(str3), capacity(str3));

    typename Size<TMe>::Type len = length(str3);
    SEQAN_ASSERT_EQ(len, 7u);

    //test begin and end
    //SEQAN_ASSERT_EQ(end(str3), begin(str3) + 7);
    SEQAN_ASSERT(end(str3) == begin(str3) + 7);
    typename Iterator<TMe, Rooted>::Type str3_begin = begin(str3, Rooted());
    *str3_begin = 'X';
    SEQAN_ASSERT_EQ(str3, "Xamster");

    //test at
    value(str3, 1) = 'Y';
    str3[2] = 'Z';
    int i1 = 3;
    value(str3, i1) = 'A';
    i1 = 4;
    str3[i1] = 'B';
    SEQAN_ASSERT_EQ(str3, "XYZABer");
    SEQAN_ASSERT_EQ(getValue(str3, 5), 'e');
    SEQAN_ASSERT_EQ(str3[6], 'r');

    //test clear and empty
    SEQAN_ASSERT_NOT(empty(str3));
    clear(str3);
    SEQAN_ASSERT(empty(str3));
    SEQAN_ASSERT_EQ(length(str3), 0u);
    SEQAN_ASSERT_EQ(length(str3), 0u);
}

//////////////////////////////////////////////////////////////////////////////
//test some basic string features for strings that can change length
//note: capacity of TMe strings must be >= 200

template <typename TMe>
void TestStringResize()
{
    TMe str1;
    resize(str1, 5, 1);
    SEQAN_ASSERT_EQ(length(str1), 5u);
    SEQAN_ASSERT_EQ((int)ordValue(getValue(str1, 1)), 1);
    SEQAN_ASSERT_EQ((int)ordValue(getValue(str1, 3)), 1);

    resize(str1, 50, 2);
    SEQAN_ASSERT_EQ(length(str1), 50u);
    SEQAN_ASSERT_EQ((int)ordValue(getValue(str1, 4)), 1);
    SEQAN_ASSERT_EQ((int)ordValue(getValue(str1, 5)), 2);

    resize(str1, 100, 3);
    SEQAN_ASSERT_EQ(length(str1), 100u);
    SEQAN_ASSERT_EQ((int)ordValue(getValue(str1, 49)), 2);
    SEQAN_ASSERT_EQ((int)ordValue(getValue(str1, 50)), 3);
    SEQAN_ASSERT_EQ((int)ordValue(getValue(str1, 51)), 3);
    SEQAN_ASSERT_EQ((int)ordValue(getValue(str1, 99)), 3);

}

//////////////////////////////////////////////////////////////////////////////

SEQAN_DEFINE_TEST(String_Alloc)
{
    TestStringBasics<String<char> >();
    TestStringResize<String<char> >();

    String<char, Alloc<> > str1 = "hello";
    SEQAN_ASSERT_EQ(str1[1], 'e');

    String<char, Alloc<> > const str2 = "hello";
    SEQAN_ASSERT_EQ(str2[4], 'o');
    SEQAN_ASSERT_GEQ(capacity(str2), length(str2));

    String<char, Alloc<> > str3;
    move(str3, str1);
    SEQAN_ASSERT_EQ(str3, "hello");
    SEQAN_ASSERT_EQ(length(str1), 0u);

    String<char, Alloc<> > left = "left";
    String<char, Alloc<> > right = "right";
    seqan::swap(left, right);
    SEQAN_ASSERT_EQ(left, "right");
    SEQAN_ASSERT_EQ(right, "left");
}

//////////////////////////////////////////////////////////////////////////////

SEQAN_DEFINE_TEST(String_Array)
{
    TestStringBasics<String<char, Array<100> > >();

    String<char, Array<100> > str1 = "hello";
    SEQAN_ASSERT_EQ(str1[1], 'e');

    String<char, Array<100> > const str2 = "hello";
    SEQAN_ASSERT_EQ(str2[4], 'o');
    SEQAN_ASSERT_GEQ(capacity(str2), length(str2));
}

//////////////////////////////////////////////////////////////////////////////

SEQAN_DEFINE_TEST(String_Stack)
{
    TestStringBasics<String<char, Block<3> > >();

    String<char, Block<3> > str1 = "hello";
    SEQAN_ASSERT_EQ(str1[1], 'e');

    String<char, Block<3> > const str2 = "hello";
    SEQAN_ASSERT_EQ(str2[4], 'o');
    SEQAN_ASSERT_GEQ(capacity(str2), length(str2));

    resize(str1, 30);
}

//////////////////////////////////////////////////////////////////////////////

SEQAN_DEFINE_TEST(String_Packed)
{
    TestStringBasics<String<char, Packed<> > >();
    TestStringBasics<String<Dna, Packed<> > >();
    TestStringBasics<String<Dna5, Packed<> > >();

    TestStringResize<String<char, Packed<> > >();
    TestStringResize<String<Dna, Packed<> > >();
    TestStringResize<String<Dna5, Packed<> > >();

    typedef String<Dna5, Packed<> > TPackedString;
    Dna5String str1;
    TPackedString str2;
    Dna5String motif = "AGAGCCCTGACACACTGATTATGACACGTAAGACTCTGAGAGTCCGTGTCTCTATAG";

    for (int i = 0; i < 5; ++i)
    {
        append(str1, "GATTACATATA", Exact());
        append(str2, "GATTACATATA", Exact());
        SEQAN_ASSERT_EQ(str1, str2);
    }

    for (int i = 0; i < PackedTraits_<TPackedString>::VALUES_PER_HOST_VALUE; ++i)
    {
        appendValue(str1, (Dna5)(i % 5), Exact());
        appendValue(str2, (Dna5)(i % 5), Exact());
        SEQAN_ASSERT_EQ(str1, str2);
    }

    resize(str1, 131, 'G', Exact());
    resize(str2, 131, 'G', Exact());
    SEQAN_ASSERT_EQ(str1, str2);

    insert(str1, 19, "ACGTACGT", Exact());
    insert(str2, 19, "ACGTACGT", Exact());
    SEQAN_ASSERT_EQ(str1, str2);

    for (int i = 0; i < PackedTraits_<TPackedString>::VALUES_PER_HOST_VALUE; ++i)
    {
        Dna5String str1_ = str1;
        TPackedString str2_ = str2;

        insertValue(str1_, 64-i, (Dna5)(i % 5), Exact());
        insertValue(str2_, 64-i, (Dna5)(i % 5), Exact());
        SEQAN_ASSERT_EQ(str1_, str2_);
    }

    for (int i = 0; i < PackedTraits_<TPackedString>::VALUES_PER_HOST_VALUE; ++i)
    {
        Dna5String str1_ = str1;
        TPackedString str2_ = str2;

        insert(str1_, 64-i, "ACGTACGT", Exact());
        insert(str2_, 64-i, "ACGTACGT", Exact());
        SEQAN_ASSERT_EQ(str1_, str2_);
    }

    for (int i = 0; i < PackedTraits_<TPackedString>::VALUES_PER_HOST_VALUE; ++i)
        for (int j = 0; j < PackedTraits_<TPackedString>::VALUES_PER_HOST_VALUE; ++j)
        {
            Dna5String str1_ = str1;
            TPackedString str2_ = str2;

            insert(str1_, 64-i, suffix(motif,j), Exact());
            insert(str2_, 64-i, suffix(motif,j), Exact());
            SEQAN_ASSERT_EQ(str1_, str2_);
        }

    for (int i = 0; i < PackedTraits_<TPackedString>::VALUES_PER_HOST_VALUE; ++i)
    {
        Dna5String str1_ = str1;
        TPackedString str2_ = str2;

        erase(str1_, 64-i);
        erase(str2_, 64-i);
        SEQAN_ASSERT_EQ(str1_, str2_);
    }

    for (int i = 0; i < PackedTraits_<TPackedString>::VALUES_PER_HOST_VALUE; ++i)
        for (int j = 0; j < PackedTraits_<TPackedString>::VALUES_PER_HOST_VALUE; ++j)
        {
            Dna5String str1_ = str1;
            TPackedString str2_ = str2;

            erase(str1_, 64-i, 65-i+j);
            erase(str2_, 64-i, 65-i+j);
            SEQAN_ASSERT_EQ(str1_, str2_);
        }
}

//////////////////////////////////////////////////////////////////////////////

SEQAN_DEFINE_TEST(String_Pointer)
{
    char str1[200] = "hello";

    SEQAN_ASSERT_EQ(getValue(str1, 0), 'h');
    SEQAN_ASSERT_EQ(getValue(str1, 1), 'e');
    SEQAN_ASSERT_EQ(getValue("hello", 1), 'e');

    int str2[100] = { 10, 20, 30, 0 };
    int const str3[100] = { 10, 20, 30, 0 };

    SEQAN_ASSERT_EQ(length(str2), 3u);
    SEQAN_ASSERT_EQ(length(str3), 3u);

    SEQAN_ASSERT_EQ(length(str1), 5u);
    SEQAN_ASSERT_EQ(length("hello"), 5u);

    clear(str1);
    SEQAN_ASSERT_EQ(length(str1), 0u);
    SEQAN_ASSERT(empty(str1));

    SEQAN_ASSERT_EQ(reserve(str1, 100, Insist()), 100u);
    SEQAN_ASSERT_EQ(reserve(str1, 100, Limit()), capacity(str1));

    resize(str1, 20, 'A');

    SEQAN_ASSERT(isEqual(str1, "AAAAAAAAAAAAAAAAAAAA"));

    resize(str1, 10);
    SEQAN_ASSERT(isEqual(str1, "AAAAAAAAAA"));

//____________________________________________________________________________
// compare operators

//operators disabled due to ambiguity pointer/iterator vs. c-style string
/*
    assign(str1, "hello");
    String<char> str4 = str1;
    SEQAN_ASSERT(str1 == str4);

    str1 += str4;
    SEQAN_ASSERT(isEqual(str1, "hellohello"));

    SEQAN_ASSERT(str1 != str4);
    SEQAN_ASSERT(!isNotEqual(str1, "hellohello"));

    SEQAN_ASSERT(!(str1 < str4));
    SEQAN_ASSERT(!isLess(str1, str4));

    SEQAN_ASSERT(!(str1 <= str4));
    SEQAN_ASSERT(!isLessOrEqual(str1, str4));

    SEQAN_ASSERT(str1 > str4);
    SEQAN_ASSERT(isGreater(str1, str4));

    SEQAN_ASSERT(str1 >= str4);
    SEQAN_ASSERT(isGreaterOrEqual(str1, str4));
*/
}

//////////////////////////////////////////////////////////////////////////////

SEQAN_DEFINE_TEST(String_CStyle)
{
    String<char, CStyle> str1;
    char strq [200] = "hello seqan";

    String<char, CStyle> str2(strq);
    SEQAN_ASSERT_EQ(str2, strq);
    SEQAN_ASSERT_EQ(strlen(str2), strlen(strq));
    SEQAN_ASSERT_EQ(begin(str2), begin(strq));
    SEQAN_ASSERT_EQ(end(str2), end(strq));

    SEQAN_ASSERT_EQ(capacity(str2), capacity(strq));

    String<char, CStyle> str3(str2);

    String<char const, CStyle> str4("a const string");

    String<char> stra("alloc string");
    String<char, CStyle>str5(stra);
    SEQAN_ASSERT_EQ(str5, stra);

    String<char> const strac("const alloc string");
    String<char, CStyle> str7(strac);
    SEQAN_ASSERT_EQ(str7, strac);

    str2 = stra;
    SEQAN_ASSERT_EQ(str2, stra);

    str2 = strac;
    SEQAN_ASSERT_EQ(str2, strac);

    str1 = str5;
    SEQAN_ASSERT_EQ(str1, str5);

    char * cp1 = str1;
    void* const nullPtr = 0;
    SEQAN_ASSERT_NEQ(cp1, nullPtr);

    String<char, CStyle> const str6(strq);
    char const * cp2 = str6;
    SEQAN_ASSERT_NEQ(cp2, nullPtr);

    str1 = str6;
    SEQAN_ASSERT_EQ(str1, str6);

    String<char, CStyle> str8(str6);
    SEQAN_ASSERT_EQ(str8, str6);

    str2 = strac;
    SEQAN_ASSERT_EQ(str2, strac);
    SEQAN_ASSERT_NEQ(getObjectId(str2), getObjectId(strac));

    clear(str2);
    SEQAN_ASSERT_EQ(length(str2), 0u);

    assign(str2, stra);
    SEQAN_ASSERT_EQ(str2, stra);
    SEQAN_ASSERT_EQ(getObjectId(str2), getObjectId(stra));
    SEQAN_ASSERT_EQ(str2, toCString(stra));

//  weese: Is there a reason why CStyle strings don't support append?
//    append(str2, " whatever");
//    SEQAN_ASSERT_EQ(str2, "const alloc string whatever");

    str2 = strac;
    SEQAN_ASSERT_EQ(str2, strac);
    SEQAN_ASSERT_NEQ(getObjectId(str2), getObjectId(strac));

    String<Dna> str_dna("acgt");
    str2 = str_dna;
    SEQAN_ASSERT_EQ(str2, str_dna);
    SEQAN_ASSERT_NEQ(getObjectId(str2), getObjectId(strac));

    String<Dna, CStyle> str9(str_dna);
    SEQAN_ASSERT_EQ(str2, str_dna);

    char * strp = (char *) "this is a long array of chars";
    create(str2, strp);
    SEQAN_ASSERT_EQ(str2, strp);
    SEQAN_ASSERT_NEQ(getObjectId(str2), getObjectId(strp));

    assign(str2, strp);
    SEQAN_ASSERT_EQ(str2, strp);
    SEQAN_ASSERT_EQ(getObjectId(str2), getObjectId(strp));

    str2 = "hello";
    String<char, CStyle > str10;
    move(str10, str2);
    SEQAN_ASSERT_EQ(str10, "hello");
    SEQAN_ASSERT_EQ(length(str2), 0u);
}

//////////////////////////////////////////////////////////////////////////////

SEQAN_DEFINE_TEST(Segment)
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

    SEQAN_ASSERT_EQ(begin(infix_5), begin(str_1));
    SEQAN_ASSERT_EQ(beginPosition(infix_5), 0u);
    SEQAN_ASSERT_EQ(end(infix_5), end(str_1));
    SEQAN_ASSERT_EQ(endPosition(infix_5), length(str_1));

    SEQAN_ASSERT(begin(infix(str_1, 0, length(str_1))) == begin(str_1));
    SEQAN_ASSERT(end(infix(str_1, 0, length(str_1))) == end(str_1));
    SEQAN_ASSERT(length(infix(str_1, 0, length(str_1))) == length(str_1));

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
// suffix

    str_1 = "this is a string";

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
    SEQAN_ASSERT(endPosition(suffix_5) == length(str_1));

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

//____________________________________________________________________________
// compare operators

    str_1 = "hellohello";
    Infix<String<char> >::Type infix_6(str_1, 0, 10);

//     SEQAN_ASSERT_EQ(infix_6, str_1);
//
//     infix_6 += str_1;
    SEQAN_ASSERT(isEqual(infix_6, "hellohello"));

    SEQAN_ASSERT_NEQ(infix_6, "bla");
    SEQAN_ASSERT(!isNotEqual(infix_6, "hellohello"));

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

SEQAN_DEFINE_TEST(Std_String)
{
//____________________________________________________________________________

    std::string str_1("hamster");
    SEQAN_ASSERT_EQ(0, 0);
    size_t len1 = seqan::end(str_1) - seqan::begin(str_1);
    size_t len2 = length(str_1);
    SEQAN_ASSERT_EQ(len1, len2);

    std::string const str_2("goldfish");
    len1 = seqan::end(str_2) - seqan::begin(str_2);
    len2 = length(str_2);
    SEQAN_ASSERT_EQ(len1, len2);

    SEQAN_ASSERT_EQ(getValue(str_1, 0), 'h');
    SEQAN_ASSERT_EQ(getValue(str_1, 1), 'a');
    SEQAN_ASSERT_EQ(getValue(str_2, 1), 'o');

    SEQAN_ASSERT_LEQ(length(str_1), capacity(str_1));

    clear(str_1);
    SEQAN_ASSERT(empty(str_1));

    reserve(str_1, 200);
    SEQAN_ASSERT_GEQ(capacity(str_1), 200u);

    resize(str_1, 100);
    SEQAN_ASSERT_EQ(length(str_1), 100u);
    resize(str_1, 150, 'x');
    SEQAN_ASSERT_EQ(length(str_1), 150u);

//____________________________________________________________________________

}


//////////////////////////////////////////////////////////////////////////////

SEQAN_DEFINE_TEST(Lexical)
{
    Lexical<> lex1;
    compare(lex1, "abc", "abcd");

    SEQAN_ASSERT(isLess(lex1));
    SEQAN_ASSERT(isLess(lex1, TagPrefixLess()));
    SEQAN_ASSERT_NOT(isLess(lex1, TagPrefixGreater()));

    SEQAN_ASSERT(isLessOrEqual(lex1));
    SEQAN_ASSERT(isLessOrEqual(lex1, TagPrefixLess()));
    SEQAN_ASSERT_NOT(isLessOrEqual(lex1, TagPrefixGreater()));

    Lexical<> lex2(lex1);

    SEQAN_ASSERT(isPrefix(lex2));
    SEQAN_ASSERT_NOT(hasPrefix(lex2));

    Lexical<> lex3("abcd", "abc");

    SEQAN_ASSERT(isGreater(lex3));
    SEQAN_ASSERT(isGreater(lex3, TagPrefixLess()));
    SEQAN_ASSERT_NOT(isGreater(lex3, TagPrefixGreater()));

    SEQAN_ASSERT(isGreaterOrEqual(lex3));
    SEQAN_ASSERT(isGreaterOrEqual(lex3, TagPrefixLess()));
    SEQAN_ASSERT_NOT(isGreaterOrEqual(lex3, TagPrefixGreater()));

    lex2 = lex3;
    SEQAN_ASSERT_NOT(isPrefix(lex2));
    SEQAN_ASSERT(hasPrefix(lex2));

    String<char> str1 = "alpha";

    SEQAN_ASSERT(isEqual(str1, "alpha"));
    SEQAN_ASSERT_EQ(str1, "alpha");

    SEQAN_ASSERT(isNotEqual(str1, "beta"));
    SEQAN_ASSERT_NEQ(str1, "beta");

    SEQAN_ASSERT(isLess(str1, "beta"));
    SEQAN_ASSERT_LT(str1, "beta");
    SEQAN_ASSERT_NOT(isLess(str1, "aaa"));

    SEQAN_ASSERT(isLessOrEqual(str1, "beta"));
    SEQAN_ASSERT_LEQ(str1, "beta");
    SEQAN_ASSERT_NOT(isLessOrEqual(str1, "aaa"));
    SEQAN_ASSERT(isLessOrEqual(str1, "alpha"));

    SEQAN_ASSERT(isGreater(str1, "aaa"));
    SEQAN_ASSERT_GT(str1, "aaa");
    SEQAN_ASSERT_NOT(isGreater(str1, "beta"));

    SEQAN_ASSERT(isGreaterOrEqual(str1, "aaa"));
    SEQAN_ASSERT_GEQ(str1, "aaa");
    SEQAN_ASSERT_NOT(isGreaterOrEqual(str1, "beta"));
    SEQAN_ASSERT(isGreaterOrEqual(str1, "alpha"));

    SEQAN_ASSERT(isPrefix(str1, "alpha romeo"));
    SEQAN_ASSERT(isPrefix(str1, "alpha"));
    SEQAN_ASSERT_NOT(isPrefix(str1, ""));
    SEQAN_ASSERT_NOT(isPrefix(str1, "alp"));
    SEQAN_ASSERT_NOT(isPrefix(str1, "b"));

    SEQAN_ASSERT(hasPrefix(str1, "alp"));
    SEQAN_ASSERT(hasPrefix(str1, "alpha"));
    SEQAN_ASSERT(hasPrefix(str1, ""));
    SEQAN_ASSERT_NOT(hasPrefix(str1, "alphas"));
    SEQAN_ASSERT_NOT(hasPrefix(str1, "b"));

    SEQAN_ASSERT_EQ(lcpLength("hello", "hellmaker"), 4u);
    SEQAN_ASSERT_EQ(lcpLength("", "not empty"), 0u);
    SEQAN_ASSERT_EQ(lcpLength("hello", "hello"), 5u);
    SEQAN_ASSERT_EQ(lcpLength("hello", "good evening"), 0u);

    SEQAN_ASSERT_EQ(lcpLength("hello", 'h'), 1u);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TTarget, typename TSource, typename TExpand>
void Test_Assignments_Combinatoric(TTarget & target, TSource source, Tag<TExpand> tag, size_t limit = ~0)
{
    // assign(target, source, tag);
    // SEQAN_ASSERT_EQ(infix(source, 0, length(target)), target);

    // assign(target, source, limit, tag);
    // SEQAN_ASSERT_EQ(infix(source, 0, length(target)), target);

    TSource const source_const(source);
    // assign(target, source_const, tag);
    // SEQAN_ASSERT_EQ(infix(source, 0, length(target)), target);

    // assign(target, source_const, limit, tag);
    // SEQAN_ASSERT_EQ(infix(source, 0, length(target)), target);

    typename Size<TTarget>::Type len = length(target);

    append(target, source, tag);
    if (len < length(target))
    {
        SEQAN_ASSERT_EQ(infix(source, 0, length(target) - len), infix(target, len, length(target)) );
    }

    len = length(target);
    append(target, source, limit, tag);

    if (len < length(target))
    {
        SEQAN_ASSERT_EQ(infix(source, 0, length(target) - len), infix(target, len, length(target)) );
    }

    append(target, source_const, tag);
    append(target, source_const, limit, tag); //p

    assign(target, "my miss is a hippi");
    replace(target, 9, 11, source, tag);
    SEQAN_ASSERT_EQ(infix(target, 0, 9), "my miss i");

    len = length(target);
    if (len >= length(source) + 9)
    {
        SEQAN_ASSERT_EQ(infix(target, 9, 9 + length(source)), source);
    }

    assign(target, "my miss is a hippi");
    replace(target, 9, 11, source, limit, tag);

    assign(target, "my miss is a hippi");
    replace(target, 9, 11, source, tag);
    SEQAN_ASSERT_EQ(infix(target, 0, 9), "my miss i");

    len = length(target);
    if (len >= length(source) + 9)
    {
        SEQAN_ASSERT_EQ(infix(target, 9, 9 + length(source)), source);
    }

    assign(target, "my miss is a hippi");
    replace(target, 9, 11, source, limit, tag);

    assign(target, "my miss is a hippi");
    replace(target, 9, 11, source_const, tag);

    assign(target, "my miss is a hippi");
    replace(target, 9, 11, source_const, limit, tag);
}

SEQAN_DEFINE_TEST(Combinatoric)
{
    String<char> str1("hello");
    String<char> str2("this is test");
    String<char> const str3("this is const string");

    Test_Assignments_Combinatoric(str1, str2, Exact());
    Test_Assignments_Combinatoric(str1, str2, Generous());
    Test_Assignments_Combinatoric(str1, str2, Limit());

    Test_Assignments_Combinatoric(str1, str2, Exact(), 3);
    Test_Assignments_Combinatoric(str1, str2, Generous(), 3);
    Test_Assignments_Combinatoric(str1, str2, Limit(), 3);

    Test_Assignments_Combinatoric(str1, str1, Exact());
    Test_Assignments_Combinatoric(str1, str1, Generous());
    Test_Assignments_Combinatoric(str1, str1, Limit());

    Test_Assignments_Combinatoric(str1, str1, Exact(), 3);
    Test_Assignments_Combinatoric(str1, str1, Generous(), 3);
    Test_Assignments_Combinatoric(str1, str1, Limit(), 3);

    Test_Assignments_Combinatoric(str1, str3, Exact());
    Test_Assignments_Combinatoric(str1, str3, Generous());
    Test_Assignments_Combinatoric(str1, str3, Limit());

    Test_Assignments_Combinatoric(str1, str3, Exact(), 3);
    Test_Assignments_Combinatoric(str1, str3, Generous(), 3);
    Test_Assignments_Combinatoric(str1, str3, Limit(), 3);

//____________________________________________________________________________

    String<char> str4("hello");
    Segment<String<char> > str5(str4, 1, 2);

    Test_Assignments_Combinatoric(str1, str5, Exact());
    Test_Assignments_Combinatoric(str1, str5, Generous());
    Test_Assignments_Combinatoric(str1, str5, Limit());

    Test_Assignments_Combinatoric(str1, str5, Exact(), 3);
    Test_Assignments_Combinatoric(str1, str5, Generous(), 3);
    Test_Assignments_Combinatoric(str1, str5, Limit(), 3);

    char str6 = 'x';
    Test_Assignments_Combinatoric(str1, str6, Exact());
    Test_Assignments_Combinatoric(str1, str6, Generous());
    Test_Assignments_Combinatoric(str1, str6, Limit());

    char str7[800] = "hello again";
    reserve(str1, 10000);
    Test_Assignments_Combinatoric(str1, str7, Insist());
    Test_Assignments_Combinatoric(str7, str4, Insist());
    Test_Assignments_Combinatoric(str7, str7, Insist());
    Test_Assignments_Combinatoric(str7, "sisyphos", Insist());

    Test_Assignments_Combinatoric(str7, 'c', Insist());

    Test_Assignments_Combinatoric(str1, str7, Insist(), 3);
    Test_Assignments_Combinatoric(str7, str4, Insist(), 3);
    Test_Assignments_Combinatoric(str7, str7, Insist(), 3);
//____________________________________________________________________________

    // assign(str7, "begin middle end");
    // Segment<char *> infix_1(str7, 6, 12);
    // SEQAN_ASSERT_EQ(infix_1, "middle");

    // Test_Assignments_Combinatoric(infix_1, str1, Insist());
    // SEQAN_ASSERT_GEQ(lcpLength(end(infix_1, Standard()), " end"), 4u);
    // SEQAN_ASSERT_EQ(beginPosition(infix_1), 6u);
    // SEQAN_ASSERT_EQ(infix(str7, 0, 6), "begin ");

    // Test_Assignments_Combinatoric(infix_1, str1, Insist(), 10);
    // str4 = "begin middle end";
    // Infix<String<char> >::Type infix_2(str4, 6, 12);
    // SEQAN_ASSERT_EQ(infix_2, "middle");

    // Test_Assignments_Combinatoric(infix_2, str1, Exact());
    // Test_Assignments_Combinatoric(infix_2, str1, Generous());
    // Test_Assignments_Combinatoric(infix_2, str1, Limit());

    // SEQAN_ASSERT_EQ(beginPosition(infix_2), 6u);
    // SEQAN_ASSERT_EQ(infix(str4, 0, 6), "begin ");

    // Test_Assignments_Combinatoric(infix_2, str1, Exact(), 10);
    // Test_Assignments_Combinatoric(infix_2, str1, Generous(), 10);
    // Test_Assignments_Combinatoric(infix_2, str1, Limit(), 10);

//____________________________________________________________________________

    str1 = "seqan string";
    std::string str8("i am the standard");
    Test_Assignments_Combinatoric(str8, str1, Generous());
    Test_Assignments_Combinatoric(str8, str1, Limit());

    Test_Assignments_Combinatoric(str8, str1, Generous(), 10);
    Test_Assignments_Combinatoric(str8, str1, Limit(), 10);

    str8 = "standard string";
    Test_Assignments_Combinatoric(str1, str8, Generous());
    Test_Assignments_Combinatoric(str1, str8, Generous(), 10);

    str8 = "standard string";
    Test_Assignments_Combinatoric(str7, str8, Insist());
    Test_Assignments_Combinatoric(str7, str8, Insist(), 10);

//____________________________________________________________________________

    str1 = "this is a test string";
    String<char, Array<100> > str9;
    Test_Assignments_Combinatoric(str9, str1, Limit());
    Test_Assignments_Combinatoric(str9, str1, Limit(), 10);
}

SEQAN_DEFINE_TEST(ticket901)
{
    using namespace seqan;

    String<char> source = "acgtgcat";
    String<Dna> target;
    // The in-place move conversion.
    move(target, source);

    SEQAN_ASSERT_EQ(length(source), 0u);
    SEQAN_ASSERT_EQ(target, DnaString("acgtgcat"));
}

SEQAN_DEFINE_TEST(ticket1108)
{
    static const unsigned int ARRAY_SIZE = 3u;
    typedef String<unsigned int, Array<ARRAY_SIZE> > TArray;

    {
        typedef String<TArray, MMap<> > TString;

        TString str;
        for (unsigned int i = 0u; i < 64u; ++i)
        {
//            std::cerr << i << ", " << std::endl;
            TArray elem;
            resize(elem, ARRAY_SIZE);
            appendValue(str, elem);
            SEQAN_ASSERT_LEQ(length(elem), ARRAY_SIZE);
            SEQAN_ASSERT_LEQ(length(str[0u]), ARRAY_SIZE); // Is violated as soon as i reaches 48.
            if (length(str) > 32u) SEQAN_ASSERT_LEQ(length(str[32u]), ARRAY_SIZE);
        }
        shrinkToFit(str);
    }

    {
        typedef ExternalConfig<ExternalConfig<>::TFile, 16u, ExternalConfig<>::FRAMES> TExternalConfig;
        typedef String<TArray, External<TExternalConfig> > TString;

        TString str;
        for (unsigned int i = 0u; i < 64u; ++i)
        {
//            std::cerr << i << ", " << std::endl;
            TArray elem;
            resize(elem, ARRAY_SIZE);
            appendValue(str, elem);
            SEQAN_ASSERT_LEQ(length(elem), ARRAY_SIZE);
            SEQAN_ASSERT_LEQ(length(str[0u]), ARRAY_SIZE); // Is violated as soon as i reaches 48.
            if (length(str) > 32u) SEQAN_ASSERT_LEQ(length(str[32u]), ARRAY_SIZE);
        }
    }
}
