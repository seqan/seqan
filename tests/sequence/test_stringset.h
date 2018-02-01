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
// Test sequence default interface:  Non-container objects are treated like
// containers of length 1.
// ==========================================================================

template <typename TSpec>
void Test_StringSet()
{
    typedef StringSet<CharString, TSpec> TStringSet;
    TStringSet set;

    resize(set, 3);
    set[0] = "Hallo ";
    set[1] = "schlauer ";
    set[2] = "Hamster!";

    SEQAN_ASSERT_EQ(length(set), 3u);

    SEQAN_ASSERT(isEqual(set[0], "Hallo "));
    SEQAN_ASSERT(isEqual(set[1], "schlauer "));
    SEQAN_ASSERT(isEqual(set[2], "Hamster!"));
/*
    // currently, this won't work for Owner<ConcatDirect<..> > StringSets
    // to fix it, we need to introduce Modifiers for Segments
    // which propagate their resize events to their StringSets
    resize(set[0], 9);
    infix(set[0], 6, 9) = "du ";
    SEQAN_ASSERT(isEqual(set[0], "Hallo du "));
*/

    //StringSet iterators
    typedef typename Iterator<TStringSet>::Type TIterator;
    int i = 0;
    for (TIterator it = begin(set); it != end(set); goNext(it))
    {
        SEQAN_ASSERT_EQ(*it, set[i]);
        ++i;
    }

    TIterator itBegin = begin(set);
    SEQAN_ASSERT(atBegin(itBegin));
    SEQAN_ASSERT_NOT(atEnd(itBegin));
    TIterator itEnd = end(set);
    SEQAN_ASSERT_NOT(atBegin(itEnd));
    SEQAN_ASSERT(atEnd(itEnd));
    SEQAN_ASSERT_EQ(i, 3);

    // Test erase().
    {
        TStringSet set;
        appendValue(set, "one");
        appendValue(set, "two");
        appendValue(set, "three");
        appendValue(set, "four");
        appendValue(set, "five");
        appendValue(set, "six");

        erase(set, 0);
        SEQAN_ASSERT_EQ(length(set), 5u);
        SEQAN_ASSERT_EQ(set[0], CharString("two"));

        erase(set, 2, 4);
        SEQAN_ASSERT_EQ(length(set), 3u);
        SEQAN_ASSERT_EQ(set[0], CharString("two"));
        SEQAN_ASSERT_EQ(set[1], CharString("three"));
        SEQAN_ASSERT_EQ(set[2], CharString("six"));
    }
}

//____________________________________________________________________________


template <typename TSpec>
void Test_StringSet_Concat()
{
    StringSet<CharString, TSpec> set;

    CharString s1 = "Hallo ";

    appendValue(set, s1);
    appendValue(set, "schlauer ");
    appendValue(set, "");
    appendValue(set, "Hamster!");

    SEQAN_ASSERT_EQ(length(set), 4u);

    CharString all = concat(set);

    SEQAN_ASSERT_EQ(concat(set)[10], 'a');
    SEQAN_ASSERT(isEqual(set[0], "Hallo "));
    SEQAN_ASSERT(isEqual(set[1], "schlauer "));
    SEQAN_ASSERT(isEqual(set[2], ""));
    SEQAN_ASSERT(isEqual(set[3], "Hamster!"));
    SEQAN_ASSERT(isEqual(all, "Hallo schlauer Hamster!"));

    SEQAN_ASSERT_EQ(stringSetLimits(set)[0], 0u);
    SEQAN_ASSERT_EQ(stringSetLimits(set)[1], 6u);
    SEQAN_ASSERT_EQ(stringSetLimits(set)[2], 15u);
    SEQAN_ASSERT_EQ(stringSetLimits(set)[3], 15u);
    SEQAN_ASSERT_EQ(stringSetLimits(set)[4], 23u);

    StringSet<CharString, TSpec> const &cset = set;

    all = concat(cset);
    SEQAN_ASSERT_EQ(concat(cset)[10], 'a');
    SEQAN_ASSERT(isEqual(all, "Hallo schlauer Hamster!"));
}


//////////////////////////////////////////////////////////////////////////////

template <typename TStringSet>
void Test_StringSetIdHolder() {
    typedef    typename Id<TStringSet>::Type TId;
    //typedef StringSet<String<char>, Dependent<Tight> > TSetTight;
    //typedef StringSet<String<char>, Dependent<Generous> > TSetGenerous;

    TStringSet str;
    String<char> bla("a");
    TId id0 = assignValueById(str, bla);
    SEQAN_ASSERT_EQ(id0, 0u);
    SEQAN_ASSERT_EQ(idToPosition(str, id0), 0u);
    SEQAN_ASSERT_EQ(positionToId(str, 0), id0);
    SEQAN_ASSERT_EQ(length(str), 1u);
    SEQAN_ASSERT_EQ(str[0], "a");
    SEQAN_ASSERT_EQ(getValueById(str, id0), "a");
    String<char> bla1("b");
    TId id1 = assignValueById(str, bla1);
    SEQAN_ASSERT_EQ(id1, 1u);
    SEQAN_ASSERT_EQ(idToPosition(str, id1), 1u);
    SEQAN_ASSERT_EQ(positionToId(str, 1), id1);
    SEQAN_ASSERT_EQ(str[1], "b");
    SEQAN_ASSERT_EQ(length(str), 2u);
    SEQAN_ASSERT_EQ(getValueById(str, id1), "b");
    String<char> bla2("c");
    TId id2 = assignValueById(str, bla2);
    SEQAN_ASSERT_EQ(id2, 2U);
    SEQAN_ASSERT_EQ(str[2], "c");
    SEQAN_ASSERT_EQ(length(str), 3u);
    SEQAN_ASSERT_EQ(getValueById(str, id2), "c");
    String<char> bla3("d");
    TId id3 = assignValueById(str, bla3);
    SEQAN_ASSERT_EQ(id3, 3u);
    SEQAN_ASSERT_EQ(str[3], "d");
    SEQAN_ASSERT_EQ(length(str), 4u);
    SEQAN_ASSERT_EQ(getValueById(str, id3), "d");
    removeValueById(str,id1);
    SEQAN_ASSERT_EQ(getValueById(str, id0), "a");
    SEQAN_ASSERT_EQ(getValueById(str, id2), "c");
    SEQAN_ASSERT_EQ(getValueById(str, id3), "d");
    SEQAN_ASSERT_EQ(length(str), 3u);
    removeValueById(str,id2);
    SEQAN_ASSERT_EQ(getValueById(str, id0), "a");
    SEQAN_ASSERT_EQ(getValueById(str, id3), "d");
    SEQAN_ASSERT_EQ(length(str), 2u);

    String<char> bla4("e");
    TId id4 = assignValueById(str, bla4, 100);
    SEQAN_ASSERT_EQ(id4, 100u);
    SEQAN_ASSERT_EQ(getValueById(str, id4), "e");
    SEQAN_ASSERT_EQ(getValueById(str, id0), "a");
    SEQAN_ASSERT_EQ(getValueById(str, id3), "d");
    removeValueById(str,id3);
    SEQAN_ASSERT_EQ(getValueById(str, id0), "a");
    SEQAN_ASSERT_EQ(getValueById(str, id4), "e");
    SEQAN_ASSERT_EQ(length(str), 2u);
    String<char> bla5("f");
    TId id5 = assignValueById(str, bla5);
    SEQAN_ASSERT_EQ(getValueById(str, id0), "a");
    SEQAN_ASSERT_EQ(getValueById(str, id4), "e");
    SEQAN_ASSERT_EQ(getValueById(str, id5), "f");
    assignValueById(str, bla5, id4);
    SEQAN_ASSERT_EQ(getValueById(str, id0), "a");
    SEQAN_ASSERT_EQ(getValueById(str, id4), "f");
    SEQAN_ASSERT_EQ(getValueById(str, id5), "f");
    removeValueById(str,id4);
    SEQAN_ASSERT_EQ(getValueById(str, id0), "a");
    SEQAN_ASSERT_EQ(getValueById(str, id5), "f");
    SEQAN_ASSERT_EQ(length(str), 2u);
    clear(str);
    id1 = assignValueById(str, bla1);
    id2 = assignValueById(str, bla2);
    id3 = assignValueById(str, bla3);
    SEQAN_ASSERT_EQ(getValueById(str, id1), "b");
    SEQAN_ASSERT_EQ(getValueById(str, id2), "c");
    SEQAN_ASSERT_EQ(getValueById(str, id3), "d");
    SEQAN_ASSERT_EQ(length(str), 3u);
}

//____________________________________________________________________________


template <typename TSpec>
void Test_StringSet_Id()
{
    StringSet<CharString, Owner<Default> > origin;
    StringSet<CharString, TSpec> set;

    resize(origin, 3);
    origin[0] = "Hallo ";
    origin[1] = "schlauer ";
    origin[2] = "Hamster!";

    appendValue(set, origin[0]);
    appendValue(set, origin[1]);
    appendValue(set, origin[2]);

    SEQAN_ASSERT_EQ(length(set), 3u);

    CharString all = concat(set);

    SEQAN_ASSERT_EQ(concat(set)[10], 'a');
    SEQAN_ASSERT(isEqual(set[0], "Hallo "));
    SEQAN_ASSERT(isEqual(set[1], "schlauer "));
    SEQAN_ASSERT(isEqual(set[2], "Hamster!"));
    SEQAN_ASSERT(isEqual(all, "Hallo schlauer Hamster!"));

    SEQAN_ASSERT_EQ(stringSetLimits(set)[0], 0u);
    SEQAN_ASSERT_EQ(stringSetLimits(set)[1], 6u);
    SEQAN_ASSERT_EQ(stringSetLimits(set)[2], 15u);
    SEQAN_ASSERT_EQ(stringSetLimits(set)[3], 23u);

    StringSet<CharString, TSpec> const &cset = set;

    all = concat(cset);
    SEQAN_ASSERT_EQ(concat(cset)[10], 'a');
    SEQAN_ASSERT(isEqual(all, "Hallo schlauer Hamster!"));
}

SEQAN_DEFINE_TEST(StringSet_Owner_Default) {
    Test_StringSet< Owner<Default> >();
}

SEQAN_DEFINE_TEST(StringSet_Concat_Owner_Default) {
    Test_StringSet_Concat< Owner<Default> >();
}

SEQAN_DEFINE_TEST(StringSet_Concat_Owner_ConcatDirect) {
    Test_StringSet_Concat< Owner<ConcatDirect<> > >();
}

SEQAN_DEFINE_TEST(StringSet_Id_Dependent_Tight) {
    Test_StringSet_Id< Dependent<Tight> >();
}

SEQAN_DEFINE_TEST(StringSet_Id_Dependent_Generous) {
    Test_StringSet_Id< Dependent<Generous> >();
}

SEQAN_DEFINE_TEST(StringSetIdHolder_Char_Dependent_Tight) {
    Test_StringSetIdHolder<StringSet<String<char>, Dependent<Tight> > >();
}

SEQAN_DEFINE_TEST(StringSetIdHolder_Char_Dependent_Generous) {
    Test_StringSetIdHolder<StringSet<String<char>, Dependent<Generous> > >();
}

struct TestContainer
{
    String<int> string;

    TestContainer()
    {
        resize(string, 10, 0);
        std::cerr << "DEFAULT CONSTRUCTING TestContainer " << (void*)(this) << std::endl;
        std::cerr << "  string.data_begin " << (void*)(string.data_begin) << std::endl;
    }

    TestContainer(TestContainer const & other)
        : string(other.string)
    {
        std::cerr << "COPY CONSTRUCTING TestContainer " << (void*)(this) << std::endl;
        std::cerr << "  other " << (void*)(&other) << std::endl;
        std::cerr << "  string.data_begin " << (void*)(string.data_begin) << std::endl;
    }

    ~TestContainer()
    {
        std::cerr << "DECONSTRUCTING TestContainer " << (void*)(this) << std::endl;
        std::cerr << "  string.data_begin " << (void*)(string.data_begin) << std::endl;
    }
};

SEQAN_DEFINE_TEST(test_find_motif_memory_leak_ticket_364)
{
    // This shows the bug described in #364, also occurs with String
    // instead of StringSet.
//    {
//        StringSet<String<TestContainer> > x;
//        TestContainer c;
//        appendValue(x, c);
//    }
    {
        String<TestContainer> x;
        TestContainer c;
        appendValue(x, c);
    }
}
