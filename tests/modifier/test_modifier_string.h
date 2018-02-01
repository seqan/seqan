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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de
// ==========================================================================
// Various tests for modifier.
// ==========================================================================

// TODO(holtgrew): We need more systematic testing of modified string and iterator.

#ifndef SEQAN_TESTS_MODIFIER_TEST_MODIFIER_STRING_H_
#define SEQAN_TESTS_MODIFIER_TEST_MODIFIER_STRING_H_

#include <seqan/basic.h>
#include <seqan/modifier.h>


SEQAN_DEFINE_TEST(test_modifier_modified_string_metafunctions)
{
    typedef seqan::Dna5String                           TString;
    typedef seqan::ModifiedString<TString>              TInnerModifiedString;
    typedef seqan::ModifiedString<TInnerModifiedString> TOuterModifiedString;

    SEQAN_ASSERT(+(seqan::IsSameType<typename seqan::Pointer_<TString>::Type, TString *>::VALUE));
    SEQAN_ASSERT(+(seqan::IsSameType<typename seqan::Pointer_<TInnerModifiedString>::Type, TInnerModifiedString>::VALUE));
    SEQAN_ASSERT(+(seqan::IsSameType<typename seqan::Pointer_<TOuterModifiedString>::Type, TOuterModifiedString>::VALUE));
}

SEQAN_DEFINE_TEST(test_modifier_modified_string_construct)
{
    typedef seqan::Dna5String                                       TString;
    typedef seqan::ModifiedString<TString>                          TInnerModifiedString;
    typedef seqan::ModifiedString<TString const>                    TInnerConstModifiedString;
    typedef seqan::ModifiedString<TInnerModifiedString>             TOuterModifiedString;
    typedef seqan::ModifiedString<TInnerModifiedString const>       TOuterConstModifiedString;
    typedef seqan::ModifiedString<TInnerConstModifiedString>        TOuterModifiedConstString;
    typedef seqan::ModifiedString<TInnerConstModifiedString const>  TOuterConstModifiedConstString;

    // Default Construction with one level.
    {
        TInnerModifiedString modified;
    }

    // Construct with underlying string and one level.
    {
        TString original;
        TInnerModifiedString        inner(original);
        TInnerConstModifiedString   innerC(original);
        TOuterModifiedString        outer(original);
        TOuterModifiedConstString   outerC(original);
        TOuterModifiedString        outer2(inner);
        TOuterModifiedConstString   outer2C(innerC);
//        TOuterModifiedString        outer3(innerC);     // should fail
    }

    // Construct with underlying string and two levels.
    {
        // standard way
        TString original;
        TInnerModifiedString inner(original);
        TOuterModifiedString outer(inner);

        TOuterConstModifiedString outer1(inner);
        TOuterModifiedConstString outerC(original);
        TOuterConstModifiedConstString outerCC(original);
    }

    // Construct with underlying string and two levels, direct construction.
    {
        TString original;
        TOuterModifiedString outer((TInnerModifiedString(original)));
    }

    {
        typedef seqan::ModifiedString<seqan::Infix<TString>::Type> TModifiedString;

        TString original;
        TModifiedString modified(original);
//NOTE(h-2): ModStrings cannot be constructed from temporaries
//         TModifiedString modified2(suffix(original, 0));
    }

    {
        TString original;
        typedef seqan::Infix<TString const>::Type TFragment;
        TFragment frag = infix(original, 0, 0);

        typedef seqan::ModifiedString<TFragment>            TInnerModifiedString;
        typedef seqan::ModifiedString<TInnerModifiedString> TOuterModifiedString;

        TInnerModifiedString modified(frag);
        TOuterModifiedString modified2(frag);
//NOTE(h-2): ModStrings cannot be constructed from temporaries
//         TOuterModifiedString modified3(suffix(frag, 0));
    }
}

SEQAN_DEFINE_TEST(test_modifier_modified_string_assignment)
{
    typedef seqan::Dna5String                           TString;
    typedef seqan::ModifiedString<TString>              TInnerModifiedString;

    // Copy with one level.
    {
        TString original = "CGAT";
        TInnerModifiedString  inner(original);

        TInnerModifiedString inner2;
        inner2 = inner;
    }

    // Copy with two levels.
    {
        TString original = "CGAT";
        TInnerModifiedString  inner(original);

        TInnerModifiedString inner2;
        inner2 = inner;
    }
}

SEQAN_DEFINE_TEST(test_modifier_modified_string_length)
{
    typedef seqan::Dna5String                           TString;
    typedef seqan::ModifiedString<TString>              TInnerModifiedString;
    typedef seqan::ModifiedString<TInnerModifiedString> TOuterModifiedString;

    // Default Construction with one level.
    {
        TInnerModifiedString modified;
    }

    // Construct with underlying string and one level.
    {
        TString original = "CGAT";
        TInnerModifiedString inner(original);

        SEQAN_ASSERT_EQ(length(original), 4u);
        SEQAN_ASSERT_EQ(length(inner), 4u);
    }

    // Construct with underlying string and two levels.
    {
        TString original = "CGAT";
        TInnerModifiedString inner(original);
        TOuterModifiedString outer(inner);

        SEQAN_ASSERT_EQ(length(original), 4u);
        SEQAN_ASSERT_EQ(length(inner), 4u);
        SEQAN_ASSERT_EQ(length(outer), 4u);
    }

    // Construct with underlying string and two levels, direct construction.
    {
        TString original = "CGAT";
        TOuterModifiedString outer((TInnerModifiedString(original)));

        SEQAN_ASSERT_EQ(length(original), 4u);
        SEQAN_ASSERT_EQ(length(host(outer)), 4u);
        SEQAN_ASSERT_EQ(length(outer), 4u);
    }
}

// Construct modified string cascade from innermost host.
SEQAN_DEFINE_TEST(test_modifier_modified_string_cascade)
{
    // Host is a string.
    {
        typedef seqan::Dna5String                           TString;
        typedef seqan::ModifiedString<TString>              TInnerModifiedString;
        typedef seqan::ModifiedString<TInnerModifiedString> TOuterModifiedString;

        TString original = "CGAT";
        TInnerModifiedString inner(original);
        TOuterModifiedString outer(original);
    }

    // Host is a segment.
    {
        typedef seqan::Dna5String                           TString;
        typedef typename seqan::Infix<TString>::Type        TInfix;
        typedef seqan::ModifiedString<TInfix>               TInnerModifiedString;
        typedef seqan::ModifiedString<TInnerModifiedString> TOuterModifiedString;

        TString original = "CGAT";
        TInfix origInfix = infix(original, 1, 3);
        TInnerModifiedString inner(origInfix);
        TOuterModifiedString outer(origInfix);
    }
}

SEQAN_DEFINE_TEST(test_modifier_modified_iterator_construct)
{
    using namespace seqan;

    typedef Dna5String                          TString;
    typedef Iterator<TString, Standard>::Type   TIterator;
    typedef Iterator<TString, Rooted>::Type     TRootedIterator;
    typedef ModifiedIterator<TIterator>         TModifiedIterator;
    typedef ModifiedIterator<TRootedIterator>   TModifiedRootedIterator;

    // Default construction.
    {
        TModifiedIterator itM;
    }

    // Construction from host iterator.
    {
        TString seq = "CGAT";
        TIterator it = begin(seq);
        TModifiedIterator itM(it);
    }
    {
        TString seq = "CGAT";
        TModifiedIterator it(begin(seq, Standard()));
        TModifiedIterator it2(begin(seq, Rooted()));
        TModifiedRootedIterator itR(begin(seq, Rooted()));
        TModifiedIterator it3(itR);
    }
    {
        TModifiedIterator it;
        TModifiedIterator it2(it);
        ignoreUnusedVariableWarning(it2);
    }
}

struct LowerFunctor : std::unary_function<char, char>
{
    char operator()(char c) const
    {
        return tolower(c);
    }
};

struct CaesarFunctor : std::unary_function<char, char>
{
    char operator()(char c) const
    {
        int x = c;
        x += 3;
        return x % 256;
    }
};

SEQAN_DEFINE_TEST(test_modifier_modified_string_mod_view)
{
    typedef seqan::CharString TString;
    typedef seqan::ModView<LowerFunctor> TModViewLower;
    typedef seqan::ModifiedString<TString, TModViewLower> TInnerModifiedString;
    typedef seqan::ModView<CaesarFunctor> TModViewCaesar;
    typedef seqan::ModifiedString<TInnerModifiedString, TModViewCaesar> TOuterModifiedString;

    // One level only.
    {
        TString original = "CGAT";
        TInnerModifiedString modified(original);

        SEQAN_ASSERT_EQ(modified, "cgat");
    }
    // Two levels.
    {
        TString original = "CGAT";
        TInnerModifiedString modified(original);
        TOuterModifiedString outer(modified);

        SEQAN_ASSERT_EQ(outer, "fjdw");
    }
    // Two levels directly constructed with original.
    {
        TString original = "CGAT";
        TOuterModifiedString outer(original);

        SEQAN_ASSERT_EQ(outer, "fjdw");
    }
}

SEQAN_DEFINE_TEST(test_modifier_modified_string_mod_view_segment)
{
    typedef seqan::CharString TString;
    typedef seqan::Infix<TString>::Type TInfix;
    typedef seqan::ModView<LowerFunctor> TModViewLower;
    typedef seqan::ModifiedString<TInfix, TModViewLower> TInnerModifiedString;
    typedef seqan::ModView<CaesarFunctor> TModViewCaesar;
    typedef seqan::ModifiedString<TInnerModifiedString, TModViewCaesar> TOuterModifiedString;

    // One level only.
    {
        TString original = "CGAT";
        TInfix origInfix = infix(original, 1, 3);
        TInnerModifiedString modified(origInfix);

        SEQAN_ASSERT_EQ(modified, "ga");
    }
    // Two levels.
    {
        TString original = "CGAT";
        TInfix origInfix = infix(original, 1, 3);
        TInnerModifiedString modified(origInfix);
        TOuterModifiedString outer(modified);

        SEQAN_ASSERT_EQ(outer, "jd");
    }
    // Two levels directly constructed with original.
    {
        TString original = "CGAT";
        TInfix origInfix = infix(original, 1, 3);
        TOuterModifiedString outer(origInfix);

        SEQAN_ASSERT_EQ(outer, "jd");
    }
}

SEQAN_DEFINE_TEST(test_modifier_modified_string_mod_pos)
{
    typedef seqan::CharString                       TString;
    typedef seqan::Position<TString>::Type          TPos;
    typedef seqan::String<TPos>                     TPositions;
    typedef seqan::ModPos<TPositions>               TModPos;
    typedef seqan::ModifiedString<TString, TModPos> TModString;

    // One level only.
    {
        TString original = "CGAT";

        TModString modified(original);
        TPositions positions;
        setCargo(modified, positions);
        SEQAN_ASSERT(empty(modified));

        resize(cargo(modified), length(original), 0, seqan::Exact());
        SEQAN_ASSERT_EQ(modified, "CCCC");

        assign(cargo(modified), seqan::Range<TPos>(0, length(original)));
        SEQAN_ASSERT_EQ(original, modified);

        reverse(cargo(modified));
        SEQAN_ASSERT_EQ(modified, "TAGC");

        seqan::sort(modified);
        SEQAN_ASSERT_EQ(modified, "ACGT");
        SEQAN_ASSERT_EQ(original, "CGAT");

        SEQAN_ASSERT_EQ(infix(modified, 1, 3), "CG");
    }
}


SEQAN_DEFINE_TEST(test_modifier_modified_string_literal)
{
    // Host a literal.
    {
        typedef seqan::ModifiedString<char[5]> TModifiedString;

        char str[] = "CGAT";
        TModifiedString modStr(str);
        SEQAN_ASSERT_EQ(modStr, "CGAT");
    }
    // Reverse a literal.
    {
        typedef seqan::ModifiedString<char[5], seqan::ModReverse> TModifiedString;

        char str[] = "CGAT";
        TModifiedString modStr(str);
        SEQAN_ASSERT_EQ(modStr, "TAGC");
    }
    // Lowercase a literal.
    {
        typedef seqan::ModifiedString<char[5], seqan::ModView<LowerFunctor> > TModifiedString;

        char str[] = "CGAT";
        TModifiedString modStr(str);
        SEQAN_ASSERT_EQ(modStr, "cgat");
    }
}

SEQAN_DEFINE_TEST(test_modifier_modified_string_const_literal)
{
    // Host a literal.
    {
        typedef seqan::ModifiedString<const char[5]> TModifiedString;

        const char str[] = "CGAT";
        TModifiedString modStr(str);
        SEQAN_ASSERT_EQ(modStr, "CGAT");
    }
    // Reverse a literal.
    {
        typedef seqan::ModifiedString<const char[5], seqan::ModReverse> TModifiedString;

        const char str[] = "CGAT";
        TModifiedString modStr(str);
        SEQAN_ASSERT_EQ(modStr, "TAGC");
    }
    // Lowercase a literal.
    {
        typedef seqan::ModifiedString<const char[5], seqan::ModView<LowerFunctor> > TModifiedString;

        const char str[] = "CGAT";
        TModifiedString modStr(str);
        SEQAN_ASSERT_EQ(modStr, "cgat");
    }
}

SEQAN_DEFINE_TEST(test_modifier_modified_string_reverse_segment)
{
    // Inner is lower, outer is reverse.
    {
        typedef seqan::CharString TString;
        typedef seqan::Infix<TString>::Type TInfix;

        typedef seqan::ModView<LowerFunctor> TModViewLower;
        typedef seqan::ModifiedString<TInfix, TModViewLower> TInnerModifiedString;

        typedef seqan::ModifiedString<TInnerModifiedString, seqan::ModReverse> TOuterModifiedString;

        TString original = "CGAT";
        TInfix origInfix = infix(original, 1, 3);
        TOuterModifiedString outer(origInfix);
        SEQAN_ASSERT_EQ(outer, "ag");
    }
    // Inner is reverse, outer is lower.
    {
        typedef seqan::CharString TString;
        typedef seqan::Infix<TString>::Type TInfix;

        typedef seqan::ModifiedString<TInfix, seqan::ModReverse> TInnerModifiedString;

        typedef seqan::ModView<LowerFunctor> TModViewLower;
        typedef seqan::ModifiedString<TInnerModifiedString, TModViewLower> TOuterModifiedString;

        TString original = "CGAT";
        TInfix origInfix = infix(original, 1, 3);
        TOuterModifiedString outer(origInfix);
        SEQAN_ASSERT_EQ(outer, "ag");
    }
}

SEQAN_DEFINE_TEST(test_modifier_minimal)
{
    typedef seqan::CharString TString;
    typedef seqan::ModifiedString<TString, seqan::ModReverse> TInnerModifiedString;

    TString original = "The QUICK brown fox.";
    TInnerModifiedString inner(original);

    // static_cast<typename seqan::Parameter_<TString>::Type *>(seqan::Nothing());  // #=> non-const reference
    // static_cast<typename TInnerModifiedString::THostPointer_ *>(seqan::Nothing());  // #=> const pointer
    // static_cast<typename seqan::Pointer_<TString>::Type *>(seqan::Nothing());
}

SEQAN_DEFINE_TEST(test_modifier_reverse_back_front)
{
    using namespace seqan;

    {
        typedef CharString TString;
        typedef ModifiedString<TString, ModReverse> TInnerModifiedString;
        typedef typename Iterator<TInnerModifiedString, Standard>::Type TIt;

        TString original = "The QUICK brown fox.";
        TInnerModifiedString inner(original);

        TIt it = begin(inner, Standard());
        SEQAN_ASSERT_EQ(*it, '.');
        SEQAN_ASSERT_EQ((*begin(inner, Standard())), '.');
        it = end(inner, Standard());
        --it;
        SEQAN_ASSERT_EQ(*it, 'T');

        SEQAN_ASSERT_EQ(front(inner), '.');
        SEQAN_ASSERT_EQ(back(inner), 'T');
    }
    {
        typedef CharString const TString;
        typedef ModifiedString<TString, ModReverse> TInnerModifiedString;
        typedef typename Iterator<TInnerModifiedString, Standard>::Type TIt;

        TString original = "The QUICK brown fox.";
        TInnerModifiedString inner(original);

        TIt it = begin(inner, Standard());
        SEQAN_ASSERT_EQ(*it, '.');
        SEQAN_ASSERT_EQ((*begin(inner, Standard())), '.');
        it = end(inner, Standard());
        --it;
        SEQAN_ASSERT_EQ(*it, 'T');

        SEQAN_ASSERT_EQ(front(inner), '.');
        SEQAN_ASSERT_EQ(back(inner), 'T');
    }
}

SEQAN_DEFINE_TEST(test_modifier_reverse_iterator_metafunctions)
{
    using namespace seqan;

    typedef ModifiedIterator<CharString, ModReverse> TModifiedIterator;

    {
        typedef char TExpected;
        typedef Value<TModifiedIterator>::Type TResult;
        bool res = IsSameType<TExpected, TResult>::VALUE;
        SEQAN_ASSERT(res);
        static_assert(std::is_same<TExpected, TResult>::value, "Different type expected.");
    }
    {
        typedef char const & TExpected;
        typedef GetValue<TModifiedIterator>::Type TResult;
        bool res = IsSameType<TExpected, TResult>::VALUE;
        SEQAN_ASSERT(res);
        static_assert(std::is_same<TExpected, TResult>::value, "Different type expected.");
    }
    {
        typedef char & TExpected;
        typedef Reference<TModifiedIterator>::Type TResult;
        bool res = IsSameType<TExpected, TResult>::VALUE;
        static_assert(std::is_same<TExpected, TResult>::value, "Different type expected.");
        SEQAN_ASSERT(res);
    }
}


#endif // #ifndef SEQAN_TESTS_MODIFIER_TEST_MODIFIER_STRING_H_
