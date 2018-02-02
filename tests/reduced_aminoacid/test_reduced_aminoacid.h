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
// Author: Hannes Hauswedell <hannes.hauswedell@fu-berlin.de>
// ==========================================================================
// Tests for reduced_aminoacid module
// ==========================================================================

#ifndef SEQAN_TESTS_REDUCED_ALPHABET_H_
#define SEQAN_TESTS_REDUCED_ALPHABET_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>

#include <seqan/reduced_aminoacid.h>
#include <seqan/modifier.h>
#include <seqan/index.h>

using namespace seqan;

#if 0
SEQAN_DEFINE_TEST(test_reduced_aminoacid_cluster_red)
{
    typedef SimpleType<unsigned char, ReducedAminoAcid_<ClusterReduction<8> > >
            ReducedAminoAcid24to8;
    typedef SimpleType<unsigned char, ReducedAminoAcid_<ClusterReduction<10> > >
            ReducedAminoAcid24to10;
    typedef SimpleType<unsigned char, ReducedAminoAcid_<ClusterReduction<12> > >
            ReducedAminoAcid24to12;

    CharString str = "AaBbCcDdEeFfGgHhIiJjKkLlMmNnOoPpQqRrSsTtUuVvWwXxYyZz*+#";
    String<AminoAcid> aas = "ARNDCQEGHILKMFPSTWYVBZX*";

    // N = 12
    {
        String<ReducedAminoAcid24to12> conv = str;
        SEQAN_ASSERT_EQ(
            CharString(conv),
            "AANNCCNNRRFFGGHHIISSRRIIIINNSSPPRRRRSSSSSSIIWWSSFFRR*SS");
        conv = aas;
        SEQAN_ASSERT_EQ(CharString(conv), "ARNNCRRGHIIRIFPSSWFINRS*");
    }

    // N = 10
    {
        String<ReducedAminoAcid24to10> conv = str;
        SEQAN_ASSERT_EQ(
            CharString(conv),
            "AANNCCNNRRFFGGHHIIAARRIIIINNAAPPRRRRAAAAAAIIFFAAFFRR*AA");
        conv = aas;
        SEQAN_ASSERT_EQ(CharString(conv), "ARNNCRRGHIIRIFPAAFFINRA*");
    }

    // N = 8
    {
        String<ReducedAminoAcid24to8> conv = str;
        SEQAN_ASSERT_EQ(
            CharString(conv),
            "AARRCCRRRRFFGGRRIIAARRIIIIRRAAPPRRRRAAAAAAIIFFAAFFRR*AA");
        conv = aas;
        SEQAN_ASSERT_EQ(CharString(conv), "ARRRCRRGRIIRIFPAAFFIRRA*");
    }
}
#endif

SEQAN_DEFINE_TEST(test_reduced_aminoacid_buchfink11)
{
    typedef SimpleType<unsigned char, ReducedAminoAcid_<Buchfink11> >
            ReducedAminoAcidBuchfink11;

    CharString str = "AaBbCcDdEeFfGgHhIiJjKkLlMmNnOoPpQqRrSsTtUuVvWwXxYyZz*+#";
    String<AminoAcid> aas = "ABCDEFGHIJKLMNOPQRSTUVWYZX*";

    // N = 11
    {
        String<ReducedAminoAcidBuchfink11> conv = str;
        SEQAN_ASSERT_EQ(
            CharString(conv),
            "AABBCCBBBBFFGGHHIIIIBBIIMMBBBBPPBBBBAAAACCIIWWAAYYBBFAA");
        conv = aas;
        SEQAN_ASSERT_EQ(CharString(conv), "ABCBBFGHIIBIMBBPBBAACIWYBAF");
    }
}

SEQAN_DEFINE_TEST(test_reduced_aminoacid_cannata10)
{
    typedef SimpleType<unsigned char, ReducedAminoAcid_<Cannata10> >
            ReducedAminoAcidCannata10;

    CharString str = "AaBbCcDdEeFfGgHhIiJjKkLlMmNnOoPpQqRrSsTtUuVvWwXxYyZz*+#";
    String<AminoAcid> aas = "ABCDEFGHIJKLMNOPQRSTUVWYZX*";

    // N = 10
    {
        String<ReducedAminoAcidCannata10> conv = str;
        SEQAN_ASSERT_EQ(
            CharString(conv),
            "AABBCCBBEEFFAAHHIIIIKKIIIIBBKKPPEEKKAAAACCIIWWAAFFEEFAA");
        conv = aas;
        SEQAN_ASSERT_EQ(CharString(conv), "ABCBEFAHIIKIIBKPEKAACIWFEAF");
    }
}

SEQAN_DEFINE_TEST(test_reduced_aminoacid_li10)
{
    typedef SimpleType<unsigned char, ReducedAminoAcid_<Li10> >
            ReducedAminoAcidLi10;

    CharString str = "AaBbCcDdEeFfGgHhIiJjKkLlMmNnOoPpQqRrSsTtUuVvWwXxYyZz*+#";
    String<AminoAcid> aas = "ABCDEFGHIJKLMNOPQRSTUVWYZX*";

    // N = 10
    {
        String<ReducedAminoAcidLi10> conv = str;
        SEQAN_ASSERT_EQ(
            CharString(conv),
            "AABBCCBBBBFFGGHHIIJJKKJJJJHHKKPPBBKKAAAACCIIFFAAFFBBFAA");
        conv = aas;
        SEQAN_ASSERT_EQ(CharString(conv), "ABCBBFGHIJKJJHKPBKAACIFFBAF");
    }
}

SEQAN_DEFINE_TEST(test_reduced_aminoacid_solis10)
{
    typedef SimpleType<unsigned char, ReducedAminoAcid_<Solis10> >
            ReducedAminoAcidSolis10;

    CharString str = "AaBbCcDdEeFfGgHhIiJjKkLlMmNnOoPpQqRrSsTtUuVvWwXxYyZz*+#";
    String<AminoAcid> aas = "ABCDEFGHIJKLMNOPQRSTUVWYZX*";

    // N = 10
    {
        String<ReducedAminoAcidSolis10> conv = str;
        SEQAN_ASSERT_EQ(
            CharString(conv),
            "AABBCCBBBBFFGGHHIIIIKKIIIIGGHHPPGGHHGGPPCCIIWWAAWWBBFAA");
        conv = aas;
        SEQAN_ASSERT_EQ(CharString(conv), "ABCBBFGHIIKIIGHPGHGPCIWWBAF");
    }
}

SEQAN_DEFINE_TEST(test_reduced_aminoacid_murphy5)
{
    typedef SimpleType<unsigned char, ReducedAminoAcid_<Murphy5> >
            ReducedAminoAcidMurphy5;

    CharString str = "AaBbCcDdEeFfGgHhIiJjKkLlMmNnOoPpQqRrSsTtUuVvWwXxYyZz*+#";
    String<AminoAcid> aas = "ABCDEFGHIJKLMNOPQRSTUVWYZX*";

    // N = 5
    {
        String<ReducedAminoAcidMurphy5> conv = str;
        SEQAN_ASSERT_EQ(
            CharString(conv),
            "AABBCCBBBBFFAAHHCCCCHHCCCCBBHHAABBHHAAAACCCCFFAAFFBBFAA");
        conv = aas;
        SEQAN_ASSERT_EQ(CharString(conv), "ABCBBFAHCCHCCBHABHAACCFFBAF");
    }
}

SEQAN_DEFINE_TEST(test_reduced_aminoacid_murphy10)
{
    typedef SimpleType<unsigned char, ReducedAminoAcid_<Murphy10> >
            ReducedAminoAcidMurphy10;

    CharString str = "AaBbCcDdEeFfGgHhIiJjKkLlMmNnOoPpQqRrSsTtUuVvWwXxYyZz*+#";
    String<AminoAcid> aas = "ABCDEFGHIJKLMNOPQRSTUVWYZX*";

    // N = 10
    {
        String<ReducedAminoAcidMurphy10> conv = str;
        SEQAN_ASSERT_EQ(
            CharString(conv),
            "AABBCCBBBBFFGGHHIIIIKKIIIIBBAAPPBBKKSSSSAAIIFFAAFFBBFAA");
        conv = aas;
        SEQAN_ASSERT_EQ(CharString(conv), "ABCBBFGHIIKIIBAPBKSSAIFFBAF");
    }
}

template <typename TModString>
void _testReducedAminoAcidMurphy10ModIteratorsImpl(TModString & conv)
{
    typedef typename Iterator<TModString, Standard>::Type TIt;
    typedef typename Iterator<TModString, Rooted>::Type TItR;

    CharString toCharString = conv;
    SEQAN_ASSERT_EQ(toCharString,
                    "ABCBBFGHIIKIIBAPBKSSAIFFBAF");

    // iterating
    {
        unsigned c = 0;
        for (TIt it = begin(conv, Standard()), itEnd = end(conv, Standard());
             it != itEnd;
             ++it, ++c)
            SEQAN_ASSERT_EQ(char(*it), toCharString[c]);
    }

    // atBegin, atEnd, position (Standard)
    {
        TIt it = begin(conv, Standard());
        SEQAN_ASSERT(atBegin(it, conv));
        SEQAN_ASSERT_EQ(position(it, conv), 0u);

        it = end(conv, Standard());
        SEQAN_ASSERT(atEnd(it, conv));
        SEQAN_ASSERT_EQ(position(it, conv), length(conv));
    }

    // atBegin, atEnd, position (Rooted)
    {
        TItR it = begin(conv, Rooted());
        SEQAN_ASSERT(atBegin(it));
        SEQAN_ASSERT(atBegin(it, conv));
        SEQAN_ASSERT_EQ(position(it), 0u);
        SEQAN_ASSERT_EQ(position(it, conv), 0u);

        it = end(conv, Rooted());
        SEQAN_ASSERT(atEnd(it));
        SEQAN_ASSERT(atEnd(it, conv));
        SEQAN_ASSERT_EQ(position(it), length(conv));
        SEQAN_ASSERT_EQ(position(it, conv), length(conv));
    }
}

SEQAN_DEFINE_TEST(test_reduced_aminoacid_murphy10_moditerators)
{
    typedef SimpleType<unsigned char, ReducedAminoAcid_<Murphy10> >
            ReducedAminoAcidMurphy10;
    typedef ModifiedString<String<AminoAcid>,
                           ModView<FunctorConvert<AminoAcid, ReducedAminoAcidMurphy10>>> TModString;
    String<AminoAcid> aas = "ABCDEFGHIJKLMNOPQRSTUVWYZX*";

    TModString conv(aas);
    _testReducedAminoAcidMurphy10ModIteratorsImpl(conv);

    TModString const conv2(aas);
    _testReducedAminoAcidMurphy10ModIteratorsImpl(conv2);

    Segment<TModString, InfixSegment> convinf = infix(conv, 0, length(conv));
    _testReducedAminoAcidMurphy10ModIteratorsImpl(convinf);

    Segment<TModString const, InfixSegment> conv2inf = infix(conv2, 0, length(conv));
    _testReducedAminoAcidMurphy10ModIteratorsImpl(conv2inf);
}

struct ReducedFMIndexConfig_
{
    typedef size_t                                                          LengthSum;
    typedef WaveletTree<void, WTRDConfig<LengthSum, Alloc<>, 1> >           Bwt;
    typedef Levels<void, LevelsRDConfig<LengthSum, Alloc<>, 1> >            Sentinels;

    static const unsigned SAMPLING = 10;
};

SEQAN_DEFINE_TEST(test_reduced_aminoacid_murphy10_modview_fmindex)
{
    typedef String<AminoAcid>                                               TOrigString;
    typedef StringSet<TOrigString, Owner<ConcatDirect<> > >                 TOrigSet;

    typedef SimpleType<unsigned char, ReducedAminoAcid_<Murphy10> >         ReducedAminoAcidMurphy10;
    typedef ModView<FunctorConvert<AminoAcid, ReducedAminoAcidMurphy10> >   TModView;
    typedef ModifiedString<TOrigString, TModView>                           TModString;
    typedef StringSet<TModString, Owner<ConcatDirect<> > >                  TModSet;

    typedef FMIndex<void, ReducedFMIndexConfig_>                            TFMIndex;

    TOrigSet origSet;
    appendValue(origSet, "ABCDEFGHIJKLMNOPQRSTUVWYZX*");
    appendValue(origSet, "ABABABABABABILMVILMVILMVABABABABAB");
    appendValue(origSet, "ABCDEFGHIJKLMNOPQRSTUVWYZX*LLLLL");
    reverse(origSet); // FM-Index is reversed o_O

    TModSet modSet(origSet);
    SEQAN_ASSERT_EQ(modSet[0], "FABFFIASSKBPABIIKIIHGFBBCBA");
    SEQAN_ASSERT_EQ(modSet[1], "BABABABABAIIIIIIIIIIIIBABABABABABA");
    SEQAN_ASSERT_EQ(modSet[2], "IIIIIFABFFIASSKBPABIIKIIHGFBBCBA");

    TOrigString query = "VVVVV";
    TModString modQuery(query);
    SEQAN_ASSERT_EQ(modQuery, "IIIII");

    Index<TModSet, TFMIndex> index(modSet);
    indexRequire(index, FibreSALF());         // instantiate

    // actual search is only done if lambdas are available
    typedef typename Iterator<Index<TModSet, TFMIndex>, TopDown<>>::Type TIndexIt;

    std::vector<std::pair<uint64_t, uint64_t>> hits;
    auto callback = [&] (TIndexIt & indexIt, int)
    {
        auto const & occurrences = getOccurrences(indexIt);
        for (auto it = begin(occurrences), itEnd = end(occurrences); it != itEnd; ++it)
        {
            auto subjOcc = *it;
            // reverse positions again
            setSeqOffset(subjOcc,
                         length(origSet[getSeqNo(subjOcc)])
                         - getSeqOffset(subjOcc)
                         - length(query));
            hits.emplace_back(getSeqNo(subjOcc), getSeqOffset(subjOcc));
        }
    };

    Nothing nothing;
    _findImpl(nothing, index, modQuery, int(0), callback, Backtracking<Exact>());

    SEQAN_ASSERT_EQ(length(hits), 9u);
    SEQAN_ASSERT_EQ(std::get<0>(hits[0]), 1u); SEQAN_ASSERT_EQ(std::get<1>(hits[0]), 12u);
    SEQAN_ASSERT_EQ(std::get<0>(hits[1]), 2u); SEQAN_ASSERT_EQ(std::get<1>(hits[1]), 27u);
    for (unsigned i = 1; i < 7; ++i)
    {
        SEQAN_ASSERT_EQ(std::get<0>(hits[1+i]), 1u); SEQAN_ASSERT_EQ(std::get<1>(hits[1+i]), 12u +i);
    }
}

#endif  // SEQAN_TESTS_REDUCED_ALPHABET_H_
