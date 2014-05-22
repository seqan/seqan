// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2014, Knut Reinert, FU Berlin
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
// Author: Hannes Hauswedell <hauswedell@mi.fu-berlin.de>
// ==========================================================================
// Tests for reduced_aminoacid module
// ==========================================================================

#ifndef SEQAN_EXTRAS_TESTS_REDUCED_ALPHABET_H_
#define SEQAN_EXTRAS_TESTS_REDUCED_ALPHABET_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>

#include <seqan/reduced_aminoacid.h>

using namespace seqan;

SEQAN_DEFINE_TEST(test_reduced_aminoacid_cluster_red)
{
    typedef ReducedAminoAcid<ClusterReduction<8> >  ReducedAminoAcid24to8;
    typedef ReducedAminoAcid<ClusterReduction<10> > ReducedAminoAcid24to10;
    typedef ReducedAminoAcid<ClusterReduction<12> > ReducedAminoAcid24to12;

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

SEQAN_DEFINE_TEST(test_reduced_aminoacid_murphy10)
{
    typedef ReducedAminoAcid<Murphy10> ReducedAminoAcidMurphy10;

    CharString str = "AaBbCcDdEeFfGgHhIiJjKkLlMmNnOoPpQqRrSsTtUuVvWwXxYyZz*+#";
    String<AminoAcid> aas = "ARNDCQEGHILKMFPSTWYVBZX*";

    // N = 10
    {
        String<ReducedAminoAcidMurphy10> conv = str;
        SEQAN_ASSERT_EQ(
            CharString(conv),
            "AAAACCNNNNFFGGHHIIAARRIIIINNAAPPNNRRSSSSAAIIFFAAFFAAAAA");
        conv = aas;
        SEQAN_ASSERT_EQ(CharString(conv), "ARNNCNNGHIIRIFPSSFFIAAAA");
    }

}


#endif  // SEQAN_EXTRAS_TESTS_REDUCED_ALPHABET_H_
