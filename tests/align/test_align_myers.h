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

#ifndef TESTS_ALIGN_TEST_ALIGN_MYERS_H_
#define TESTS_ALIGN_TEST_ALIGN_MYERS_H_

//#define SEQAN_TEST

#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/align.h>

using namespace std;
using namespace seqan;

template <typename TAlphabet>
String<TAlphabet> generate_random(int length_of_sequence)
{
    // init string
    String<TAlphabet> ret;
    resize(ret,length_of_sequence);
    int alphabet_size = ValueSize<TAlphabet>::VALUE;
    // generate random sequence of length "length_of_sequence"
    for (int i = 0; i < length_of_sequence; ++i)
    ret[i] = static_cast<TAlphabet>((rand() >> 4) % alphabet_size);

    return ret;
}

template <typename TAlphabet>
String<TAlphabet> generate_second_sequence(int error_count,String<TAlphabet> copy_of_org)
{
    int length_of_org = length(copy_of_org);

    int alphabet_size = ValueSize<TAlphabet>::VALUE;

    for(int i = 0;i < error_count;++i)
    {
        // introduce errors into sequence
        // 1. choose position
        int pos = static_cast<int>((int)length_of_org * rand() / (RAND_MAX +  1.0));

        // generate new char
        TAlphabet new_char = static_cast<TAlphabet>((rand() >> 4) % alphabet_size);

        // replace char
        copy_of_org[pos] = new_char;
    }
    return copy_of_org;
}


template <typename TAlphabet>
void erase_sequence_parts(int erase_count,String<TAlphabet> & sequence)
{
    // erase single characters
    int len = length(sequence) - 1;
    for(int i = 0;i < erase_count;++i)
    {
        // calc position
        int pos = static_cast<int>((int)(len - i) * rand() / (RAND_MAX +  1.0));
        erase(sequence,pos,pos+1);
    }
}

#define ALPHABET Dna

int edit_distance(Align<String<ALPHABET>, ArrayGaps> & ali)
{
    int len_ali = length(row(ali,0));

    Iterator<Row<Align<String<ALPHABET>, ArrayGaps> >::Type >::Type ali_row_0 = iter(row(ali, 0), 0);
    Iterator<Row<Align<String<ALPHABET>, ArrayGaps> >::Type >::Type ali_row_1 = iter(row(ali, 1), 0);

    int score = 0;
    int i;
    // iteration ueber das alignment
    for(i = 0;i < len_ali;++i)
    {
        if(isGap(ali_row_0))
        {
            --score;
        }
        else if(isGap(ali_row_1))
        {
            --score;
        }
        else if(value(ali_row_0) != value(ali_row_1))
        {
            --score;
        }

        goNext(ali_row_0);
        goNext(ali_row_1);
    }

    return score;
}


SEQAN_DEFINE_TEST(test_align_myers_test_short) {
    int nw_score,m_score,hm_score;
    int test_repeat = 1;
    int test_count = 0;

    while(test_count < test_repeat)    {
        // create random sequences
        String<ALPHABET> s_str0 = generate_random<ALPHABET>(20);
        String<ALPHABET> s_str1 = generate_second_sequence<ALPHABET>(3,s_str0);
        erase_sequence_parts(3,s_str1);

        // test alignment with random sequences
        // use needleman wunsch as reference
        Align<String<ALPHABET>, ArrayGaps> s_nw_ali;
        resize(rows(s_nw_ali), 2);
        assignSource(row(s_nw_ali, 0), s_str0);
        assignSource(row(s_nw_ali, 1), s_str1);

        nw_score = globalAlignment(s_nw_ali,SimpleScore());

        // compute only score
        Align<String<ALPHABET>, ArrayGaps> s_m_ali;
        resize(rows(s_m_ali), 2);
        assignSource(row(s_m_ali, 0), s_str0);
        assignSource(row(s_m_ali, 1), s_str1);

        m_score = globalAlignment(s_m_ali,SimpleScore(), MyersBitVector());
        SEQAN_ASSERT_EQ(nw_score, m_score);

        // compute complete alignments
        Align<String<ALPHABET>, ArrayGaps> s_hm_ali;
        resize(rows(s_hm_ali), 2);
        assignSource(row(s_hm_ali, 0), s_str0);
        assignSource(row(s_hm_ali, 1), s_str1);

        hm_score = globalAlignment(s_hm_ali,SimpleScore(), MyersHirschberg());
        SEQAN_ASSERT_EQ(nw_score, hm_score);
        SEQAN_ASSERT_EQ(edit_distance(s_hm_ali), hm_score);

        ++test_count;
    }
}


SEQAN_DEFINE_TEST(test_align_myers_test_long) {
    int nw_score,m_score,hm_score;
    int test_repeat = 1;
    int test_count = 0;

    while(test_count < test_repeat) {
        // create random sequences
        String<ALPHABET> l_str0 = generate_random<ALPHABET>(200);
        String<ALPHABET> l_str1 = generate_second_sequence<ALPHABET>(10,l_str0);
        erase_sequence_parts(5,l_str1);

        // test alignment with random sequences
        // use needleman wunsch as reference
        Align<String<ALPHABET>, ArrayGaps> l_nw_ali;
        resize(rows(l_nw_ali), 2);
        assignSource(row(l_nw_ali, 0), l_str0);
        assignSource(row(l_nw_ali, 1), l_str1);

        nw_score = globalAlignment(l_nw_ali,SimpleScore());

        // compute only score
        Align<String<ALPHABET>, ArrayGaps> l_m_ali;
        resize(rows(l_m_ali), 2);
        assignSource(row(l_m_ali, 0), l_str0);
        assignSource(row(l_m_ali, 1), l_str1);

        m_score = globalAlignment(l_m_ali,SimpleScore(), MyersBitVector());

        SEQAN_ASSERT_EQ(nw_score, m_score);

        // compute complete alignments
        Align<String<ALPHABET>, ArrayGaps> l_hm_ali;
        resize(rows(l_hm_ali), 2);
        assignSource(row(l_hm_ali, 0), l_str0);
        assignSource(row(l_hm_ali, 1), l_str1);

        hm_score = globalAlignment(l_hm_ali, SimpleScore(), MyersHirschberg());

        SEQAN_ASSERT_EQ(nw_score, hm_score);
        SEQAN_ASSERT_EQ(edit_distance(l_hm_ali), hm_score);

        ++test_count;
    }
}


SEQAN_DEFINE_TEST(test_align_hirschberger) {
    int nw_score, hm_score;
    int test_repeat = 1;
    int test_count = 0;

    while(test_count < test_repeat) {
        // create random sequences
        String<ALPHABET> str0 = generate_random<ALPHABET>(20);
        String<ALPHABET> str1 = generate_second_sequence<ALPHABET>(2,str0);

        erase_sequence_parts(5,str1);

        // test alignment with random sequences
        // use needleman wunsch as reference
        Align<String<ALPHABET>, ArrayGaps> nw_ali;
        resize(rows(nw_ali), 2);
        assignSource(row(nw_ali, 0), str0);
        assignSource(row(nw_ali, 1), str1);

        nw_score = globalAlignment(nw_ali,SimpleScore());

        // compute complete alignments with hirschberg algorithm
        Align<String<ALPHABET>, ArrayGaps> hirsch_ali;
        resize(rows(hirsch_ali), 2);
        assignSource(row(hirsch_ali, 0), str0);
        assignSource(row(hirsch_ali, 1), str1);

        hm_score = globalAlignment(hirsch_ali,SimpleScore(), Hirschberg());

        SEQAN_ASSERT_EQ(nw_score,  hm_score);
        SEQAN_ASSERT_EQ(edit_distance(hirsch_ali), hm_score);

        ++test_count;
    }
}

#endif
