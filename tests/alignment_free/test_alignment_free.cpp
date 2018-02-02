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
// Author: Jonathan Göke <goeke@molgen.mpg.de>
// ==========================================================================
// Tests for the alignment_free module.
// ==========================================================================

#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/sequence.h>
#include <seqan/alignment_free.h>

template <typename TStringSet>
void alfTestHelperGetSequences(TStringSet & sequences)
{
    using namespace seqan;
    CharString seqIID1 =
        "TAGGTTTTCCGAAAAGGTAGCAACTTTACGTGATCAAACCTCTGACGGGGTTTTCCCCGTCGAAATTGGGTG"
        "TTTCTTGTCTTGTTCTCACTTGGGGCATCTCCGTCAAGCCAAGAAAGTGCTCCCTGGATTCTGTTGCTAACG"
        "AGTCTCCTCTGCATTCCTGCTTGACTGATTGGGCGGACGGGGTGTCCACCTGACGCTGAGTATCGCCGTCAC"
        "GGTGCCACATGTCTTATCTATTCAGGGATCAGAATTTATTCAGGAAATCAGGAGATGCTACACTTGGGTTAT"
        "CGAAGCTCCTTCCAAGGCGTAGCAAGGGCGACTGAGCGCGTAAGCTCTAGATCTCCTCGTGTTGCAACTACA"
        "CGCGCGGGTCACTCGAAACACATAGTATGAACTTAACGACTGCTCGTACTGAACAATGCTGAGGCAGAAGAT"
        "CGCAGACCAGGCATCCCACTGCTTGAAAAAACTATNNNNCTACCCGCCTTTTTATTATCTCATCAGATCAAG";
    CharString seqIID2 =
        "ACCGACGATTAGCTTTGTCCGAGTTACAACGGTTCAATAATACAAAGGATGGCATAAACCCATTTGTGTGAA"
        "AGTGCCCATCACATTATGATTCTGTCTACTATGGTTAATTCCCAATATACTCTCGAAAAGAGGGTATGCTCC"
        "CACGGCCATTTACGTCACTAAAAGATAAGATTGCTCAAANNNNNNNNNACTGCCAACTTGCTGGTAGCTTCA"
        "GGGGTTGTCCACAGCGGGGGGTCGTATGCCTTTGTGGTATACCTTACTAGCCGCGCCATGGTGCCTAAGAAT"
        "GAAGTAAAACAATTGATGTGAGACTCGACAGCCAGGCTTCGCGCTAAGGACGCAAAGAAATTCCCTACATCA"
        "GACGGCCGCGNNNAACGATGCTATCGGTTAGGACATTGTGCCCTAGTATGTACATGCCTAATACAATTGGAT"
        "CAAACGTTATTCCCACACACGGGTAGAAGAACNNNNATTACCCGTAGGCACTCCCCGATTCAAGTAGCCGCG";

    clear(sequences);
    appendValue(sequences, seqIID1);
    appendValue(sequences, seqIID2);
}

SEQAN_DEFINE_TEST(test_alignment_free_alignment_free_comparison)
{
    // This test is the example for the function alignmentFreeComparison
    using namespace seqan;
    StringSet<Dna5String> sequences;
    Dna5String seq1 =
        "TAGGTTTTCCGAAAAGGTAGCAACTTTACGTGATCAAACCTCTGACGGGGTTTTCCCCGTCGAAATTGGGTG"
        "TTTCTTGTCTTGTTCTCACTTGGGGCATCTCCGTCAAGCCAAGAAAGTGCTCCCTGGATTCTGTTGCTAACG"
        "AGTCTCCTCTGCATTCCTGCTTGACTGATTGGGCGGACGGGGTGTCCACCTGACGCTGAGTATCGCCGTCAC"
        "GGTGCCACATGTCTTATCTATTCAGGGATCAGAATTTATTCAGGAAATCAGGAGATGCTACACTTGGGTTAT"
        "CGAAGCTCCTTCCAAGGCGTAGCAAGGGCGACTGAGCGCGTAAGCTCTAGATCTCCTCGTGTTGCAACTACA"
        "CGCGCGGGTCACTCGAAACACATAGTATGAACTTAACGACTGCTCGTACTGAACAATGCTGAGGCAGAAGAT"
        "CGCAGACCAGGCATCCCACTGCTTGAAAAAACTATNNNNCTACCCGCCTTTTTATTATCTCATCAGATCAAG";
    Dna5String seq2 =
        "ACCGACGATTAGCTTTGTCCGAGTTACAACGGTTCAATAATACAAAGGATGGCATAAACCCATTTGTGTGAA"
        "AGTGCCCATCACATTATGATTCTGTCTACTATGGTTAATTCCCAATATACTCTCGAAAAGAGGGTATGCTCC"
        "CACGGCCATTTACGTCACTAAAAGATAAGATTGCTCAAANNNNNNNNNACTGCCAACTTGCTGGTAGCTTCA"
        "GGGGTTGTCCACAGCGGGGGGTCGTATGCCTTTGTGGTATACCTTACTAGCCGCGCCATGGTGCCTAAGAAT"
        "GAAGTAAAACAATTGATGTGAGACTCGACAGCCAGGCTTCGCGCTAAGGACGCAAAGAAATTCCCTACATCA"
        "GACGGCCGCGNNNAACGATGCTATCGGTTAGGACATTGTGCCCTAGTATGTACATGCCTAATACAATTGGAT"
        "CAAACGTTATTCCCACACACGGGTAGAAGAACNNNNATTACCCGTAGGCACTCCCCGATTCAAGTAGCCGCG";

    clear(sequences);
    appendValue(sequences, seq1);
    appendValue(sequences, seq2);

    Matrix<double, 2> myMatrix;

    unsigned kmerSize = 5;
    unsigned bgModelOrder = 1;
    String<char>  revCom = "both_strands";
    unsigned mismatches = 1;
    double mismatchWeight = 0.5;
    AFScore<N2> myScoreN2(kmerSize, bgModelOrder, revCom, mismatches, mismatchWeight);

    alignmentFreeComparison(myMatrix, sequences, myScoreN2);

    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 0, 0), 1.0, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 0, 1), -0.129431, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 1, 0), -0.129431, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 1, 1), 1.0, 0.01);
}

SEQAN_DEFINE_TEST(test_alignment_free_calculate_probability)
{
    using namespace seqan;
    double p = 0.0;
    DnaString word = "CCCAAGTTT";
    String<double> model;
    resize(model, 4);
    model[0] = 0.3;  // p(A)
    model[1] = 0.2;  // p(C)
    model[2] = 0.2;  // p(G)
    model[3] = 0.3;  // p(T)
    calculateProbability(p, word, model);
    SEQAN_ASSERT_IN_DELTA(p, 0.00000387, 0.000001);
}

SEQAN_DEFINE_TEST(test_alignment_free_calculate_variance)
{
    using namespace seqan;
    double var = 0.0;
    int n = 10000;
    DnaString word = "CAAGTC";
    String<double> model;
    resize(model, 4);
    model[0] = 0.3;  // p(A)
    model[1] = 0.2;  // p(C)
    model[2] = 0.2;  // p(G)
    model[3] = 0.3;  // p(T)
    calculateVariance(var, word, model, n);  // var = 2.16
    SEQAN_ASSERT_IN_DELTA(var, 2.15845, 0.001);
    StringSet<DnaString> sequences;
    appendValue(sequences, "CAGAAAAAAACACTGATTAACAGGAATAAGCAGTTTACTTATTTTGGGCCTGGGACCCGTGTCTCTAATTTAATTAGGTGATCCCTGCGAAGTTTCTCCA");
    MarkovModel<Dna, double> modelMM0(0);  // Bernoulli model
    modelMM0.build(sequences);
    calculateVariance(var, word, modelMM0, n);  // var = 2.16
    SEQAN_ASSERT_IN_DELTA(var, 2.15845, 0.001);
    MarkovModel<Dna, double> modelMM1(1);  // First order Markov model
    modelMM1.build(sequences);
    calculateVariance(var, word, modelMM1, n);  // var = 1.69716
    SEQAN_ASSERT_IN_DELTA(var, 1.69716, 0.001);
}

SEQAN_DEFINE_TEST(test_alignment_free_calculate_covariance)
{
    using namespace seqan;
    double covar = 0.0;
    int n = 10000;
    DnaString word1 = "ATATAT";
    DnaString word2 = "TATATA";
    String<double> model;
    resize(model, 4);
    model[0] = 0.3;  // p(A)
    model[1] = 0.2;  // p(C)
    model[2] = 0.2;  // p(G)
    model[3] = 0.3;  // p(T)
    calculateCovariance(covar, word1, word2, model, n);  // covar = 4.74
    SEQAN_ASSERT_IN_DELTA(covar, 4.741, 0.001);
    StringSet<DnaString> sequences;
    appendValue(sequences, "CAGCACTGATTAACAGGAATAAGCAGTTTACTTCTGTCAGAATATTGGGCATATATACTGGGACCCGTGTAATACTCTAATTTAATTAGGTGATCCCTGCGAAGTCTCCA");
    MarkovModel<Dna, double> modelMM0(0);  // Bernoulli model
    modelMM0.build(sequences);
    calculateCovariance(covar, word1, word2, modelMM0, n);  // covar = 4.74
    SEQAN_ASSERT_IN_DELTA(covar, 4.741, 0.001);
    MarkovModel<Dna, double> modelMM1(1);  // First order Markov model
    modelMM1.build(sequences);
    calculateCovariance(covar, word1, word2, modelMM1, n);  // covar = 4.74
    SEQAN_ASSERT_IN_DELTA(covar, 13.1541, 0.001);
}

SEQAN_DEFINE_TEST(test_alignment_free_d2_dna)
{
    using namespace seqan;
    StringSet<DnaString> sequences;
    alfTestHelperGetSequences(sequences);

    typedef Matrix<double, 2> TMatrix;
    TMatrix myMatrix;

    // kmerSize = 3
    unsigned kmerSize = 3;
    AFScore<D2> myScoreD2(kmerSize);
    alignmentFreeComparison(myMatrix, sequences, myScoreD2);

    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 0, 0), 4424.0, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 0, 1), 3965.0, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 1, 0), 3965.0, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 1, 1), 4814.0, 0.01);
}

SEQAN_DEFINE_TEST(test_alignment_free_d2_dna5)
{
    using namespace seqan;
    StringSet<Dna5String> sequences;
    alfTestHelperGetSequences(sequences);

    typedef Matrix<double, 2> TMatrix;
    TMatrix myMatrix;

    // kmerSize = 3
    unsigned kmerSize = 3;
    AFScore<D2> myScoreD2(kmerSize);
    alignmentFreeComparison(myMatrix, sequences, myScoreD2);

    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 0, 0), 4322.0, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 0, 1), 3652.0, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 1, 0), 3652.0, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 1, 1), 4030.0, 0.01);

    kmerSize=5;
    myScoreD2.kmerSize=kmerSize;
    alignmentFreeComparison(myMatrix, sequences, myScoreD2);

    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 0, 0), 762.0, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 0, 1), 216.0, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 1, 0), 216.0, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 1, 1), 688.0, 0.01);
}

SEQAN_DEFINE_TEST(test_alignment_free_d2star_dna5)
{
    using namespace seqan;
    StringSet<Dna5String> sequences;
    alfTestHelperGetSequences(sequences);

    typedef Matrix<double, 2> TMatrix;
    TMatrix myMatrix;

    // kmerSize = 3
    unsigned kmerSize = 3;
    unsigned bgModelOrder = 0;
    AFScore<D2Star> myScoreD2Star(kmerSize,bgModelOrder);
    alignmentFreeComparison(myMatrix, sequences, myScoreD2Star);

    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 0, 0), 58.2374, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 0, 1), -10.949, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 1, 0), -10.949, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 1, 1), 48.358, 0.01);

    bgModelOrder = 1;
    myScoreD2Star.bgModelOrder=bgModelOrder;
    alignmentFreeComparison(myMatrix, sequences,myScoreD2Star);

    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 0, 0), 30.336, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 0, 1), -15.104, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 1, 0), -15.104, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 1, 1), 36.3913, 0.01);
}

SEQAN_DEFINE_TEST(test_alignment_free_d2z_dna5)
{
    using namespace seqan;
    StringSet<Dna5String> sequences;
    alfTestHelperGetSequences(sequences);

    typedef Matrix<double, 2> TMatrix;
    TMatrix myMatrix;

    //kmerSize = 3
    unsigned kmerSize = 3;
    unsigned bgModelOrder = 0;
    AFScore<D2z> myScoreD2z(kmerSize,bgModelOrder);
    alignmentFreeComparison(myMatrix, sequences, myScoreD2z);

    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 0, 0), 5.13022, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 0, 1), -0.614828, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 1, 0), -0.614828, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 1, 1), 3.40064, 0.01);

    bgModelOrder=1;
    myScoreD2z.bgModelOrder=bgModelOrder;
    alignmentFreeComparison(myMatrix, sequences,myScoreD2z);

    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 0, 0), 1.61939, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 0, 1), 0.218295, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 1, 0), 0.218295, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 1, 1), 2.47939, 0.01);
}

SEQAN_DEFINE_TEST(test_alignment_free_n2_dna5)
{
    using namespace seqan;
    StringSet<Dna5String> sequences;
    alfTestHelperGetSequences(sequences);

    typedef Matrix<double, 2> TMatrix;
    TMatrix myMatrix;

    unsigned kmerSize = 3;
    unsigned bgModelOrder = 0;
    String<char>  revCom = "";
    unsigned mismatches = 0;
    double mismatchWeight = 0.5;
    AFScore<N2> myScoreN2(kmerSize, bgModelOrder, revCom, mismatches, mismatchWeight);

    alignmentFreeComparison(myMatrix, sequences, myScoreN2);

    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 0, 0), 1.0, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 0, 1), -0.161242, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 1, 0), -0.161242, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 1, 1), 1.0, 0.01);

    bgModelOrder = 1;
    myScoreN2.bgModelOrder = bgModelOrder;
    alignmentFreeComparison(myMatrix, sequences, myScoreN2);

    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 0, 0), 1.0, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 0, 1), 0.143021, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 1, 0), 0.143021, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 1, 1), 1.0, 0.01);

    myScoreN2.revCom = "both_strands";
    alignmentFreeComparison(myMatrix, sequences, myScoreN2);

    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 0, 0), 1.0, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 0, 1), -0.236594, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 1, 0), -0.236594, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 1, 1), 1.0, 0.01);

    mismatches = 1;
    mismatchWeight = 0.5;
    myScoreN2.mismatches = mismatches;
    myScoreN2.mismatchWeight = mismatchWeight;
    alignmentFreeComparison(myMatrix, sequences, myScoreN2);

    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 0, 0), 1.0, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 0, 1), -0.382932, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 1, 0), -0.382932, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 1, 1), 1.0, 0.01);

    bgModelOrder=0;
    myScoreN2.bgModelOrder = bgModelOrder;
    alignmentFreeComparison(myMatrix, sequences, myScoreN2);

    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 0, 0), 1.0, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 0, 1), -0.675514, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 1, 0), -0.675514, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 1, 1), 1.0, 0.01);

    kmerSize = 4;
    myScoreN2.kmerSize = kmerSize;
    alignmentFreeComparison(myMatrix, sequences, myScoreN2);

    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 0, 0), 1.0, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 0, 1), -0.479099, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 1, 0), -0.479099, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 1, 1), 1.0, 0.01);

    bgModelOrder = 2;
    myScoreN2.bgModelOrder = bgModelOrder;
    alignmentFreeComparison(myMatrix, sequences, myScoreN2);

    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 0, 0), 1.0, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 0, 1), -0.0858019, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 1, 0), -0.0858019, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 1, 1), 1.0, 0.01);

}

SEQAN_DEFINE_TEST(test_alignment_free_count_kmers)
{
    using namespace seqan;
    StringSet<Dna5String> sequencesDna5;
    alfTestHelperGetSequences(sequencesDna5);
    String<unsigned> kmerCounts;
    unsigned k = 2;
    countKmers(kmerCounts, sequencesDna5[0], k);
    SEQAN_ASSERT_EQ(kmerCounts[0], 34u);  // AA
    SEQAN_ASSERT_EQ(kmerCounts[1], 30u);  // AC
    SEQAN_ASSERT_EQ(kmerCounts[2], 28u);  // AG
    SEQAN_ASSERT_EQ(kmerCounts[3], 27u);  // AT
    SEQAN_ASSERT_EQ(kmerCounts[4], 32u);  // CA
    SEQAN_ASSERT_EQ(kmerCounts[5], 24u);  // ...
    SEQAN_ASSERT_EQ(kmerCounts[6], 27u);
    SEQAN_ASSERT_EQ(kmerCounts[7], 44u);
    SEQAN_ASSERT_EQ(kmerCounts[8], 31u);
    SEQAN_ASSERT_EQ(kmerCounts[9], 29u);
    SEQAN_ASSERT_EQ(kmerCounts[10], 31u);
    SEQAN_ASSERT_EQ(kmerCounts[11], 27u);
    SEQAN_ASSERT_EQ(kmerCounts[12], 22u);
    SEQAN_ASSERT_EQ(kmerCounts[13], 43u);
    SEQAN_ASSERT_EQ(kmerCounts[14], 33u);  // ...
    SEQAN_ASSERT_EQ(kmerCounts[15], 36u);  // TT

    String<double> nucleotideFrequencies;
    countKmers(kmerCounts, nucleotideFrequencies, sequencesDna5[0], k);

    SEQAN_ASSERT_EQ(kmerCounts[0], 34u);  // AA
    SEQAN_ASSERT_EQ(kmerCounts[1], 30u);  // AC
    SEQAN_ASSERT_EQ(kmerCounts[2], 28u);  // AG
    SEQAN_ASSERT_EQ(kmerCounts[3], 27u);  // AT
    SEQAN_ASSERT_EQ(kmerCounts[4], 32u);  // CA
    SEQAN_ASSERT_EQ(kmerCounts[5], 24u);  // ...
    SEQAN_ASSERT_EQ(kmerCounts[6], 27u);
    SEQAN_ASSERT_EQ(kmerCounts[7], 44u);
    SEQAN_ASSERT_EQ(kmerCounts[8], 31u);
    SEQAN_ASSERT_EQ(kmerCounts[9], 29u);
    SEQAN_ASSERT_EQ(kmerCounts[10], 31u);
    SEQAN_ASSERT_EQ(kmerCounts[11], 27u);
    SEQAN_ASSERT_EQ(kmerCounts[12], 22u);
    SEQAN_ASSERT_EQ(kmerCounts[13], 43u);
    SEQAN_ASSERT_EQ(kmerCounts[14], 33u);  // ...
    SEQAN_ASSERT_EQ(kmerCounts[15], 36u);  // TT

    SEQAN_ASSERT_EQ(nucleotideFrequencies[0],0.238);  // p(A)
    SEQAN_ASSERT_EQ(nucleotideFrequencies[1],0.254);  // p(C)
    SEQAN_ASSERT_EQ(nucleotideFrequencies[2],0.238);  // p(G)
    SEQAN_ASSERT_EQ(nucleotideFrequencies[3],0.27);   // p(T)

    MarkovModel<Dna, double>  backgroundModel(1);
    countKmers(kmerCounts, backgroundModel, sequencesDna5[0], k);
    SEQAN_ASSERT_EQ(kmerCounts[0], 34u);  // AA
    SEQAN_ASSERT_EQ(kmerCounts[1], 30u);  // AC
    SEQAN_ASSERT_EQ(kmerCounts[2], 28u);  // AG
    SEQAN_ASSERT_EQ(kmerCounts[3], 27u);  // AT
    SEQAN_ASSERT_EQ(kmerCounts[4], 32u);  // CA
    SEQAN_ASSERT_EQ(kmerCounts[5], 24u);  // ...
    SEQAN_ASSERT_EQ(kmerCounts[6], 27u);
    SEQAN_ASSERT_EQ(kmerCounts[7], 44u);
    SEQAN_ASSERT_EQ(kmerCounts[8], 31u);
    SEQAN_ASSERT_EQ(kmerCounts[9], 29u);
    SEQAN_ASSERT_EQ(kmerCounts[10], 31u);
    SEQAN_ASSERT_EQ(kmerCounts[11], 27u);
    SEQAN_ASSERT_EQ(kmerCounts[12], 22u);
    SEQAN_ASSERT_EQ(kmerCounts[13], 43u);
    SEQAN_ASSERT_EQ(kmerCounts[14], 33u);  // ...
    SEQAN_ASSERT_EQ(kmerCounts[15], 36u);  // TT
    SEQAN_ASSERT_IN_DELTA(value(backgroundModel.transition, 0, 0), 0.2857143, 0.0001);  // p(A->A)
    SEQAN_ASSERT_IN_DELTA(value(backgroundModel.transition, 1, 2), 0.2125984, 0.0001);  // p(C->G)
    SEQAN_ASSERT_IN_DELTA(value(backgroundModel.transition, 3, 0), 0.1641791, 0.0001);  // p(T->A)
}

SEQAN_DEFINE_TEST(test_alignment_free_calculate_periodicity)
{
    using namespace seqan;
    DnaString word1 = "ATATA";
    DnaString word2 = "TATAT";
    String<int> periodicity;
    calculatePeriodicity(periodicity, word1, word2);
    // periodocity[0] = 1:
    // 01234
    // ATATA
    // -TATAT
    SEQAN_ASSERT_EQ(periodicity[0], 1);
    // periodocity[1] = 3:
    // 01234
    // ATATA
    // ---TATAT
    SEQAN_ASSERT_EQ(periodicity[1], 3);
}

SEQAN_DEFINE_TEST(test_alignment_free_calculate_overlap_indicator)
{
    using namespace seqan;
    DnaString word1 = "ATATA";
    DnaString word2 = "TATAT";
    String<int> epsilon;
    calculateOverlapIndicator(epsilon, word1, word2);
    // epsilon =         01010:
    // word1             ATATA
    // word2 overlap 1:  -TATAT
    // word2 overlap 2:  ---TATAT
    SEQAN_ASSERT_EQ(epsilon[0], 0);
    SEQAN_ASSERT_EQ(epsilon[1], 1);
    SEQAN_ASSERT_EQ(epsilon[2], 0);
    SEQAN_ASSERT_EQ(epsilon[3], 1);
    SEQAN_ASSERT_EQ(epsilon[4], 0);
}


SEQAN_DEFINE_TEST(test_alignment_free_string_to_string_set)
{
    using namespace seqan;
    Dna5String sequenceDna5 =
        "NNNNNNTTTCCGAAAAGGTANNNNNGCAACTTTANNNCGTGATCAAAGTTTTCCCCGTCGAAATTGGGNNTG";
    StringSet<DnaString> sequencesDna;
    stringToStringSet(sequencesDna, sequenceDna5);
    SEQAN_ASSERT_EQ(sequencesDna[0], "TTTCCGAAAAGGTA");
    SEQAN_ASSERT_EQ(sequencesDna[1], "GCAACTTTA");
    SEQAN_ASSERT_EQ(sequencesDna[2], "CGTGATCAAAGTTTTCCCCGTCGAAATTGGG");
    SEQAN_ASSERT_EQ(sequencesDna[3], "TG");
}

SEQAN_DEFINE_TEST(test_alignment_free_cut_ns)
{
    using namespace seqan;
    Dna5String sequenceMasked =
        "NNNNNNTTTCCGAAAAGGTANNNNNGCAACTTTANNNCGTGATCAAAGTTTTCCCCGTCGAAATTGGGNNTG";
    Dna5String sequenceMaskedPartsRemoved;
    cutNs(sequenceMaskedPartsRemoved, sequenceMasked);
    SEQAN_ASSERT_EQ(sequenceMaskedPartsRemoved, "TTTCCGAAAAGGTAGCAACTTTACGTGATCAAAGTTTTCCCCGTCGAAATTGGGTG");
}

SEQAN_BEGIN_TESTSUITE(test_alignment_free)
{
    // Call tests.
    SEQAN_CALL_TEST(test_alignment_free_d2_dna);
    SEQAN_CALL_TEST(test_alignment_free_d2_dna5);
    SEQAN_CALL_TEST(test_alignment_free_d2star_dna5);
    SEQAN_CALL_TEST(test_alignment_free_d2z_dna5);
    SEQAN_CALL_TEST(test_alignment_free_n2_dna5);
    SEQAN_CALL_TEST(test_alignment_free_calculate_probability);
    SEQAN_CALL_TEST(test_alignment_free_calculate_variance);
    SEQAN_CALL_TEST(test_alignment_free_calculate_covariance);
    SEQAN_CALL_TEST(test_alignment_free_alignment_free_comparison);
    SEQAN_CALL_TEST(test_alignment_free_count_kmers);
    SEQAN_CALL_TEST(test_alignment_free_calculate_periodicity);
    SEQAN_CALL_TEST(test_alignment_free_calculate_overlap_indicator);
    SEQAN_CALL_TEST(test_alignment_free_string_to_string_set);
    SEQAN_CALL_TEST(test_alignment_free_cut_ns);
}
SEQAN_END_TESTSUITE
