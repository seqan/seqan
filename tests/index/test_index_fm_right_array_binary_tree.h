// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
// Copyright (c) 2013 NVIDIA Corporation
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
// Author: Jochen Singer <jochen.singer@fu-berlin.de>
// ==========================================================================

#ifndef TESTS_WAVELT_TREE_STRUCTURE_BETA_H_
#define TESTS_WAVELT_TREE_STRUCTURE_BETA_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/index.h>

using namespace seqan;

template <typename TIter, typename TBorderString, typename TChar>
void determineBorders(TIter & it, TBorderString & borderString, TChar left, TChar right)
{
    TIter it2 = it;
    borderString[getPosition(it)] = Pair<TChar>(left, right);

    TChar pivot = getCharacter(it);
    if (goLeftChild(it))
    {
        determineBorders(it, borderString, left, (TChar)(ordValue(pivot) - 1));
    }
    if (goRightChild(it2))
    {
        determineBorders(it2, borderString, pivot, right);
    }
}


template <typename TRightArrayBinaryTree>
void waveletTreeStructureConstructor(TRightArrayBinaryTree & /*tag*/)
{
    typedef typename Fibre<TRightArrayBinaryTree, FibreTreeStructureEncoding>::Type TWaveletTreeVertices;
    typedef typename Value<TWaveletTreeVertices>::Type TWaveletTreeVertex;
    typedef typename Value<TWaveletTreeVertex, 1>::Type TChar;
    //typedef typename Value<TWaveletTreeVertex, 2>::Type TPos;

    {
        TRightArrayBinaryTree waveletTreeStructure;
        SEQAN_ASSERT_EQ(length(waveletTreeStructure.treeVertices), 0u);
    }
    {
        String<TChar> text;
        generateText(text);

        String<Pair<TChar> > borderString;
        resize(borderString, ValueSize<TChar>::VALUE);

        String<unsigned> freq;
        resize(freq, ValueSize<TChar>::VALUE, 0);
        _getFrequencies(freq, text);

        TRightArrayBinaryTree waveletTreeStructure(text);
        typename Iterator<TRightArrayBinaryTree, TopDown<> >::Type it(waveletTreeStructure, 0u);
        determineBorders(it, borderString, (TChar)0, (TChar)(ValueSize<TChar>::VALUE - 1));

        for (unsigned i = 0; i < length(getFibre(waveletTreeStructure, FibreTreeStructureEncoding())); ++i)
        {
            // determine the current difference
            TChar left = borderString[i].i1;
            int leftCounter = 0;
            for (; left < getFibre(waveletTreeStructure, FibreTreeStructureEncoding())[i].i1; ++left)
            {
                leftCounter += freq[ordValue(left)];
            }
            int rightCounter = 0;

            for (; left <= borderString[i].i2 && left > borderString[i].i1; ++left)
            {
                rightCounter+= freq[ordValue(left)];
            }
            int currentSum = leftCounter - rightCounter;
            if (currentSum < 0)
                currentSum *= -1;

            for (TChar newPivot = (unsigned)borderString[i].i1 + 1; newPivot <= borderString[i].i2 && newPivot > borderString[i].i1; ++newPivot)
            {
                TChar newLeft = borderString[i].i1;
                int newLeftCounter = 0;
                for (; newLeft < newPivot; ++newLeft)
                {
                    newLeftCounter += freq[ordValue(newLeft)];
                }
                int newRightCounter = 0;
                for (; newLeft <= borderString[i].i2 && newLeft > borderString[i].i1; ++newLeft)
                {
                    newRightCounter += freq[ordValue(newLeft)];
                }
                int newCurrentSum = newLeftCounter - newRightCounter;
                if (newCurrentSum < 0)
                    newCurrentSum *= -1;
                SEQAN_ASSERT_LEQ(currentSum, newCurrentSum);
            }
        }
    }
}

template <typename TRightArrayBinaryTree>
void waveletTreeStructureClear(TRightArrayBinaryTree & /*tag*/)
{
    typedef typename Fibre<TRightArrayBinaryTree, FibreTreeStructureEncoding>::Type TWaveletTreeVertices;
    typedef typename Value<TWaveletTreeVertices>::Type TWaveletTreeVertex;
    typedef typename Value<TWaveletTreeVertex, 1>::Type TChar;
    //typedef typename Value<TWaveletTreeVertex, 2>::Type TPos;

    {
        TRightArrayBinaryTree waveletTreeStructure;
        SEQAN_ASSERT_EQ(length(waveletTreeStructure.treeVertices), 0u);
        clear(waveletTreeStructure);
        SEQAN_ASSERT_EQ(length(waveletTreeStructure.treeVertices), 0u);
    }
    {
        String<TChar> text;
        generateText(text);

        TRightArrayBinaryTree waveletTreeStructure(text);
        SEQAN_ASSERT_NEQ(length(waveletTreeStructure.treeVertices), 0u);
        clear(waveletTreeStructure);
        SEQAN_ASSERT_EQ(length(waveletTreeStructure.treeVertices), 0u);
    }
}

template <typename TRightArrayBinaryTree>
void waveletTreeStructureEmpty(TRightArrayBinaryTree & /*tag*/)
{
    typedef typename Fibre<TRightArrayBinaryTree, FibreTreeStructureEncoding>::Type TWaveletTreeVertices;
    typedef typename Value<TWaveletTreeVertices>::Type TWaveletTreeVertex;
    typedef typename Value<TWaveletTreeVertex, 1>::Type TChar;
    //typedef typename Value<TWaveletTreeVertex, 2>::Type TPos;

    {
        TRightArrayBinaryTree waveletTreeStructure;
        SEQAN_ASSERT_EQ(empty(waveletTreeStructure), true);
        clear(waveletTreeStructure);
        SEQAN_ASSERT_EQ(empty(waveletTreeStructure), true);
    }
    {
        String<TChar> text;
        generateText(text);

        TRightArrayBinaryTree waveletTreeStructure(text);
        SEQAN_ASSERT_EQ(empty(waveletTreeStructure), false);
        clear(waveletTreeStructure);
        SEQAN_ASSERT_EQ(empty(waveletTreeStructure), true);
    }
}

template <typename TRightArrayBinaryTree>
void waveletTreeStructureGetFibre(TRightArrayBinaryTree & /*tag*/)
{
    typedef typename Fibre<TRightArrayBinaryTree, FibreTreeStructureEncoding>::Type TWaveletTreeVertices;
    //typedef typename Value<TWaveletTreeVertices>::Type TWaveletTreeVertex;
    //typedef typename Value<TWaveletTreeVertex, 1>::Type TChar;
    //typedef typename Value<TWaveletTreeVertex, 2>::Type TPos;

    {
        TRightArrayBinaryTree waveletTreeStructure;

        TWaveletTreeVertices & temp = getFibre(waveletTreeStructure, FibreTreeStructureEncoding());
        SEQAN_ASSERT_EQ(length(waveletTreeStructure.treeVertices), 0u);

        resize(temp, 10);

        SEQAN_ASSERT_EQ(length(waveletTreeStructure.treeVertices), 10u);
    }
}

template <typename TRightArrayBinaryTree>
void waveletTreeStructureLength(TRightArrayBinaryTree & /*tag*/)
{
    typedef typename Fibre<TRightArrayBinaryTree, FibreTreeStructureEncoding>::Type TWaveletTreeVertices;
    //typedef typename Value<TWaveletTreeVertices>::Type TWaveletTreeVertex;
    //typedef typename Value<TWaveletTreeVertex, 1>::Type TChar;
    //typedef typename Value<TWaveletTreeVertex, 2>::Type TPos;

    {
        TRightArrayBinaryTree waveletTreeStructure;

        TWaveletTreeVertices & temp = getFibre(waveletTreeStructure, FibreTreeStructureEncoding());
        SEQAN_ASSERT_EQ(length(waveletTreeStructure), 0u);

        resize(temp, 10);

        SEQAN_ASSERT_EQ(length(waveletTreeStructure), 10u);
    }
}

template <typename TRightArrayBinaryTree>
void waveletTreeStructureResize(TRightArrayBinaryTree & /*tag*/)
{
    typedef typename Fibre<TRightArrayBinaryTree, FibreTreeStructureEncoding>::Type TWaveletTreeVertices;
    typedef typename Value<TWaveletTreeVertices>::Type TWaveletTreeVertex;
    //typedef typename Value<TWaveletTreeVertex, 1>::Type TChar;
    //typedef typename Value<TWaveletTreeVertex, 2>::Type TPos;

    {
        TRightArrayBinaryTree waveletTreeStructure;

        TWaveletTreeVertices & temp = getFibre(waveletTreeStructure, FibreTreeStructureEncoding());
        SEQAN_ASSERT_EQ(length(waveletTreeStructure), 0u);

        resize(temp, 10);

        SEQAN_ASSERT_EQ(length(waveletTreeStructure), 10u);
    }
    {
        TRightArrayBinaryTree waveletTreeStructure;

        TWaveletTreeVertices & temp = getFibre(waveletTreeStructure, FibreTreeStructureEncoding());
        SEQAN_ASSERT_EQ(length(waveletTreeStructure), 0u);

        resize(temp, 10, TWaveletTreeVertex(0, 10));

        SEQAN_ASSERT_EQ(length(waveletTreeStructure), 10u);
        SEQAN_ASSERT_EQ(getFibre(waveletTreeStructure, FibreTreeStructureEncoding())[9],  TWaveletTreeVertex(0, 10));
    }

}

// NOTE(esiragusa): OpenSave test should not use operator==()
//template <typename TRightArrayBinaryTree>
//void waveletTreeStructureOpenSave(TRightArrayBinaryTree & /*tag*/)
//{
//    //typedef typename Fibre<TRightArrayBinaryTree, FibreTreeStructureEncoding>::Type TWaveletTreeVertices;
//    //typedef typename Value<TWaveletTreeVertices>::Type TWaveletTreeVertex;
//    //typedef typename Value<TWaveletTreeVertex, 1>::Type TChar;
//    //typedef typename Value<TWaveletTreeVertex, 2>::Type TPos;
//
//
//    TRightArrayBinaryTree waveletTreeStructure;
//    _resize(waveletTreeStructure, 10, Exact());
//
//    CharString tempFilename = SEQAN_TEMP_FILENAME();
//
//    save(waveletTreeStructure, toCString(tempFilename));
//
//    TRightArrayBinaryTree openWaveletTreeStructure;
//    open(openWaveletTreeStructure, toCString(tempFilename));
//
//    SEQAN_ASSERT(waveletTreeStructure == openWaveletTreeStructure);
//
//    waveletTreeStructure = openWaveletTreeStructure;
//
//    SEQAN_ASSERT(waveletTreeStructure == openWaveletTreeStructure);
//}

SEQAN_DEFINE_TEST(wavelet_tree_structure_constructor)
{
    using namespace seqan;

    RightArrayBinaryTree<Dna, void> dnaTag;
    RightArrayBinaryTree<Dna5, void> dna5Tag;
    RightArrayBinaryTree<AminoAcid, void> asTag;
    RightArrayBinaryTree<signed char, void> charTag;
    RightArrayBinaryTree<unsigned char, void> uCharTag;
    waveletTreeStructureConstructor(dnaTag);
    waveletTreeStructureConstructor(dna5Tag);
    waveletTreeStructureConstructor(asTag);
    waveletTreeStructureConstructor(charTag);
    waveletTreeStructureConstructor(uCharTag);
}

SEQAN_DEFINE_TEST(wavelet_tree_structure_clear)
{
    using namespace seqan;

    RightArrayBinaryTree<Dna, void> dnaTag;
    RightArrayBinaryTree<Dna5, void> dna5Tag;
    RightArrayBinaryTree<AminoAcid, void> asTag;
    RightArrayBinaryTree<char, void> charTag;
    RightArrayBinaryTree<unsigned char, void> uCharTag;
    waveletTreeStructureClear(dnaTag);
    waveletTreeStructureClear(dna5Tag);
    waveletTreeStructureClear(asTag);
    waveletTreeStructureClear(charTag);
    waveletTreeStructureClear(uCharTag);
}

SEQAN_DEFINE_TEST(wavelet_tree_structure_empty)
{
    using namespace seqan;

    RightArrayBinaryTree<Dna, void> dnaTag;
    RightArrayBinaryTree<Dna5, void> dna5Tag;
    RightArrayBinaryTree<AminoAcid, void> asTag;
    RightArrayBinaryTree<char, void> charTag;
    RightArrayBinaryTree<unsigned char, void> uCharTag;
    waveletTreeStructureEmpty(dnaTag);
    waveletTreeStructureEmpty(dna5Tag);
    waveletTreeStructureEmpty(asTag);
    waveletTreeStructureEmpty(charTag);
    waveletTreeStructureEmpty(uCharTag);
}

SEQAN_DEFINE_TEST(wavelet_tree_structure_get_fibre)
{
    using namespace seqan;

    RightArrayBinaryTree<Dna, void> dnaTag;
    RightArrayBinaryTree<Dna5, void> dna5Tag;
    RightArrayBinaryTree<AminoAcid, void> asTag;
    RightArrayBinaryTree<char, void> charTag;
    RightArrayBinaryTree<unsigned char, void> uCharTag;
    waveletTreeStructureGetFibre(dnaTag);
    waveletTreeStructureGetFibre(dna5Tag);
    waveletTreeStructureGetFibre(asTag);
    waveletTreeStructureGetFibre(charTag);
    waveletTreeStructureGetFibre(uCharTag);
}

SEQAN_DEFINE_TEST(wavelet_tree_structure_length)
{
    using namespace seqan;

    RightArrayBinaryTree<Dna, void> dnaTag;
    RightArrayBinaryTree<Dna5, void> dna5Tag;
    RightArrayBinaryTree<AminoAcid, void> asTag;
    RightArrayBinaryTree<char, void> charTag;
    RightArrayBinaryTree<unsigned char, void> uCharTag;
    waveletTreeStructureLength(dnaTag);
    waveletTreeStructureLength(dna5Tag);
    waveletTreeStructureLength(asTag);
    waveletTreeStructureLength(charTag);
    waveletTreeStructureLength(uCharTag);
}

SEQAN_DEFINE_TEST(wavelet_tree_structure_resize)
{
    using namespace seqan;

    RightArrayBinaryTree<Dna, void> dnaTag;
    RightArrayBinaryTree<Dna5, void> dna5Tag;
    RightArrayBinaryTree<AminoAcid, void> asTag;
    RightArrayBinaryTree<char, void> charTag;
    RightArrayBinaryTree<unsigned char, void> uCharTag;
    waveletTreeStructureResize(dnaTag);
    waveletTreeStructureResize(dna5Tag);
    waveletTreeStructureResize(asTag);
    waveletTreeStructureResize(charTag);
    waveletTreeStructureResize(uCharTag);
}

// NOTE(esiragusa): OpenSave test should not use operator==()
//SEQAN_DEFINE_TEST(wavelet_tree_structure_open_save)
//{
//    using namespace seqan;
//
//    RightArrayBinaryTree<AminoAcid, void> asTag;
//    waveletTreeStructureOpenSave(asTag);
//}


#endif  // TESTS_WAVELT_TREE_STRUCTURE_BETA_H_
