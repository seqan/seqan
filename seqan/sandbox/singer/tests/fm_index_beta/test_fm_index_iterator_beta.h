// ==========================================================================
//                               fm_index_beta
// ==========================================================================
// Copyright (c) 2006-2011, Knut Reinert, FU Berlin
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
// Author: Your Name <your.email@example.net>
// ==========================================================================

#ifndef TEST_FM_INDEX_ITERATOR_BETA_H_
#define TEST_FM_INDEX_ITERATOR_BETA_H_

#include <seqan/basic.h>
#include <seqan/index.h>
#include <seqan/sequence.h>
#include <seqan/fm_index.h>

#include <seqan/random.h>

using namespace seqan;

template <typename TIter>
void fmIndexIteratorConstuctor(TIter & /*tag*/)
{
	typedef typename Container<TIter>::Type TIndex;
    typedef typename Value<TIndex>::Type TText;

	TText text;
	generateText(text);

	Index<TText> esa(text);
	typename Iterator<Index<TText>, TopDown<> >::Type esaIt(esa);
	goDown(esaIt);

	TIndex fmIndex(text);

	TIter it(fmIndex);

	SEQAN_ASSERT_EQ(isRoot(it), true);
	SEQAN_ASSERT_EQ(repLength(it), 0u);
	SEQAN_ASSERT_EQ(goDown(it), true);
	SEQAN_ASSERT_EQ(isRoot(it), false);
	SEQAN_ASSERT_EQ(repLength(it), 1u);
	SEQAN_ASSERT_EQ(representative(it), 'A');
	SEQAN_ASSERT_EQ(goRight(it), true);
	SEQAN_ASSERT_EQ(representative(it), 'C');
	SEQAN_ASSERT_EQ(repLength(it), 1u);
	SEQAN_ASSERT_EQ(goRight(it), true);
	SEQAN_ASSERT_EQ(representative(it), 'G');
	SEQAN_ASSERT_EQ(repLength(it), 1u);
	SEQAN_ASSERT_EQ(goRight(it), true);
	SEQAN_ASSERT_EQ(representative(it), 'T');
	SEQAN_ASSERT_EQ(repLength(it), 1u);
	SEQAN_ASSERT_EQ(goRight(it), false);
	SEQAN_ASSERT_EQ(representative(it), 'T');
	SEQAN_ASSERT_EQ(repLength(it), 1u);
	SEQAN_ASSERT_EQ(goDown(it), true);
	SEQAN_ASSERT_EQ(goDown(it), true);
	SEQAN_ASSERT_EQ(goUp(it), true);
	SEQAN_ASSERT_EQ(representative(it), "AT");
	SEQAN_ASSERT_EQ(goUp(it), true);
	SEQAN_ASSERT_EQ(representative(it), "T");
	
}

template <typename TIter>
void fmIndexIteratorGoDown(TIter & /*tag*/)
{
    typedef typename Container<TIter>::Type TIndex;
    typedef typename Value<TIndex>::Type TText;

	TText text = "ACGACG";
	TIndex fmIndex(text);
	
    {
        TIter it(fmIndex);
        SEQAN_ASSERT_EQ(goDown(it), true);
        SEQAN_ASSERT_EQ(goDown(it), true);
        SEQAN_ASSERT_EQ(goDown(it), true);
        SEQAN_ASSERT_EQ(goDown(it), true);
        SEQAN_ASSERT_EQ(representative(it), "ACGA");
        SEQAN_ASSERT_EQ(goDown(it), false);
        SEQAN_ASSERT_EQ(representative(it), "ACGA");
    }
    {
        TIter it(fmIndex);
        SEQAN_ASSERT_EQ(goDown(it, 'G'), true);
        SEQAN_ASSERT_EQ(representative(it), "G");
        SEQAN_ASSERT_EQ(goDown(it, 'G'), false);
        SEQAN_ASSERT_EQ(representative(it), "G");
    }
	{
        TIter it(fmIndex);
        SEQAN_ASSERT_EQ(goDown(it, 'G'), true);
        SEQAN_ASSERT_EQ(goDown(it, 'C'), true);
        SEQAN_ASSERT_EQ(goDown(it, 'A'), true);
        SEQAN_ASSERT_EQ(goDown(it, 'G'), true);
        SEQAN_ASSERT_EQ(goDown(it, 'C'), true);
        SEQAN_ASSERT_EQ(goDown(it, 'A'), true);
        SEQAN_ASSERT_EQ(representative(it), "ACGACG");
        SEQAN_ASSERT_EQ(goDown(it, 'G'), false);
        SEQAN_ASSERT_EQ(representative(it), "ACGACG");
    }
	{
        TIter it(fmIndex);
        std::cerr << "START" << std::endl;
        SEQAN_ASSERT_EQ(goDown(it, "ACGACG"), true);
        std::cerr << "END" << std::endl;
        SEQAN_ASSERT_EQ(representative(it), "ACGACG");
        SEQAN_ASSERT_EQ(goDown(it, 'G'), false);
        SEQAN_ASSERT_EQ(representative(it), "ACGACG");
    }
}

template <typename TIter>
void fmIndexIteratorIsLeaf(TIter & /*tag*/)
{
    typedef typename Container<TIter>::Type TIndex;
    typedef typename Value<TIndex>::Type TText;

	TText text = "ACGACG";
	TIndex fmIndex(text);
	
    {
        TIter it(fmIndex);
        SEQAN_ASSERT_EQ(isLeaf(it), false);
        SEQAN_ASSERT_EQ(goDown(it), true);
        SEQAN_ASSERT_EQ(isLeaf(it), false);
        SEQAN_ASSERT_EQ(goDown(it), true);
        SEQAN_ASSERT_EQ(isLeaf(it), false);
        SEQAN_ASSERT_EQ(goDown(it), true);
        SEQAN_ASSERT_EQ(isLeaf(it), false);
        SEQAN_ASSERT_EQ(goDown(it), true);
        SEQAN_ASSERT_EQ(isLeaf(it), true);
        SEQAN_ASSERT_EQ(representative(it), "ACGA");
        SEQAN_ASSERT_EQ(goDown(it), false);
        SEQAN_ASSERT_EQ(representative(it), "ACGA");
    }
}

template <typename TIter>
void fmIndexIteratorGoRight(TIter & /*tag*/)
{
    typedef typename Container<TIter>::Type TIndex;
    typedef typename Value<TIndex>::Type TText;

	TText text = "ACGACG";
	TIndex fmIndex(text);
	
    {
        TIter it(fmIndex);
        SEQAN_ASSERT_EQ(goDown(it), true);
        SEQAN_ASSERT_EQ(representative(it), "A");
        SEQAN_ASSERT_EQ(goRight(it), true);
        SEQAN_ASSERT_EQ(representative(it), "C");
        SEQAN_ASSERT_EQ(goRight(it), true);
        SEQAN_ASSERT_EQ(representative(it), "G");
        SEQAN_ASSERT_EQ(goRight(it), false);
        SEQAN_ASSERT_EQ(representative(it), "G");
    }
}

template <typename TIter>
void fmIndexIteratorGoUp(TIter & /*tag*/)
{
    typedef typename Container<TIter>::Type TIndex;
    typedef typename Value<TIndex>::Type TText;

	TText text = "ACGACG";
	TIndex fmIndex(text);
	
    {
        TIter it(fmIndex);
        SEQAN_ASSERT_EQ(goDown(it, "ACGACG"), true);
        SEQAN_ASSERT_EQ(representative(it), "ACGACG");
        SEQAN_ASSERT_EQ(goUp(it), true);
        SEQAN_ASSERT_EQ(representative(it), "CGACG");
        SEQAN_ASSERT_EQ(goUp(it), true);
        SEQAN_ASSERT_EQ(representative(it), "GACG");
        SEQAN_ASSERT_EQ(goUp(it), true);
        SEQAN_ASSERT_EQ(representative(it), "ACG");
        SEQAN_ASSERT_EQ(goUp(it), true);
        SEQAN_ASSERT_EQ(representative(it), "CG");
        SEQAN_ASSERT_EQ(goUp(it), true);
        SEQAN_ASSERT_EQ(representative(it), "G");
        SEQAN_ASSERT_EQ(goUp(it), true);
        SEQAN_ASSERT_EQ(representative(it), "");
        SEQAN_ASSERT_EQ(goUp(it), false);
        SEQAN_ASSERT_EQ(representative(it), "");
    }
}

template <typename TIter>
void fmIndexIteratorIsRoot(TIter & /*tag*/)
{
    typedef typename Container<TIter>::Type TIndex;
    typedef typename Value<TIndex>::Type TText;

	TText text = "ACGACG";
	TIndex fmIndex(text);
	
    {
        TIter it(fmIndex);

        SEQAN_ASSERT_EQ(isRoot(it), true);
        SEQAN_ASSERT_EQ(goDown(it, "ACGACG"), true);
        SEQAN_ASSERT_EQ(isRoot(it), false);
    }
}

template <typename TIter>
void fmIndexIteratorCountOccurrences(TIter & /*tag*/)
{
    typedef typename Container<TIter>::Type TIndex;
    typedef typename Value<TIndex>::Type TText;
    typedef Index<TText, IndexEsa<> > TEsaIndex;
    typedef typename Iterator<TEsaIndex, TopDown<ParentLinks<> > >::Type TEsaIter;
    
    TText text;
	generateText(text);

	StringSet<String<typename Alphabet<TText>::Type> > pattern;
	generatePattern(pattern, text);

	TIndex fmIndex(text);
	TEsaIndex esaIndex(text);
        
    TIter it(fmIndex);
    TEsaIter esaIt(esaIndex);

    for (unsigned i = 0; i < length(pattern); ++i)
    {
        
        bool _goDown = goDown(it, pattern[i]);
        SEQAN_ASSERT_EQ(_goDown, goDown(esaIt, pattern[i]));
        if (_goDown)
            SEQAN_ASSERT_EQ(countOccurrences(it), countOccurrences(esaIt));
        goRoot(it);
        goRoot(esaIt);
    }
}

template <typename TIter>
void fmIndexIteratorRange(TIter & /*tag*/)
{
    typedef typename Container<TIter>::Type TIndex;
    typedef typename Value<TIndex>::Type TText;
    typedef Index<TText, IndexEsa<> > TEsaIndex;
    typedef typename Iterator<TEsaIndex, TopDown<ParentLinks<> > >::Type TEsaIter;
    
    TText text;
	generateText(text);

	StringSet<String<typename Alphabet<TText>::Type> > pattern;
	generatePattern(pattern, text);

	TIndex fmIndex(text);
    TEsaIndex esaIndex(text);
        
    TIter it(fmIndex);
    TEsaIter esaIt(esaIndex);

    for (unsigned i = 0; i < length(pattern); ++i)
    {
        bool _goDown = goDown(it, pattern[i]);
        SEQAN_ASSERT_EQ(_goDown, goDown(esaIt, pattern[i]));
        if (_goDown)
        {
            SEQAN_ASSERT_EQ(range(it).i1 - 1, range(esaIt).i1);
            SEQAN_ASSERT_EQ(range(it).i2 - 1, range(esaIt).i2);
        }
        goRoot(it);
        goRoot(esaIt);
    }
}

SEQAN_DEFINE_TEST(fm_index_iterator_constuctor)
{
    using namespace seqan;

    typedef FmIndex<WaveletTreeBased<FmiDollarSubstituted<SingleDollar<void> > >, void> TDefaultIndex;
    typedef FmIndex<WaveletTreeBased<FmiDollarSubstituted<SingleDollar<void> > >, CompressText> TCompressedIndex;
    typedef TopDown<> TIterSpec;
    typedef TopDown<ParentLinks<> > TParentLinksIterSpec;

    DnaString genome = "AAA";
    
//     {
//         Index<DnaString,TDefaultIndex> index(genome);
//         Iterator<Index<DnaString,TDefaultIndex>, TIterSpec>::Type dnaTag(index);
//         fmIndexIteratorConstuctor(dnaTag);
//     }
//     {
//         Index<DnaString,TCompressedIndex> index(genome);
//         Iterator<Index<DnaString,TCompressedIndex>, TIterSpec>::Type dnaTag(index);
//         fmIndexIteratorConstuctor(dnaTag);
//     }
    {
        Index<DnaString,TDefaultIndex> index(genome);
        Iterator<Index<DnaString,TDefaultIndex>, TParentLinksIterSpec>::Type dnaTag(index);
        fmIndexIteratorConstuctor(dnaTag);
    }
    {
        Index<DnaString,TCompressedIndex> index(genome);
        Iterator<Index<DnaString,TCompressedIndex>, TParentLinksIterSpec>::Type dnaTag(index);
        fmIndexIteratorConstuctor(dnaTag);
    }
}

SEQAN_DEFINE_TEST(fm_index_iterator_go_down)
{
    using namespace seqan;

    typedef FmIndex<WaveletTreeBased<FmiDollarSubstituted<SingleDollar<void> > >, void> TDefaultIndex;
    typedef FmIndex<WaveletTreeBased<FmiDollarSubstituted<SingleDollar<void> > >, CompressText> TCompressedIndex;
    typedef TopDown<> TIterSpec;
    typedef TopDown<ParentLinks<> > TParentLinksIterSpec;

    DnaString genome = "A";
    
    {
        Index<DnaString,TDefaultIndex> index(genome);
        Iterator<Index<DnaString,TDefaultIndex>, TIterSpec>::Type dnaTag(index);
        fmIndexIteratorGoDown(dnaTag);
    }
    {
        Index<DnaString,TCompressedIndex> index(genome);
        Iterator<Index<DnaString,TCompressedIndex>, TIterSpec>::Type dnaTag(index);
        fmIndexIteratorGoDown(dnaTag);
    }
    {
        Index<DnaString,TDefaultIndex> index(genome);
        Iterator<Index<DnaString,TDefaultIndex>, TParentLinksIterSpec>::Type dnaTag(index);
        fmIndexIteratorGoDown(dnaTag);
    }
    {
        Index<DnaString,TCompressedIndex> index(genome);
        Iterator<Index<DnaString,TCompressedIndex>, TParentLinksIterSpec>::Type dnaTag(index);
        fmIndexIteratorGoDown(dnaTag);
    }
}

SEQAN_DEFINE_TEST(fm_index_iterator_is_leaf)
{
    using namespace seqan;

    typedef FmIndex<WaveletTreeBased<FmiDollarSubstituted<SingleDollar<void> > >, void> TDefaultIndex;
    typedef FmIndex<WaveletTreeBased<FmiDollarSubstituted<SingleDollar<void> > >, CompressText> TCompressedIndex;
    typedef TopDown<> TIterSpec;
    typedef TopDown<ParentLinks<> > TParentLinksIterSpec;

    DnaString genome = "A";
    
    {
        Index<DnaString,TDefaultIndex> index(genome);
        Iterator<Index<DnaString,TDefaultIndex>, TIterSpec>::Type dnaTag(index);
        fmIndexIteratorIsLeaf(dnaTag);
    }
    {
        Index<DnaString,TCompressedIndex> index(genome);
        Iterator<Index<DnaString,TCompressedIndex>, TIterSpec>::Type dnaTag(index);
        fmIndexIteratorIsLeaf(dnaTag);
    }
    {
        Index<DnaString,TDefaultIndex> index(genome);
        Iterator<Index<DnaString,TDefaultIndex>, TParentLinksIterSpec>::Type dnaTag(index);
        fmIndexIteratorIsLeaf(dnaTag);
    }
    {
        Index<DnaString,TCompressedIndex> index(genome);
        Iterator<Index<DnaString,TCompressedIndex>, TParentLinksIterSpec>::Type dnaTag(index);
        fmIndexIteratorIsLeaf(dnaTag);
    }
}

SEQAN_DEFINE_TEST(fm_index_iterator_go_right)
{
    using namespace seqan;

    typedef FmIndex<WaveletTreeBased<FmiDollarSubstituted<SingleDollar<void> > >, void> TDefaultIndex;
    typedef FmIndex<WaveletTreeBased<FmiDollarSubstituted<SingleDollar<void> > >, CompressText> TCompressedIndex;
    typedef TopDown<> TIterSpec;
    typedef TopDown<ParentLinks<> > TParentLinksIterSpec;

    DnaString genome = "A";
    
    {
        Index<DnaString,TDefaultIndex> index(genome);
        Iterator<Index<DnaString,TDefaultIndex>, TIterSpec>::Type dnaTag(index);
        fmIndexIteratorGoRight(dnaTag);
    }
    {
        Index<DnaString,TCompressedIndex> index(genome);
        Iterator<Index<DnaString,TCompressedIndex>, TIterSpec>::Type dnaTag(index);
        fmIndexIteratorGoRight(dnaTag);
    }
    {
        Index<DnaString,TDefaultIndex> index(genome);
        Iterator<Index<DnaString,TDefaultIndex>, TParentLinksIterSpec>::Type dnaTag(index);
        fmIndexIteratorGoRight(dnaTag);
    }
    {
        Index<DnaString,TCompressedIndex> index(genome);
        Iterator<Index<DnaString,TCompressedIndex>, TParentLinksIterSpec>::Type dnaTag(index);
        fmIndexIteratorGoRight(dnaTag);
    }
}

SEQAN_DEFINE_TEST(fm_index_iterator_go_up)
{
    using namespace seqan;

    typedef FmIndex<WaveletTreeBased<FmiDollarSubstituted<SingleDollar<void> > >, void> TDefaultIndex;
    typedef FmIndex<WaveletTreeBased<FmiDollarSubstituted<SingleDollar<void> > >, CompressText> TCompressedIndex;
    typedef TopDown<> TIterSpec;
    typedef TopDown<ParentLinks<> > TParentLinksIterSpec;

    DnaString genome = "A";
    
    {
        Index<DnaString,TDefaultIndex> index(genome);
        Iterator<Index<DnaString,TDefaultIndex>, TParentLinksIterSpec>::Type dnaTag(index);
        fmIndexIteratorGoUp(dnaTag);
    }
    {
        Index<DnaString,TCompressedIndex> index(genome);
        Iterator<Index<DnaString,TCompressedIndex>, TParentLinksIterSpec>::Type dnaTag(index);
        fmIndexIteratorGoUp(dnaTag);
    }
}

SEQAN_DEFINE_TEST(fm_index_iterator_is_root)
{
    using namespace seqan;

    typedef FmIndex<WaveletTreeBased<FmiDollarSubstituted<SingleDollar<void> > >, void> TDefaultIndex;
    typedef FmIndex<WaveletTreeBased<FmiDollarSubstituted<SingleDollar<void> > >, CompressText> TCompressedIndex;
    typedef TopDown<> TIterSpec;
    typedef TopDown<ParentLinks<> > TParentLinksIterSpec;

    DnaString genome = "A";
    
    {
        Index<DnaString,TDefaultIndex> index(genome);
        Iterator<Index<DnaString,TDefaultIndex>, TIterSpec>::Type dnaTag(index);
        fmIndexIteratorIsRoot(dnaTag);
    }
    {
        Index<DnaString,TCompressedIndex> index(genome);
        Iterator<Index<DnaString,TCompressedIndex>, TIterSpec>::Type dnaTag(index);
        fmIndexIteratorIsRoot(dnaTag);
    }
    {
        Index<DnaString,TDefaultIndex> index(genome);
        Iterator<Index<DnaString,TDefaultIndex>, TParentLinksIterSpec>::Type dnaTag(index);
        fmIndexIteratorIsRoot(dnaTag);
    }
    {
        Index<DnaString,TCompressedIndex> index(genome);
        Iterator<Index<DnaString,TCompressedIndex>, TParentLinksIterSpec>::Type dnaTag(index);
        fmIndexIteratorIsRoot(dnaTag);
    }
}

SEQAN_DEFINE_TEST(fm_index_iterator_count_occurrences)
{
    using namespace seqan;

    typedef FmIndex<WaveletTreeBased<FmiDollarSubstituted<SingleDollar<void> > >, void> TDefaultIndex;
    typedef FmIndex<WaveletTreeBased<FmiDollarSubstituted<SingleDollar<void> > >, CompressText> TCompressedIndex;
    typedef TopDown<> TIterSpec;
    typedef TopDown<ParentLinks<> > TParentLinksIterSpec;

    DnaString genome = "A";
    
    {
        Index<DnaString,TDefaultIndex> index(genome);
        Iterator<Index<DnaString,TDefaultIndex>, TIterSpec>::Type dnaTag(index);
        fmIndexIteratorCountOccurrences(dnaTag);
    }
    {
        Index<DnaString,TCompressedIndex> index(genome);
        Iterator<Index<DnaString,TCompressedIndex>, TIterSpec>::Type dnaTag(index);
        fmIndexIteratorCountOccurrences(dnaTag);
    }
    {
        Index<DnaString,TDefaultIndex> index(genome);
        Iterator<Index<DnaString,TDefaultIndex>, TParentLinksIterSpec>::Type dnaTag(index);
        fmIndexIteratorCountOccurrences(dnaTag);
    }
    {
        Index<DnaString,TCompressedIndex> index(genome);
        Iterator<Index<DnaString,TCompressedIndex>, TParentLinksIterSpec>::Type dnaTag(index);
        fmIndexIteratorCountOccurrences(dnaTag);
    }
}

SEQAN_DEFINE_TEST(fm_index_iterator_range)
{
    using namespace seqan;

    typedef FmIndex<WaveletTreeBased<FmiDollarSubstituted<SingleDollar<void> > >, void> TDefaultIndex;
    typedef FmIndex<WaveletTreeBased<FmiDollarSubstituted<SingleDollar<void> > >, CompressText> TCompressedIndex;
    typedef TopDown<> TIterSpec;
    typedef TopDown<ParentLinks<> > TParentLinksIterSpec;

    DnaString genome = "A";
    
    {
        Index<DnaString,TDefaultIndex> index(genome);
        Iterator<Index<DnaString,TDefaultIndex>, TIterSpec>::Type dnaTag(index);
        fmIndexIteratorRange(dnaTag);
    }
    {
        Index<DnaString,TCompressedIndex> index(genome);
        Iterator<Index<DnaString,TCompressedIndex>, TIterSpec>::Type dnaTag(index);
        fmIndexIteratorRange(dnaTag);
    }
    {
        Index<DnaString,TDefaultIndex> index(genome);
        Iterator<Index<DnaString,TDefaultIndex>, TParentLinksIterSpec>::Type dnaTag(index);
        fmIndexIteratorRange(dnaTag);
    }
    {
        Index<DnaString,TCompressedIndex> index(genome);
        Iterator<Index<DnaString,TCompressedIndex>, TParentLinksIterSpec>::Type dnaTag(index);
        fmIndexIteratorRange(dnaTag);
    }
}



#endif // TEST_FM_INDEX_ITERATOR_BETA_H_

