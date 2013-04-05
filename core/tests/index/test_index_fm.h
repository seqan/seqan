// ==========================================================================
//                               fm_index_beta
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
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

#ifndef TEST_FM_INDEX_BETA_H_
#define TEST_FM_INDEX_BETA_H_

#include <seqan/basic.h>
#include <seqan/index.h>
#include <seqan/sequence.h>
#include <seqan/random.h>

using namespace seqan;


unsigned const SEED = 41;

template< typename TText>
void generateText(TText & text, unsigned num)
{

    typedef typename Value<TText>::Type TChar;
    unsigned textLength = num;

    int minChar = MinValue<TChar>::VALUE;
    unsigned alphabetSize = ValueSize<TChar>::VALUE;

    Rng<MersenneTwister> rng(SEED);


    resize(text, textLength);

    for (unsigned i = 0; i < textLength; ++i)
        text[i] = pickRandomNumber(rng) % alphabetSize - minChar;
}

template< typename TText>
void generateText(TText & text)
{
    generateText(text, 100000);
}

void generateText(String<char> & text, unsigned num)
{

    typedef char TChar;
    unsigned textLength = num;

    int minChar = -128;
    unsigned alphabetSize = ValueSize<TChar>::VALUE;

    Rng<MersenneTwister> rng(SEED);

    resize(text, textLength);

    for (unsigned i = 0; i < textLength; ++i)
        text[i] = pickRandomNumber(rng) % alphabetSize - minChar;
}

void generateText(String<char> & text)
{
    generateText(text, 100000);
}

template< typename TText>
void generateText(StringSet<TText> & text, unsigned num)
{

    typedef typename Value<TText>::Type TChar;
    unsigned numSeq = num;
    unsigned seqLength = 2000;

    int minChar = MinValue<TChar>::VALUE;
    unsigned alphabetSize = ValueSize<TChar>::VALUE;

    Rng<MersenneTwister> rng(SEED);

    resize(text, numSeq);

    for (unsigned i = 0; i < numSeq; ++i)
    {
        resize(text[i], pickRandomNumber(rng) % (seqLength-1) + 1);
        for (unsigned j = 0; j < length(text[i]); ++j)
            text[i][j] = pickRandomNumber(rng) % alphabetSize - minChar;
    }
}

template< typename TText>
void generateText(StringSet<TText> & text)
{
    generateText(text, 1000);
}

template< typename TText>
void generatePattern(StringSet<TText> & pattern, TText const & text)
{

    typedef typename Value<TText>::Type TChar;
    unsigned patternLength = 1000;

    int minChar = MinValue<TChar>::VALUE;
    unsigned alphabetSize = ValueSize<TChar>::VALUE;

    Rng<MersenneTwister> rng(SEED);


    resize(pattern, patternLength);

    for (unsigned i = 0; i < patternLength; i = i + 2)
    {
        TText localPattern;
        for (unsigned j = 0; j <= i; ++j)
            appendValue(pattern[i], (TChar)(pickRandomNumber(rng) % alphabetSize - minChar));
        unsigned readPos = pickRandomNumber(rng) % (length(text) - i - 1);
        pattern[i + 1] = infix(text, readPos, readPos + i + 1);
    }
}

template< typename TText>
void generatePattern(StringSet<TText> & pattern, StringSet<TText> const & text)
{

    //typedef typename Value<TText>::Type TChar;
    unsigned patternLength = 1000;

    Rng<MersenneTwister> rng(SEED);

    resize(pattern, patternLength);

    for (unsigned i = 0; i < patternLength; ++i)
    {
        unsigned readLocation;
        do
        {
            readLocation = pickRandomNumber(rng) % length(text);
        } while (length(text[readLocation]) < (i + 10));

        TText localPattern;
        //for (unsigned j = 0; j <= i; ++j)
          //  appendValue(pattern[i], (TChar)(pickRandomNumber(rng) % alphabetSize - minChar));
        unsigned readPos = pickRandomNumber(rng) % (length(text[readLocation]) - i - 1);
        pattern[i] = infix(text[readLocation], readPos, readPos + i + 1);
    }
}

template <typename TText, typename TIndexSpec, typename TOptimization>
void fmIndexConstructor(Index<TText, FMIndex<TIndexSpec, TOptimization> > /*tag*/)
{
	typedef Index<TText, FMIndex<TIndexSpec, TOptimization> > TIndex;
	//typedef Index<TText, IndexEsa<> > TIndexEsa;
	typedef typename Value<TText>::Type TChar;

    {
	    TText text;
	    TIndex fmiIndex(text);
    }
    {
	    TText text;
        generateText(text);
        TIndex fmiIndex(text);
        indexCreate(fmiIndex);

        // checking the BWT
        unsigned pos = 0;
        for (unsigned i = 0; i < length(text); ++ i)
        {
//             static_cast<Nothing>(getFibre(getFibre(fmiIndex, FibreLfTable()), FibreOccTable()));
            TChar character = getValue(getFibre(getFibre(fmiIndex, FibreLfTable()), FibreOccTable()), pos);
            SEQAN_ASSERT_EQ(character, text[length(text) - 1 - i]);
            pos = lfMapping(getFibre(fmiIndex, FibreLfTable()), pos);
        }

        // checking the occurrence table
        StringSet<String<unsigned> > occ;
        resize(occ, ValueSize<TChar>::VALUE);
        for (unsigned i = 0; i <  ValueSize<TChar>::VALUE; ++i)
            resize(occ[i], length(text) + 1, 0);

        ++occ[ordValue(getValue(getFibre(getFibre(fmiIndex, FibreLfTable()), FibreOccTable()), 0u))][0];
        for (unsigned i = 1; i < length(getFibre(getFibre(fmiIndex, FibreLfTable()), FibreOccTable())); ++i)
        {
            for (unsigned j = 0; j <  ValueSize<TChar>::VALUE; ++j)
                occ[j][i] = occ[j][i - 1];

            TChar temp = getValue(getFibre(getFibre(fmiIndex, FibreLfTable()), FibreOccTable()), i);
            ++occ[ordValue(temp)][i];
            if (temp == getSentinelSubstitute(getFibre(getFibre(fmiIndex, FibreLfTable()), FibreOccTable())) &&
                i == _getSentinelPosition(getFibre(getFibre(fmiIndex, FibreLfTable()), FibreOccTable())))
                --occ[ordValue(temp)][i];
        }

        unsigned counter = 0;
        for (TChar character = MinValue<TChar>::VALUE; counter < ValueSize<TChar>::VALUE; ++character)
        {
            ++counter;
            for (unsigned i = 0; i < length(getFibre(getFibre(fmiIndex, FibreLfTable()), FibreOccTable())); ++i)
                SEQAN_ASSERT_EQ(countOccurrences(getFibre(getFibre(fmiIndex, FibreLfTable()), FibreOccTable()), character, i), occ[ordValue(character)][i]);
        }

        // check the compressed suffix array
        String<unsigned> sa;
        resize(sa, length(text));
        createSuffixArray(sa, text, Skew7());

        for (unsigned i = 0; i < length(sa); ++i)
            SEQAN_ASSERT_EQ(getFibre(fmiIndex, FibreSA())[i + 1], sa[i]);

    }
}

template <typename TText, typename TIndexSpec, typename TOptimization>
void fmIndexClear(Index<TText, FMIndex<TIndexSpec, TOptimization> > /*tag*/)
{
	typedef Index<TText, FMIndex<TIndexSpec, TOptimization> > TIndex;
	//typedef typename Value<TText>::Type TChar;
	//typedef typename Fibre<TIndex, FibreLfTable>::Type TLfTable;
	//typedef typename Fibre<TLfTable, FibrePrefixSumTable>::Type TPrefixSumTable;

	//typedef Index<TText, IndexEsa<> > TIndexEsa;

    {
        TIndex fmIndex;
	    SEQAN_ASSERT_LEQ(empty(fmIndex), true);
    }
    {
	    TText text;
	    generateText(text);

	    TIndex fmIndex(text);

        clear(fmIndex);

        TIndex controlIndex;
	    //SEQAN_ASSERT(getFibre(fmIndex, FibreLfTable()) == getFibre(controlIndex, FibreLfTable()));
	    //SEQAN_ASSERT(getFibre(fmIndex, FibreLfTable()) == getFibre(controlIndex, FibreLfTable()));
	    SEQAN_ASSERT(empty(fmIndex) == empty(controlIndex));
	    SEQAN_ASSERT(empty(fmIndex) == empty(controlIndex));
    }
}

template <typename TText, typename TIndexSpec, typename TOptimization>
void _fmIndexDetermineSentinelSubstitute(Index<TText, FMIndex<TIndexSpec, TOptimization> > /*tag*/)
{
	typedef Index<TText, FMIndex<TIndexSpec, TOptimization> > TIndex;
	typedef typename Value<TText>::Type TChar;
	typedef typename Fibre<TIndex, FibreLfTable>::Type TLfTable;
	typedef typename Fibre<TLfTable, FibrePrefixSumTable>::Type TPrefixSumTable;

	//typedef Index<TText, IndexEsa<> > TIndexEsa;

	TText text;
	generateText(text);

	String<unsigned> freq;
	resize(freq, ValueSize<TChar>::VALUE, 0);

	for (unsigned i = 0; i < length(text); ++i)
	    ++freq[ordValue(text[i])];
	
	TPrefixSumTable pst(text);
	TChar sentinelSub;

	_determineSentinelSubstitute(pst, sentinelSub);

	for (unsigned i = 0; i < length(freq); ++i)
    {
        if(freq[i])
    	    SEQAN_ASSERT_LEQ(freq[ordValue(sentinelSub)], freq[i]);
    }
}

template <typename TText, typename TIndexSpec, typename TOptimization>
void fmIndexEmpty(Index<TText, FMIndex<TIndexSpec, TOptimization> > /*tag*/)
{
	typedef Index<TText, FMIndex<TIndexSpec, TOptimization> > TIndex;
	//typedef typename Value<TText>::Type TChar;
	//typedef typename Fibre<TIndex, FibreLfTable>::Type TLfTable;
	//typedef typename Fibre<TLfTable, FibrePrefixSumTable>::Type TPrefixSumTable;

	//typedef Index<TText, IndexEsa<> > TIndexEsa;

    {
        TIndex fmIndex;
	    SEQAN_ASSERT_LEQ(empty(fmIndex), true);
    }
    {
	    TText text;
	    generateText(text);

	    TIndex fmIndex(text);
        SEQAN_ASSERT_LEQ(empty(fmIndex), false);

        clear(fmIndex);
	    SEQAN_ASSERT_LEQ(empty(fmIndex), true);
    }
}

template <typename TText, typename TIndexSpec, typename TOptimization>
void fmIndexFindFirstIndex_(Index<TText, FMIndex<TIndexSpec, TOptimization> > /*tag*/)
{
	typedef Index<TText, FMIndex<TIndexSpec, TOptimization> > TIndex;
	//typedef typename Value<TText>::Type TChar;
	//typedef typename Fibre<TIndex, FibreLfTable>::Type TLfTable;
	//typedef typename Fibre<TLfTable, FibrePrefixSumTable>::Type TPrefixSumTable;

	typedef Index<TText, IndexEsa<> > TIndexEsa;

    TText text;
    generateText(text);

    TIndex fmIndex(text);
	Finder<TIndex> fmiFinder(fmIndex);

    TText pattern = infix(text, 5u, 10u);

    _findFirstIndex(fmiFinder, pattern, FinderFMIndex());

	TIndexEsa esaIndex(text);
	Finder<TIndexEsa> esaFinder(esaIndex);
	find(esaFinder, pattern);

    SEQAN_ASSERT_EQ(value(fmiFinder.range.i1), position(esaFinder)); 
}

template <typename TText, typename TIndexSpec, typename TOptimization>
void fmIndexGetFibre(Index<TText, FMIndex<TIndexSpec, TOptimization> > /*tag*/)
{
	typedef Index<TText, FMIndex<TIndexSpec, TOptimization> > TIndex;
	//typedef Index<TText, IndexEsa<> > TIndexEsa;

	TText text;
	generateText(text);

	TIndex fmIndex(text, 100000u);

    typename Fibre<TIndex, FibreLfTable>::Type & lfTable = getFibre(fmIndex, FibreLfTable());
    resize(getFibre(lfTable, FibrePrefixSumTable()), 1000u);
    SEQAN_ASSERT_EQ(length(getFibre(getFibre(fmIndex, FibreLfTable()), FibrePrefixSumTable())), 1000u);    
}

template <typename TText, typename TLfTable>
void lfTableLfMapping(StringSet<TText> /*tag*/, TLfTable &/*tag*/)
{
    typedef typename Value<TText>::Type TChar;

    {
        StringSet<TText> text;
        generateText(text);

        TLfTable lfTable;
        createLfTable(lfTable, text);

        for (unsigned i = 0; i < length(text); ++i)
        {
            unsigned pos = i;
            for (unsigned j = 0; j < length(text[length(text) - 1 - i]); ++j)
            {
                TChar character = getValue(getFibre(lfTable, FibreOccTable()), pos);
                SEQAN_ASSERT_EQ(character, text[length(text) - 1 - i][length(text[length(text) - 1 - i]) - 1 - j]);
                pos = lfMapping(lfTable, pos);
            }
        }
    }
}

template <typename TText, typename TLfTable>
void lfTableLfMapping(TText /*tag*/, TLfTable &/*tag*/)
{
    TText text;
    generateText(text);

    TLfTable lfTable;
    createLfTable(lfTable, text);

    unsigned pos = 0;
    for (unsigned i = 0; i < length(text); ++ i)
    {
        SEQAN_ASSERT_EQ(getValue(getFibre(lfTable, FibreOccTable()), pos), text[length(text) - 1 - i]);
        pos = lfMapping(lfTable, pos);
    }
}

template <typename TText, typename TIndexSpec, typename TOptimization>
void fmIndexSearch(Index<TText, FMIndex<TIndexSpec, TOptimization> > /*tag*/)
{
	typedef Index<TText, FMIndex<TIndexSpec, TOptimization> > TIndex;
	typedef Index<TText, IndexEsa<> > TIndexEsa;

	TText text;
	generateText(text);

	TIndex fmiIndex(text);
	Finder<TIndex> fmiFinder(fmiIndex);
	indexRequire(fmiIndex, FibreSaLfTable());

	TIndexEsa esaIndex(text);
	Finder<TIndexEsa> esaFinder(esaIndex);

	StringSet<String<typename Value<TIndex>::Type> > pattern;
	generatePattern(pattern, text);

    for(unsigned i = 36; i < length(pattern); ++i)
    {  
        clear(fmiFinder);
        clear(esaFinder);
        String<typename Value<TIndex>::Type> localPattern = pattern[i];

        while(find(fmiFinder, localPattern))
        {
            SEQAN_ASSERT(find(esaFinder, localPattern));
            SEQAN_ASSERT_EQ(position(fmiFinder), position(esaFinder));
        }
        SEQAN_ASSERT_NOT(find(esaFinder, localPattern));
    }
}

template <typename TText, typename TIndexSpec, typename TOptimization>
void fmIndexOpenSave(Index<TText, FMIndex<TIndexSpec, TOptimization> > /*tag*/)
{
	typedef Index<TText, FMIndex<TIndexSpec, TOptimization> > TIndex;
	//typedef Index<TText, IndexEsa<> > TIndexEsa;
	typedef typename Value<TIndex>::Type TAlphabet;
    typedef String<TAlphabet> TString;

	TText text;
	generateText(text, 1000);

	CharString tempFilename = SEQAN_TEMP_FILENAME();

    TIndex indexSave(text);
    indexCreate(indexSave);
    save(indexSave, toCString(tempFilename));

    TIndex indexOpen;
    open(indexOpen, toCString(tempFilename));

    Finder<TIndex> saveFinder(indexSave);
    Finder<TIndex> openFinder(indexOpen);

    TString pattern = "C";

    for(unsigned i = 0; i < length(pattern); ++i)
    {  
        clear(saveFinder);
        clear(openFinder);

        while(find(saveFinder, pattern))
        {
            SEQAN_ASSERT(find(openFinder, pattern));
            SEQAN_ASSERT_EQ(position(openFinder), position(saveFinder));
        }
        SEQAN_ASSERT_NOT(find(openFinder, pattern));
    }
}

// A test for strings.
SEQAN_DEFINE_TEST(test_fm_index_constructor)
{
    using namespace seqan;
    {
        Index<DnaString, FMIndex<> > dnaTag;
        Index<String<Dna5>, FMIndex<WT<>, void > > dna5Tag;
        Index<String<AminoAcid>, FMIndex<WT<>, void > > asTag;
        Index<String<signed char>, FMIndex<WT<>, void > > charTag;
        Index<String<unsigned char>, FMIndex<WT<>, void > > uCharTag;
        fmIndexConstructor(dnaTag);
        fmIndexConstructor(dna5Tag);
        fmIndexConstructor(asTag);
        fmIndexConstructor(uCharTag);
        fmIndexConstructor(charTag);
    }
    {
        Index<String<Dna>, FMIndex<SBM<>, void > > dnaTag;
        Index<String<Dna5>, FMIndex<SBM<>, void > > dna5Tag;
//         Index<String<AminoAcid>, FMIndex<SBM<>, void > > asTag;
//         Index<String<signed char>, FMIndex<SBM<>, void > > charTag;
//         Index<String<unsigned char>, FMIndex<SBM<>, void > > uCharTag;
        fmIndexConstructor(dnaTag);
        fmIndexConstructor(dna5Tag);
//         fmIndexConstructor(asTag);
//         fmIndexConstructor(uCharTag);
//         fmIndexConstructor(charTag);
    }
}

SEQAN_DEFINE_TEST(test_fm_index_clear)
{
    using namespace seqan;

    Index<DnaString, FMIndex<WT<>, void > > dnaTag;
    Index<String<Dna5>, FMIndex<WT<>, void > > dna5Tag;
    Index<String<AminoAcid>, FMIndex<WT<>, void > > asTag;
    Index<String<signed char>, FMIndex<WT<>, void > > charTag;
    Index<String<unsigned char>, FMIndex<WT<>, void > > uCharTag;
    fmIndexClear(dnaTag);
    fmIndexClear(dna5Tag);
    fmIndexClear(asTag);
    fmIndexClear(uCharTag);
    fmIndexClear(charTag);
}

SEQAN_DEFINE_TEST(test_fm_index_determine_sentinel_substitute_)
{
    using namespace seqan;

    Index<DnaString, FMIndex<WT<>, void > > dnaTag;
    Index<String<Dna5>, FMIndex<WT<>, void > > dna5Tag;
    Index<String<AminoAcid>, FMIndex<WT<>, void > > asTag;
    Index<String<signed char>, FMIndex<WT<>, void > > charTag;
    Index<String<unsigned char>, FMIndex<WT<>, void > > uCharTag;
    _fmIndexDetermineSentinelSubstitute(dnaTag);
    _fmIndexDetermineSentinelSubstitute(dna5Tag);
    _fmIndexDetermineSentinelSubstitute(asTag);
    _fmIndexDetermineSentinelSubstitute(uCharTag);
    _fmIndexDetermineSentinelSubstitute(charTag);
}

SEQAN_DEFINE_TEST(test_fm_index_empty)
{
    using namespace seqan;

    Index<DnaString, FMIndex<WT<>, void > > dnaTag;
    Index<String<Dna5>, FMIndex<WT<>, void > > dna5Tag;
    Index<String<AminoAcid>, FMIndex<WT<>, void > > asTag;
    Index<String<signed char>, FMIndex<WT<>, void > > charTag;
    Index<String<unsigned char>, FMIndex<WT<>, void > > uCharTag;
    fmIndexEmpty(dnaTag);
    fmIndexEmpty(dna5Tag);
    fmIndexEmpty(asTag);
    fmIndexEmpty(uCharTag);
    fmIndexEmpty(charTag);
}

SEQAN_DEFINE_TEST(test_fm_index_find_first_index_)
{
    using namespace seqan;

    Index<DnaString, FMIndex<WT<>, void > > dnaTag;
    Index<String<Dna5>, FMIndex<WT<>, void > > dna5Tag;
    Index<String<AminoAcid>, FMIndex<WT<>, void > > asTag;
    Index<String<signed char>, FMIndex<WT<>, void > > charTag;
    Index<String<unsigned char>, FMIndex<WT<>, void > > uCharTag;
    fmIndexFindFirstIndex_(dnaTag);
    fmIndexFindFirstIndex_(dna5Tag);
    fmIndexFindFirstIndex_(asTag);
    fmIndexFindFirstIndex_(uCharTag);
    fmIndexFindFirstIndex_(charTag);
}

SEQAN_DEFINE_TEST(test_fm_index_get_fibre)
{
    using namespace seqan;

    Index<DnaString, FMIndex<WT<>, void > > dnaTag;
    Index<String<Dna5>, FMIndex<WT<>, void > > dna5Tag;
    Index<String<AminoAcid>, FMIndex<WT<>, void > > asTag;
    Index<String<signed char>, FMIndex<WT<>, void > > charTag;
    Index<String<unsigned char>, FMIndex<WT<>, void > > uCharTag;    
    fmIndexGetFibre(dnaTag);
    fmIndexGetFibre(dna5Tag);
    fmIndexGetFibre(asTag);
    fmIndexGetFibre(uCharTag);
    fmIndexGetFibre(charTag);
}

// A test for strings.
SEQAN_DEFINE_TEST(test_fm_index_search)
{
    using namespace seqan;
    {
        Index<DnaString, FMIndex<WT<>, void > > dnaTag;
        Index<String<Dna5>, FMIndex<WT<>, void > > dna5Tag;
        Index<String<AminoAcid>, FMIndex<WT<>, void > > asTag;
        Index<String<signed char>, FMIndex<WT<>, void > > sCharTag;
        Index<String<unsigned char>, FMIndex<WT<>, void > > uCharTag;
        Index<String<char>, FMIndex<WT<>, void > > charTag;
        fmIndexSearch(dnaTag);
        fmIndexSearch(dna5Tag);
        fmIndexSearch(asTag);
        fmIndexSearch(uCharTag);
        fmIndexSearch(sCharTag);
        fmIndexSearch(charTag);
    }
    {
        Index<StringSet<DnaString>, FMIndex<WT<>, void > > dnaTag;
        Index<StringSet<Dna5String>, FMIndex<WT<>, void > > dna5Tag;
        Index<StringSet<String<AminoAcid> >, FMIndex<WT<>, void > > asTag;
        Index<StringSet<String<unsigned char> >, FMIndex<WT<>, void > > uCharTag;
        Index<StringSet<String<signed char> >, FMIndex<WT<>, void > > sCharTag;
        Index<StringSet<String<char> >, FMIndex<WT<>, void > > charTag;
        fmIndexSearch(dnaTag);
        fmIndexSearch(dna5Tag);
        fmIndexSearch(asTag);
        fmIndexSearch(uCharTag);
        fmIndexSearch(sCharTag);
        fmIndexSearch(charTag);
    }    
    {
        Index<DnaString, FMIndex<WT<>, CompressText> > dnaTag;
        Index<String<Dna5>, FMIndex<WT<>, CompressText > > dna5Tag;
        Index<String<AminoAcid>, FMIndex<WT<>,CompressText> > asTag;
        Index<String<signed char>, FMIndex<WT<>, CompressText > > sCharTag;
        Index<String<unsigned char>, FMIndex<WT<>, CompressText > > uCharTag;
        Index<String<char>, FMIndex<WT<>, CompressText > > charTag;
        fmIndexSearch(dnaTag);
        fmIndexSearch(dna5Tag);
        fmIndexSearch(asTag);
        fmIndexSearch(uCharTag);
        fmIndexSearch(sCharTag);
        fmIndexSearch(charTag); 
    }
    {
        Index<StringSet<DnaString>, FMIndex<WT<>, CompressText> > dnaTag;
        Index<StringSet<Dna5String>, FMIndex<WT<>, CompressText> > dna5Tag;
        Index<StringSet<String<AminoAcid> >, FMIndex<WT<>, CompressText> > asTag;
        Index<StringSet<String<unsigned char> >, FMIndex<WT<>, CompressText> > uCharTag;
        Index<StringSet<String<signed char> >, FMIndex<WT<>, CompressText> > sCharTag;
        Index<StringSet<String<char> >, FMIndex<WT<>, CompressText> > charTag;
        fmIndexSearch(dnaTag);
        fmIndexSearch(dna5Tag);
        fmIndexSearch(asTag);
        fmIndexSearch(uCharTag);
        fmIndexSearch(sCharTag);
        fmIndexSearch(charTag);
    }  
}

SEQAN_DEFINE_TEST(test_fm_index_open_save)
{
    using namespace seqan;
    {
        Index<DnaString, FMIndex<WT<>, void > > dnaTag;
        fmIndexOpenSave(dnaTag);
    }
    {
        Index<StringSet<DnaString>, FMIndex<WT<>, void > > dnaTag;
        fmIndexOpenSave(dnaTag);
    }    
}


SEQAN_DEFINE_TEST(test_lf_table_lf_mapping)
{
    using namespace seqan;

//     {
//         typedef Dna TChar;
//         typedef String<TChar> TString;
//         typedef WaveletTree<TString, FmiSentinelSubstituted<SingleSentinel<void> > > TOccTable;
//         typedef PrefixSumTable<TChar> TPrefixSumTable;
//         typedef LfTable<TOccTable, TPrefixSumTable> TLfTable;
// 
//         TLfTable tag;
//         lfTableLfMapping(TString(), tag);
//     }
//     {
//         typedef Dna TChar;
//         typedef StringSet<String<TChar> > TString;
//         typedef WaveletTree<String<TChar>, FmiSentinelSubstituted<MultiSentinel<void> > > TOccTable;
//         typedef PrefixSumTable<TChar> TPrefixSumTable;
//         typedef LfTable<TOccTable, TPrefixSumTable> TLfTable;
// 
//         TLfTable tag;
//         lfTableLfMapping(TString(), tag);
//     }
    
}



#endif  // TEST_FM_INDEX_BETA_H_
