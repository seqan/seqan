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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================

#ifndef TESTS_INDEX_TEST_STREE_ITERATORS_H
#define TESTS_INDEX_TEST_STREE_ITERATORS_H


//////////////////////////////////////////////////////////////////////////////

namespace seqan
{

SEQAN_DEFINE_TEST(testEmptyIndex)
{
    typedef String<Dna5> TString;
    typedef Index<TString> TIndex;

    TIndex index;

    Iterator<TIndex, TopDown<> >::Type iterator(index);

}

SEQAN_DEFINE_TEST(testIssue509)
{
    // goDown causes assertion on small index #509
    // Originally by John on Trac
    typedef String<Dna5> TString;
    typedef StringSet<TString> TStrings;
    typedef Index<TStrings> TIndex;

    TStrings seqs;
    TString seq("A");
    appendValue(seqs, seq);
    TIndex index(seqs);

    Iterator<TIndex, TopDown<> >::Type iterator(index);
    SEQAN_ASSERT(isRoot(iterator));
    SEQAN_ASSERT_NOT(isLeaf(iterator));

    SEQAN_ASSERT(goDown(iterator));
    SEQAN_ASSERT_NOT(isRoot(iterator));
    SEQAN_ASSERT(isLeaf(iterator));

    goRoot(iterator);
    SEQAN_ASSERT(isRoot(iterator));
    SEQAN_ASSERT_NOT(isLeaf(iterator));
    SEQAN_ASSERT(goDown(iterator, seq));
}

SEQAN_DEFINE_TEST(testIssue509b)
{
    // goDown causes assertion on small index #509
    // Originally by John on Trac
    typedef String<Dna5> TString;
    typedef Index<TString> TIndex;

    TString seq("A");
    TIndex index(seq);

    Iterator<TIndex, TopDown<> >::Type iterator(index);
    SEQAN_ASSERT(isRoot(iterator));
    SEQAN_ASSERT_NOT(isLeaf(iterator));

    SEQAN_ASSERT(goDown(iterator));
    SEQAN_ASSERT_NOT(isRoot(iterator));
    SEQAN_ASSERT(isLeaf(iterator));

    goRoot(iterator);
    SEQAN_ASSERT(isRoot(iterator));
    SEQAN_ASSERT_NOT(isLeaf(iterator));
    SEQAN_ASSERT(goDown(iterator, seq));
}

SEQAN_DEFINE_TEST(goDownOnEmptyString)
{
    typedef Index<DnaString, FMIndex<> > TIndex;
    typedef Iter<TIndex, VSTree<TopDown<> > > TIter;

    DnaString text("GCCTC");
    TIndex index(text);
    TIter it(index);

    goDown(it, "C");
    goDown(it);
    goDown(it, "");
    goRight(it);
    SEQAN_ASSERT_EQ(representative(it), DnaString("GC"));
    goRight(it);
    SEQAN_ASSERT_EQ(representative(it), DnaString("TC"));
}

SEQAN_DEFINE_TEST(testBuild)
{
        typedef String<char> TText;
        typedef StringSet< TText, Owner<ConcatDirect<> > > TMulti;

        String<char> gen1, gen2, gen3;
        std::cout << open(gen1, "corpus/NC_000117.txt");
        std::cout << open(gen2, "corpus/NC_002620.txt");
        std::cout << open(gen3, "corpus/NC_007429.txt");

        Index<TMulti> esa;
        appendValue(indexText(esa), gen1);
        appendValue(indexText(esa), gen2);
        appendValue(indexText(esa), gen3);


        indexRequire(esa, EsaSA());
        indexRequire(esa, EsaLcp());
        indexRequire(esa, EsaBwt());
        indexRequire(esa, EsaChildtab());

        save(esa, "corpus/chlamydia");
}

template <typename TIter>
inline void _printNode(TIter const it)
{
        std::cout << countOccurrences(it) << "\t";
        std::cout << representative(it) << "\t";
//        std::cout << "parentEdgeLabel:" << parentEdgeLabel(it);
        std::cout << std::endl;
}

SEQAN_DEFINE_TEST(testMultiIndex)
{
        typedef String<Dna5> TText;
        typedef StringSet< TText, Owner<> > TMulti;

        String<Dna5> t[6];
        //t[0] = "caterpillar";
        //t[1] = "catwoman";
        //t[2] = "pillow";
        //t[3] = "willow";
        //t[4] = "ill";
        //t[5] = "wow";

        // test empty stringsets
        {
            typedef Index<TMulti> TIndex;
            TMulti emptyStringSet;

            TIndex index(emptyStringSet);
            Iterator<
                TIndex,
                TopDown<>
            >::Type iterator(index);
        }

        t[0] = "caggctcgcgt";
        t[1] = "caggaacg";
        t[2] = "tcgttg";
        t[3] = "tggtcg";
        t[4] = "agg";
        t[5] = "ctg";

        t[0] = "ac";
        t[1] = "ac";
//        t[2] = "aatt";

        Index<TMulti> esa;
        for(unsigned i=0; i<2; ++i)
            appendValue(indexText(esa), t[i]);

        // efficient dfs iterator (hiding edges with empty labels)
        {
            std::cout << "BottomUp without empty edges" << std::endl;
            Iter<Index<TMulti>, VSTree< BottomUp<> > > it(esa);
            while (!atEnd(it)) {
                _printNode(it);
                goNext(it);
            }
        }

        // efficient dfs iterator
        {
            std::cout << std::endl << "BottomUp with empty edges" << std::endl;
            Iter<Index<TMulti>, VSTree< BottomUp<PostorderEmptyEdges> > > it(esa);
            while (!atEnd(it)) {
                _printNode(it);
                goNext(it);
            }
        }

        // topdown dfs iterator (hiding edges with empty labels)
        {
            std::cout << std::endl << "TopDown postorder without empty edges" << std::endl;
            Iter<Index<TMulti>, VSTree< TopDown<ParentLinks<Postorder> > > > it(esa);
            while (goDown(it)) ;
            while (!atEnd(it)) {
                _printNode(it);
                goNext(it);
            }
        }

        // topdown dfs iterator
        {
            std::cout << std::endl << "TopDown postorder with empty edges" << std::endl;
            Iter<Index<TMulti>, VSTree< TopDown<ParentLinks<PostorderEmptyEdges> > > > it(esa);
            while (goDown(it)) ;
            while (!atEnd(it)) {
                _printNode(it);
                goNext(it);
            }
        }

        // topdown dfs iterator (hiding edges with empty labels)
        {
            std::cout << std::endl << "TopDown preorder without empty edges" << std::endl;
            Iter<Index<TMulti>, VSTree< TopDown<ParentLinks<Preorder> > > > it(esa);
            goDown(it,'c');
            while (!atEnd(it)) {
                _printNode(it);
                goNext(it);
            }
        }

        // topdown dfs iterator
        {
            std::cout << std::endl << "TopDown preorder with empty edges" << std::endl;
            Iter<Index<TMulti>, VSTree< TopDown<ParentLinks<PreorderEmptyEdges> > > > it(esa);
            while (!atEnd(it)) {
                _printNode(it);
                goNext(it);
            }
        }

        // topdown iterator w/o parent links (hiding edges with empty labels)
        {
            std::cout << std::endl << "TopDown with empty edges" << std::endl;
            Iter<Index<TMulti>, VSTree< TopDown<HideEmptyEdges> > > it(esa);
            _printNode(it);
            while (goDown(it))
                _printNode(it);
        }

        // topdown iterator w/o parent links
        {
            std::cout << std::endl << "TopDown with empty edges" << std::endl;
            Iter<Index<TMulti>, VSTree< TopDown<EmptyEdges> > > it(esa);
            _printNode(it);
            while (goDown(it))
                _printNode(it);
        }

//        indexRequire(esa, EsaSA());
//        indexRequire(esa, EsaBwt());
//        for(int i=0; i<length(indexRawSA(esa)); ++i)
//            std::cout << saAt(i,esa) << " " << bwtAt(i,esa) << "    " << suffix(t[getValueI1(saAt(i,esa))], getValueI2(saAt(i,esa))) << std::endl;
//
////        resize(indexLcp(esa), length(indexRawText(esa)));
////        createLcpTableExt(indexLcp(esa), indexText(esa), indexSA(esa), Kasai());
//        indexRequire(esa, EsaLcp());
//        for(int i=0; i<length(indexRawSA(esa)); ++i)
//            std::cout << lcpAt(i,esa) << "    " << suffix(t[getValueI1(saAt(i,esa))], getValueI2(saAt(i,esa))) << std::endl;
//
//        for(int i=0; i<length(indexRawSA(esa)); ++i)
//            std::cout << saAt(i,esa) << " = " << indexRawSA(esa)[i] << "    " << std::endl;
//        for(int i=0; i<length(indexRawSA(esa)); ++i)
//            std::cout << bwtAt(i,esa) << " = " << indexBwt(esa).tab[i] << "    " << std::endl;
//        for(int i=0; i<length(indexRawSA(esa)); ++i)
//            std::cout << lcpAt(i,esa) << " = " << indexLcp(esa)[i] << "    " << std::endl;


//        resize(sa, length(indexRawText(esa)));
//        createSuffixArrayExt(sa, indexText(esa), Skew7());
//
//        for(int i=0; i<length(indexRawText(esa)); ++i)
//            std::cout << indexRawText(esa)[i] << "    ";
//
//        String<unsigned> lcp;
//        resize(lcp, length(indexRawText(esa)));
//        createLcpTableExt(lcp, indexText(esa), sa, Kasai());
}

template <typename TIndex1, typename TIndex2>
void compareTreeIterators(TIndex1 &index1, TIndex2 &index2)
{
    Iter<TIndex1, VSTree< TopDown< ParentLinks<Preorder> > > > it1(index1);
    Iter<TIndex2, VSTree< TopDown< ParentLinks<Preorder> > > > it2(index2);

    while (!atEnd(it1) && !atEnd(it2))
    {
        SEQAN_ASSERT_EQ(representative(it1), representative(it2));
        SEQAN_ASSERT_EQ(parentEdgeLabel(it1), parentEdgeLabel(it2));
        SEQAN_ASSERT_EQ(countOccurrences(it1), countOccurrences(it2));
        SEQAN_ASSERT_EQ(isRoot(it1), isRoot(it2));
//        SEQAN_ASSERT_EQ(isLeaf(it1), isLeaf(it2));
        goNext(it1);
        goNext(it2);
    }

    SEQAN_ASSERT_EQ(atEnd(it1), atEnd(it1));
}

template <typename TIndexSpec1, typename TIndexSpec2>
void compareIndices()
{
    {
        CharString text("mississippi");
        Index<CharString, TIndexSpec1> index1(text);
        Index<CharString, TIndexSpec2> index2(text);
        compareTreeIterators(index1, index2);
    }
    {
        DnaString text("acaaacatat");
        Index<DnaString, TIndexSpec1> index1(text);
        Index<DnaString, TIndexSpec2> index2(text);
        compareTreeIterators(index1, index2);
    }
    {
        StringSet<CharString> t;
        resize(t, 6);
        t[0] = "caterpillar";
        t[1] = "catwoman";
        t[2] = "pillow";
        t[3] = "willow";
        t[4] = "ill";
        t[5] = "wow";
        Index<StringSet<CharString>, TIndexSpec1> index1(t);
        Index<StringSet<CharString>, TIndexSpec2> index2(t);
        compareTreeIterators(index1, index2);
    }
    {
        StringSet<DnaString> t;
        resize(t, 6);
        t[0] = "caggctcgcgt";
        t[1] = "caggaacg";
        t[2] = "tcgttg";
        t[3] = "tggtcg";
        t[4] = "agg";
        t[5] = "ctg";
        Index<StringSet<DnaString>, TIndexSpec1> index1(t);
        Index<StringSet<DnaString>, TIndexSpec2> index2(t);
        compareTreeIterators(index1, index2);
    }
}

SEQAN_DEFINE_TEST(testCompareIndices_Esa_Wotd)
{
    compareIndices<IndexEsa<>, IndexWotd<> >();
}


template <typename TIndexSpec>
void testSTreeIterators()
{
    typedef Index<String<char>, TIndexSpec> TIndex;
    typedef typename Iterator<TIndex, TopDown<> >::Type TIterator;
    typedef typename Iterator<TIndex, TopDown<ParentLinks<> > >::Type TParentLinkIterator;

    // test empty trees
    {
        TIndex index("");
        TIterator iter(index);
        TParentLinkIterator piter(index);
        SEQAN_ASSERT_NOT(goDown(iter));
        SEQAN_ASSERT_NOT(goRight(iter));
        SEQAN_ASSERT_NOT(goUp(piter));
    }

    {
        String<char> text("acaaacatatz");
//        String<char> text("AAAAAGGGGG");
        TIndex index(text);
        Iter<TIndex, VSTree< TopDown< ParentLinks<Preorder> > > > it(index);
        Iter<TIndex, VSTree< TopDown<> > > itNoLinks(it);    // test conversion
        //Iter<TIndex, VSTree< BottomUp<> > > it(index);

//        while (goDown(it));
        while (!atEnd(it)) {
//            std::cout << countOccurrences(it) << "\t";
            std::cout << representative(it) << "\t";
            std::cout << "parentEdgeLabel: " << parentEdgeLabel(it); // << " " << value(it).node << "  " << value(it).range;
            std::cout << std::endl;
            goNext(it);
        }
            std::cout << std::endl;
        _dump(index);
//        goBegin(it);
//        while (!atEnd(it)) {
//            std::cout << countOccurrences(it) << "\t";
//            std::cout << representative(it) << "\t";
//            std::cout << "parentEdgeLabel: " << parentEdgeLabel(it) << " " << value(it).node << "  " << value(it).range;
//            std::cout << std::endl;
//            goNext(it);
//        }
//        _dump(index);
    }
}

SEQAN_DEFINE_TEST(testSTreeIterators_Wotd)
{
    testSTreeIterators<IndexWotd<> >();
}

SEQAN_DEFINE_TEST(testSTreeIterators_WotdOriginal)
{
    testSTreeIterators<IndexWotd<WotdOriginal> >();
}

SEQAN_DEFINE_TEST(testSTreeIterators_Esa)
{
    testSTreeIterators<IndexEsa<> >();
}


template <typename TPair>
struct PairLess_ :
        public std::binary_function<TPair, TPair, bool>
{
        inline bool
        operator() (TPair const& a1, TPair const& a2) const
        {
                if (a1.i1 == a2.i1) return (a1.i2 < a2.i2);
                else return (a1.i1 < a2.i1);
        }
};

SEQAN_DEFINE_TEST(testMaxRepeats)
{
//        typedef String<char, External<> > TText;
        typedef String<char> TText;

        Index<TText> esa;
//        open(esa, "corpus/NC_000117.txt");
        //                01234567890123456789
        indexText(esa) = "HALLOBALLOHALLEBALLO";

//        FILE* dotFile = fopen("stree.dot","w");
//        writeRecords(dotFile, esa, DotDrawing());
//        fclose(dotFile);

        Iterator< Index<TText>, MaxRepeats >::Type it(esa, 3);
        typedef MaxRepeat< Index<TText> > TRepeat;
        typedef Value<TRepeat>::Type TPair;

        int found = 0;
        while (!atEnd(it))
        {
//            std::cout << representative(it) << ":";
            Iterator<TRepeat, MaxRepeatOccurrences>::Type mit(it);
            String<TPair> occs;
            while (!atEnd(mit)) {
//                std::cout << "\t" << *mit << std::flush;
                appendValue(occs, *mit);
                if (back(occs).i1 > back(occs).i2)
                {
                    Value<TPair, 1>::Type tmp = back(occs).i1;
                    back(occs).i1 = back(occs).i2;
                    back(occs).i2 = tmp;
                }
                ++mit;
            }

            std::sort(begin(occs, Standard()), end(occs, Standard()), PairLess_<TPair>());
            if (representative(it) == "ALL")
            {
                SEQAN_ASSERT_EQ(length(occs), 2u);
                SEQAN_ASSERT_EQ(occs[0], TPair(6,11));
                SEQAN_ASSERT_EQ(occs[1], TPair(11,16));
            } else
            if (representative(it) == "HALL")
            {
                SEQAN_ASSERT_EQ(length(occs), 1u);
                SEQAN_ASSERT_EQ(occs[0], TPair(0,10));
            } else
            if (representative(it) == "BALLO")
            {
                SEQAN_ASSERT_EQ(length(occs), 1u);
                SEQAN_ASSERT_EQ(occs[0], TPair(5,15));
            } else
            if (representative(it) == "ALLO")
            {
                SEQAN_ASSERT_EQ(length(occs), 1u);
                if (occs[0].i2 == 6)
                    SEQAN_ASSERT_EQ(occs[0], TPair(1,6));
                else
                    SEQAN_ASSERT_EQ(occs[0], TPair(1,16));
            } else
            {
                SEQAN_ASSERT_FAIL("Unknown maximal repeat found!");
            }

            ++it;
            ++found;
        }
        SEQAN_ASSERT_EQ(found, 5);
}


SEQAN_DEFINE_TEST(testMultiMEMs)
{
        typedef String<char> TText;
        typedef StringSet< TText, Owner<ConcatDirect<> > > TMulti;

        Index<TMulti> esa;

        String<Dna5> t[6];
        t[0] = "caterpillar";
        t[1] = "catwoman";
        t[2] = "pillow";
        t[3] = "willow";
        t[4] = "ill";
        t[5] = "wow";

        std::ofstream dotFile("stree.dot");
        writeRecords(dotFile, esa, DotDrawing());
        dotFile.close();

        Iterator< Index<TMulti>, MultiMems >::Type it(esa, 3);
        typedef MultiMem< Index<TMulti> > TMultiMEM;
        while (!atEnd(it)) {
            std::cout << representative(it) << ":";
            Iterator<TMultiMEM>::Type mit(it);
            while (!atEnd(mit)) {
                std::cout << "\t" << *mit;
                ++mit;
            }
            std::cout << std::endl;
            ++it;
        }
}

template <typename TIteratorSpec>
void _testSuperMaxRepeats()
{
//        typedef String<char, External<> > TText;
        typedef String<char> TText;

        Index<TText> esa;
//        open(esa, "corpus/NC_000117.txt");
        indexText(esa) = "HALLOBALLOHALLEBALLO";

        typedef Index<TText> TIndex;
        typedef SAValue<TIndex>::Type TSAValue;
        typename Iterator<TIndex, TIteratorSpec >::Type it(esa);

        int found = 0;
        while (!atEnd(it))
        {
            String<TSAValue> occs = getOccurrences(it);
            std::sort(begin(occs, Standard()), end(occs, Standard()));

//            std::cout << representative(it) << ":";
//            for(typename Size<Index<TText> >::Type i = 0; i < countOccurrences(it); ++i)
//                std::cout << "\t" << getOccurrences(it)[i];
//            std::cout << std::endl;

            if (representative(it) == "BALLO")
            {
                SEQAN_ASSERT_EQ(length(occs), 2u);
                SEQAN_ASSERT_EQ(occs[0], 5u);
                SEQAN_ASSERT_EQ(occs[1], 15u);
            } else
            if (representative(it) == "HALL")
            {
                SEQAN_ASSERT_EQ(length(occs), 2u);
                SEQAN_ASSERT_EQ(occs[0], 0u);
                SEQAN_ASSERT_EQ(occs[1], 10u);
            } else
            {
                SEQAN_ASSERT_FAIL("Unknown supermaximal repeat found!");
            }

            ++it;
            ++found;
        }
        SEQAN_ASSERT_EQ(found, 2);
}

SEQAN_DEFINE_TEST(testSuperMaxRepeats)
{
    _testSuperMaxRepeats<SuperMaxRepeats>();
}

SEQAN_DEFINE_TEST(testSuperMaxRepeatsFast)
{
    _testSuperMaxRepeats<SuperMaxRepeatsFast>();
}


SEQAN_DEFINE_TEST(testMUMs)
{
        typedef String<char> TText;
        typedef StringSet< TText, Owner<ConcatDirect<> > > TMulti;
        typedef Index<TMulti, IndexEsa<> > TIndex;

        String<char> t[3];

        t[0] = "fefhalloballo";
        t[1] = "halloballefser";
        t[2] = "grballoballo";

        TIndex esa;
        for(int i = 0; i < 3; ++i)
            appendValue(indexText(esa), t[i]);            // add sequences to multiple index

        Iterator<TIndex, Mums>::Type  it(esa, 3);        // set minimum MUM length to 3
        typedef SAValue<TIndex>::Type TPair;
        String<TPair> occs;                                // temp. string storing the hit positions

        int found = 0;
//        std::cout << std::resetiosflags(std::ios::left);
        while (!atEnd(it))
        {
            occs = getOccurrences(it);                    // gives hit positions (seqNo,seqOfs)
            orderOccurrences(occs);                        // order them by seqNo

//            std::cout << representative(it) << ":";
//            for(unsigned i = 0; i < length(occs); ++i)
//                std::cout << "\t" << getValueI2(occs[i]);
//            std::cout << std::endl;
//            std::cout << alignment(it) << std::endl;

            if (representative(it) == "alloball")
            {
                SEQAN_ASSERT_EQ(length(occs), 3u);
                SEQAN_ASSERT_EQ(occs[0], TPair(0,4));
                SEQAN_ASSERT_EQ(occs[1], TPair(1,1));
                SEQAN_ASSERT_EQ(occs[2], TPair(2,3));
            } else {
                SEQAN_ASSERT_FAIL("Unknown MUM found!");
            }


            ++it;
            ++found;
        }
        SEQAN_ASSERT_EQ(found, 1);
}


template <typename TAlgorithmSpec>
void testFind()
{
        String<unsigned int> pos;

    //____________________________________________________________________________
    // Test1 - small needle

        String<char> haystack("Dies ist ein Haystack. Ja, das ist wirklich einer!");
        Index<String<char> > index(haystack);

        Finder<Index< String<char> >, TAlgorithmSpec> finder(index);

        String<char> needle1("ist");
        seqan::Pattern<String<char> > pattern(needle1);

        while (find(finder, pattern))
            appendValue(pos,position(finder));

        SEQAN_ASSERT_EQ(length(pos), 2u);
        SEQAN_ASSERT(pos[0] == 5u);
        SEQAN_ASSERT(pos[1] == 31u);

    //____________________________________________________________________________
    // Test2 - large needle

        haystack = "abcdefghijklmnopqrstuvwxyzabcdefghijklmnopqrstuvwxyzabcdefgaabcdef";
        clear(index);
        clear(finder);

        needle1 = "abcdefghijklmnopqrstuvwxyzabcdefg";
        setNeedle(pattern, needle1);

        clear(pos);
        while (find(finder, pattern))
            appendValue(pos,position(finder));

        SEQAN_ASSERT_EQ(length(pos), 2u);
        SEQAN_ASSERT(pos[1] == 0u);
        SEQAN_ASSERT(pos[0] == 26u);
}

SEQAN_DEFINE_TEST(testFind_Esa_Mlr)
{
    testFind<FinderMlr>();
}

SEQAN_DEFINE_TEST(testMultipleStrings_Ticket1109)
{
    StringSet<String<char> > text;
    appendValue(text, "How many");
    appendValue(text, " wood would");
    appendValue(text, " a woodchuck chuck?");

    typedef Index<StringSet<CharString> > TIndex;
    TIndex index(text);

    Iterator< TIndex, TopDown<> >::Type it(index);
}

SEQAN_DEFINE_TEST(testTrieIterator)
{
    typedef StringSet<CharString>               TText;
    typedef Index<TText, IndexSa<> >            TIndex;
    typedef Iterator<TIndex, TopDown<> >::Type  TIter;

    TText text;
    appendValue(text, "bananamama");
    appendValue(text, "bananajoe");

    TIndex index(text);
    indexCreate(index, FibreSA(), Trie());
    TIter it(index);
    goDown(it);
    SEQAN_ASSERT_EQ(parentEdgeLabel(it), 'b');
    SEQAN_ASSERT_EQ(countOccurrences(it), 2u);
    goDown(it, "anana");
    goDown(it);
    SEQAN_ASSERT_EQ(parentEdgeLabel(it), 'j');
}

SEQAN_DEFINE_TEST(testRadixTreeIterator)
{
    typedef StringSet<CharString>               TText;
    typedef Index<TText, IndexWotd<> >          TIndex;
    typedef Iterator<TIndex, TopDown<> >::Type  TIter;

    TText text;
    appendValue(text, "bananamama");
    appendValue(text, "bananajoe");

    TIndex index(text);
    indexCreate(index, WotdDir(), Trie());
    TIter it(index);
    goDown(it);
    SEQAN_ASSERT_EQ(parentEdgeLabel(it), "banana");
    SEQAN_ASSERT_EQ(countOccurrences(it), 2u);
    goDown(it);
    SEQAN_ASSERT_EQ(parentEdgeLabel(it), "joe");
}

//////////////////////////////////////////////////////////////////////////////


} //namespace seqan

#endif //#ifndef SEQAN_HEADER_...
