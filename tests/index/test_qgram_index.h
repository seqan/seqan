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

#ifndef TESTS_INDEX_TEST_QGRAM_INDEX_H
#define TESTS_INDEX_TEST_QGRAM_INDEX_H


//////////////////////////////////////////////////////////////////////////////

namespace seqan
{

SEQAN_DEFINE_TEST(testGappedShapes)
{
    String<char> shape_string = "0010011011101";
    Shape<Dna,GenericShape> shape1;
    stringToShape(shape1, shape_string);
    Shape<Dna,GenericShape> shape2 = Shape<Dna,GenericShape>(shape1);

    SEQAN_ASSERT(shape1.weight == shape2.weight);
    SEQAN_ASSERT(shape1.span == shape2.span);
    SEQAN_ASSERT(shape1.diffs == shape2.diffs);
    SEQAN_ASSERT(length(shape1) == length(shape2));
    SEQAN_ASSERT(weight(shape1) == weight(shape2));
/*
    Shape<Dna,GenericShape> shape3 = Shape<Dna,GenericShape>(5, 13);
    shape3[0]=2;
    shape3[1]=2;
    shape3[2]=1;
    shape3[3]=1;
    shape3[4]=2;
    shape3[5]=1;
    shape3[6]=1;
    shape3[7]=2;
    for(int i = 0; i < 8; ++i)
        SEQAN_ASSERT(shape1[i] == shape3[i]);
*/
}


SEQAN_DEFINE_TEST(testUngappedShapes)
{
    Shape<Dna,SimpleShape> shape1;
    resize(shape1, 4);
    Shape<Dna,SimpleShape> shape2 = Shape<Dna,SimpleShape>(shape1);

    SEQAN_ASSERT(shape1.span == shape2.span);
    SEQAN_ASSERT(shape1.leftFactor == shape2.leftFactor);
    SEQAN_ASSERT(length(shape1) == length(shape2));
    SEQAN_ASSERT(weight(shape1) == weight(shape2));


    Shape<Dna,SimpleShape> shape3 = Shape<Dna,SimpleShape>(4);
    SEQAN_ASSERT(shape3.leftFactor == 64);
    SEQAN_ASSERT(shape1.leftFactor == shape3.leftFactor);


}

template <typename TIndex>
void testStepSize()
{
    TIndex index("CATGATTACATA");
    setStepSize(index, 2);
    hash(indexShape(index), "CAT");
    String<typename Position<DnaString>::Type> occs;
    occs = getOccurrences(index, indexShape(index));
    SEQAN_ASSERT_EQ(length(occs), 2u);
    SEQAN_ASSERT_EQ(occs[0], 0u);
    SEQAN_ASSERT_EQ(occs[1], 8u);
}

SEQAN_DEFINE_TEST(testStepSize)
{
    typedef Index<DnaString, IndexQGram< UngappedShape<3> > > TIndex1;
    typedef Index<DnaString, IndexQGram< UngappedShape<3>, OpenAddressing > > TIndex2;

    testStepSize<TIndex1>();
    testStepSize<TIndex2>();
}

/*
void testQGramIndexSchnell()
{
    clock_t start, finish;
    double duration;

    String<Dna> text;
    fstream strm_t;
    strm_t.open(TEST_PATH "fasta.txt", ios_base::in);
    read(strm_t, text, Fasta());
    String<Dna> next;
    resize(next,length(text));
    for(int i = 0; i < 1; ++i)                            // datei zu ende?
    {
        arrayCopyForward(begin(text),end(text),begin(next));
        append(text, next);
    }
    strm_t.close();

    String<char> shape_string = "1000100100001110011";
    int q = length(shape_string);
    Shape<Dna,GenericShape> shape;
    stringToShape(shape, shape_string);

    typedef Position<String<Dna> >::Type TPosition;
    String<TPosition> pos;
    resize(pos, length(text) - q + 2);

    String<TPosition> dir;
    int pos_size = _intPow((unsigned)ValueSize<Dna>::VALUE, weight(shape));
    pos_size += 1;
    resize(dir, pos_size);

    start = clock();
    Nothing nothing;
    createQGramIndex(pos, dir, nothing, text, shape, 1);
    finish = clock();
    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    //std::cout << "\nQGramIndex bauen dauert: " << duration << " Sekunden.\n\n";


}
*/
/*
void testGappedQGramIndex()
{
    String<Dna> text = "CTGAACCCTAAACCCT";
    String<char> shape_string = "101";
    int q = length(shape_string);
    Shape<Dna,GenericShape> shape;
    stringToShape(shape, shape_string);

    typedef Position<String<Dna> >::Type TPosition;
    String<TPosition> pos;
    resize(pos, length(text) - q + 2);

    String<TPosition> dir;
    int pos_size = _intPow((unsigned)ValueSize<Dna>::VALUE, weight(shape));
    pos_size += 1;
    resize(dir, pos_size);

    Nothing nothing;
    createQGramIndex(pos, dir, nothing, text, shape, 1);

    SEQAN_ASSERT(dir[0] == 0);
    SEQAN_ASSERT(dir[1] == 1);
    SEQAN_ASSERT(dir[2] == 5);
    SEQAN_ASSERT(dir[3] == 5);
    SEQAN_ASSERT(dir[4] == 5);
    SEQAN_ASSERT(dir[5] == 6);
    SEQAN_ASSERT(dir[6] == 8);
    SEQAN_ASSERT(dir[7] == 9);
    SEQAN_ASSERT(dir[8] == 11);
    SEQAN_ASSERT(dir[9] == 12);
    SEQAN_ASSERT(dir[10] == 12);
    SEQAN_ASSERT(dir[11] == 12);
    SEQAN_ASSERT(dir[12] == 12);
    SEQAN_ASSERT(dir[13] == 14);
    SEQAN_ASSERT(dir[14] == 14);
    SEQAN_ASSERT(dir[15] == 14);
    SEQAN_ASSERT(dir[16] == 14);

    SEQAN_ASSERT(pos[0] == 9);
    SEQAN_ASSERT(pos[1] == 11);
    SEQAN_ASSERT(pos[2] == 10);
    SEQAN_ASSERT(pos[3] == 4);
    SEQAN_ASSERT(pos[4] == 3);
    SEQAN_ASSERT(pos[5] == 7);
    SEQAN_ASSERT(pos[6] == 12);
    SEQAN_ASSERT(pos[7] == 5);
    SEQAN_ASSERT(pos[8] == 0);
    SEQAN_ASSERT(pos[9] == 13);
    SEQAN_ASSERT(pos[10] == 6);
    SEQAN_ASSERT(pos[11] == 2);
    SEQAN_ASSERT(pos[12] == 8);
    SEQAN_ASSERT(pos[13] == 1);

}
*/
SEQAN_DEFINE_TEST(testUngappedQGramIndex)
{
    typedef String<Dna> TString;
    typedef Shape<Dna, UngappedShape<2> > TShape;
    typedef Index<TString, IndexQGram<TShape> > TIndex;
    typedef Position<TString>::Type TPosition;

    TString text("CTGAACCCTAAACCCT");

    TIndex index(text);
    indexCreate(index, QGramSADir());

    String<TPosition> pos(getFibre(index, QGramSA()));
    String<TPosition> dir(getFibre(index, QGramDir()));

    SEQAN_ASSERT(dir[0] == 0);
    SEQAN_ASSERT(dir[1] == 3);
    SEQAN_ASSERT(dir[2] == 5);
    SEQAN_ASSERT(dir[3] == 5);
    SEQAN_ASSERT(dir[4] == 5);
    SEQAN_ASSERT(dir[5] == 5);
    SEQAN_ASSERT(dir[6] == 9);
    SEQAN_ASSERT(dir[7] == 9);
    SEQAN_ASSERT(dir[8] == 12);
    SEQAN_ASSERT(dir[9] == 13);
    SEQAN_ASSERT(dir[10] == 13);
    SEQAN_ASSERT(dir[11] == 13);
    SEQAN_ASSERT(dir[12] == 13);
    SEQAN_ASSERT(dir[13] == 14);
    SEQAN_ASSERT(dir[14] == 14);
    SEQAN_ASSERT(dir[15] == 15);

    SEQAN_ASSERT(pos[0] == 3);
    SEQAN_ASSERT(pos[1] == 9);
    SEQAN_ASSERT(pos[2] == 10);
    SEQAN_ASSERT(pos[3] == 4);
    SEQAN_ASSERT(pos[4] == 11);
    SEQAN_ASSERT(pos[5] == 5);
    SEQAN_ASSERT(pos[6] == 6);
    SEQAN_ASSERT(pos[7] == 12);
    SEQAN_ASSERT(pos[8] == 13);
    SEQAN_ASSERT(pos[9] == 0);
    SEQAN_ASSERT(pos[10] == 7);
    SEQAN_ASSERT(pos[11] == 14);
    SEQAN_ASSERT(pos[12] == 2);
    SEQAN_ASSERT(pos[13] == 8);
    SEQAN_ASSERT(pos[14] == 1);
}

inline bool
_qgramDisableBuckets(Index<StringSet<DnaString>, IndexQGram<Shape<Dna, UngappedShape<3> > > > &index)
{
    indexDir(index)[1] = -1;
    return true;
}

SEQAN_DEFINE_TEST(testUngappedQGramIndexMulti)
{
    typedef StringSet<DnaString>                    TStrings;
    //typedef SAValue<TStrings>::Type                 TSAValue;
    //typedef String<TSAValue>                        TPos;
    typedef Shape<Dna, UngappedShape<3> >           TShape;
    typedef Index<TStrings, IndexQGram<TShape> >    TIndex;

    TStrings strings;
    TIndex refIndex(strings);
    TIndex testIndex(strings);

                       //           111111
                       // 0123456789012345
    appendValue(strings, "CTGAACCCTAAACCCT");
//    appendValue(strings, "");                   // TODO: fix tupler to cope with strings smaller than q
    appendValue(strings, "GAAGGAGTGTGTGT");     //       requires to adapt pipe interface to return limitsString
//    appendValue(strings, "GT");
    appendValue(strings, "AAAACCCCAAACCCC");

    TShape &shape = indexShape(refIndex);
    indexCreate(refIndex, QGramSADir());
    resize(indexSA(refIndex), lengthSum(strings) - length(strings) * (length(shape) - 1));
    resize(indexDir(refIndex), _intPow((unsigned)ValueSize<Dna>::VALUE, length(shape)) + 1);
    createQGramIndex(refIndex);

    // test classic external index construction
    resize(indexSA(testIndex), length(indexSA(refIndex)));
    resize(indexDir(testIndex), length(indexDir(refIndex)));
    createQGramIndexExt(testIndex);

    for (unsigned i = 0; i < length(indexDir(refIndex)); ++i)
        SEQAN_ASSERT_EQ_MSG(dirAt(i, refIndex), dirAt(i, testIndex), "i is %d", i);
    for (unsigned i = 0; i < length(indexSA(refIndex)); ++i)
        SEQAN_ASSERT_EQ_MSG(saAt(i, refIndex), saAt(i, testIndex), "i is %d", i);

    clear(indexSA(testIndex));
    clear(indexDir(testIndex));

    // test new external index construction
    resize(indexSA(testIndex), length(indexSA(refIndex)));
    resize(indexDir(testIndex), length(indexDir(refIndex)));
    createQGramIndexExtSA(testIndex);

    for (unsigned i = 0; i < length(indexDir(refIndex)); ++i)
        SEQAN_ASSERT_EQ_MSG(dirAt(i, refIndex), dirAt(i, testIndex), "i is %d", i);
    for (unsigned i = 0; i < length(indexSA(refIndex)); ++i)
        SEQAN_ASSERT_EQ_MSG(saAt(i, refIndex), saAt(i, testIndex), "i is %d", i);
}


//////////////////////////////////////////////////////////////////////////////

SEQAN_DEFINE_TEST(testQGramFind)
{
    typedef Index<String<char>, IndexQGram<UngappedShape<2> > > TQGramIndex;
    TQGramIndex idx("to be or not to be");
    Finder<TQGramIndex> finder(idx);

    SEQAN_ASSERT(find(finder, "be"));
    SEQAN_ASSERT(position(finder) == 3);
    SEQAN_ASSERT(find(finder, "be"));
    SEQAN_ASSERT(position(finder) == 16);
    SEQAN_ASSERT(!find(finder, "be"));
/*
    while (find(finder, "be"))
    {
        std::cout << position(finder) << "\n";
    }
*/
}

//////////////////////////////////////////////////////////////////////////////


} //namespace seqan

#endif //#ifndef SEQAN_HEADER_...
