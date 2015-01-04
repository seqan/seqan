// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================

#ifndef TESTS_INDEX_TEST_INDEX_CREATION_H
#define TESTS_INDEX_TEST_INDEX_CREATION_H

#include <seqan/random.h>

#define SEQAN_PROFILE
//#define SEQAN_DEBUG
//#define SEQAN_DEBUG_INDEX

//#define SEQAN_TEST
//#define SEQAN_TEST_SKEW3
//#define SEQAN_TEST_SKEW7


//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN
{

SEQAN_DEFINE_TEST(testIndexModifiedStringReverseEsa)
{
    typedef String<AminoAcid> TString;
    typedef ModifiedString<TString, ModReverse> TReverse;
    typedef Index<TReverse> TIndex;

    TString org("ACGXN");
    TReverse rev(org);

    TIndex index(rev);
    Iterator<TIndex, TopDown<> >::Type iter(index);
}

SEQAN_DEFINE_TEST(testIndexModifiedStringReverseFM)
{
    typedef String<AminoAcid> TString;
    typedef ModifiedString<TString, ModReverse> TReverse;
    typedef Index<TReverse, FMIndex<> > TIndex;

    TString org("ACGXN");
    TReverse rev(org);

    TIndex index(rev);
    Iterator<TIndex, TopDown<> >::Type iter(index);
}

SEQAN_DEFINE_TEST(testIndexModifiedStringViewEsa)
{
    typedef String<AminoAcid> TString;
    typedef ModifiedString<TString, FunctorConvert<AminoAcid, char> > TConvert;
    typedef Index<TConvert> TIndex;

    TString org("ACGXN");
    TConvert conv(org);

    TIndex index(conv);
    Iterator<TIndex, TopDown<> >::Type iter(index);
}

SEQAN_DEFINE_TEST(testIndexModifiedStringViewFM)
{
    typedef String<AminoAcid> TString;
    typedef ModifiedString<TString, FunctorConvert<AminoAcid, char> > TConvert;
    typedef Index<TConvert, FMIndex<> > TIndex;

    TString org("ACGXN");
    TConvert conv(org);

    TIndex index(conv);
    Iterator<TIndex, TopDown<> >::Type iter(index);
}

SEQAN_DEFINE_TEST(testIssue519)
{
    // Originally from Sascha on Trac
    // Bug in SAQSort: For StringSets, the sorting method gives a wrong order.
    CharString text = "bananamama";
    CharString text2 = "bananajoe";
    CharString text3 = "joesmama";

    StringSet<CharString> strSet;
    appendValue(strSet, text); appendValue(strSet, text2); appendValue(strSet, text3);
    Index<StringSet<CharString>, IndexEsa<> > index1(strSet);
    Index<StringSet<CharString>, IndexEsa<> > index2(strSet);

    indexCreate(index1, EsaSA(), Skew7());
    indexCreate(index2, EsaSA(), SAQSort());

    SEQAN_ASSERT_EQ(indexSA(index1), indexSA(index2));

//    Iterator<String<SAValue<StringSet<CharString> >::Type> >::Type iterSet = begin(indexSA(index2));
//    std::cout << "Suffix Array: " << std::endl;
//    for(; iterSet != end(indexSA(index2)); ++iterSet)
//        std::cout << getSeqNo(*iterSet) << "," << getSeqOffset(*iterSet) << "\t"
//                  << suffix(getValue(strSet, getSeqNo(*iterSet)), getSeqOffset(*iterSet)) << std::endl;
}

SEQAN_DEFINE_TEST(testIndexCreation)
{
    Rng<> rng(/*seed=*/1);

        typedef String<char> TText;
        typedef String<unsigned> TArray;

        TText   text;
        TArray  sa;
        TArray  isa;
        TArray  lcp;
        TArray  child, childExt;
        TText   bwt;

        const int runs = 2;                     // conduct 10 test runs
        const int maxSize = 20 * 1024 * 1024;	// max text size is 20 megabyte
        bool result = true;
        (void)result;  // Is never read...

        _proFloat timeDelta[14];
        _proFloat timeSum[14];
        for(int i = 0; i < 10; ++i)
            timeSum[i] = 0;
        __int64 textSum = 0;

        static const char* algNames[] = {
            "Skew3         ",
            "Skew7         ",
            "ManberMyers   ",
            "LarssonSadake ",
            "SAQSort       ",
            "Skew3Ext      ",
            "Skew7Ext      ",
            "InvSA         ",
            "InvSA Parallel",
            "Kasai         ",
            "KasaiInPlace  ",
            "KasaiExt      ",
            "Childtab      ",
            "ChildTabExt   "
        };

        int TI;
        for(int i = 0; i < runs; ++i) {

            std::cout << "*** RUN " << i << " ***";

            Pdf<Uniform<int> > pdf(0, maxSize);
            int size = pickRandomNumber(rng, pdf);
            TI = 0;

//___randomize_text___________________________________________________________

            resize(text,size);
/*          if (i < runs/2)
                randomize(text);
            else
*/          textRandomize(text);
/*          String<char,External<> > errorText;	// read in text causing an error
            open(errorText,"error.txt");
            text = errorText;
*/
/*          text = "MISSISSIPPI";
            size = length(text);
            std::cout << "text created (n=" << size << ")" << std::endl;
*/
            std::cout << "   textSize: " << length(text) << std::endl;

//___create_suffix_array______________________________________________________

            resize(sa, size);
            timeDelta[TI] = -SEQAN_PROGETTIME;
            createSuffixArray(sa, text, Skew3());
            timeDelta[TI++] += SEQAN_PROGETTIME;
            if (!isSuffixArray(sa, text)) {
                std::cout << "suffix array creation (internal Skew3) failed" << std::endl;
                result = false;
            }
            std::cout << "."; std::cout.flush();

            blank(sa);
            timeDelta[TI] = -SEQAN_PROGETTIME;
            createSuffixArray(sa, text, Skew7());
            timeDelta[TI++] += SEQAN_PROGETTIME;
            if (!isSuffixArray(sa, text)) {
                std::cout << "suffix array creation (internal Skew7) failed" << std::endl;
                result = false;
            }
            std::cout << "."; std::cout.flush();

            blank(sa);
            timeDelta[TI] = -SEQAN_PROGETTIME;
            createSuffixArray(sa, text, ManberMyers());
            timeDelta[TI++] += SEQAN_PROGETTIME;
            if (!isSuffixArray(sa, text)) {
                std::cout << "suffix array creation (internal ManberMyers) failed" << std::endl;
                result = false;
            }
            std::cout << "."; std::cout.flush();

            blank(sa);
            timeDelta[TI] = -SEQAN_PROGETTIME;
            _createSuffixArrayPipelining(sa, text, LarssonSadakane());
            timeDelta[TI++] += SEQAN_PROGETTIME;
            if (!isSuffixArray(sa, text)) {
                std::cout << "suffix array creation (external LarssonSadakane) failed" << std::endl;
                result = false;
            }
            std::cout << "."; std::cout.flush();

            blank(sa);
            timeDelta[TI] = -SEQAN_PROGETTIME;
            createSuffixArray(sa, text, SAQSort());
            timeDelta[TI++] += SEQAN_PROGETTIME;
            if (!isSuffixArray(sa, text)) {
                std::cout << "suffix array creation (internal SAQSort) failed" << std::endl;
                result = false;
            }
            std::cout << "."; std::cout.flush();

/*            blank(sa);
            timeDelta[TI] = -SEQAN_PROGETTIME;
            createSuffixArray(sa, text, QSQGSR(), 3);
            timeDelta[TI++] += SEQAN_PROGETTIME;
            if (!isSuffixArray(sa, text)) {
                std::cout << "suffix array creation (internal QSQGSR) failed" << std::endl;
                result = false;
            }
            std::cout << "."; std::cout.flush();
*/
            blank(sa);
            timeDelta[TI] = -SEQAN_PROGETTIME;
            _createSuffixArrayPipelining(sa, text, Skew3());
            timeDelta[TI++] += SEQAN_PROGETTIME;
            if (!isSuffixArray(sa, text)) {
                std::cout << "suffix array creation (external Skew3) failed" << std::endl;
                result = false;
            }
            std::cout << "."; std::cout.flush();

            blank(sa);
            timeDelta[TI] = -SEQAN_PROGETTIME;
            _createSuffixArrayPipelining(sa, text, Skew7());
            timeDelta[TI++] += SEQAN_PROGETTIME;
            if (!isSuffixArray(sa, text)) {
                std::cout << "suffix array creation (external Skew7) failed" << std::endl;
                result = false;
            }
            std::cout << "."; std::cout.flush();

//___create_inverse_suffix_array______________________________________________

            resize(isa, size);
            timeDelta[TI] = -SEQAN_PROGETTIME;
            createInvSuffixArray(isa, sa, FromSortedSa<Serial>());
            timeDelta[TI++] += SEQAN_PROGETTIME;
            if (!isInvSuffixArray(isa, sa, text)) {
                std::cout << "inverse suffix array creation (in-memory) failed" << std::endl;
                result = false;
            }
            std::cout << "."; std::cout.flush();
            blank(isa);

            timeDelta[TI] = -SEQAN_PROGETTIME;
            createInvSuffixArray(isa, sa, FromSortedSa<Parallel>());
            timeDelta[TI++] += SEQAN_PROGETTIME;
            if (!isInvSuffixArray(isa, sa, text)) {
                std::cout << "parallel inverse suffix array creation (in-memory) failed" << std::endl;
                result = false;
            }
            std::cout << "."; std::cout.flush();
            blank(isa);

//___create_lcp_table_________________________________________________________

            resize(lcp, size);
            timeDelta[TI] = -SEQAN_PROGETTIME;
            createLcpTable(lcp, text, sa, KasaiOriginal());
            timeDelta[TI++] += SEQAN_PROGETTIME;
            if (!isLCPTable(lcp, sa, text)) {
                std::cout << "suffix array creation (internal Kasai) failed" << std::endl;
                result = false;
            }
            std::cout << "."; std::cout.flush();

            blank(lcp);
            timeDelta[TI] = -SEQAN_PROGETTIME;
            createLcpTable(lcp, text, sa, Kasai());
            timeDelta[TI++] += SEQAN_PROGETTIME;
            if (!isLCPTable(lcp, sa, text)) {
                std::cout << "suffix array creation (internal in-place Kasai) failed" << std::endl;
                result = false;
            }
            std::cout << "."; std::cout.flush();

            blank(lcp);
            timeDelta[TI] = -SEQAN_PROGETTIME;
            _createLCPTablePipelining(lcp, text, sa, Kasai());
            timeDelta[TI++] += SEQAN_PROGETTIME;
            if (!isLCPTable(lcp, sa, text)) {
                std::cout << "suffix array creation (external Kasai) failed" << std::endl;
                result = false;
            }
            std::cout << "."; std::cout.flush();

//___create_child_table_______________________________________________________

            resize(child, size);
            for(int i=0; i<size; ++i)
                child[i] = maxValue<unsigned>();
            timeDelta[TI] = -SEQAN_PROGETTIME;
            createChildtab(child, lcp);
            timeDelta[TI++] += SEQAN_PROGETTIME;
            std::cout << "."; std::cout.flush();

            unsigned undefs=0;
            for(int i=0; i<size; ++i)
                if (child[i] == maxValue<unsigned>()) ++undefs;
            if (undefs) std::cout << undefs << " undefined values";

            resize(childExt, size);
            timeDelta[TI] = -SEQAN_PROGETTIME;
            createChildtabExt(childExt, lcp);
            timeDelta[TI++] += SEQAN_PROGETTIME;
            std::cout << "."; std::cout.flush();

            if (!isEqual(child, childExt)) {
                std::cout << "child table creation failed" << std::endl;
                result = false;
            }

//___update_performance_table_________________________________________________

            for(int i=0; i<TI; ++i) {
                timeSum[i] += timeDelta[i];
                textSum += length(text);
            }

            std::cout << " OK!" << std::endl;

        }
        std::cout << "*** TIME RESULTS (sec/MB) ***" << std::endl;
        for(int i=0; i<TI; ++i)
            std::cout << algNames[i] << " " << 1024.0*1024.0 * timeSum[i] / textSum << std::endl;
}

//////////////////////////////////////////////////////////////////////////////


} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
