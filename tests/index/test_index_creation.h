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

#ifndef TESTS_INDEX_TEST_INDEX_CREATION_H
#define TESTS_INDEX_TEST_INDEX_CREATION_H

#define SEQAN_PROFILE
//#define SEQAN_DEBUG
//#define SEQAN_DEBUG_INDEX

//#define SEQAN_TEST
//#define SEQAN_TEST_SKEW3
//#define SEQAN_TEST_SKEW7


//////////////////////////////////////////////////////////////////////////////

namespace seqan
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
    typedef String<char>        TText;
    typedef String<unsigned>    TArray;

    TText   text;
    TArray  sa;
    TArray  isa;
    TArray  lcp;
    TArray  child, childExt;
    TText   bwt;

    std::string path = getAbsolutePath("/tests/index/m_tuberculosis_h37rv.fa");

    SeqFileIn inputFile(path.c_str());
    CharString id;
    readRecord(id, text, inputFile);

    unsigned size = 10000;
    resize(text, size);

//___create_suffix_array______________________________________________________

    resize(sa, size);

    blank(sa);
    createSuffixArray(sa, text, Skew3());
    if (!isSuffixArray(sa, text)) {
        std::cout << "suffix array creation (internal Skew3) failed." << std::endl;
    }

    blank(sa);
    createSuffixArray(sa, text, Skew7());
    if (!isSuffixArray(sa, text)) {
        std::cout << "suffix array creation (internal Skew7) failed." << std::endl;
    }

    blank(sa);
    createSuffixArray(sa, text, ManberMyers());
    if (!isSuffixArray(sa, text)) {
        std::cout << "suffix array creation (internal ManberMyers) failed." << std::endl;
    }

    blank(sa);
    _createSuffixArrayPipelining(sa, text, LarssonSadakane());
    if (!isSuffixArray(sa, text)) {
        std::cout << "suffix array creation (external LarssonSadakane) failed." << std::endl;
    }

    blank(sa);
    createSuffixArray(sa, text, SAQSort());
    if (!isSuffixArray(sa, text)) {
        std::cout << "suffix array creation (internal SAQSort) failed." << std::endl;
    }

//    blank(sa);
//    createSuffixArray(sa, text, QSQGSR(), 3);
//    if (!isSuffixArray(sa, text)) {
//        std::cout << "suffix array creation (internal QSQGSR) failed." << std::endl;
//    }
//
    blank(sa);
    _createSuffixArrayPipelining(sa, text, Skew3());
    if (!isSuffixArray(sa, text)) {
        std::cout << "suffix array creation (external Skew3) failed." << std::endl;
    }

    blank(sa);
    _createSuffixArrayPipelining(sa, text, Skew7());
    if (!isSuffixArray(sa, text)) {
        std::cout << "suffix array creation (external Skew7) failed." << std::endl;
    }

//___create_inverse_suffix_array______________________________________________

    resize(isa, size);

    blank(isa);
    createInvSuffixArray(isa, sa, FromSortedSa<Serial>());
    if (!isInvSuffixArray(isa, sa, text)) {
        std::cout << "inverse suffix array creation (in-memory) failed." << std::endl;
    }

    blank(isa);
    createInvSuffixArray(isa, sa, FromSortedSa<Parallel>());
    if (!isInvSuffixArray(isa, sa, text)) {
        std::cout << "parallel inverse suffix array creation (in-memory) failed." << std::endl;
    }

//___create_lcp_table_________________________________________________________

    resize(lcp, size);

    blank(lcp);
    createLcpTable(lcp, text, sa, KasaiOriginal());
    if (!isLCPTable(lcp, sa, text)) {
        std::cout << "suffix array creation (internal Kasai) failed." << std::endl;
    }

    blank(lcp);
    createLcpTable(lcp, text, sa, Kasai());
    if (!isLCPTable(lcp, sa, text)) {
        std::cout << "suffix array creation (internal in-place Kasai) failed." << std::endl;
    }

    blank(lcp);
    _createLCPTablePipelining(lcp, text, sa, Kasai());
    if (!isLCPTable(lcp, sa, text)) {
        std::cout << "suffix array creation (external Kasai) failed." << std::endl;
    }

//___create_child_table_______________________________________________________

    resize(child, size);
    for(unsigned i=0; i<size; ++i)
        child[i] = std::numeric_limits<unsigned>::max();
    createChildtab(child, lcp);

    unsigned undefs=0;
    for(unsigned i=0; i<size; ++i)
        if (child[i] == std::numeric_limits<unsigned>::max()) ++undefs;
    if (undefs) std::cout << undefs << " undefined values";

    resize(childExt, size);
    createChildtabExt(childExt, lcp);

    if (!isEqual(child, childExt)) {
        std::cout << "child table creation failed." << std::endl;
    }

}

//////////////////////////////////////////////////////////////////////////////


} //namespace seqan

#endif //#ifndef SEQAN_HEADER_...
