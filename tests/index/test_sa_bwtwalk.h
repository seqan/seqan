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
// Author: Marcel Martin <marcel.martin@tu-dortmund.de>
// Author: Tobias Marschall <tobias.marschall@tu-dortmund.de>
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================

#ifndef TESTS_INDEX_TEST_SA_BWTWALK_H
#define TESTS_INDEX_TEST_SA_BWTWALK_H


//////////////////////////////////////////////////////////////////////////////

namespace seqan
{

template < typename TSA, typename TText>
void print_sa(TSA &sa, TText &s) {
    typedef typename Iterator<TSA>::Type TIter;
    int i = 0;
    for (TIter it = begin(sa); it!=end(sa); ++it, ++i) {
//         std::cout << i << " " << *it << " " << suffix(s,*it) << std::endl;
        std::cout << i << " " << *it << " " << infix(s,*it, min(length(s),*it+100)) << std::endl;
    }
    std::cout << std::endl;
}

template < typename TSA1, typename TSA2 >
bool compare_sa(TSA1 &sa1, TSA2 &sa2) {
    if (length(sa1) != length(sa2)) return false;
    typedef typename Iterator<TSA1>::Type TIter1;
    typedef typename Iterator<TSA2>::Type TIter2;
    TIter1 it1 = begin(sa1);
    TIter2 it2 = begin(sa2);
    int i = 0;
    for (; it1!=end(sa1); ++it1, ++it2, ++i) {
        if (*it1!=*it2) {
            std::cout << "Mismatch at position " << i << std::endl;
            return false;
        }
    }
    return true;
}

//template < typename TText, typename TTag, typename TValue, typename TAllowsFastRandomAccess >
//bool check_sa_algorithm(TText& text, TAllowsFastRandomAccess&);

// a suffix array of 'text' is computed using the algorithm indicated by 'tag'
// If 'saReference' is not empty, it is assumed to be the correct suffix array
// and compared against.
template < typename TText, typename TTag, typename TValue>
bool check_sa_algorithm(TText& text, const False&)
{
    String<TValue> sa;
    resize(sa, length(text));
//    std::cout << "check_sa_algorithm<" << typeid(TTag).name() << "," << typeid(TValue).name() << ",False> ... " << std::flush;
    createSuffixArray(sa, text, TTag());
//    std::cout << "done" << std::endl;
    return isSuffixArray(sa, text);
}

// a suffix array of 'text' is computed using the algorithm indicated by 'tag'
// If 'saReference' is not empty, it is assumed to be the correct suffix array
// and compared against.
template < typename TText, typename TTag, typename TValue>
bool check_sa_algorithm(TText& text, const True&)
{
    String<TValue, External<> > sa;
    resize(sa, length(text));
//    std::cout << "check_sa_algorithm<" << typeid(TTag).name() << "," << typeid(TValue).name() << ",True> ... " << std::flush;
    createSuffixArray(sa, text, TTag());
//    std::cout << "done" << std::endl;
    bool ok = isSuffixArray(sa, text);
    return ok;
}

#define MYASSERT(tag, type, use64) if (!check_sa_algorithm<CharString, BwtWalk<tag>, type>(text, use64())) { std::cerr << "Assertion failed in line " << __LINE__ << ": " << #tag << " " << #type << " " << #use64 << std::endl; exit(1); }

SEQAN_DEFINE_TEST(testBWTWalk)
{
//#if defined(__GNUC__)
//#  if defined(__OPTIMIZE__)
//    std::cout << "Optimized build: Yes\n";
//#  else
//    std::cout << "Optimized build: No\n";
//#  endif
//#endif

    std::string path = getAbsolutePath("/tests/index/m_tuberculosis_h37rv.fa");

    SeqFileIn inputFile(path.c_str());
    CharString text, id;
    readRecord(id, text, inputFile);
    resize(text, 10000);
//    std::cout << "textsize: " << length(text) << std::endl;

    MYASSERT(BwtWalkFast, unsigned, False);
    MYASSERT(BwtWalkFast, unsigned, True);
    MYASSERT(BwtWalkFast, uint64_t, False);
    MYASSERT(BwtWalkFast, uint64_t, True);

    MYASSERT(BwtWalkInPlace, unsigned, False);
    MYASSERT(BwtWalkInPlace, unsigned, True);
    MYASSERT(BwtWalkInPlace, uint64_t, False);
    MYASSERT(BwtWalkInPlace, uint64_t, True);
}

//////////////////////////////////////////////////////////////////////////////


} //namespace seqan

#endif //#ifndef SEQAN_HEADER_...
