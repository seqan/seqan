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

#ifndef TESTS_INDEX_TEST_INDEX_HELPERS_H
#define TESTS_INDEX_TEST_INDEX_HELPERS_H


//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN
{


template < typename TBuffer >
void permute(TBuffer &buf) {
	typename Size<TBuffer>::Type i, j, s = length(buf);
//        srand( (unsigned)time( NULL ) );
	for(i = 0; i < s; i++)
		buf[i] = s-i-1;
	for(i = 0; i < s; i++) {
		if (i > 0) {
			j = i - (rand() % i) - 1;
			assert(0 <= j && j < s);
		} else
			j = 0;
		unsigned tmp = buf[i];
		buf[i] = buf[j];
		buf[j] = tmp;
	}
}

template < typename TBuffer >
void blank(TBuffer &buf) {
	typename Size<TBuffer>::Type i, s = length(buf);
	typename Value<TBuffer>::Type c = typename Value<TBuffer>::Type();
	for(i = 0; i < s; i++)
		buf[i] = c;
}

template < typename TBuffer >
void randomize(TBuffer &buf) {
	mtRandInit(false);
	typename Size<TBuffer>::Type i, s = length(buf);
	for(i = 0; i < s; i++)
		buf[i] = mtRand() % s;
}

template < typename TBuffer >
void textRandomize(TBuffer &buf) {
	mtRandInit(false);
	typename Size<TBuffer>::Type i, s = length(buf);
	for(i = 0; i < s; i++)
		buf[i] = '@' + mtRand() % 2;//('z'-'@');
}

template < typename TValue >
struct IdentityMap : public ::std::unary_function< TValue, TValue > {
	inline TValue operator() (TValue const i) { return i; }
};

template < typename TValue >
struct SimpleCompare : public ::std::binary_function< TValue const, TValue const, int > {
	inline int operator() (TValue const a, TValue const b) const {
		if (a < b) return -1;
		if (a > b) return 1;
		return 0;
	}
};

template <typename TSequence>
void printArray(TSequence const &SA) {
	typedef typename Iterator<TSequence const, Standard>::Type TIter;
	TIter it = begin(SA);
	TIter itEnd = end(SA);
	for(; it != itEnd; ++it)
		::std::cout << *it << "  ";
	::std::cout << ::std::endl;
}

template <typename TSequence>
bool isPermutation(TSequence const &SA) {
	typedef typename Value<TSequence>::Type		TSize;
	typedef typename MakeSigned_<TSize>::Type	TSigned;
	TSize n = length(SA);
	bool *seen = new bool[n];
	TSize i;
	for (i = 0;  i < n;  i++) seen[i] = 0;
	for (i = 0;  i < n;  i++)
		if ((TSigned)SA[i] >= 0 && SA[i] < n) {
			if (seen[SA[i]])
				printf("isPermutation: not unique %d->%d\n", (int)i, (int)SA[i]);
			seen[SA[i]] = true;
		} else
			printf("isPermutation: SA index out of range (n=%d) SA[%d]=%d\n", (int)n, (int)i, (int)SA[i]);

	for (i = 0;  i < n;  i++)
		if (!seen[i]) {
			printf("isPermutation: not surjective %d empty\n", (int)i);
			delete[] seen;
			return false;
		}
	delete[] seen;
	return true;
}

template <typename TInput, typename TSpec>
bool isPermutation(Pipe<TInput, TSpec> &SA) {
	typedef typename Value< Pipe<TInput, TSpec> >::Type TSize;
	TSize n = length(SA);
	bool *seen = new bool[n];
	TSize i;
	for (i = 0;  i < n;  i++) seen[i] = 0;
	beginRead(SA);
	for (i = 0;  i < n;  i++, ++SA)
		if ((TSize)*SA < n) {
			if (seen[*SA])
				printf("isPermutation: not unique %d->%d\n", (int)i, (int)*SA);
			seen[*SA] = true;
		} else
			printf("isPermutation: SA index out of range (n=%d) SA[%d]=%d\n", (int)n, (int)i, (int)*SA);
	endRead(SA);

	for (i = 0;  i < n;  i++)
		if (!seen[i]) {
			printf("isPermutation: not surjective %d empty\n", (int)i);
			delete[] seen;
			return false;
		}
	delete[] seen;
	return true;
}

template <typename TInput, typename TSpec>
bool isPermutation(Pipe<TInput, TSpec> const &SA) {
	return isPermutation(const_cast<Pipe<TInput, TSpec> &>(SA));
}


template <typename TIt, typename ST>
bool sleq__(TIt s1, TIt s2, ST n1, ST n2) {
	ST n = _min(n1, n2);
	for(ST i = 0; i < n; i++, ++s1, ++s2) {
		if (ordLess(*s1,*s2)) return 1;
		if (ordLess(*s2,*s1)) {
			::std::cerr<<(ordLess(*s2,*s1));
			::std::cerr<<*s1;
			::std::cerr<<*s2;
			::std::cerr << "after " << i << " compares not " << (unsigned)*s1 << " leq " << (unsigned)*s2 << ::std::endl;
			return 0;
		}
	}
	return (n1 < n2);
} 

// is SA a sorted suffix array for s?
template <typename TSequence, typename TText>
bool isSorted(TSequence const &SA, TText const &s) {
	typedef typename Value<TSequence>::Type TSize;
	TSize n = length(s);
	for(TSize i = 1; i < n; ++i) {
		if (!sleq__(begin(s) + SA[i-1], begin(s) + SA[i], n-SA[i-1], n-SA[i])) {
			printf("isSorted: sort error s_%d(SA[%d]) >= s_%d(SA[%d])\n",(int)SA[i-1],(int)(i-1),(int)SA[i],(int)i);
/*
		String<unsigned, External<> > safile;
		if (!open(safile,"error.sa")) printf("could not open ERROR.SA\n");
		safile = SA;
*/			    return false;
		}
	}
	return true;  
}

template <typename TInput, typename TSpec, typename TText>
bool isSorted(Pipe<TInput, TSpec> &SA, TText const &s) {
	typedef typename Value< Pipe<TInput, TSpec> >::Type TSize;
	TSize n = length(s);
	beginRead(SA);
	TSize prev = *SA; ++SA;
	for(TSize i = 1; !eof(SA); ++SA, ++i) {
		if (!sleq__(begin(s) + prev, begin(s) + *SA, n-prev, n-*SA)) {
			printf("isSorted: sort error s_%d(SA[%d]) >= s_%d(SA[%d])\n",(int)prev,(int)(i-1),(int)*SA,(int)i);
			endRead(SA);

/*			String<unsigned, External<> > safile;
		if (!open(safile,"error.sa")) printf("could not open ERROR.SA\n");
		safile << SA;
*/
			return false;
		}
		prev = *SA;
	}
	endRead(SA);
	return true;  
}

template <typename TInput, typename TSpec, typename TText>
bool isSorted(Pipe<TInput, TSpec> const &SA, TText const &s) {
	return isSorted(const_cast<Pipe<TInput, TSpec> &>(SA), s);
}


template <typename TIt, typename ST>
bool sleqLcp__(TIt s1, TIt s2, ST n1, ST n2, ST lcp) {
	ST n = _min(n1, n2);
	for(ST i = 0; i < n; i++, ++s1, ++s2) {
		if (ordLess(*s1,*s2)) return (i == lcp);
		if (ordLess(*s2,*s1)) {
			::std::cerr << "after " << i << " compares not " << (unsigned)*s1 << " leq " << (unsigned)*s2 << ::std::endl;
			return false;
		}
	}
	return (n1 < n2) && (n == lcp);
} 

// is SA a sorted suffix array and LCP the correct LCP-Table for s?
/*    template <
	typename TSize1, typename TSpec1,
	typename TSize2, typename TSpec2,
	typename TText >
bool isSortedLCP(String<TSize1, TSpec1> &LCP, String<TSize2, TSpec2> &SA, TText const &s) {
	TSize2 n = length(s);
	for(TSize2 i = 1; i < n; ++i) {
		if (!sleq__(begin(s) + SA[i-1], begin(s) + SA[i], n-SA[i-1], n-SA[i], LCP[i-1])) {
			printf("isSorted: sort error s_%d(%d) >= s_%d(%d)\n",i-1,SA[i-1],i,SA[i]);
			return false;
		}
	}
	return true;  
}
*/
template < typename TLCP, typename TSA, typename TText >
bool isSortedLCP(TLCP &LCP, TSA &SA, TText const &s) {
	typedef typename Value<TSA>::Type TSize;
	typedef typename Iterator<TSA, Rooted>::Type ISA;
	typedef typename Iterator<TLCP, Standard>::Type ILCP;

	TSize n  = length(s);
	ISA  sa  = begin(SA);
	ILCP lcp = begin(LCP);

	TSize prev = *sa; ++sa;
	for(TSize i = 1; sa != end(SA); ++sa, ++lcp, ++i) {
		if (!sleqLcp__(begin(s) + prev, begin(s) + *sa, n-prev, n-*sa, *lcp)) {
			printf("isLCP: sort error s_%d(%d) >= s_%d(%d)\n",(int)(i-1),(int)prev,(int)i,(int)*sa);
			return false;
		}
		prev = *sa;
	}
	return true;  
}

template <typename TSufArray, typename TText>
bool isSuffixArray(TSufArray &SA, TText const &s) {
	if (length(SA) != length(s)) {
		::std::cerr<<"isSuffixArray: length is bad: SA="<<length(SA)<<", s="<<length(s)<<::std::endl;
		return false;
	}
	
	if (!isPermutation(SA)) {
		::std::cerr<<"isSuffixArray: SA is not a permutation!"<<::std::endl;
		return false;
	}

	if (!isSorted(SA, s)) {
		::std::cerr<<"isSuffixArray: SA is not sorted!"<<::std::endl;
/*			String<unsigned char, External<> > textfile;
		if (!open(textfile,"error.txt")) printf("could not open ERROR.TXT\n");
		textfile=s;
*/			return false;
	}

//        ::std::cerr<<"SATest OK! n="<<length(s)<<std::endl;
	return true;
}

template <typename TLCP, typename TSufArray, typename TText>
bool isLCPTable(TLCP &LCP, TSufArray &SA, TText const &s) {
	if (length(SA) != length(s)) {
		printf("isLCPTable: length is bad: SA=%d, s=%d\n", (int)length(SA), (int)length(s));
		return false;
	}
	
	if (length(LCP) != length(s)) {
		printf("isLCPTable: length is bad: LCP=%d, s=%d\n", (int)length(LCP), (int)length(s));
		return false;
	}
	
	if (!isPermutation(SA)) {
		::std::cerr<<"isLCPTable: SA is not a permutation!\n";
		return false;
	}

	if (!isSortedLCP(LCP, SA, s)) {
		::std::cerr<<"isLCPTable: SA is not sorted!\n";
		return false;
	}

//        ::std::cerr<<"LCPTest OK! n="<<length(s)<<std::endl;
	return true;
}

template <typename TA, typename TB>
bool isEqual(TA &_a, TB &_b) {
	typedef typename Iterator<TA, Standard>::Type IA;
	typedef typename Iterator<TB, Standard>::Type IB;

	IA a = begin(_a), e = end(_a);
	IB b = begin(_b);
	while (a!=e) {
		if (!(*a == *b)) {
			::std::cerr << "isEqual: difference at " << (e-a) << " a=" << *a << "  b=" << *b << ::std::endl;
			return false;
		}
		++a;
		++b;
	}

//        ::std::cerr<<"EQUALTest OK! n="<<length(_a)<<std::endl;
	return true;
}


//////////////////////////////////////////////////////////////////////////////


} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
