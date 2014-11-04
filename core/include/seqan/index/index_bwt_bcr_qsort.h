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
// Author: David Iwanowitsch <iwanowit@inf.fu-berlin.de>
// ==========================================================================

#ifndef SANDBOX_IWANOWIT_INCLUDE_SEQAN_INDEX2_QUICKSORT_H_
#define SANDBOX_IWANOWIT_INCLUDE_SEQAN_INDEX2_QUICKSORT_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct qsortParallel_;
typedef Tag<qsortParallel_> qsortParallel;

struct qsortSequential_;
typedef Tag<qsortSequential_> qsortSequential;

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

//// sort SA according to sortKeys
//// choose between quickSort and insertionSort
////
//template<typename TToSort, typename TSortFunctor, typename TIndex>
//void doQuickSortSeq(TToSort &toSort, TSortFunctor &sortFunctor, const TIndex start, const TIndex end) {
//	if (end - start > 15)
//		quickSort(qsortSequential(), toSort, sortFunctor, start, end);
//	else
//		insertionSort(toSort, sortFunctor, start, end);
//}
//
//// sort SA according to sortKeys
//// choose between quickSort and insertionSort
////
//template<typename TToSort, typename TSortFunctor, typename TIndex>
//void doQuickSortPar(TToSort &toSort, TSortFunctor &sortFunctor, const TIndex start, const TIndex end) {
//	if (end - start > 15)
//		quickSort(qsortParallel(), toSort, sortFunctor, start, end);
//	else
//		insertionSort(toSort, sortFunctor, start, end);
//}

// sort SA according to sortKeys
// choose between quickSort and insertionSort
//
template<typename TToSort, typename TSortFunctor, typename TIndex>
void doQuickSort(const qsortSequential &tag, TToSort &toSort, TSortFunctor &sortFunctor, const TIndex start, const TIndex end) {
	if (end - start > 15)
		quickSort(tag, toSort, sortFunctor, start, end);
	else
		insertionSort(toSort, sortFunctor, start, end);
}

template<typename TToSort, typename TSortFunctor, typename TIndex>
void doQuickSort(const qsortParallel, TToSort &toSort, TSortFunctor &sortFunctor, const TIndex start, const TIndex end) {
	if (end - start > 15)
		quickSort(qsortParallel(), toSort, sortFunctor, start, end, 0);
	else
		insertionSort(toSort, sortFunctor, start, end);
}

template<typename TToSort, typename TSortFunctor, typename TIndex>
void doQuickSort(const qsortParallel, TToSort &toSort, TSortFunctor &sortFunctor, const TIndex start, const TIndex end, const long depth) {
	if (end - start > 15)
		quickSort(qsortParallel(), toSort, sortFunctor, start, end, depth);
	else
		insertionSort(toSort, sortFunctor, start, end);
}

template<typename Val>
int medianOfThree(const Val a, const Val b, const Val c) {
	if (a == b)
		return 0;
	if (a == c || b == c)
		return 2;
	if (a < b) {
		if (b < c)
			return 1;
		else if (a < c)
			return 2;
		else
			return 0;
	} else if (b > c)
		return 1;
	else if (a < c)
		return 0;
	else
		return 2;
}

template<typename TResult, typename TSA, typename TSortFunctor, typename TIndex>
TResult _findPivot(TSA &SA, TSortFunctor &sortFunctor, const TIndex left, const TIndex right) {

	typedef typename TSortFunctor::result_type TSortKey;

	const TIndex halfBucketSize = (right - left) / 2;
	TSortKey pivot = sortFunctor(right); // getBptrVal(SA, bptr, limits, bptrExtPerString, offset, right);
	const TSortKey pivotB = sortFunctor(left); //getBptrVal(SA, bptr, limits, bptrExtPerString, offset, left);
	const TSortKey pivotC = sortFunctor(left + halfBucketSize); // getBptrVal(SA, bptr, limits, bptrExtPerString, offset, left + halfBucketSize);
	int medianNumber = medianOfThree(pivot, pivotB, pivotC);

	if (medianNumber != 0) {
		pivot = medianNumber == 1 ? pivotB : pivotC;
		const TIndex swapIndex =
				medianNumber == 1 ? left : left + halfBucketSize;
		std::swap(SA[swapIndex], SA[right]);
	}
	return pivot;
}

template<typename TResult, typename TPivot, typename TSA, typename TSortFunctor, typename TIndex>
TResult _partition(TSA &SA, const TPivot pivot, TSortFunctor &sortFunctor, const TIndex left, const TIndex right) {

		//lomutos partitioning scheme
		TIndex leftTmp = left;
		TIndex rightTmp = left;
		while (rightTmp < right) {
			if (sortFunctor(rightTmp) <= pivot) {
				std::swap(SA[leftTmp], SA[rightTmp]);

				leftTmp++;
			}
			rightTmp++;
		}

		std::swap(SA[leftTmp], SA[rightTmp]);

		return leftTmp;
}

template<typename TSA, typename TSortFunctor, typename TIndex>
void quickSort(const qsortSequential &tag, TSA &SA, TSortFunctor &sortFunctor, const TIndex left, const TIndex right) {

	typedef typename TSortFunctor::result_type TSortKey;

	TSortKey pivot = _findPivot<TSortKey>(SA, sortFunctor, left, right);

	TIndex leftTmp = _partition<TIndex>(SA, pivot, sortFunctor, left, right);

	if (right - leftTmp > 1) {
		doQuickSort(tag, SA, sortFunctor, leftTmp + 1, right);
	}
	do {
		--leftTmp;
	} while (leftTmp > left && sortFunctor(leftTmp) == pivot);
	if (leftTmp - left > 0) {
		doQuickSort(tag, SA, sortFunctor, left, leftTmp);
	}
}

template<typename TSA, typename TSortFunctor, typename TIndex>
void quickSort(const qsortParallel &tag, TSA &SA, TSortFunctor &sortFunctor, const TIndex left, const TIndex right, long depth) {

	if(depth > 10){
		//Avoid the creation of too many tasks. Their overhead would make it very slow
		quickSort(qsortSequential(), SA, sortFunctor, left, right);
		return;
	}

	typedef typename TSortFunctor::result_type TSortKey;

	TSortKey pivot = _findPivot<TSortKey>(SA, sortFunctor, left, right);

	TIndex leftTmp = _partition<TIndex>(SA, pivot, sortFunctor, left, right);

	const TIndex leftTmpCopy = leftTmp;

	SEQAN_OMP_PRAGMA(task shared(SA, tag, sortFunctor)) {
		if (right - leftTmp > 1) {
			doQuickSort(tag, SA, sortFunctor, leftTmpCopy + 1, right, depth + 1);
		}
	}
//	SEQAN_OMP_PRAGMA(task shared(SA, tag, sortFunctor)) {
		do {
			--leftTmp;
		} while (leftTmp > left && sortFunctor(leftTmp) == pivot);
		if (leftTmp - left > 0) {
			doQuickSort(tag, SA, sortFunctor, left, leftTmp, depth + 1);
		}
//	}
}

template<typename TSA, typename TSortFunctor, typename TIndex>
void insertionSort(TSA &SA, TSortFunctor &sortFunctor, const TIndex left, const TIndex right) {

	typedef typename Value<TSA>::Type TSAVal;
	typedef typename TSortFunctor::result_type TSortKey;

	TIndex j = 0;
	TSortKey sortkey;
	TSAVal tmp;

	for (TIndex i = left + 1; i <= right; ++i) {
		sortkey = sortFunctor(i);
		tmp = SA[i];
		j = i - 1;

		int c = 1;
		while ((j >= left) && (sortFunctor(j) > sortkey)) {
			SA[j + 1] = SA[j];
			if (j == 0) {
				c = 0;
				break; //underflow with unsigned
			}
			--j;
		}
		SA[j + c] = tmp;
	}
}

}  // namespace seqan

#endif  // #ifndef SANDBOX_IWANOWIT_INCLUDE_SEQAN_INDEX2_QUICKSORT_H_
