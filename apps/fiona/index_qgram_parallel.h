// ==========================================================================
//                           index_qgram_parallel.h
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

#ifndef INCLUDE_SEQAN_INDEX_INDEX_QGRAM_PARALLEL_H_
#define INCLUDE_SEQAN_INDEX_INDEX_QGRAM_PARALLEL_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// TODO(weese): should go to basic/metaprogramming...
template <bool VALUE_>
struct Bool: public Eval<VALUE_> {};

template <unsigned VALUE_>
struct Unsigned
{
    enum { VALUE = VALUE_ };
};

template <int VALUE_>
struct Integer
{
    enum { VALUE = VALUE_ };
};


// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

template < typename TIndex, typename TParallelTag >
inline bool _qgramDisableBuckets(TIndex &index, Tag<TParallelTag>)
{
	// use the serial version (if no parallel overload is available)
    return _qgramDisableBuckets(index);
}

//////////////////////////////////////////////////////////////////////////////
// Counting sort - Step 3: Cumulative sum
//
// a disabled bucket 3,2,x,4   (x = disabled)               y_i=sum_{j<=i} x_j
// results in        0,0,3,x,5 (_qgramCummulativeSum)       z_i=y_{i-2}
// or in             0,3,5,5,9 (_qgramCummulativeSumAlt)    z_i=y_{i-1}

//    3 2 1 | x x 5 |
// 0 0 3 | 5 x x | 6 11
//   0 3 5 6 6 6 11

// 3 2 1 | 3*  2  5 |
// 3 5 6 | 9 11 16 |        
// 0 3 5 | 6  9 11 | 16     SHIFT 0
// 0 0 3 | 5  6*  9 | 11 16  SHIFT 1


// First two entries are 0.
// Step 4 increments the entries hash(qgram)+1 on-the-fly while filling the SA table.
// After step 4 each entry (0..n-1) is the beginning of a qgram bucket.

template < typename TSequence, typename TParallelTag >
inline typename Value<TSequence>::Type
_sumIgnoreDisabled(TSequence const &seq, Tag<TParallelTag>)
{
    typedef typename Value<TSequence>::Type TValue;
    typename Iterator<TSequence const>::Type it = begin(seq, Standard());
    typename Iterator<TSequence const>::Type itEnd = end(seq, Standard());
    TValue sum = 0;
    for (; it != itEnd; ++it)
        if (*it != (TValue)-1)
            sum += *it;
    return sum;
}

//////////////////////////////////////////////////////////////////////////////
// Counting sort - Step 2: Count q-grams
template < typename TDir, typename TBucketMap, typename TText, typename TShape, typename TStepSize, typename TParallelTag >
inline void
_qgramCountQGrams(TDir &dir, TBucketMap &bucketMap, TText const &text, TShape shape, TStepSize stepSize, Tag<TParallelTag> parallelTag)
{
    typedef typename Iterator<TText const, Standard>::Type  TIterator;
    typedef typename Iterator<TDir, Standard>::Type         TDirIterator;
    typedef typename Value<TDir>::Type                      TSize;

    if (empty(shape) || length(text) < length(shape))
        return;

    TSize num_qgrams = (length(text) - length(shape)) / stepSize + 1;
    TDirIterator dirBegin = begin(dir, Standard());
    Splitter<TSize> splitter(0, num_qgrams, parallelTag);

    if (stepSize == 1)
    {
        SEQAN_OMP_PRAGMA(parallel for firstprivate(shape))
        for (int job = 0; job < (int)length(splitter); ++job)
        {
            TIterator itText = begin(text, Standard()) + splitter[job];
            TIterator itTextEnd = begin(text, Standard()) + splitter[job + 1];

            atomicInc(*(dirBegin + requestBucket(bucketMap, hash(shape, itText), parallelTag)), parallelTag);
            for (++itText; itText != itTextEnd; ++itText)
                atomicInc(*(dirBegin + requestBucket(bucketMap, hashNext(shape, itText), parallelTag)), parallelTag);
        }
    }
    else
    {
        SEQAN_OMP_PRAGMA(parallel for firstprivate(shape))
        for (int job = 0; job < (int)length(splitter); ++job)
        {
            TIterator itText = begin(text, Standard()) + (splitter[job] * stepSize);
            TIterator itTextEnd = begin(text, Standard()) + (splitter[job + 1] * stepSize);

            for (; itText != itTextEnd; itText += stepSize)
                atomicInc(*(dirBegin + requestBucket(bucketMap, hash(shape, itText), parallelTag)), parallelTag);
        }
    }
}

template < typename TDir, typename TBucketMap, typename TString, typename TSpec, typename TShape, typename TStepSize, typename TParallelTag >
inline void
_qgramCountQGrams(TDir &dir, TBucketMap &bucketMap, StringSet<TString, TSpec> const &stringSet, TShape shape, TStepSize stepSize, Tag<TParallelTag> parallelTag)
{
    typedef typename Iterator<TString const, Standard>::Type    TIterator;
    typedef typename Iterator<TDir, Standard>::Type             TDirIterator;
    typedef typename Value<TDir>::Type                          TSize;

    if (empty(shape) || empty(stringSet))
        return;

    TDirIterator dirBegin = begin(dir, Standard());
    Splitter<TSize> seqSplitter(0, countSequences(stringSet), parallelTag);

    if (stepSize == 1)
    {
        SEQAN_OMP_PRAGMA(parallel for firstprivate(shape))
        for (int job = 0; job < (int)length(seqSplitter); ++job)
        {
			for(unsigned seqNo = seqSplitter[job]; seqNo < seqSplitter[job + 1]; ++seqNo)
			{
				TString const &sequence = value(stringSet, seqNo);
				if (length(sequence) < length(shape)) continue;

                TIterator itText = begin(stringSet[seqNo], Standard());
                TIterator itTextEnd = itText + (length(sequence) - length(shape) + 1);
                atomicInc(*(dirBegin + requestBucket(bucketMap, hash(shape, itText), parallelTag)), parallelTag);
                for (++itText; itText != itTextEnd; ++itText)
                    atomicInc(*(dirBegin + requestBucket(bucketMap, hashNext(shape, itText), parallelTag)), parallelTag);
            }
        }
    }
    else
    {
        SEQAN_OMP_PRAGMA(parallel for firstprivate(shape))
        for (int job = 0; job < (int)length(seqSplitter); ++job)
        {
			for(unsigned seqNo = seqSplitter[job]; seqNo < seqSplitter[job + 1]; ++seqNo)
			{
				TString const &sequence = value(stringSet, seqNo);
				if (length(sequence) < length(shape)) continue;

                TIterator itText = begin(stringSet[seqNo], Standard());
                TIterator itTextEnd = itText + ((length(sequence) - length(shape)) / stepSize + 1) * stepSize;
                for (; itText != itTextEnd; itText += stepSize)
                    atomicInc(*(dirBegin + requestBucket(bucketMap, hash(shape, itText), parallelTag)), parallelTag);
            }
        }
    }
}

template < typename TDir, typename TWithConstraints, typename TKeepDisabledBuckets, unsigned SHIFT, typename TParallelTag >
inline typename Value<TDir>::Type
_qgramCummulativeSum(TDir &dir, TWithConstraints, TKeepDisabledBuckets, Unsigned<SHIFT>, Tag<TParallelTag> parallelTag)
{
    typedef typename Value<TDir>::Type TValue;
    typedef typename Size<TDir>::Type TSize;
    typedef String<TValue> TBuffer;
    typedef typename Iterator<TDir const, Standard>::Type TConstIterator;
    typedef typename Iterator<TDir, Standard>::Type TIterator;
    typedef typename Iterator<TBuffer, Standard>::Type TConstBufferIterator;

    if (empty(dir))
        return 0;

    Splitter<TSize> splitter(0, length(dir), parallelTag);
    String<TValue> localSums;
    TBuffer prevCounts;
    resize(localSums, length(splitter), Exact());
    resize(prevCounts, length(splitter) * SHIFT, 0, Exact());
    localSums[0] = 0;

    // STEP 1: compute sums of all subintervals (in parallel)
    //
    SEQAN_OMP_PRAGMA(parallel for)
    for (int job = 0; job < (int)length(splitter); ++job)
    {
        unsigned len = _min(SHIFT, splitter[job]);
        replace(prevCounts, (job + 1) * SHIFT - len, (job + 1) * SHIFT, infix(dir, splitter[job] - len, splitter[job]));

        typename Infix<TDir>::Type dirInfix = infix(
            dir,
            _max((int64_t)0, (int64_t)splitter[job] - (int64_t)SHIFT),
            _max((int64_t)0, (int64_t)splitter[job + 1] - (int64_t)SHIFT));

        if (TWithConstraints::VALUE)
            localSums[job] = _sumIgnoreDisabled(dirInfix, Serial());
        else
            localSums[job] = sum(dirInfix, Serial());        
    }

    // STEP 2: compute partial sums (of subinterval sums) from position 0 to the end of each subinterval
    //
    for (int job = 1; job < (int)length(splitter); ++job)
        localSums[job] += localSums[job - 1];

    // STEP 3: compute partial sums of each subinterval starting from offset (in parallel)
    //
    SEQAN_OMP_PRAGMA(parallel for)
    for (int job = 0; job < (int)length(splitter); ++job)
    {
        TConstIterator itBegin = begin(dir, Standard()) + splitter[job];
        TIterator dstIt = begin(dir, Standard()) + splitter[job + 1];
        TValue sum = localSums[job];

        // read over our subinterval
        {
            TConstIterator it = begin(dir, Standard()) + splitter[job + 1] - SHIFT;
            while (it != itBegin)
            {
                TValue counter = *(--it);
                if (!TWithConstraints::VALUE || counter != (TValue)-1)
                    sum -= counter;
                else
                    if (TKeepDisabledBuckets::VALUE)
                    {
                        *(--dstIt) = (TValue)-1;
                        continue;
                    }
                *(--dstIt) = sum;
            }
        }

        // read suffix of the previous subinterval
        if (SHIFT != 0u)
        {
            TConstBufferIterator it = begin(prevCounts, Standard()) + (job + 1) * SHIFT;
            while (dstIt != itBegin)
            {
                TValue counter = *(--it);
                if (!TWithConstraints::VALUE || counter != (TValue)-1)
                    sum -= counter;
                else
                    if (TKeepDisabledBuckets::VALUE)
                    {
                        *(--dstIt) = (TValue)-1;
                        continue;
                    }
                *(--dstIt) = sum;
            }
        }
    }
    
    return back(localSums);
}

//////////////////////////////////////////////////////////////////////////////
// Counting sort - Step 4: Fill suffix array
// w/o constraints
template <
    typename TSA,
    typename TText, 
    typename TShape, 
    typename TDir, 
    typename TBucketMap, 
    typename TWithConstraints, 
    typename TStepSize,
    typename TParallelTag >
inline void
_qgramFillSuffixArray(
    TSA &sa, 
    TText const &text, 
    TShape shape, 
    TDir &dir, 
    TBucketMap &bucketMap, 
    TStepSize stepSize,
    TWithConstraints const,
    Tag<TParallelTag> parallelTag)
{
    typedef typename Iterator<TText const, Standard>::Type  TIterator;
    typedef typename Iterator<TDir, Standard>::Type         TDirIterator;
    typedef typename Value<TDir>::Type                      TSize;

    if (empty(shape) || length(text) < length(shape))
        return;

    TSize num_qgrams = (length(text) - length(shape)) / stepSize + 1;
    TDirIterator dirBegin = begin(dir, Standard());
    TDirIterator dirBegin1 = dirBegin + 1;
    Splitter<TSize> splitter(0, num_qgrams, parallelTag);

    if (stepSize == 1)
    {
        SEQAN_OMP_PRAGMA(parallel for firstprivate(shape))
        for (int job = 0; job < (int)length(splitter); ++job)
        {
            TSize pos = splitter[job];
            TSize posEnd = splitter[job + 1];
            if (pos == posEnd) continue;
            TIterator itText = begin(text, Standard()) + pos;

            // first hash
            TDirIterator const bktPtr = dirBegin1 + getBucket(bucketMap, hash(shape, itText));
            if (!TWithConstraints::VALUE || *bktPtr != (TSize)-1)       // ignore disabled buckets
                sa[atomicPostInc(*bktPtr, parallelTag)] = pos;          

            for (++pos; pos != posEnd; ++pos)
            {
                TDirIterator const bktPtr = dirBegin1 + getBucket(bucketMap, hashNext(shape, ++itText));
                if (!TWithConstraints::VALUE || *bktPtr != (TSize)-1)   // ignore disabled buckets
                    sa[atomicPostInc(*bktPtr, parallelTag)] = pos;
            }
        }
    }
    else
    {
        SEQAN_OMP_PRAGMA(parallel for firstprivate(shape))
        for (int job = 0; job < (int)length(splitter); ++job)
        {
            TSize pos = splitter[job] * stepSize;
            TSize posEnd = splitter[job + 1] * stepSize;
            TIterator itText = begin(text, Standard()) + pos;

            for (; pos != posEnd; pos += stepSize, itText += stepSize)
            {
                TDirIterator const bktPtr = dirBegin1 + getBucket(bucketMap, hash(shape, itText));
                if (!TWithConstraints::VALUE || *bktPtr != (TSize)-1)   // ignore disabled buckets
                    sa[atomicPostInc(*bktPtr, parallelTag)] = pos;
            }
        }
    }
}

template <
    typename TSA,
    typename TString,
    typename TSpec,
    typename TShape, 
    typename TDir, 
    typename TBucketMap, 
    typename TWithConstraints, 
    typename TStepSize,
    typename TParallelTag >
inline void
_qgramFillSuffixArray(
    TSA &sa, 
    StringSet<TString, TSpec> const &stringSet,
    TShape shape, 
    TDir &dir, 
    TBucketMap &bucketMap, 
    TStepSize stepSize,
    TWithConstraints const,
    Tag<TParallelTag> parallelTag)
{
    typedef typename Iterator<TString const, Standard>::Type  TIterator;
    typedef typename Iterator<TDir, Standard>::Type         TDirIterator;
    typedef typename Value<TDir>::Type                      TSize;

    if (empty(shape) || empty(stringSet))
        return;

    TDirIterator dirBegin = begin(dir, Standard());
    TDirIterator dirBegin1 = dirBegin + 1;
    Splitter<TSize> seqSplitter(0, countSequences(stringSet), parallelTag);

    if (stepSize == 1)
    {
        SEQAN_OMP_PRAGMA(parallel for firstprivate(shape))
        for (int job = 0; job < (int)length(seqSplitter); ++job)
        {
			for(unsigned seqNo = seqSplitter[job]; seqNo < seqSplitter[job + 1]; ++seqNo)
			{
				TString const &sequence = value(stringSet, seqNo);
				if (length(sequence) < length(shape)) continue;

                TIterator itText = begin(stringSet[seqNo], Standard());
                TIterator itTextEnd = itText + (length(sequence) - length(shape) + 1);

				typename Value<TSA>::Type localPos;
				assignValueI1(localPos, seqNo);
				assignValueI2(localPos, 0);

                // first hash
                TDirIterator const bktPtr = dirBegin1 + getBucket(bucketMap, hash(shape, itText));
                if (!TWithConstraints::VALUE || *bktPtr != (TSize)-1)       // ignore disabled buckets
                    sa[atomicPostInc(*bktPtr, parallelTag)] = localPos;          

                for (++itText; itText != itTextEnd; ++itText)
                {
                    posInc(localPos);
                    TDirIterator const bktPtr = dirBegin1 + getBucket(bucketMap, hashNext(shape, itText));
                    if (!TWithConstraints::VALUE || *bktPtr != (TSize)-1)   // ignore disabled buckets
                        sa[atomicPostInc(*bktPtr, parallelTag)] = localPos;
                }
            }
        }
    }
    else
    {
        SEQAN_OMP_PRAGMA(parallel for firstprivate(shape))
        for (int job = 0; job < (int)length(seqSplitter); ++job)
        {
			for(unsigned seqNo = seqSplitter[job]; seqNo < seqSplitter[job + 1]; ++seqNo)
			{
				TString const &sequence = value(stringSet, seqNo);
				if (length(sequence) < length(shape)) continue;

                TIterator itText = begin(stringSet[seqNo], Standard());
                TIterator itTextEnd = itText + ((length(sequence) - length(shape)) / stepSize + 1) * stepSize;

				typename Value<TSA>::Type localPos;
				assignValueI1(localPos, seqNo);
				assignValueI2(localPos, 0);

                for (; itText != itTextEnd; ++itText)
                {
                    posInc(localPos, stepSize);
                    TDirIterator const bktPtr = dirBegin1 + getBucket(bucketMap, hash(shape, itText));
                    if (!TWithConstraints::VALUE || *bktPtr != (TSize)-1)   // ignore disabled buckets
                        sa[atomicPostInc(*bktPtr, parallelTag)] = localPos;
                }
            }
        }
    }
}

//////////////////////////////////////////////////////////////////////////////
// Step 5: Correct disabled buckets
template < typename TDir, typename TParallelTag >
inline void
_qgramPostprocessBuckets(TDir &dir, Tag<TParallelTag> parallelTag)
{
    typedef typename Iterator<TDir, Standard>::Type TDirIterator;
    typedef typename Size<TDir>::Type               TSize;
    typedef typename Value<TDir>::Type              TValue;

    Splitter<TDirIterator> splitter(begin(dir, Standard()), end(dir, Standard()), parallelTag);
    String<TValue> last;
    resize(last, length(splitter), Exact());

    SEQAN_OMP_PRAGMA(parallel for)
    for (int job = 0; job < (int)length(splitter); ++job)
    {
        TDirIterator it = splitter[job];
        TDirIterator itEnd = splitter[job + 1];
        TSize prev = (job == 0)? 0: (TSize)-1;
        for (; it != itEnd; ++it)
            if (*it == (TSize)-1)
                *it = prev;
            else
                prev = *it;
        last[job] = prev;
    }

    for (int job = 1; job < (int)length(splitter); ++job)
        if (last[job] == (TSize)-1)
            last[job] = last[job - 1];

    SEQAN_OMP_PRAGMA(parallel for)
    for (int job = 1; job < (int)length(splitter); ++job)
    {
        TDirIterator it = splitter[job];
        TDirIterator itEnd = splitter[job + 1];
        TSize prev = last[job - 1];
        for (; it != itEnd && *it == (TSize)-1; ++it)
            *it = prev;
    }
}

template < typename TIndex, typename TParallelTag >
void createQGramIndex(TIndex &index, Tag<TParallelTag> parallelTag)
{
    typename Fibre<TIndex, QGramText>::Type const &text      = indexText(index);
    typename Fibre<TIndex, QGramSA>::Type         &sa        = indexSA(index);
    typename Fibre<TIndex, QGramDir>::Type        &dir       = indexDir(index);
    typename Fibre<TIndex, QGramShape>::Type      &shape     = indexShape(index);
    typename Fibre<TIndex, QGramBucketMap>::Type  &bucketMap = index.bucketMap;
    
    // 1. clear counters
    _qgramClearDir(dir, bucketMap, parallelTag);

    // 2. count q-grams
    _qgramCountQGrams(dir, bucketMap, text, shape, getStepSize(index), parallelTag);

    if (_qgramDisableBuckets(index, parallelTag))
    {
        // 3. cumulative sum

        // with disabled buckets
        // disabled buckets should still be marked in the partial sum
        // shift all entries by one towards the end (will be corrected by _qgramFillSuffixArray)
        _qgramCummulativeSum(dir, True(), True(), Unsigned<1>(), parallelTag);

        // 4. fill suffix array
        _qgramFillSuffixArray(sa, text, shape, dir, bucketMap, getStepSize(index), True(), parallelTag);

        // 5. correct disabled buckets
        _qgramPostprocessBuckets(dir, parallelTag);
    }
    else
    {
        // 3. cumulative sum

        // without disabled buckets
        // shift all entries by one towards the end (will be corrected by _qgramFillSuffixArray)
        _qgramCummulativeSum(dir, False(), False(), Unsigned<1>(), parallelTag);
        
        // 4. fill suffix array
        _qgramFillSuffixArray(sa, text, shape, dir, bucketMap, getStepSize(index), False(), parallelTag);
    }
}

template < typename TIndex, typename TParallelTag >
void createQGramIndexDirOnly(TIndex &index, Tag<TParallelTag> parallelTag)
{
    typename Fibre<TIndex, QGramText>::Type const &text      = indexText(index);
    typename Fibre<TIndex, QGramSA>::Type         &sa        = indexSA(index);
    typename Fibre<TIndex, QGramDir>::Type        &dir       = indexDir(index);
    typename Fibre<TIndex, QGramShape>::Type      &shape     = indexShape(index);
    typename Fibre<TIndex, QGramBucketMap>::Type  &bucketMap = index.bucketMap;
    
    // 1. clear counters
    _qgramClearDir(dir, bucketMap, parallelTag);

    // 2. count q-grams
    _qgramCountQGrams(dir, bucketMap, text, shape, getStepSize(index), parallelTag);

    // 3. cumulative sum (Step 4 is ommited)
    if (_qgramDisableBuckets(index, parallelTag))
    {
        // with disabled buckets
        // disabled buckets should not be marked in the partial sum
        // don't shift entries as there is no _qgramFillSuffixArray call
        _qgramCummulativeSum(dir, True(), False(), Unsigned<0>(), parallelTag); // no shift
    }
    else
    {
        // without disabled buckets
        // don't shift entries as there is no _qgramFillSuffixArray call
        _qgramCummulativeSum(dir, False(), False(), Unsigned<0>(), parallelTag);
    }
}

}

#endif  // #ifndef INCLUDE_SEQAN_INDEX_INDEX_QGRAM_PARALLEL_H_
