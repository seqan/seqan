// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2012, Knut Reinert, FU Berlin
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

#ifndef SANDBOX_IWANOWIT_INCLUDE_SEQAN_INDEX2_INDEX_SA_BPR_H_
#define SANDBOX_IWANOWIT_INCLUDE_SEQAN_INDEX2_INDEX_SA_BPR_H_

#include <seqan/basic.h>
#include <seqan/parallel.h>

namespace SEQAN_NAMESPACE_MAIN {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct Bpr_;
typedef Tag<Bpr_> Bpr;

template <typename InType, typename TSa, typename TBptr, typename TLimits, typename TOffset, typename OutType = typename Value<TBptr>::Type>
struct _bprComparator :
    public std::unary_function<InType, OutType>
{

    const TSa & SA;
    const TBptr & bptr;
    const TLimits & limits;
    const long bptrExtPerString;
    const TOffset offset;

    _bprComparator(const TSa & _sa, const TBptr & _bptr, const TLimits & _limits, const long _bptrExtPerString, const TOffset _offset) :
        SA(_sa), bptr(_bptr), limits(_limits), bptrExtPerString(_bptrExtPerString), offset(_offset)
    {}

    OutType operator()(InType index)
    {
        return getBptrVal(SA, bptr, limits, bptrExtPerString, offset, index);
    }

};

// ============================================================================
// Metafunctions
// ============================================================================

// TODO: Alphabet should refer to Value by default.

template <typename TObject>
struct AlphabetType
{
    typedef TObject Type;
};

template <typename TString, typename TSpec>
struct AlphabetType<StringSet<TString, TSpec> >
{
    typedef typename Value<TString>::Type Type;
};

template <typename TString, typename TSpec>
struct AlphabetType<String<TString, TSpec> >
{
    typedef typename Value<String<TString, TSpec> >::Type Type;
};

// ============================================================================
// Functions
// ============================================================================

//////////////////////////////////////////////////////////////////////////////
// Bucket Pointer Refinement Implementation
// find the suffix array SA of s[0..n-1] in {0..K}^n
//
// creates suffix array SA of s
// chars have to be in the range [0,K)
//
// d defines the suffixlength for the initial bucket sort
template <typename TSA, typename TText>
void createSuffixArray(TSA & SA, TText & s, Bpr const &, const unsigned short d)
{
    typedef typename Value<TSA>::Type TSAValue;
    typedef typename AlphabetType<TText>::Type TAlphabet;
    typedef typename ValueSize<TAlphabet>::Type TAlphabetSize;
    typedef typename Size<TText>::Type TTextSize;
    typedef typename StringSetLimits<TText>::Type TStringSetLimits;

    if (lengthSum(s) < 1)
        return;

    const TStringSetLimits limits = stringSetLimits(s);
    const TAlphabetSize ALPHABETSIZE = ValueSize<TAlphabet>::VALUE;
    const TAlphabetSize ALPHSIZEWITHDOLLAR = ValueSize<TAlphabet>::VALUE + 1;

    //TODO: remove primitive types

/////////////////////////
// Phase 1:
// sort suffixes with length d and create bucketPointer bptr
//
/////////////////////////

    //set contains the count of the alphabetsize^d buckets
    String<TTextSize> bkt;

    //contains for each entry in text, the index of the last entry of the corresponding Bucket in SA
    //remark: contains negative values to get the right order of the $signs
    String<long> bptr;

#if (SEQAN_ENABLE_DEBUG || SEQAN_ENABLE_TESTING) && SEQAN_ENABLE_PARALLELISM
    const double begin_time = omp_get_wtime();
#endif

    fillBucketsPhase1(SA, s, bptr, bkt, limits, d, ALPHSIZEWITHDOLLAR);

#if (SEQAN_ENABLE_DEBUG || SEQAN_ENABLE_TESTING) && SEQAN_ENABLE_PARALLELISM
    std::cout << "\t Phase 1: " << double(omp_get_wtime() - begin_time)
              << "s. ";
#endif

/////////////////////////
// Phase 2:
// recursively refine buckets with size > 1 by sorting buckets with corresponding
// offset and refresh bptr
/////////////////////////

#if (SEQAN_ENABLE_DEBUG || SEQAN_ENABLE_TESTING) && SEQAN_ENABLE_PARALLELISM
    const double begin_time_2 = omp_get_wtime();
#endif

    sewardCopyPhase2(SA, s, bptr, bkt, limits, d, ALPHSIZEWITHDOLLAR);

    clear(bkt);
    clear(bptr);

#if (SEQAN_ENABLE_DEBUG || SEQAN_ENABLE_TESTING) && SEQAN_ENABLE_PARALLELISM
    std::cout << "\t Phase 2: "
              << double(omp_get_wtime() - begin_time_2) << "s. ";
#endif

}

/*
 * sort suffixes with length d and create Bucketpointer
 *
 * Set bucketpointer to buckets last position in SA
 */
template <typename TSA, typename TText, typename TBptr, typename TBkt,
          typename TLimit, typename TD, typename TAlpha>
inline void fillBucketsPhase1(TSA & SA, const TText & s, TBptr & bptr, TBkt & bkt,
                              const TLimit & limits, const TD d, const TAlpha alphabeSizeWithDollar)
{

    typedef typename SAValue<TText>::Type TSAValue;
    typedef typename Size<TText>::Type TTextSize;

    //after each sequence we insert 2*d+1 values to sort the $ correctly
    const unsigned bptrExtPerString = (2 * d + 1);

    const long tmpModulo = pow(alphabeSizeWithDollar, (d - 1));
    const long bucketCount = tmpModulo * alphabeSizeWithDollar;

    const TTextSize n = lengthSum(s);
    const TTextSize stringCount = countSequences(s);

    resize(bptr, n + stringCount * bptrExtPerString, Exact());
    resize(bkt, bucketCount + 1, Exact());

    Splitter<TTextSize> splitter(0, n);
    String<TBkt> bktPerThread;
    resize(bktPerThread, length(splitter), Exact());

    //parallel: Split text in parts of equal size for each thread
    SEQAN_OMP_PRAGMA(parallel num_threads(length(splitter)))
    {
        const unsigned threadNum = omp_get_thread_num();
        const TTextSize start = splitter[threadNum];
        const TTextSize end = splitter[threadNum + 1];
        TSAValue saValue;
        long bptrOffset;

        TBkt & threadBucket = bktPerThread[threadNum];
        resize(threadBucket, bucketCount + 1, 0, Exact());

        //first iteration: compute rank and bucketsizes
        posLocalize(saValue, start, limits);
        bptrOffset = getSeqNo(saValue) * bptrExtPerString;

        long codeD = code_d(
            getSequenceByNo(getSeqNo(saValue), s), d, saValue,
            alphabeSizeWithDollar);
        bptr[start + bptrOffset] = codeD;
        threadBucket[codeD]++;

        for (TTextSize i = start + 1; i < end; ++i)
        {
            posLocalize(saValue, i, limits);
            bptrOffset = getSeqNo(saValue) * bptrExtPerString;

            codeD = code_d(
                getSequenceByNo(getSeqNo(saValue), s), d, saValue,
                alphabeSizeWithDollar, tmpModulo,
                codeD);
            bptr[i + bptrOffset] = codeD;

            threadBucket[codeD]++;
        }

        SEQAN_OMP_PRAGMA(barrier)

        //summarize buckets: each thread needs to know how many entries the previous threads are going to insert into the bucket
        //bkt contains the starting positions (=number of entrys in smaller buckets)
        SEQAN_OMP_PRAGMA(single)
        {
            bkt[0] = 0;
            int last = 0;
            int tmp;
            for (unsigned i = 0; i < length(bktPerThread[0]); ++i)
            {

                //TODO: kommentieren, hier wird gleichzeit bkt und bktperTrhead kummuliert
                int perThreadLast = 0;

                int bktOffset = 1;
                if (i == bucketCount)
                    bktOffset = 0;
                //TODO: wofür brauche ich den bktOffset?

                bkt[i + bktOffset] = last;
                for (int j = 0; j < length(splitter); ++j)
                {
                    bkt[i + bktOffset] += bktPerThread[j][i];
                    tmp = bktPerThread[j][i];
                    bktPerThread[j][i] = perThreadLast;
                    perThreadLast += tmp;
                }
                last = bkt[i + bktOffset];
            }
            bkt[bucketCount] = n;
        }

        //now fill the buckets
        for (TTextSize i = start; i < end; ++i)
        {
            posLocalize(saValue, i, limits);
            bptrOffset = getSeqNo(saValue) * bptrExtPerString;
            long index = (bkt[bptr[i + bptrOffset] + 1]
                          - threadBucket[bptr[i + bptrOffset]]) - 1;

            threadBucket[bptr[i + bptrOffset]]++;
            SA[index] = saValue;
        }

        //up to now we misused bptr to save the rank-value of a suffix
        //now bptr should be a pointer from suffix to bucket
        for (TTextSize i = start; i < end; ++i)
        {
            posLocalize(saValue, i, limits);
            bptrOffset = getSeqNo(saValue) * bptrExtPerString;
            bptr[i + bptrOffset] = bkt[bptr[i + bptrOffset] + 1] - 1;
        }

        SEQAN_OMP_PRAGMA(barrier)

        //the last d-suffixes get special handling corresponding to the sorting of the $-signs
        SEQAN_OMP_PRAGMA(single)
        {

            for (int seq = stringCount - 1; seq >= 0; --seq)
            {
                bptrOffset = seq * bptrExtPerString;

                unsigned seqEndIndex = sequenceLength(seq, s);
                unsigned seqStartIndex = posGlobalize(
                    Pair<unsigned, unsigned>(seq, 0), limits);
                const int tmpValue = code_d(getSequenceByNo(seq, s), d,
                                            seqEndIndex - d - 1, alphabeSizeWithDollar);
                int lastValue = tmpValue;

                for (long i = seqEndIndex - d; i < seqEndIndex; ++i)
                {
                    lastValue = code_d(getSequenceByNo(seq, s), d, i,
                                       alphabeSizeWithDollar, tmpModulo, lastValue);

                    bptr[bptrOffset + seqStartIndex + i] = bkt[lastValue];
                    bkt[lastValue] = bkt[lastValue] + 1;
                }

                //insert negative values after each sequence
                long j = seq * 2 * d + seq;
                j *= -1;
                for (long i = seqEndIndex; i <= seqEndIndex + 2 * d; ++i)
                {
                    bptr[seqStartIndex + i + bptrOffset] = --j;
                }
            }

            //reinsert correct bkt values (needed for seward copy)
            for (int seq = stringCount - 1; seq >= 0; --seq)
            {
                bptrOffset = seq * bptrExtPerString;

                unsigned seqEndIndex = sequenceLength(seq, s);
                long lastValue = code_d(getSequenceByNo(seq, s), d,
                                        seqEndIndex - d - 1, alphabeSizeWithDollar);
                for (long i = seqEndIndex - d; i < seqEndIndex; ++i)
                {
                    lastValue = code_d(getSequenceByNo(seq, s), d, i,
                                       alphabeSizeWithDollar, tmpModulo, lastValue);
                    bkt[lastValue] = bkt[lastValue] - 1;
                }
            }
        }
    }

    clear(bktPerThread);
}

// computes the "multiple character encoding" (the rank)
//
//&s: SuffixArray
//d: länge
//i: pos in s
template <typename TText, typename TPos, typename TAlphabetSize>
unsigned int code_d(TText & s, const unsigned short & d, TPos pos,
                    const TAlphabetSize & alphabetSize)
{
    unsigned int result = 0;
    unsigned local = getSeqOffset(pos);

    for (unsigned short k = 1; k <= d; ++k)
    {
        if (length(s) <= (local + k - 1))
            result += 0;
        else
            result += pow(alphabetSize, (d - k))
                      * (unsigned) (1 + ordValue(s[local + k - 1]));
    }
    return result;
}

// computes the "multiple character encoding" using the previous result
//
//&s: SuffixArray
//d: länge
//i: pos in s
//code_d_i: wert vom vorherigen suffix
template <typename TText, typename TPos, typename TAlphabetSize>
unsigned code_d(TText & s, const unsigned short & d, TPos pos,
                const TAlphabetSize & alphabetSize, long modulo, unsigned code_d_i)
{
    unsigned local = getSeqOffset(pos);
    if (local == 0)
        return code_d(s, d, pos, alphabetSize);

    if (length(s) <= local + d - 1)
    {
        return alphabetSize * (code_d_i % modulo) + 0;
    }
    else
        return alphabetSize * (code_d_i % modulo) + 1
               + ordValue(s[local + d - 1]);
}

template <typename TSA, typename TText, typename TBptr, typename TBkt,
          typename TLimit, typename TD, typename TAlpha>
inline void sewardCopyPhase2(TSA & SA, const TText & s, TBptr & bptr,
                             const TBkt & bkt, const TLimit & limits, const TD d,
                             const TAlpha alphabeSizeWithDollar)
{

    typedef typename Value<TSA>::Type TSAValue;
    typedef typename AlphabetType<TText>::Type TAlphabet;
    typedef typename ValueSize<TAlphabet>::Type TAlphabetSize;
    typedef typename Size<TText>::Type TTextSize;
    typedef typename Value<TBptr>::Type TBptrValue;

    const TAlphabetSize ALPHABETSIZE = ValueSize<TAlphabet>::VALUE;
    const unsigned bptrExtPerString = (2 * d + 1);

    //Find L1 and L2 Buckets to Sort the characters according to their count
    String<TAlphabet> orderedAlphabet = getOrdererdAlphabet(bkt, s, d);

    const long bucketsInL2Bucket = pow(alphabeSizeWithDollar, d - 2); //d minus 2, since these chars are equal in L2
    const long bucketsInL1Bucket = bucketsInL2Bucket * alphabeSizeWithDollar;

    StringSet<String<Triple<long, long, long> > > bucketIndicesSet;
    resize(bucketIndicesSet, 1);

    /*
     * Collect all Buckets starting with two different characters, starting with the smallest
     */
    String<Triple<long, long, long> > & bucketIndicesCur = bucketIndicesSet[0];
    for (unsigned i = 0; i < ALPHABETSIZE; ++i)
    {
        const TAlphabetSize firstChar = 1
                                        + ordValue(orderedAlphabet[ALPHABETSIZE - 1 - i]);

        for (TAlphabetSize j = 0; j < ALPHABETSIZE; ++j)
        {
            if (j == i)
                continue;
            const TAlphabetSize secondChar = 1
                                             + ordValue(orderedAlphabet[ALPHABETSIZE - 1 - j]);

            const long bktStartIndex = firstChar * bucketsInL1Bucket
                                       + secondChar * bucketsInL2Bucket; //startbucket for suffixes starting with firstchar, secondchar
            const long bktEndIndex = bktStartIndex + bucketsInL2Bucket;


            for (long k = bktStartIndex; k < bktEndIndex; ++k)   //all buckets starting with both characters
            {
                const long left = bkt[k]; //Start
                const long right = bkt[k + 1] - 1; //end: next bucket start -1

                long bucketSize = right - left;
                if (bucketSize > 0)
                {
                    appendValue(bucketIndicesCur, Triple<long, long, long>(left, right, d));
                }
            }
        }
    }

#if (SEQAN_ENABLE_DEBUG || SEQAN_ENABLE_TESTING) && SEQAN_ENABLE_PARALLELISM
    const clock_t begin_time = omp_get_wtime();
#endif

    /*
     * Now sort these Buckets
     */

    StringSet<String<Triple<long, long, long> > > nextBucketIndices;
    String<Triple<long, long, long> > bucketIndices = concat(bucketIndicesSet);
    String<long> split;
    String<Pair<long, long> > resultValues;

    String<String<Pair<long, long> > > nextBptrValuesPerThread;
    resize(nextBptrValuesPerThread, omp_get_max_threads());

    while (length(bucketIndices) > 0)
    {

        //Split buckets in reasonable chunks, after each chunk sync refined bptr
        clear(split);
        appendValue(split, 0);
        long currentLength = 0;
        long lastSplit = 0;
        const long minBuckets = 100 * omp_get_max_threads();
        for (unsigned k = 0; k < length(bucketIndices); ++k)
        {
            if (k - lastSplit > minBuckets)
            {
                appendValue(split, k);
                currentLength = 0;
                lastSplit = k;
            }
            currentLength += bucketIndices[k].i2 - bucketIndices[k].i1;
        }
        appendValue(split, length(bucketIndices));
        resize(resultValues, length(bucketIndices));

        for (unsigned i = 1; i < length(split); ++i)
        {

            long startIndex = split[i - 1];
            long endIndex = split[i];

            SEQAN_OMP_PRAGMA(parallel shared(bucketIndices, nextBucketIndices, SA, bptr, limits))
            {

                //parallely sort all collected buckets
                SEQAN_OMP_PRAGMA(for schedule(guided))
                for (unsigned k = startIndex; k < endIndex; ++k)
                {

                    long offset = bucketIndices[k].i3;
                    long start = bucketIndices[k].i1;
                    long end = bucketIndices[k].i2;

                    if (end - start == 1)
                    {
                        sortSizeTwo(SA, bptr, limits, d, offset, start, end);
                        bucketIndices[k].i3 = 0; //marker to ignore this bucket after refinement (no recursion needed)
                    }
                    else
                    {
                        typedef _bprComparator<TTextSize, TSA, TBptr, TLimit, TTextSize> TSortFunctor;
                        TSortFunctor sortFunctor = TSortFunctor(SA, bptr, limits, bptrExtPerString, offset);
                        doQuickSort(qsortSequential(), SA, sortFunctor, start, end);
                    }
                }

                //collect all changed bptr values per thread:
                String<Pair<long, long> > & nextBptrValues = nextBptrValuesPerThread[omp_get_thread_num()];

                //parallely refine the previously sorted buckets
                SEQAN_OMP_PRAGMA(for schedule(guided))
                for (unsigned k = startIndex; k < endIndex; ++k)
                {

                    long offset = bucketIndices[k].i3;
                    long start = bucketIndices[k].i1;
                    long end = bucketIndices[k].i2;

                    //TODO: erkläre resultValues
                    resultValues[k] = refineBucket(SA, bptr, nextBptrValues, limits, start, end, offset, d);
                }

                // note, the implicit barrier after the for-loop above

                //sync bptr, to use refined pointers in next chunk
                for (unsigned k = 0; k < length(nextBptrValues); ++k)
                {
                    bptr[nextBptrValues[k].i1] = nextBptrValues[k].i2;
                }

                clear(nextBptrValues);
            }
        }

        resize(nextBucketIndices, omp_get_max_threads());

        //parallely search next buckets for next iteration
        SEQAN_OMP_PRAGMA(parallel for schedule(guided))
        for (unsigned k = 0; k < length(bucketIndices); ++k)
        {
            long offset = bucketIndices[k].i3;
            if (offset == 0)
                continue;
            long start = bucketIndices[k].i1;
            long end = bucketIndices[k].i2;

            Pair<long, long> middleValues = resultValues[k];

            String<Triple<long, long, long> > & nextBucketIndicesThread = nextBucketIndices[omp_get_thread_num()];

            unsigned int newOffset = offset + d;
            long bptrOffset = getSeqNo(SA[start]) * bptrExtPerString;
            //increase offset
            if (bptr[bptrOffset + posGlobalize(SA[start], limits)] == end)
            {
                newOffset = computeLCPAndOffset(SA, bptr, limits, start, end,
                                                newOffset, d);
            }

            long leftTmp = start;
            while (leftTmp < middleValues.i1)
            {
                bptrOffset = getSeqNo(SA[leftTmp]) * bptrExtPerString;
                long rightTmp = bptr[bptrOffset + posGlobalize(SA[leftTmp], limits)];

                long tmp = rightTmp - leftTmp;
                if (tmp > 0)
                {
                    appendValue(nextBucketIndicesThread, Triple<long, long, long>(leftTmp, rightTmp, newOffset));
                }
                leftTmp = rightTmp + 1;
            }

            leftTmp = middleValues.i2 + 1;
            while (leftTmp < end)
            {
                bptrOffset = getSeqNo(SA[leftTmp]) * bptrExtPerString;
                long rightTmp = bptr[bptrOffset + posGlobalize(SA[leftTmp], limits)];
                long tmp = rightTmp - leftTmp;
                if (tmp > 0)
                {
                    appendValue(nextBucketIndicesThread, Triple<long, long, long>(leftTmp, rightTmp, newOffset));
                }
                leftTmp = rightTmp + 1;
            }

            if (middleValues.i2 > middleValues.i1 + 1)
            {
                appendValue(nextBucketIndicesThread, Triple<long, long, long>(middleValues.i1 + 1, middleValues.i2, 2 * offset));
            }
        }

        clear(resultValues);
        clear(bucketIndicesSet);
        swap(bucketIndicesSet, nextBucketIndices);
        bucketIndices = concat(bucketIndicesSet);
    }

#if (SEQAN_ENABLE_DEBUG || SEQAN_ENABLE_TESTING) && SEQAN_ENABLE_PARALLELISM
    std::cout << " bucketSort: " << float(omp_get_wtime() - begin_time) << std::endl;

    const clock_t begin_time2 = omp_get_wtime();
#endif

    //Now use seward copy to generate the order for all buckets starting with the same two characters
    SEQAN_OMP_PRAGMA(parallel for)
    for (unsigned i = 0; i < ALPHABETSIZE; ++i)
    {
        const TAlphabetSize firstChar = 1 + ordValue(orderedAlphabet[i]);

        long leftVal = bkt[firstChar * bucketsInL1Bucket
                           + firstChar * bucketsInL2Bucket];
        long rightVal = bkt[firstChar * bucketsInL1Bucket
                            + (firstChar + 1) * bucketsInL2Bucket];

        long left = bkt[firstChar * bucketsInL1Bucket];
        long right = bkt[(firstChar + 1) * bucketsInL1Bucket];

        while (left < leftVal)
        {
            TSAValue tmp = SA[left]; //firstChar-bucket ist in der tmp-ten Stelle des Textes //TODO
            TTextSize seqOffset = getSeqOffset(tmp);
            const long seqNum = getSeqNo(tmp);
            if (seqOffset <= 0)
            {
                ++left;
                continue;
            }
            TAlphabetSize character = 1
                                      + ordValue(getSequenceByNo(seqNum, s)[--seqOffset]);
            if (firstChar == character)
            {
                const long bptrOffset = seqNum * bptrExtPerString;
                setSeqOffset(tmp, seqOffset);
                bptr[posGlobalize(tmp, limits) + bptrOffset] = leftVal;
                SA[leftVal] = tmp;
                ++leftVal;
            }
            ++left;
        }

        while (left < right)
        {
            --right;
            TSAValue tmp = SA[right];
            TTextSize seqOffset = getSeqOffset(tmp);
            const long seqNum = getSeqNo(tmp);
            TAlphabetSize character;
            if (seqOffset > 0 && firstChar == (character =
                                                   1 + ordValue(getSequenceByNo(seqNum, s)[--seqOffset])))
            {

                --rightVal;
                const long bptrOffset = bptrExtPerString * seqNum;
                setSeqOffset(tmp, seqOffset);
                bptr[posGlobalize(tmp, limits) + bptrOffset] = rightVal;
                SA[rightVal] = tmp;
            }
        }
    }

#if (SEQAN_ENABLE_DEBUG || SEQAN_ENABLE_TESTING) && SEQAN_ENABLE_PARALLELISM
    std::cout << " seward: " << float(omp_get_wtime() - begin_time2) << std::endl;
#endif
}

template <typename TText>
TText getOrdererdAlphabet(const String<long> & bkt, const StringSet<TText> & s,
                          const unsigned short & d)
{
    return getOrdererdAlphabet(bkt, s[0], d);
}

// Sorts all characters according to their count within the text
// Theirfore it uses the sizes of the different Buckets
//
template <typename TText>
TText getOrdererdAlphabet(const String<long> & bkt, const TText & /*s*/,
                          const unsigned short & d)
{
    typedef typename AlphabetType<TText>::Type TAlphabet;
    typedef typename ValueSize<TAlphabet>::Type TAlphabetSize;
    const TAlphabetSize ALPHABETSIZE = ValueSize<TAlphabet>::VALUE;

    String<TAlphabet> orderedAlphabet;
    resize(orderedAlphabet, ALPHABETSIZE, 0);
    String<long> alphabetWeights;
    resize(alphabetWeights, ALPHABETSIZE, 0);

    String<TAlphabet> t1;
    resize(t1, d, (TAlphabet) 0);    //ie: "aaa"
    String<TAlphabet> t2;
    resize(t2, d, (TAlphabet) (ALPHABETSIZE - 1));    //ie "zzz"
    String<TAlphabet> t3;
    resize(t3, d, (TAlphabet) 0);    //ie: "aaa"
    String<TAlphabet> t4;
    resize(t4, d, (TAlphabet) (ALPHABETSIZE - 1));    //ie "zzz"

    for (TAlphabetSize i = 0; i < ALPHABETSIZE; ++i)
    {
        t1[0] = (TAlphabet) i;    //ie: "caa"
        t2[0] = (TAlphabet) i;    //ie: "czz"

        t3[0] = (TAlphabet) i;
        t4[0] = (TAlphabet) i;
        t3[1] = (TAlphabet) i;    //ie: "cca"
        t4[1] = (TAlphabet) i;    //ie: "ccz"
        unsigned int codeD_t1 = code_d(t1, d, 0, ALPHABETSIZE + 1);
        unsigned int codeD_t2 = code_d(t2, d, 0, ALPHABETSIZE + 1);

        unsigned int codeD_t3 = code_d(t3, d, 0, ALPHABETSIZE + 1);
        unsigned int codeD_t4 = code_d(t4, d, 0, ALPHABETSIZE + 1);

        if (codeD_t2 < length(bkt) - 1)
            ++codeD_t2;
        if (codeD_t4 < length(bkt) - 1)
            ++codeD_t4;

        alphabetWeights[i] = (bkt[codeD_t2] - bkt[codeD_t1])
                             - (bkt[codeD_t4] - bkt[codeD_t3]);
        orderedAlphabet[i] = (TAlphabet) i;
        //alphabetWeights[i] = Number of buckets starting with i-th character - number of Buckets starting with 2x i-th character
    }

    //now sort characters according to their weight
    for (TAlphabetSize i = 1; i < ALPHABETSIZE; ++i)
    {
        const long tmpWeight = alphabetWeights[i];
        TAlphabetSize j = i;
        while (j > 0
              && tmpWeight < alphabetWeights[ordValue(orderedAlphabet[j - 1])])
        {
            orderedAlphabet[j] = orderedAlphabet[j - 1];
            --j;
        }
        orderedAlphabet[j] = (TAlphabet) i;
    }

    return orderedAlphabet;
}

//update bptr after sort
//this part has its origin in the original algorithm by schürmann and stoye: updatePtrAndRefineBuckets_SaBucket
template <typename TSA, typename TBptr, typename TSize, typename TLimits>
Pair<TSize, TSize> refineBucket(TSA & SA, TBptr & bptr, String<Pair<TSize, TSize> > & nextBptr, TLimits & limits, TSize left, TSize right,
                                unsigned int offset, unsigned short d)
{

    unsigned bptrExtPerString = (2 * d + 1);

    if (right - left == 1)
    {
        long bptrOffset1 = getSeqNo(SA[left]) * bptrExtPerString;
        appendValue(nextBptr, Pair<TSize, TSize>(bptrOffset1 + posGlobalize(SA[left], limits), left));
        return Pair<TSize, TSize>(0, 0);
    }

    typedef typename Value<TSA>::Type TSAVal;

    /*
     * from right to left compare sort keys, as long as they are equal they are in the same bucket.
     * set the bucket to the rightmost position
     */
    TSize leftInterval = right;
    TSize rightInterval = right;
    TSize tmp;

    TSize middleRight;
    TSize middleLeft;

    while (left <= leftInterval
          && right < (tmp = getBptrVal(SA, bptr, limits, bptrExtPerString, offset, leftInterval)))     // bptr[SA[leftInterval] + offset]))
    {
        do
        {
            long bptrOffset = getSeqNo(SA[leftInterval]) * bptrExtPerString;
            // nextBptr[bptrOffset + posGlobalize(SA[leftInterval], limits)] =	rightInterval;
            appendValue(nextBptr, Pair<TSize, TSize>(bptrOffset + posGlobalize(SA[leftInterval], limits), rightInterval));
            --leftInterval;
        }
        while (left <= leftInterval
              && getBptrVal(SA, bptr, limits, bptrExtPerString, offset, leftInterval)  /*bptr[SA[leftInterval] + offset]*/
               == tmp);
        rightInterval = leftInterval;
    }

    rightInterval = leftInterval;
    while (left <= leftInterval
          && left <= getBptrVal(SA, bptr, limits, bptrExtPerString, offset, leftInterval)   //bptr[SA[leftInterval] + offset]
          && getBptrVal(SA, bptr, limits, bptrExtPerString, offset, leftInterval) <= right)
    {
        long bptrOffset = getSeqNo(SA[leftInterval]) * bptrExtPerString;
        //nextBptr[bptrOffset + posGlobalize(SA[leftInterval], limits)] =	rightInterval;
        appendValue(nextBptr, Pair<TSize, TSize>(bptrOffset + posGlobalize(SA[leftInterval], limits), rightInterval));
        --leftInterval;
    }

    middleRight = rightInterval;
    middleLeft = leftInterval;

    rightInterval = leftInterval;
    while (left <= leftInterval)
    {
        const TSize tmp2 = getBptrVal(SA, bptr, limits, bptrExtPerString, offset, leftInterval); //bptr[SA[leftInterval] + offset];
        do
        {
            long bptrOffset = getSeqNo(SA[leftInterval]) * bptrExtPerString;
            //nextBptr[bptrOffset + posGlobalize(SA[leftInterval], limits)] =	rightInterval;
            appendValue(nextBptr, Pair<TSize, TSize>(bptrOffset + posGlobalize(SA[leftInterval], limits), rightInterval));
            --leftInterval;
        }
        while (left <= leftInterval
              && getBptrVal(SA, bptr, limits, bptrExtPerString, offset, leftInterval)  /*bptr[SA[leftInterval] + offset] */
               == tmp2);
        rightInterval = leftInterval;
    }

    return Pair<TSize, TSize>(middleLeft, middleRight);
}

//finds the longest common prefix and increases the current offset
template <typename TSa, typename TBptr, typename TSize, typename TLimits>
unsigned int computeLCPAndOffset(TSa & SA, TBptr & bptr,
                                 const TLimits & limits, const TSize left, const TSize right,
                                 const unsigned int offset, const unsigned short d)
{

    typedef typename Value<TBptr>::Type TBptrVal;

    unsigned bptrExtPerString = (2 * d + 1);

    unsigned int lcp = offset;
    while (true)
    {
        TSize index = left;
        long bptrOffset = getSeqNo(SA[right]) * bptrExtPerString;
        TBptrVal tmp = bptr[bptrOffset + posGlobalize(SA[right], limits) + lcp];
        while (index < right)
        {
            bptrOffset = getSeqNo(SA[index]) * bptrExtPerString;
            if (bptr[bptrOffset + posGlobalize(SA[index], limits) + lcp] != tmp)
            {
                return lcp;
            }
            ++index;
        }
        lcp += d;
    }
}

//Sort and Refine a Bucket with Size 2 without recursion
template <typename TSA, typename TBptr, typename TIndex, typename TLimits>
void sortSizeTwo(TSA & SA, TBptr & bptr, const TLimits & limits,
                 const unsigned short d, const unsigned int offset, const TIndex left,
                 const TIndex right)
{

    typedef typename Value<TSA>::Type TSAVal;

    unsigned bptrExtPerString = (2 * d + 1);
    long bptrOffset1 = getSeqNo(SA[left]) * bptrExtPerString;
    long bptrOffset2 = getSeqNo(SA[right]) * bptrExtPerString;

    unsigned suffix1 = posGlobalize(SA[left], limits) + offset;
    unsigned suffix2 = posGlobalize(SA[right], limits) + offset;

    while (bptr[bptrOffset1 + suffix1] == bptr[bptrOffset2 + suffix2])
    {
        suffix1 += d;
        suffix2 += d;
    }

    if (bptr[bptrOffset1 + suffix1] > bptr[bptrOffset2 + suffix2])
    {
        const TSAVal tmp = SA[left];
        SA[left] = SA[right];
        SA[right] = tmp;
    }
}

// gets the correct value from bptr, respects alls necesary offsets
//
template <typename TSA, typename TBptrVal, typename TLimits, typename TExtString, typename TOffset, typename TIndex>
TBptrVal getBptrVal(const TSA & SA, const String<TBptrVal> & bptr, const TLimits & limits, const TExtString bptrExtPerString, TOffset offset, TIndex index)
{

    const long bptrOffset = getSeqNo(SA[index]) * bptrExtPerString;
    return bptr[posGlobalize(SA[index], limits) + offset + bptrOffset];
}

}
// namespace seqan

#endif  // #ifndef SANDBOX_IWANOWIT_INCLUDE_SEQAN_INDEX2_INDEX_SA_BPR_H_
