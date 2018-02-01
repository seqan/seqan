// ==========================================================================
//                      Yara - Yet Another Read Aligner
// ==========================================================================
// Copyright (c) 2011-2018, Enrico Siragusa, FU Berlin
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
//     * Neither the name of Enrico Siragusa or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ENRICO SIRAGUSA OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================
// This file contains functions to map read sequences to read/pair ids.
// ==========================================================================

#ifndef APP_YARA_BITS_READS_H_
#define APP_YARA_BITS_READS_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

using namespace seqan;

// ============================================================================
// Functions
// ============================================================================

template <typename TReadSeqs>
inline typename Size<TReadSeqs>::Type
getReadSeqsCount(TReadSeqs const & readSeqs)
{
    return length(readSeqs);
}

template <typename TReadSeqs>
inline typename Size<TReadSeqs>::Type
getReadsCount(TReadSeqs const & readSeqs)
{
    return length(readSeqs) / 2;
}

template <typename TReadSeqs>
inline typename Size<TReadSeqs>::Type
getPairsCount(TReadSeqs const & readSeqs)
{
    return length(readSeqs) / 4;
}

template <typename TReadSeqs, typename TReadSeqId>
inline bool isFwdReadSeq(TReadSeqs const & readSeqs, TReadSeqId readSeqId)
{
    SEQAN_ASSERT_LT(readSeqId, getReadSeqsCount(readSeqs));
    return readSeqId < getReadsCount(readSeqs);
}

template <typename TReadSeqs, typename TReadSeqId>
inline bool isRevReadSeq(TReadSeqs const & readSeqs, TReadSeqId readSeqId)
{
    SEQAN_ASSERT_LT(readSeqId, getReadSeqsCount(readSeqs));
    return !isFwdReadSeq(readSeqs, readSeqId);
}

template <typename TReadSeqs, typename TReadSeqId>
inline bool isFirstMate(TReadSeqs const & readSeqs, TReadSeqId readSeqId)
{
    SEQAN_ASSERT_LT(readSeqId, getReadSeqsCount(readSeqs));
    return isFwdReadSeq(readSeqs, readSeqId) ? readSeqId < getPairsCount(readSeqs) : readSeqId < getPairsCount(readSeqs) + getReadsCount(readSeqs);
}

template <typename TReadSeqs, typename TReadSeqId>
inline bool isSecondMate(TReadSeqs const & readSeqs, TReadSeqId readSeqId)
{
    SEQAN_ASSERT_LT(readSeqId, getReadSeqsCount(readSeqs));
    return !isFirstMate(readSeqs, readSeqId);
}

template <typename TReadSeqs, typename TPairId>
inline typename Size<TReadSeqs>::Type
getFirstMateFwdSeqId(TReadSeqs const & /* readSeqs */, TPairId pairId)
{
    return pairId;
}

template <typename TReadSeqs, typename TPairId>
inline typename Size<TReadSeqs>::Type
getSecondMateFwdSeqId(TReadSeqs const & readSeqs, TPairId pairId)
{
    SEQAN_ASSERT_LT(pairId, getPairsCount(readSeqs));
    return pairId + getPairsCount(readSeqs);
}

template <typename TReadSeqs, typename TPairId>
inline typename Size<TReadSeqs>::Type
getFirstMateRevSeqId(TReadSeqs const & readSeqs, TPairId pairId)
{
//    SEQAN_ASSERT_LT(pairId, getPairsCount(readSeqs));
    return getFirstMateFwdSeqId(readSeqs, pairId) + getReadsCount(readSeqs);
}

template <typename TReadSeqs, typename TPairId>
inline typename Size<TReadSeqs>::Type
getSecondMateRevSeqId(TReadSeqs const & readSeqs, TPairId pairId)
{
    SEQAN_ASSERT_LT(pairId, getPairsCount(readSeqs));
    return getSecondMateFwdSeqId(readSeqs, pairId) + getReadsCount(readSeqs);
}

template <typename TReadSeqs, typename TReadSeqId>
inline typename Size<TReadSeqs>::Type
getReadId(TReadSeqs const & readSeqs, TReadSeqId readSeqId)
{
    return isFwdReadSeq(readSeqs, readSeqId) ? readSeqId : readSeqId - getReadsCount(readSeqs);
}

template <typename TReadSeqs, typename TReadSeqId>
inline typename Size<TReadSeqs>::Type
getPairId(TReadSeqs const & readSeqs, TReadSeqId readSeqId)
{
    SEQAN_ASSERT_LT(readSeqId, getReadSeqsCount(readSeqs));

    typename Size<TReadSeqs>::Type pairId = readSeqId;

    if (isRevReadSeq(readSeqs, readSeqId))
        pairId -= getReadsCount(readSeqs);

    if (isSecondMate(readSeqs, readSeqId))
        pairId -= getPairsCount(readSeqs);

    SEQAN_ASSERT_LT(pairId, getPairsCount(readSeqs));
    
    return pairId;
}

template <typename TReadSeqs, typename TReadId>
inline typename Size<TReadSeqs>::Type
getMateId(TReadSeqs const & readSeqs, TReadId readId)
{
    typename Size<TReadSeqs>::Type pairId = getPairId(readSeqs, readId);

    return isFirstMate(readSeqs, readId) ? getSecondMateFwdSeqId(readSeqs, pairId) : getFirstMateFwdSeqId(readSeqs, pairId);
}

template <typename TReadSeqs, typename TReadSeqId>
inline typename Size<TReadSeqs>::Type
getMateSeqId(TReadSeqs const & readSeqs, TReadSeqId readSeqId)
{
    SEQAN_ASSERT_LT(readSeqId, getReadSeqsCount(readSeqs));

    typename Size<TReadSeqs>::Type pairId = getPairId(readSeqs, readSeqId);

    if (isFirstMate(readSeqs, readSeqId))
    {
        if (isFwdReadSeq(readSeqs, readSeqId))
            return getSecondMateRevSeqId(readSeqs, pairId);
        else
            return getSecondMateFwdSeqId(readSeqs, pairId);
    }
    else
    {
        if (isFwdReadSeq(readSeqs, readSeqId))
            return getFirstMateRevSeqId(readSeqs, pairId);
        else
            return getFirstMateFwdSeqId(readSeqs, pairId);
    }
}

template <typename TReadSeqs, typename TReadId, typename TErrors>
inline float toErrorRate(TReadSeqs const & readSeqs, TReadId readId, TErrors errors)
{
    return (float)errors / length(readSeqs[readId]);
}

#endif  // #ifndef APP_YARA_BITS_READS_H_
