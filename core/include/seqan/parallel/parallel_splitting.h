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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Utility macros for parallelism.
// ==========================================================================

#ifndef SEQAN_PARALLEL_PARALLEL_SPLITTING_H_
#define SEQAN_PARALLEL_PARALLEL_SPLITTING_H_

namespace seqan {

/**
.Function.computeSplitters
..cat:Parallelism
..summary:Compute splitters for a sequence of objects.
..signature:computeSplitters(splitters, size, count)
..param.splitters:Resulting splitters, will be resized to contain $count + 1$ elements.
...type:Spec.Alloc String
..param.size:The number of objects to split.
..param.count:The number of chunks.
..remarks:The first $count - 1$ chunks will have the size $ceil(size / count)$, the last chunk will contain the rest.
..example.text:Most simple case for splitting.
..example.code:String<unsigned> splitters;
computeSplitters(splitters, 10, 5);
// splitters == {0, 5, 10}
..example.text:In this case, the last chunks will stay empty.
..example.code:computeSplitters(splitters, 3, 5);
// splitters == {0, 1, 2, 3, 3, 3}
..include:seqan/parallel.h
 */

template <typename TPosString, typename TSize, typename TCount>
void computeSplitters(TPosString & splitters, TSize size, TCount count)
{
    resize(splitters, count + 1);
    splitters[0] = 0;
    TSize blockLength = size / count;
    TCount rest = size % count;
    for (TCount i = 1; i <= count; ++i)
    {
        splitters[i] = splitters[i - 1] + blockLength;
        if (i <= rest)
            splitters[i] += 1;
    }

    typedef typename Value<TPosString>::Type TPos SEQAN_TYPEDEF_FOR_DEBUG;
    SEQAN_ASSERT_EQ(back(splitters), static_cast<TPos>(size));
}

}  // namespace seqan

#endif  // SEQAN_PARALLEL_PARALLEL_SPLITTING_H_
