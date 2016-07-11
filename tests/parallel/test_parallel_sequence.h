// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2016, Knut Reinert, FU Berlin
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
// Author: Rene Rahn <rene.rahn@fu-berlin.de>
// ==========================================================================

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/parallel.h>

template <typename TString, typename TRange, typename TTag>
void _concurrentAppend(TString & string,
                       TRange const range,
                       TTag const tag)
{
    auto val = begin(range);
    while (val != end(range))
    {
        appendValue(string, getValue(val), tag, seqan::Parallel());
        ++val;
    }
}

SEQAN_DEFINE_TEST(test_parallel_sequence_concurrent_append)
{
    using namespace seqan;
    using TString = String<unsigned, Alloc<ConcurrentAppend<> > >;
    unsigned scale = 1000;
    unsigned numThreads = std::thread::hardware_concurrency();
    TString string;

    std::vector<std::thread> vec;
    for (unsigned i = 0; i < numThreads; ++i)
        vec.push_back(std::thread(_concurrentAppend<TString, Range<unsigned>, Exact>, std::ref(string),
                                  Range<unsigned>(i * scale, (i + 1) * scale), Exact()));

    for (auto & t : vec)
        t.join();

    SEQAN_ASSERT_EQ(length(string), scale * numThreads);
    std::sort(begin(string), end(string));
    auto it = std::unique(begin(string), end(string));
    SEQAN_ASSERT_EQ(it - begin(string), scale * numThreads);
}

SEQAN_DEFINE_TEST(test_parallel_sequence_concurrent_append_insist)
{
    using namespace seqan;
    using TString = String<unsigned, Alloc<ConcurrentAppend<> > >;
    unsigned scale = 1000;
    unsigned numThreads = std::thread::hardware_concurrency();
    TString string;
    reserve(string, numThreads * scale);

    std::vector<std::thread> vec;
    for (unsigned i = 0; i < numThreads; ++i)
        vec.push_back(std::thread(_concurrentAppend<TString, Range<unsigned>, Insist>, std::ref(string),
                                  Range<unsigned>(i * scale, (i + 1) * scale), Insist()));

    for (auto & t : vec)
        t.join();

    SEQAN_ASSERT_EQ(length(string), scale * numThreads);
    std::sort(begin(string), end(string));
    auto it = std::unique(begin(string), end(string));
    SEQAN_ASSERT_EQ(it - begin(string), scale * numThreads);
}