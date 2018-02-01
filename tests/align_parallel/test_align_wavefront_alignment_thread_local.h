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
// Author: Rene Rahn <rene.rahn@fu-berlin.de>
// ==========================================================================

#include <seqan/align_parallel.h>

namespace test_align_parallel
{
struct TestConfig
{
    using TCache        = std::vector<int>;
    using TIntermediate = double;

    using TLocalHost    = std::tuple<TIntermediate, TCache>;
};
}  // namespace test_align_parallel

SEQAN_DEFINE_TEST(test_align_parallel_wavefront_alignment_thread_local_storage_construt)
{
    using namespace seqan;

    using TLocalStorage = WavefrontAlignmentThreadLocalStorage<test_align_parallel::TestConfig>;
    SEQAN_ASSERT(std::is_default_constructible<TLocalStorage>::value);
    SEQAN_ASSERT(std::is_copy_constructible<TLocalStorage>::value);
    SEQAN_ASSERT(std::is_move_constructible<TLocalStorage>::value);
    SEQAN_ASSERT(std::is_copy_assignable<TLocalStorage>::value);
    SEQAN_ASSERT(std::is_move_assignable<TLocalStorage>::value);

    {
        TLocalStorage store;
        SEQAN_ASSERT_EQ(store._multiAlignmentThreadLocal.size(), 1u);
    }

    {
        TLocalStorage store{3};
        SEQAN_ASSERT_EQ(store._multiAlignmentThreadLocal.size(), 3u);
    }
}

SEQAN_DEFINE_TEST(test_align_parallel_wavefront_alignment_thread_local_storage_intermediate)
{
    using namespace seqan;

    using TLocalStorage = WavefrontAlignmentThreadLocalStorage<test_align_parallel::TestConfig>;
    {
        TLocalStorage store{3};
        auto & interim = intermediate(store, 1);
        interim = 2.3;
        SEQAN_ASSERT_EQ(intermediate(store, 1), 2.3);
    }
}

SEQAN_DEFINE_TEST(test_align_parallel_wavefront_alignment_thread_local_storage_cache)
{
    using namespace seqan;

    using TLocalStorage = WavefrontAlignmentThreadLocalStorage<test_align_parallel::TestConfig>;
    {
        TLocalStorage store{3};
        auto & local = cache(store, 1);
        local.resize(5, 10);
        SEQAN_ASSERT_EQ(cache(store, 1).size(), 5u);
        SEQAN_ASSERT_EQ(cache(store, 1)[0], 10);
    }
}
