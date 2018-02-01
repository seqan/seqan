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

#include <algorithm>

#include <seqan/parallel.h>

namespace test_align_parallel
{
struct TestValue
{
    std::string mMsg{"default constructed"};

    //NOTE(rrahn) Bug in g++-4.9 prevents us from using as aggregate type.
    TestValue() = default;

    TestValue(std::string const _msg) : mMsg(std::move(_msg))
    {}
};

template <typename TEtl>
inline void
run(TEtl & tls, size_t const threadCount)
{
    unsigned counter{0};
    std::mutex mutexCout;
    // spawn sveral threads that call for a number of repetitions the local variable.
    auto task = [&counter, &mutexCout, &tls](unsigned const tid)
    {
        while (true)
        {
            auto & val = local(tls);
            if (val.first == "master")
            {
                std::stringstream ss;
                ss << "thread_" << tid;
                val.first = ss.str();
            }
            {
                std::lock_guard<std::mutex> lck(mutexCout);
                if (counter == 1000)
                    break;
                ++counter;
            }
            --val.second;
        }
    };

    std::vector<std::thread> pool;
    for (unsigned tid = 0; tid < threadCount; ++tid)
        pool.emplace_back(task, tid);

    std::for_each(std::begin(pool), std::end(pool),
                  [](auto & thread)
                  {
                      if (thread.joinable())
                          thread.join();
                  });
}

template <typename TEtl>
void
testEnumerate(TEtl & etl)
{
    unsigned count{0};

    for (auto it = begin(etl); it != end(etl); ++it)
    {
        count += 1000 - it->second;
    }
    SEQAN_ASSERT_EQ(count, 1000u);
}

template <typename TEtl>
void
testEnumerateConst(TEtl const & etl)
{
    unsigned count{0};

    for (auto it = begin(etl); it != end(etl); ++it)
    {
        count += 1000 - it->second;
    }
    SEQAN_ASSERT_EQ(count, 1000u);
}

}

SEQAN_DEFINE_TEST(test_parallel_enumerable_thread_local_construct)
{
    using namespace seqan;

    SEQAN_ASSERT(std::is_default_constructible<EnumerableThreadLocal<int>>::value);
    SEQAN_ASSERT(!std::is_copy_constructible<EnumerableThreadLocal<int>>::value);
    SEQAN_ASSERT(!std::is_move_constructible<EnumerableThreadLocal<int>>::value);
    SEQAN_ASSERT(!std::is_copy_assignable<EnumerableThreadLocal<int>>::value);
    SEQAN_ASSERT(!std::is_move_assignable<EnumerableThreadLocal<int>>::value);

    {  // Default construction.
        EnumerableThreadLocal<test_align_parallel::TestValue> tls;
        SEQAN_ASSERT_EQ(tls._initValue.mMsg, "default constructed");
    }

    {  // Predefined initialization value.
        EnumerableThreadLocal<test_align_parallel::TestValue> tls{test_align_parallel::TestValue{"predefined"}};
        SEQAN_ASSERT_EQ(tls._initValue.mMsg, "predefined");
    }
}

SEQAN_DEFINE_TEST(test_parallel_enumerable_thread_local_local)
{
    using namespace seqan;

    using TPair = std::pair<std::string, unsigned>;
    EnumerableThreadLocal<TPair> tls{TPair{"master", 1000}};

    size_t threadCount = std::min(std::thread::hardware_concurrency(), static_cast<unsigned>(4));
    test_align_parallel::run(tls, threadCount);

    SEQAN_ASSERT_EQ(tls._map.size(), threadCount);

    std::for_each(std::begin(tls._map), std::end(tls._map),
    [](auto const & mapValue)
    {
        auto const& val = mapValue.second;
        SEQAN_ASSERT_NEQ(val.first, "master");
        SEQAN_ASSERT_LEQ(val.second, 1000u);
        SEQAN_ASSERT_GEQ(val.second, 0u);
    });
}

SEQAN_DEFINE_TEST(test_parallel_enumerable_thread_local_enumerate)
{
    using namespace seqan;

    using TPair = std::pair<std::string, unsigned>;
    EnumerableThreadLocal<TPair> tls{TPair{"master", 1000}};

    size_t threadCount = std::min(std::thread::hardware_concurrency(), static_cast<unsigned>(4));
    test_align_parallel::run(tls, threadCount);

    SEQAN_ASSERT_EQ(tls._map.size(), threadCount);

    test_align_parallel::testEnumerate(tls);
    test_align_parallel::testEnumerateConst(tls);
}

SEQAN_DEFINE_TEST(test_parallel_enumerable_thread_local_combine_unary)
{
    using namespace seqan;

    using TPair = std::pair<std::string, unsigned>;
    EnumerableThreadLocal<TPair> tls{TPair{"master", 1000}};

    size_t threadCount = std::min(std::thread::hardware_concurrency(), static_cast<unsigned>(4));
    test_align_parallel::run(tls, threadCount);

    SEQAN_ASSERT_EQ(tls._map.size(), threadCount);

    unsigned count{0};
    combineEach(tls, [&](TPair const & p)
    {
        count += 1000 - p.second;
    });
    SEQAN_ASSERT_EQ(count, 1000u);
}

SEQAN_DEFINE_TEST(test_parallel_enumerable_thread_local_combine_binary)
{
    using namespace seqan;

    using TPair = std::pair<std::string, unsigned>;
    EnumerableThreadLocal<TPair> tls{TPair{"master", 1000}};

    size_t threadCount = std::min(std::thread::hardware_concurrency(), static_cast<unsigned>(4));
    test_align_parallel::run(tls, threadCount);

    SEQAN_ASSERT_EQ(tls._map.size(), threadCount);

    TPair count = combine(tls, [](TPair const & initial, TPair const & val)
                          {
                              return TPair{initial.first, initial.second + (1000 - val.second)};
                          });
    SEQAN_ASSERT_EQ(count.first, "");
    SEQAN_ASSERT_EQ(count.second, 1000u);
}
