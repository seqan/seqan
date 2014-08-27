// ==========================================================================
//                                  ANISE
// ==========================================================================
// Copyright (c) 2006-2014, Knut Reinert, FU Berlin
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

#include "parallel_utils.h"

#include <thread>
#include <vector>

// --------------------------------------------------------------------------
// Class JobGenerator
// --------------------------------------------------------------------------

int JobGenerator::operator()()
{
    int lastSeen = nextJob;
    if (lastSeen >= numJobs)
        return -1;
    int newValue = lastSeen + 1;
    while (!std::atomic_compare_exchange_weak(&nextJob, &lastSeen, newValue))
    {
        lastSeen = nextJob;
        if (lastSeen >= numJobs)
            return -1;
        newValue = lastSeen + 1;
    }

    return lastSeen;
}


// --------------------------------------------------------------------------
// Function forkJoin()
// --------------------------------------------------------------------------

void forkJoin(int count, std::function<void()> fn)
{
    // Fork out threads.
    std::vector<std::thread> threads;
    for (int i = 0; i < count; ++i)
        threads.push_back(std::thread(fn));
    // Join threads.
    for (auto & thread : threads)
        thread.join();
}
