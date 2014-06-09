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
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================

#include <seqan/basic.h>

using namespace seqan;

// ============================================================================
// Tests
// ============================================================================

// ----------------------------------------------------------------------------
// Test test_cuda_arch
// ----------------------------------------------------------------------------

SEQAN_GLOBAL void testCudaArch()
{
#if defined(__CUDA_ARCH__) && __CUDA_ARCH__ < 200
#error CUDA architecture 2.0 or higher is required.
#endif
}

SEQAN_DEFINE_TEST(test_cuda_arch)
{
    testCudaArch<<<1,1>>>();
    cudaDeviceSynchronize();
    SEQAN_ASSERT_EQ(cudaGetLastError(), cudaSuccess);
}

// ----------------------------------------------------------------------------
// Test test_cuda_assert
// ----------------------------------------------------------------------------

SEQAN_GLOBAL void testCudaAssert()
{
//    asm("trap;");
    SEQAN_ASSERT(false);
}

SEQAN_DEFINE_TEST(test_cuda_assert)
{
    testCudaAssert<<<1,1>>>();
    cudaDeviceSynchronize();
    SEQAN_ASSERT_NEQ(cudaGetLastError(), cudaSuccess);
}

// ============================================================================
// Register Tests
// ============================================================================

SEQAN_BEGIN_TESTSUITE(test_cuda_basic)
{
    // Call tests.
    SEQAN_CALL_TEST(test_cuda_arch);
    SEQAN_CALL_TEST(test_cuda_assert);
}
SEQAN_END_TESTSUITE
