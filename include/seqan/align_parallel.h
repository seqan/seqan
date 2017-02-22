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

#ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_H_
#define INCLUDE_SEQAN_ALIGN_PARALLEL_H_

// ============================================================================
// Prerequisites
// ============================================================================

#if SEQAN_DEBUG_ENABLED
#include <typeinfo>
#include <cxxabi.h>
#include <stdlib.h>
#endif

#include <type_traits>
#include <utility>
#include <vector>

#include <seqan/basic.h>
#include <seqan/align.h>
#include <seqan/parallel.h>

// ============================================================================
// Parallel DataStructures
// ============================================================================

//#include <seqan/align_parallel/parallel_thread_pool_base.h>
//#include <seqan/align_parallel/parallel_thread_pool_std.h>
//#include <seqan/align_parallel/parallel_task_pool_base.h>
//#include <seqan/align_parallel/parallel_task_pool_std.h>
#include <seqan/align_parallel/parallel_alignment_instance.h>
#include <seqan/align_parallel/parallel_alignment_scheduler.h>

// ============================================================================
// DP Task
// ============================================================================

#include <seqan/align_parallel/dp_parallel_base.h>
#include <seqan/align_parallel/dp_parallel_scout.h>

// Simd specific code.
#include <seqan/align_parallel/dp_parallel_scout_simd.h>
#include <seqan/align_parallel/dp_task_base_simd.h>

#include <seqan/align_parallel/dp_task_base.h>
#if defined(SEQAN_TBB)
#include <seqan/align_parallel/dp_task_tbb.h>
#endif
#if defined(_OPENMP)
#include <seqan/align_parallel/dp_task_omp.h>
#endif
#include <seqan/align_parallel/dp_task_std.h>

// ============================================================================
// Helper
// ============================================================================

#include <seqan/align_parallel/dp_trace_matrix_navigator_block_wise.h>

// ============================================================================
// Interfaces
// ============================================================================

#include <seqan/align_parallel/align_interface.h>
#include <seqan/align_parallel/align_instance.h>
#include <seqan/align_parallel/align_parallel_impl.h>

#endif  // #ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_H_