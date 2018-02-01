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

#ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_H_
#define INCLUDE_SEQAN_ALIGN_PARALLEL_H_

// ============================================================================
// STD Prerequisites
// ============================================================================

#if SEQAN_DEBUG_ENABLED
#include <typeinfo>
#include <cxxabi.h>
#include <stdlib.h>
#endif

#include <type_traits>
#include <utility>
#include <vector>

// ============================================================================
// SeqAn Prerequisites
// ============================================================================

#include <seqan/basic.h>
#include <seqan/simd.h>
#include <seqan/align.h>
#include <seqan/parallel.h>

// ============================================================================
// Parallel DP Prerequisites and Adaptors
// ============================================================================

#include <seqan/align_parallel/dp_parallel_execution_policies.h>
#include <seqan/align_parallel/dp_traits.h>
#include <seqan/align_parallel/dp_kernel_adaptor.h>
#include <seqan/align_parallel/dp_settings.h>
#include <seqan/align_parallel/dp_parallel_scout.h>
#ifdef SEQAN_SIMD_ENABLED
#include <seqan/align_parallel/dp_parallel_scout_simd.h>
#endif

// ============================================================================
// Wavefront  Task
// ============================================================================

#include <seqan/align_parallel/wavefront_task_scheduler.h>
#include <seqan/align_parallel/wavefront_task_queue.h>
#include <seqan/align_parallel/wavefront_task_event.h>
#include <seqan/align_parallel/wavefront_task_util.h>
#include <seqan/align_parallel/wavefront_task.h>
#include <seqan/align_parallel/wavefront_task_executor.h>

// ============================================================================
// Wavefront Alignment Tasks
// ============================================================================

#include <seqan/align_parallel/wavefront_alignment_scheduler.h>
#include <seqan/align_parallel/wavefront_alignment_result.h>     // TODO(rrahn): rename! refactor!
#include <seqan/align_parallel/wavefront_alignment_thread_local_storage.h>
#include <seqan/align_parallel/wavefront_alignment_executor.h>
#include <seqan/align_parallel/wavefront_alignment_task.h>

// ============================================================================
// Interfaces
// ============================================================================

#include <seqan/align_parallel/async_wave_execution_interface.h>
#include <seqan/align_parallel/parallel_align_interface.h>

#endif  // #ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_H_
