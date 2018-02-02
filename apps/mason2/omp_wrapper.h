// ==========================================================================
//                               omp_wrapper.h
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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Wrapper code for OpenMP things.  Currently only RAII stuff for locks.
// ==========================================================================

#ifdef _OPENMP
#include <omp.h>
#endif   // #ifdef _OPENMP

#ifndef APPS_MASON2_OMP_WRAPPER_H_
#define APPS_MASON2_OMP_WRAPPER_H_

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// OmpLock
// ----------------------------------------------------------------------------

// A wrapper for OpenMP locks.

class OmpLock
{
public:
    friend class OmpLockGuard;

    OmpLock()
    {
#ifdef _OPENMP
        omp_init_lock(&theLock);
#endif   // #ifdef _OPENMP
    }

    ~OmpLock()
    {
#ifdef _OPENMP
        omp_destroy_lock&(theLock);
#endif   // #ifdef _OPENMP
    }

private:
    void acquire()
    {
#ifdef _OPENMP
        omp_set_lock(&theLock);
#endif   // #ifdef _OPENMP
    }

    void release()
    {
#ifdef _OPENMP
        omp_unset_lock(&theLock);
#endif   // #ifdef _OPENMP
    }

#ifdef _OPENMP
    omp_lock_t theLock;
#endif   // #ifdef _OPENMP
};

// ----------------------------------------------------------------------------
// OmpLockGuard
// ----------------------------------------------------------------------------

// Implement RAII-style locking of OmpLock objects.

class OmpLockGuard
{
public:
    OmpLockGuard(OmpLock & lock) : lock(lock)
    {
        lock.acquire();
    }

    ~OmpLockGuard()
    {
        lock.release();
    }

private:
    OmpLock & lock;
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

#endif  // #ifndef APPS_MASON2_OMP_WRAPPER_H_
