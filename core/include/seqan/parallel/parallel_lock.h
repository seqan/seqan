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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================
// 2-level spinlock and corresponding scoped locks for each level
// ==========================================================================

#ifndef SEQAN_PARALLEL_PARALLEL_LOCK_H_
#define SEQAN_PARALLEL_PARALLEL_LOCK_H_

namespace seqan {

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class ReadWriteLock
// ----------------------------------------------------------------------------
// this lock augments a class by thread-safety as follows:
//  - supports multiple concurrent readers (possibly waiting for writer to finish)
//  - supports only a single writer at a time (possibly waiting for readers or other writers to finish)
//  - the writer has higher priority than all readers

struct ReadWriteLock
{
    volatile unsigned readers;
    volatile unsigned writers;

    ReadWriteLock() :
        readers(0),
        writers(0)
    {}
};

// ----------------------------------------------------------------------------
// Class ScopedReadLock
// ----------------------------------------------------------------------------

template <typename TLock = ReadWriteLock>
struct ScopedReadLock
{
    TLock &lock;

    ScopedReadLock(TLock &lock):
        lock(lock)
    {
        lockReading(lock);
    }

    ~ScopedReadLock()
    {
        unlockReading(lock);
    }
};

// ----------------------------------------------------------------------------
// Class ScopedWriteLock
// ----------------------------------------------------------------------------

template <typename TLock = ReadWriteLock>
struct ScopedWriteLock
{
    TLock &lock;

    ScopedWriteLock(TLock &lock):
        lock(lock)
    {
        lockWriting(lock);
    }

    ~ScopedWriteLock()
    {
        unlockWriting(lock);
    }
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function lockReading()
// ----------------------------------------------------------------------------

inline void
lockReading(ReadWriteLock &lock)
{
    do
    {
        // wait for the end of a write access
        while (lock.writers != 0) ;

        atomicInc(lock.readers);

        if (lock.writers == 0)
            break;

        // writer hasn't noticed us -> retry
        atomicDec(lock.readers);
    }
    while (true);
}

// ----------------------------------------------------------------------------
// Function unlockReading()
// ----------------------------------------------------------------------------

inline void
unlockReading(ReadWriteLock &lock)
{
    atomicDec(lock.readers);
}

// ----------------------------------------------------------------------------
// Function lockWriting()
// ----------------------------------------------------------------------------

inline void
lockWriting(ReadWriteLock &lock)
{
    // wait until we are the only writer
    while (atomicCas(lock.writers, 0u, 1u) != 0)
    {}

    // wait until all readers are done
    while (lock.readers != 0)
    {}
}

// ----------------------------------------------------------------------------
// Function unlockWriting()
// ----------------------------------------------------------------------------

inline void
unlockWriting(ReadWriteLock &lock)
{
    lock.writers = 0;
}

// ----------------------------------------------------------------------------
// Function empty()
// ----------------------------------------------------------------------------

inline bool
empty(ReadWriteLock &lock)
{
    return lock.readers == 0 && lock.writers == 0;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_PARALLEL_PARALLEL_LOCK_H_
