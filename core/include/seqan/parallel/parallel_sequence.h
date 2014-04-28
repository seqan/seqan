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
// Thread-safe / lock-free sequence operations.
// ==========================================================================

#ifndef SEQAN_PARALLEL_SEQUENCE_H_
#define SEQAN_PARALLEL_SEQUENCE_H_

namespace seqan {

// ============================================================================
// Class ReadWriteLock
// ============================================================================
// this lock augments a class by thread-safety as follows:
//  - supports multiple concurrent readers (possibly waiting for writer to finish)
//  - supports only a single writer at a time (possibly waiting for readers or other writers to finish)
//  - the writer has higher priority than all readers

// ----------------------------------------------------------------------------
// Class
// ----------------------------------------------------------------------------

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
    while (atomicCas(lock.writers, 0u, 1u) != 0) ;

    // wait until all readers are done
    while (lock.readers != 0) ;
}

// ----------------------------------------------------------------------------
// Function unlockWriting()
// ----------------------------------------------------------------------------

inline void
unlockWriting(ReadWriteLock &lock)
{
    lock.writers = 0;
}

// ============================================================================
// Class ConcurrentAppendValue
// ============================================================================

// ----------------------------------------------------------------------------
// Class
// ----------------------------------------------------------------------------

template <typename TString, typename TSpec = void>
struct ConcurrentAppendValue
{
    TString *       data_host;
    ReadWriteLock   lock;

    ConcurrentAppendValue() :
        data_host()
    {}

    ConcurrentAppendValue(TString & string)
    {
        setHost(*this, string);
    }
};

// ----------------------------------------------------------------------------
// Function host()
// ----------------------------------------------------------------------------

template <typename TString, typename TSpec>
inline TString &
host(ConcurrentAppendValue<TString, TSpec> & me)
{
    return *me.data_host;
}

template <typename TString, typename TSpec>
inline TString const &
host(ConcurrentAppendValue<TString, TSpec> const & me)
{
    return *me.data_host;
}

// ----------------------------------------------------------------------------
// Function setHost()
// ----------------------------------------------------------------------------

template <typename TString, typename TSpec>
inline void
setHost(ConcurrentAppendValue<TString, TSpec> & me, TString & string)
{
    me.data_host = &string;
}

// ----------------------------------------------------------------------------
// Function appendValue()
// ----------------------------------------------------------------------------

template <typename TString, typename TSpec, typename TValue, typename TExpand, typename TParallel>
inline void
appendValue(ConcurrentAppendValue<TString, TSpec> & me,
            TValue const & val,
            Tag<TExpand> const & expandTag,
            Tag<TParallel> const & parallelTag)
{
    appendValue(host(me), val, expandTag, parallelTag);
}

// ----------------------------------------------------------------------------
// Function appendValue(Parallel)
// ----------------------------------------------------------------------------

template <typename TString, typename TSpec, typename TValue, typename TExpand>
inline void
appendValue(ConcurrentAppendValue<TString, TSpec> & me,
            TValue const & val,
            Tag<TExpand> const & expandTag,
            Parallel)
{
    typedef typename Size<TString>::Type TSize;

    TString & string = host(me);

    while (true)
    {
        // try to append the value
        lockReading(me.lock);
        TSize newLen = _incLength(string, Parallel());
        if (newLen < capacity(string))
        {
            assignValue(string, newLen - 1, val);
            unlockReading(me.lock);
            break;
        }
        else
        {
            _decLength(string, Parallel());
        }
        unlockReading(me.lock);

        // try to extend capacity
        lockWriting(me.lock);
        TSize cap = capacity(string);
        if (cap == length(string))
            reserve(string, cap + 1, expandTag);
        unlockWriting(me.lock);
    }
}

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _incLength()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TParallel>
inline typename Size<String<TValue, Alloc<TSpec> > >::Type
_incLength(String<TValue, Alloc<TSpec> > & me, Tag<TParallel> const & tag)
{
    return atomicInc(me.data_end, tag) - begin(me, Standard());
}

// ----------------------------------------------------------------------------
// Function _decLength()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TParallel>
inline typename Size<String<TValue, Alloc<TSpec> > >::Type
_decLength(String<TValue, Alloc<TSpec> > & me, Tag<TParallel> const & tag)
{
    return atomicDec(me.data_end, tag) - begin(me, Standard());
}

// ----------------------------------------------------------------------------
// Function appendValue(Insist, Parallel); Atomic
// ----------------------------------------------------------------------------

template <typename TTargetValue, typename TTargetSpec, typename TValue>
inline void
appendValue(String<TTargetValue, TTargetSpec> & me, TValue const & _value, Insist, Parallel)
{
    valueConstruct(begin(me, Standard()) + _incLength(me, Parallel()) - 1, _value);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_PARALLEL_SEQUENCE_H_
