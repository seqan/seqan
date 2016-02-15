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

#ifndef INCLUDE_SEQAN_JOURNALED_STRING_TREE_STACK_OBSERVER_H_
#define INCLUDE_SEQAN_JOURNALED_STRING_TREE_STACK_OBSERVER_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

template <typename TObject>
struct ObservedValue
{
    using Type = Nothing;
};

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Tag PushEvent
// ----------------------------------------------------------------------------

struct PushEvent_;
typedef Tag<PushEvent_> PushEvent;

// ----------------------------------------------------------------------------
// Tag PopEvent
// ----------------------------------------------------------------------------

struct PopEvent_;
typedef Tag<PopEvent_> PopEvent;

// ----------------------------------------------------------------------------
// Class StackObserver
// ----------------------------------------------------------------------------

template <typename TObject>
class StackObserver
{
public:
    using TValue = typename ObservedValue<TObject>::Type;
    using TStack = String<TValue, Block<> >;

    TObject & _obj;
    TStack _data;

    StackObserver(TObject & obj) : _obj(obj)
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function notify(); PushEvent
// ----------------------------------------------------------------------------

template <typename TObject>
inline void
notify(StackObserver<TObject> & me,
       PushEvent const & /*tag*/)
{
    appendValue(me._data, getObservedValue(me._obj));
}

// ----------------------------------------------------------------------------
// Function notify(); PopEvent
// ----------------------------------------------------------------------------

template <typename TObject>
inline void
notify(StackObserver<TObject> & me,
       PopEvent const & /*tag*/)
{
    setObservedValue(me._obj, std::move(back(me._data)));
    eraseBack(me._data);
}

}

#endif  // #ifndef INCLUDE_SEQAN_JOURNALED_STRING_TREE_STACK_OBSERVER_H_
