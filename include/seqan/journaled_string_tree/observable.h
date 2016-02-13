// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
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

#ifndef INCLUDE_SEQAN_JOURNALED_STRING_TREE_OBSERVABLE_H_
#define INCLUDE_SEQAN_JOURNALED_STRING_TREE_OBSERVABLE_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// Type alias used to hide std::tuple in SeqAn namespace.
template <typename... TObservers>
using ObserverList = std::tuple<TObservers...>;

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction LENGTH
// ----------------------------------------------------------------------------

template <typename... TObservers>
struct LENGTH<ObserverList<TObservers...> >
{
    static constexpr unsigned VALUE = sizeof...(TObservers);
};

// ----------------------------------------------------------------------------
// Metafunction Element
// ----------------------------------------------------------------------------

template <unsigned INDEX, typename TObject>
struct Element;

template <unsigned INDEX, typename TFirst, typename... TObservers>
struct Element<INDEX, ObserverList<TFirst, TObservers...> >
{
    static_assert(INDEX < LENGTH<ObserverList<TFirst, TObservers...> >::VALUE, "Trying to access element behind the end.");
    typedef typename Element<INDEX - 1, ObserverList<TObservers...> >::Type Type;
};

template <typename TFirst, typename... TObservers>
struct Element<0, ObserverList<TFirst, TObservers...> >
{
    typedef TFirst Type;
};

// ----------------------------------------------------------------------------
// Metafunction NotifyObserver
// ----------------------------------------------------------------------------

template <unsigned INDEX>
struct NotifyObserver
{
    template <typename... TObserver, typename TTag>
    inline void operator()(ObserverList<TObserver...> & subject, TTag const & /*tag*/)
    {
//        std::get<INDEX>(subject.observers)(TTag());
        notify(std::get<INDEX>(subject), TTag());
        NotifyObserver<INDEX - 1>{}(subject, TTag());
    }
};

template <>
struct NotifyObserver<0>
{

    template <typename... TObserver, typename TTag>
    inline void operator()(ObserverList<TObserver...> & subject, TTag const & /*tag*/)
    {
        notify(std::get<0>(subject), TTag());
    }
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function length()
// ----------------------------------------------------------------------------

template <typename... TObserverTypes>
inline constexpr unsigned
length(ObserverList<TObserverTypes...> const & /*subject*/)
{
    return LENGTH<ObserverList<TObserverTypes...> >::VALUE;
}

// ----------------------------------------------------------------------------
// Function setObserver()
// ----------------------------------------------------------------------------

template <unsigned INDEX, typename TObserver>
inline void
setObserver(ObserverList<> & /*subject*/,
            TObserver & /*observer*/)
{
    // no-op;
}

template <unsigned INDEX, typename... TObserverTypes, typename TObserver>
inline void
setObserver(ObserverList<TObserverTypes...> & subject,
            TObserver && observer)
{
    static_assert(INDEX < LENGTH<ObserverList<TObserverTypes...> >::VALUE, "Trying to access element behind the end.");
    std::get<INDEX>(subject) = std::forward<TObserver>(observer);
}

// ----------------------------------------------------------------------------
// Function notify()
// ----------------------------------------------------------------------------

template <typename TEventTag>
inline void
notify(ObserverList<> & /*subject*/,
       TEventTag const & /*event tag*/)
{
    // no-op;
}

template <typename... TObserverTypes, typename TEventTag>
inline void
notify(ObserverList<TObserverTypes...> & subject,
       TEventTag const & /*event tag*/)
{
    NotifyObserver<LENGTH<ObserverList<TObserverTypes...> >::VALUE - 1>()(subject, TEventTag());
}

// ----------------------------------------------------------------------------
// Function makeObserverList
// ----------------------------------------------------------------------------

template <typename... TObservers>
inline ObserverList<TObservers...>
makeObserverList(TObservers&&... observers)
{
    return std::forward_as_tuple(observers...);
}

}

#endif  // #ifndef INCLUDE_SEQAN_JOURNALED_STRING_TREE_OBSERVABLE_H_
