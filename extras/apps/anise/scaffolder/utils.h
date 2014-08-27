// ==========================================================================
//                                  ANISE
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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

#ifndef SANDBOX_HERBARIUM_APPS_ANISE_SCAFFOLDER_UTILS_H_
#define SANDBOX_HERBARIUM_APPS_ANISE_SCAFFOLDER_UTILS_H_

#include <seqan/basic.h>

namespace scaffolder {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function er()
// ----------------------------------------------------------------------------

template <typename Iter, typename Eq = std::equal_to<typename seqan::Value<Iter>::Type> >
class eqrange
{
    Iter rangeFirst, rangeLast, last;
    Eq eq;

    void adv()
    {
        rangeLast = rangeFirst;
        while (rangeLast != last && eq(*rangeFirst, *rangeLast)) { ++rangeLast; }
    }

public:
    template <typename Iter2>
    eqrange(typename std::enable_if<std::is_default_constructible<Eq>::value, Iter2>::type first, Iter last) :
            rangeFirst(first), rangeLast(first), last(last), eq()
    {
        adv();
    }

    eqrange(Iter first, Iter last, Eq eq) :
            rangeFirst(first), rangeLast(first), last(last), eq(eq)
    {
        adv();
    }

    Iter begin() { return rangeFirst; }
    Iter end()   { return rangeLast; }

    eqrange & operator++(/*prefix*/) {
        rangeFirst = rangeLast;
        adv();
        return *this;
    }

    eqrange operator++(int /*postfix*/) {
        eqrange tmp(*this);
        ++this;
        return tmp;
    }

    bool operator==(eqrange const & rhs) const
    {
        return rangeFirst == rhs.rangeFirst && rangeLast == rhs.rangeLast && last == rhs.last;
    }

    bool operator!=(eqrange const & rhs) const
    {
        return !(*this == rhs);
    }

    eqrange make_end() const
    {
        return eqrange(last, last, eq);
    }
};

template <typename Iter>
eqrange<Iter> er(Iter first, Iter last)
{
    typedef std::equal_to<typename seqan::Value<Iter>::Type> Cmp;
    return eqrange<Iter, Cmp>(first, last, Cmp());
}

template <typename Iter, typename Cmp>
eqrange<Iter, Cmp> er(Iter first, Iter last, Cmp cmp)
{
    return eqrange<Iter, Cmp>(first, last, cmp);
}

// ----------------------------------------------------------------------------
// Function enumerate()
// ----------------------------------------------------------------------------

template <typename T>
struct iterator_extractor { typedef typename T::iterator type; };

template <typename T>
struct iterator_extractor<T const> { typedef typename T::const_iterator type; };


template <typename T>
class Indexer {
public:
    class iterator {
        typedef typename iterator_extractor<T>::type inner_iterator;

        typedef typename std::iterator_traits<inner_iterator>::reference inner_reference;
    public:
        typedef std::pair<size_t, inner_reference> reference;

        iterator(inner_iterator it): _pos(0), _it(it) {}

        reference operator*() const { return reference(_pos, *_it); }

        iterator& operator++() { ++_pos; ++_it; return *this; }
        iterator operator++(int) { iterator tmp(*this); ++*this; return tmp; }

        bool operator==(iterator const& it) const { return _it == it._it; }
        bool operator!=(iterator const& it) const { return !(*this == it); }

    private:
        size_t _pos;
        inner_iterator _it;
    };

    Indexer(T& t): _container(t) {}

    iterator begin() const { return iterator(_container.begin()); }
    iterator end() const { return iterator(_container.end()); }

private:
    T& _container;
}; // class Indexer

template <typename T>
Indexer<T> enumerate(T& t) { return Indexer<T>(t); }

}  // namespace scaffolder

#endif  // #ifndef SANDBOX_HERBARIUM_APPS_ANISE_SCAFFOLDER_UTILS_H_
