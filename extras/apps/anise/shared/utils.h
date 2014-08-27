// ==========================================================================
//                                   ANISE
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
// Header with utility algorithms and iterators, no(!) dependencies.
// ==========================================================================

#ifndef SANDBOX_HERBARIUM_APPS_ANISE_UTILS_H_
#define SANDBOX_HERBARIUM_APPS_ANISE_UTILS_H_

// ----------------------------------------------------------------------------
// Class adjacent_iterator
// ----------------------------------------------------------------------------

template <typename FwdIt> class adjacent_iterator {
public:
    adjacent_iterator(FwdIt first, FwdIt last)
        : m_first(first), m_next(first == last ? first : std::next(first)) { }

    bool operator!=(const adjacent_iterator& other) const {
        return m_next != other.m_next; // NOT m_first!
    }

    adjacent_iterator& operator++() {
        ++m_first;
        ++m_next;
        return *this;
    }

    typedef typename std::iterator_traits<FwdIt>::reference Ref;
    typedef std::pair<Ref, Ref> Pair;

    Pair operator*() const {
        return Pair(*m_first, *m_next); // NOT std::make_pair()!
    }

private:
    FwdIt m_first;
    FwdIt m_next;
};

// ----------------------------------------------------------------------------
// Class adjacent_range
// ----------------------------------------------------------------------------

template <typename FwdIt> class adjacent_range {
public:
    adjacent_range(FwdIt first, FwdIt last)
        : m_first(first), m_last(last) { }

    adjacent_iterator<FwdIt> begin() const {
        return adjacent_iterator<FwdIt>(m_first, m_last);
    }

    adjacent_iterator<FwdIt> end() const {
        return adjacent_iterator<FwdIt>(m_last, m_last);
    }

private:
    FwdIt m_first;
    FwdIt m_last;
};

// ----------------------------------------------------------------------------
// Function make_adjacent_range()
// ----------------------------------------------------------------------------

template <typename C> auto make_adjacent_range(C& c) -> adjacent_range<decltype(c.begin())> {
    return adjacent_range<decltype(c.begin())>(c.begin(), c.end());
}

// ----------------------------------------------------------------------------
// Class pairwise_iterator
// ----------------------------------------------------------------------------

template <typename FwdIt> class pairwise_iterator {
public:
    pairwise_iterator(FwdIt first, FwdIt last)
        : m_first(first), m_next(first == last ? first : std::next(first)) { }

    bool operator!=(const pairwise_iterator& other) const {
        return m_first != other.m_first; // NOT m_next!
    }

    pairwise_iterator& operator++() {
        ++m_first;
        ++m_first;
        ++m_next;
        ++m_next;
        return *this;
    }

    typedef typename std::iterator_traits<FwdIt>::reference InnerRef;
    typedef typename std::remove_reference<InnerRef>::type Value;
    typedef std::reference_wrapper<Value> Ref;
    typedef std::pair<Ref, Ref> Pair;

    Pair operator*() const {
        return Pair(*m_first, *m_next); // NOT std::make_pair()!
    }

private:
    FwdIt m_first;
    FwdIt m_next;
};

// ----------------------------------------------------------------------------
// Class pairwise_range
// ----------------------------------------------------------------------------

template <typename FwdIt> class pairwise_range {
public:
    pairwise_range(FwdIt first, FwdIt last)
        : m_first(first), m_last(last) { }

    pairwise_iterator<FwdIt> begin() const {
        return pairwise_iterator<FwdIt>(m_first, m_last);
    }

    pairwise_iterator<FwdIt> begin() {
        return pairwise_iterator<FwdIt>(m_first, m_last);
    }

    pairwise_iterator<FwdIt> end() const {
        return pairwise_iterator<FwdIt>(m_last, m_last);
    }

    pairwise_iterator<FwdIt> end() {
        return pairwise_iterator<FwdIt>(m_last, m_last);
    }

private:
    FwdIt m_first;
    FwdIt m_last;
};

// ----------------------------------------------------------------------------
// Function make_pairwise_range()
// ----------------------------------------------------------------------------

template <typename C> auto make_pairwise_range(C& c) -> pairwise_range<decltype(c.begin())> {
    return pairwise_range<decltype(c.begin())>(c.begin(), c.end());
}

// ----------------------------------------------------------------------------
// Class simultaneous_iterator
// ----------------------------------------------------------------------------

template <typename FwdIt> class simultaneous_iterator {
public:
    simultaneous_iterator(FwdIt left_first, FwdIt right_first)
        : m_left_first(left_first), m_right_first(right_first) {}

    bool operator!=(const simultaneous_iterator& other) const {
        return std::make_pair(m_left_first, m_right_first) != std::make_pair(other.m_left_first, other.m_right_first);
    }

    simultaneous_iterator& operator++() {
        ++m_left_first;
        ++m_right_first;
        return *this;
    }

    typedef typename std::iterator_traits<FwdIt>::reference Ref;
    typedef std::pair<Ref, Ref> Pair;

    Pair operator*() const {
        return Pair(*m_left_first, *m_right_first); // NOT std::make_pair()!
    }

private:
    FwdIt m_left_first;
    FwdIt m_right_first;
};

// ----------------------------------------------------------------------------
// Class simultaneous_range
// ----------------------------------------------------------------------------

template <typename FwdIt> class simultaneous_range {
public:
    simultaneous_range(FwdIt left_first, FwdIt right_first, FwdIt left_last, FwdIt right_last)
        : m_left_first(left_first), m_right_first(right_first), m_left_last(left_last), m_right_last(right_last)
    {}

    simultaneous_iterator<FwdIt> begin() const {
        return simultaneous_iterator<FwdIt>(m_left_first, m_right_first);
    }

    simultaneous_iterator<FwdIt> end() const {
        return simultaneous_iterator<FwdIt>(m_left_last, m_right_last);
    }

private:
    FwdIt m_left_first;
    FwdIt m_right_first;
    FwdIt m_left_last;
    FwdIt m_right_last;
};

// ----------------------------------------------------------------------------
// Function make_simultaneous_range()
// ----------------------------------------------------------------------------

template <typename C> auto make_simultaneous_range(C& c, C& d) -> simultaneous_range<decltype(c.begin())> {
    return simultaneous_range<decltype(c.begin())>(c.begin(), d.begin(), c.end(), d.end());
}

// ----------------------------------------------------------------------------
// Function transform_if()
// ----------------------------------------------------------------------------

// Utility function.
//
// TODO(holtgrew): Should go into utils header.

template <class InputIterator, class OutputIterator, class UnaryOperator, class Condition>
OutputIterator transform_if(InputIterator first1, InputIterator last1,
                            OutputIterator result, UnaryOperator op, Condition cond)
{
    while (first1 != last1)
    {
        if (cond(*first1))
        {
            *result = op(*first1);
            ++result;
        }
        ++first1;
    }
    return result;
}

#endif  // #ifndef SANDBOX_HERBARIUM_APPS_ANISE_UTILS_H_
