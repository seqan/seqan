// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2014, Knut Reinert, FU Berlin
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
// Counting iterator implementation.
// ==========================================================================

#ifndef SEQAN_BASIC_ITERATOR_COUNTING_H_
#define SEQAN_BASIC_ITERATOR_COUNTING_H_

namespace seqan {

// ============================================================================
// Tags
// ============================================================================

struct CountingIterator;

// ============================================================================
// Classes
// ============================================================================

template <typename TIncrementable>
class Iter<TIncrementable, CountingIterator>
{
public:
    TIncrementable data_position;

    // ------------------------------------------------------------------------
    // Constructors
    // ------------------------------------------------------------------------

    Iter(TIncrementable position = 0) :
        data_position(position)
    {}

    template <typename TOther>
    Iter(TOther position) :
        data_position(position)
    {}

    Iter(Iter const & other) :
        data_position(other.data_position)
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

template <typename TIncrementable>
struct Size<Iter<TIncrementable, CountingIterator> > : Size<TIncrementable> {};

template <typename TIncrementable>
struct Position<Iter<TIncrementable, CountingIterator> > : Position<TIncrementable> {};

template <typename TIncrementable>
struct Reference<Iter<TIncrementable, CountingIterator> > : Reference<TIncrementable> {};

template <typename TIncrementable>
struct Difference<Iter<TIncrementable, CountingIterator> > : Difference<TIncrementable> {};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function position()
// ----------------------------------------------------------------------------

///.Function.position.param.iterator.type:Spec.Position Iterator
///.Function.position.class:Spec.Position Iterator

template <typename TIncrementable>
inline typename Position<TIncrementable>::Type &
position(Iter<TIncrementable, CountingIterator> & me)
{
    return me.data_position;
}

template <typename TIncrementable>
inline typename Position<TIncrementable>::Type const &
position(Iter<TIncrementable, CountingIterator> const & me)
{
    return me.data_position;
}

// ----------------------------------------------------------------------------
// Function setPosition()
// ----------------------------------------------------------------------------

///.Function.setPosition.param.iterator.type:Spec.Position Iterator
///.Function.setPosition.class:Spec.Position Iterator

template <typename TIncrementable, typename TPosition>
inline void
setPosition(Iter<TIncrementable, CountingIterator> & me, TPosition position_)
{
    me.data_position = position_;
}

// ----------------------------------------------------------------------------
// Function value()
// ----------------------------------------------------------------------------

template <typename TIncrementable>
inline typename Reference<Iter<TIncrementable, CountingIterator> >::Type
value(Iter<TIncrementable, CountingIterator> & me)
{
    return position(me);
}

template <typename TIncrementable>
inline typename Reference<Iter<TIncrementable, CountingIterator> >::Type
value(Iter<TIncrementable, CountingIterator> const & me)
{
    return position(me);
}

// ----------------------------------------------------------------------------
// Function assignValue()
// ----------------------------------------------------------------------------

template <typename TIncrementable, typename TValue>
inline void
assignValue(Iter<TIncrementable, CountingIterator> & me, TValue _value)
{
    setPosition(me, _value);
}

// ----------------------------------------------------------------------------
// Function operator==()
// ----------------------------------------------------------------------------

template <typename TIncrementable>
inline bool
operator==(Iter<TIncrementable, CountingIterator> const & left,
           Iter<TIncrementable, CountingIterator> const & right)
{
    return position(left) == position(right);
}

// ----------------------------------------------------------------------------
// Function operator!=()
// ----------------------------------------------------------------------------

template <typename TIncrementable>
inline bool
operator!=(Iter<TIncrementable, CountingIterator> const & left,
           Iter<TIncrementable, CountingIterator> const & right)
{
    return position(left) != position(right);
}

// ----------------------------------------------------------------------------
// Function operator<()
// ----------------------------------------------------------------------------

template <typename TIncrementable>
inline bool
operator<(Iter<TIncrementable, CountingIterator> const & left,
          Iter<TIncrementable, CountingIterator> const & right)
{
    return position(left) < position(right);
}

template <typename TIncrementable>
inline bool
operator>(Iter<TIncrementable, CountingIterator> const & left,
          Iter<TIncrementable, CountingIterator> const & right)
{
    return position(left) > position(right);
}

// ----------------------------------------------------------------------------
// Function operator<=()
// ----------------------------------------------------------------------------

template <typename TIncrementable>
inline bool
operator<=(Iter<TIncrementable, CountingIterator> const & left,
           Iter<TIncrementable, CountingIterator> const & right)
{
    return position(left) <= position(right);
}

// ----------------------------------------------------------------------------
// Function operator>=()
// ----------------------------------------------------------------------------

template <typename TIncrementable>
inline bool
operator>=(Iter<TIncrementable, CountingIterator> const & left,
           Iter<TIncrementable, CountingIterator> const & right)
{
    return position(left) >= position(right);
}

// ----------------------------------------------------------------------------
// Function goNext()
// ----------------------------------------------------------------------------

template <typename TIncrementable>
inline void
goNext(Iter<TIncrementable, CountingIterator> & me)
{
    setPosition(me, position(me) + 1);
}

// ----------------------------------------------------------------------------
// Function goPrevious()
// ----------------------------------------------------------------------------

template <typename TIncrementable>
inline void
goPrevious(Iter<TIncrementable, CountingIterator> & me)
{
    setPosition(me, position(me) - 1);
}

// ----------------------------------------------------------------------------
// Function operator+()
// ----------------------------------------------------------------------------

template <typename TIncrementable, typename TIntegral>
inline Iter<TIncrementable, CountingIterator>
operator+(Iter<TIncrementable, CountingIterator> const & left,
          TIntegral right)
{
    return Iter<TIncrementable, CountingIterator>(container(left), position(left) + right);
}

// for <anonymous enum> types
template <typename TIncrementable>
inline Iter<TIncrementable, CountingIterator>
operator+(Iter<TIncrementable, CountingIterator> const & left,
          int right)
{
    return Iter<TIncrementable, CountingIterator>(container(left), position(left) + right);
}

template <typename TIncrementable, typename TIntegral>
inline Iter<TIncrementable, CountingIterator>
operator+(TIntegral left,
          Iter<TIncrementable, CountingIterator> const & right)
{
    return Iter<TIncrementable, CountingIterator>(container(right), position(right) + left);
}

// for <anonymous enum> types
template <typename TIncrementable>
inline Iter<TIncrementable, CountingIterator>
operator+(int left,
          Iter<TIncrementable, CountingIterator> const & right)
{
    return Iter<TIncrementable, CountingIterator>(container(right), position(right) + left);
}

// ----------------------------------------------------------------------------
// Function operator+=()
// ----------------------------------------------------------------------------

template <typename TIncrementable, typename TIntegral>
inline Iter<TIncrementable, CountingIterator> &
operator+=(Iter<TIncrementable, CountingIterator> & left,
           TIntegral right)
{
    setPosition(left, position(left) + right);
    return left;
}

// for <anonymous enum> types
template <typename TIncrementable>
inline Iter<TIncrementable, CountingIterator> &
operator+=(Iter<TIncrementable, CountingIterator> & left,
           int right)
{
    setPosition(left, position(left) + right);
    return left;
}

// ----------------------------------------------------------------------------
// Function operator-()
// ----------------------------------------------------------------------------

template <typename TIncrementable, typename TIntegral>
inline Iter<TIncrementable, CountingIterator>
operator-(Iter<TIncrementable, CountingIterator> const & left,
          TIntegral right)
{
    return Iter<TIncrementable, CountingIterator>(container(left), position(left) - right);
}

// for <anonymous enum> types
template <typename TIncrementable>
inline Iter<TIncrementable, CountingIterator>
operator-(Iter<TIncrementable, CountingIterator> const & left,
          int right)
{
    return Iter<TIncrementable, CountingIterator>(container(left), position(left) - right);
}

template <typename TIncrementable>
inline typename Difference<TIncrementable>::Type
operator-(Iter<TIncrementable, CountingIterator> const & left,
          Iter<TIncrementable, CountingIterator> const & right)
{
    return position(left) - position(right);
}

// ----------------------------------------------------------------------------
// Function operator-=()
// ----------------------------------------------------------------------------

template <typename TIncrementable, typename TIntegral>
inline Iter<TIncrementable, CountingIterator> &
operator-=(Iter<TIncrementable, CountingIterator> & left,
           TIntegral right)
{
    setPosition(left, position(left) - right);
    return left;
}

// for <anonymous enum> types
template <typename TIncrementable>
inline Iter<TIncrementable, CountingIterator> &
operator-=(Iter<TIncrementable, CountingIterator> & left,
           int right)
{
    setPosition(left, position(left) - right);
    return left;
}

// ----------------------------------------------------------------------------
// Function assign()
// ----------------------------------------------------------------------------

// Conversion assignment.
template <typename TIncrementable, typename TSource>
inline void
assign(Iter<TIncrementable, CountingIterator> & target,
       TSource const & source)
{
    setPosition(target, position(source));
}

}  // namespace seqan

#endif  // #ifndef SEQAN_BASIC_ITERATOR_COUNTING_H_
