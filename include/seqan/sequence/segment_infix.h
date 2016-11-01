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
// Author: Andreas Gogol-Doering <andreas.doering@mdc-berlin.de>
// ==========================================================================
// Implementation of the Infix Segment specialization.
// ==========================================================================

#ifndef SEQAN_HEADER_SEGMENT_INFIX_H
#define SEQAN_HEADER_SEGMENT_INFIX_H

namespace seqan
{

//////////////////////////////////////////////////////////////////////////////
// InfixSegment
//////////////////////////////////////////////////////////////////////////////

/*!
 * @class InfixSegment Infix Segment
 * @extends Segment
 * @headerfile <seqan/sequence.h>
 * @brief An infix of a sequence.
 *
 * @signature template <typename THost>
 *            class Segment<THost, InfixSegment>;
 *
 * @tparam THost The underlying @link ContainerConcept sequence @endlink type.
 *
 * @see SuffixSegment
 * @see PrefixSegment
 * @aka substring
 */

template <typename THost_>
class Segment<THost_, InfixSegment>
{
public:
    typedef typename Host<Segment>::Type THost;

    typename Pointer_<THost>::Type data_host;
    typename Position<THost>::Type data_begin_position;
    typename Position<THost>::Type data_end_position;


//____________________________________________________________________________

public:
    // Check member variables with assertions.  This is called in the
    // constructors.
    inline void _checkMemberVariables() const {
        SEQAN_ASSERT_LEQ(data_begin_position, data_end_position);
    }

   
    Segment():
        data_host(),
        data_begin_position(0),
        data_end_position(0)
    {
        _checkMemberVariables();
    }

   
    Segment(typename Parameter_<THost>::Type _host):
        data_host(_toPointer(_host)),
        data_begin_position(0),
        data_end_position(length(value(data_host)))
    {
        _checkMemberVariables();
    }

   
    Segment(typename Parameter_<THost>::Type _host, typename Position<THost>::Type _begin_index, typename Position<THost>::Type _end_index):
        data_host(_toPointer(_host)),
        data_begin_position(_begin_index),
        data_end_position(_end_index)
    {
        _checkMemberVariables();
    }
/*
    Segment(typename Parameter_<THost>::Type _host, typename Iterator<THost, Rooted>::Type _begin, typename Iterator<THost, Rooted>::Type _end):
        data_host(_toPointer(_host)),
        data_begin_position(position(_begin)),
        data_end_position(position(_end))
    {
    }
*/

   
    Segment(typename Parameter_<THost>::Type _host, typename Iterator<THost, Standard>::Type _begin, typename Iterator<THost, Standard>::Type _end):
        data_host(_toPointer(_host)),
        data_begin_position(position(_begin, _host)),
        data_end_position(position(_end, _host))
    {
        _checkMemberVariables();
    }

    template <typename THost2, typename TSpec2>
   
    Segment(Segment<THost2, TSpec2> const & _other):
        data_host(_toPointer(host(_other))),
        data_begin_position(beginPosition(_other)),
        data_end_position(endPosition(_other))
    {
        _checkMemberVariables();
    }

   
    ~ Segment()
    {
    }

    inline Segment &
    operator = (Segment const & source)
    {
        assign(*this, source);
        return *this;
    }

    template<typename T> explicit operator T () const
    {
        T temp_copy;
        assign(temp_copy, *this);
        return temp_copy;
    }

//____________________________________________________________________________

public:


//____________________________________________________________________________

    template <typename TPos>
    inline typename Reference<Segment>::Type
    operator [] (TPos pos)
    {
        return value(*this, pos);
    }

    template <typename TPos>
    inline typename Reference<Segment const>::Type
    operator [] (TPos pos) const
    {
        return value(*this, pos);
    }

};
//////////////////////////////////////////////////////////////////////////////


// template <typename THost>
// inline void
// clear(Segment<THost, InfixSegment> & target)
// {
//     replace(host(target), beginPosition(target), endPosition(target), "");
//     setEndPosition(target, beginPosition(target));
// }

///Function.host.param.object.type:Class.Segment

template <typename THost_>
inline typename Parameter_<THost_>::Type
host(Segment<THost_, InfixSegment> & me)
{
    return _toParameter<THost_>(me.data_host);
}

template <typename THost_>
inline typename Parameter_<THost_>::Type
host(Segment<THost_, InfixSegment> const & me)
{
    return _toParameter<THost_>(me.data_host);
}


//____________________________________________________________________________

template <typename THost_>
inline typename Iterator<Segment<THost_, InfixSegment>, Standard>::Type
begin(Segment<THost_, InfixSegment> & me,
    Standard)
{
    return begin(host(me), Standard()) + me.data_begin_position;
}
template <typename THost_>
inline typename Iterator<Segment<THost_, InfixSegment> const, Standard>::Type
begin(Segment<THost_, InfixSegment> const & me,
    Standard)
{
    return begin(host(me), Standard()) + me.data_begin_position;
}

//____________________________________________________________________________

template <typename THost_>
inline typename Position<Segment<THost_, InfixSegment> >::Type
beginPosition(Segment<THost_, InfixSegment> & me)
{
    return me.data_begin_position;
}
template <typename THost_>
inline typename Position<Segment<THost_, InfixSegment> const>::Type
beginPosition(Segment<THost_, InfixSegment> const & me)
{
    return me.data_begin_position;
}

//____________________________________________________________________________

template <typename THost_, typename TIterator>
inline void
setBegin(Segment<THost_, InfixSegment> & me, TIterator new_begin)
{
    me.data_begin_position = new_begin - begin(host(me));//, Standard());
}


//____________________________________________________________________________

template <typename THost_, typename TPosition>
inline void
setBeginPosition(Segment<THost_, InfixSegment> & me, TPosition new_begin)
{
    me.data_begin_position = new_begin;
}

//____________________________________________________________________________

template <typename THost_>
inline typename Iterator<Segment<THost_, InfixSegment>, Standard>::Type
end(Segment<THost_, InfixSegment> & me,
    Standard)
{
    return begin(host(me), Standard()) + me.data_end_position;
}
template <typename THost_>
inline typename Iterator<Segment<THost_, InfixSegment> const, Standard>::Type
end(Segment<THost_, InfixSegment> const & me,
    Standard)
{
    return begin(host(me), Standard()) + me.data_end_position;
}

//____________________________________________________________________________

template <typename THost_>
inline typename Position<Segment<THost_, InfixSegment> >::Type
endPosition(Segment<THost_, InfixSegment> & me)
{
    return me.data_end_position;
}
template <typename THost_>
inline typename Position<Segment<THost_, InfixSegment> >::Type
endPosition(Segment<THost_, InfixSegment> const & me)
{
    return me.data_end_position;
}

//____________________________________________________________________________

template <typename THost_, typename TIterator>
inline void
setEnd(Segment<THost_, InfixSegment> & me, TIterator new_end)
{
    // me.data_end_position = new_end - begin(host(me)); //, Standard());
    me.data_end_position = new_end - TIterator(begin(host(me)));
}

//____________________________________________________________________________


template <typename THost_, typename TPosition>
inline void
setEndPosition(Segment<THost_, InfixSegment> & me, TPosition new_end)
{
    me.data_end_position = new_end;
}

//____________________________________________________________________________

template <typename THost_>
inline void
_setLength(
    Segment<THost_, InfixSegment> & me,
    typename Size<THost_>::Type new_length)
{
    me.data_end_position = me.data_begin_position + new_length;
}


//____________________________________________________________________________

template <typename THost_>
inline void
setHost(Segment<THost_, InfixSegment> & me, typename Parameter_<THost_>::Type _host)
{
    me.data_host = _toPointer(_host);
}

template <typename THost_>
inline void
setHost(Segment<THost_ const, InfixSegment> & me, typename Parameter_<THost_>::Type _host)
{
    me.data_host = _toPointer(_host);
}

//////////////////////////////////////////////////////////////////////////////

template <typename THost>
struct Infix
{
    typedef Segment<THost, InfixSegment> Type;
};

template <typename THost, typename TSpec>
struct Infix< Segment<THost, TSpec> >
{
    typedef Segment<THost, InfixSegment> Type;
};

template <typename THost, typename TSpec>
struct Infix< Segment<THost, TSpec> const >:
    Infix< Segment<THost, TSpec> > {};

template <typename THost>
struct Infix<THost &>:
    Infix<THost> {};

// ----------------------------------------------------------------------------
// Metafunction InfixOnValue
// ----------------------------------------------------------------------------

// The default implementation returns Infix<T>::Type.
template <typename T>
struct InfixOnValue :
    Infix<T> {};

//////////////////////////////////////////////////////////////////////////////

template <typename THost, typename TPosition1, typename TPosition2>
inline void
set(Segment<THost, InfixSegment> & me,
    THost & host_,
    TPosition1 begin_,
    TPosition2 end_)
{
    setHost(me, host_);
    setBeginPosition(me, begin_);
    setEndPosition(me, end_);
}
//____________________________________________________________________________

template <typename THost>
inline void
set(Segment<THost, InfixSegment> & me,
    THost & host_)
{
    setHost(me, host_);
    setBegin(me, begin(host_, Standard()));
    setEnd(me, end(host_, Standard()));
}
template <typename THost>
inline void
set(Segment<THost, InfixSegment> & me,
    THost const & host_)
{
    setHost(me, host_);
    setBegin(me, begin(host_, Standard()));
    setEnd(me, end(host_, Standard()));
}

//____________________________________________________________________________

template <typename THost, typename TSpec>
inline void
set(Segment<THost, InfixSegment> & me,
    Segment<THost, TSpec> & source)
{
    setHost(me, host(source));
    setBeginPosition(me, beginPosition(source));
    setEndPosition(me, endPosition(source));
}
template <typename THost, typename TSpec>
inline void
set(Segment<THost const, InfixSegment> & me,
    Segment<THost, TSpec> & source)
{
    setHost(me, host(source));
    setBeginPosition(me, beginPosition(source));
    setEndPosition(me, endPosition(source));
}
template <typename THost, typename TSpec>
inline void
set(Segment<THost, InfixSegment> & me,
    Segment<THost, TSpec> const & source)
{
    setHost(me, host(source));
    setBeginPosition(me, beginPosition(source));
    setEndPosition(me, endPosition(source));
}
template <typename THost, typename TSpec>
inline void
set(Segment<THost const, InfixSegment> & me,
    Segment<THost, TSpec> const & source)
{
    setHost(me, host(source));
    setBeginPosition(me, beginPosition(source));
    setEndPosition(me, endPosition(source));
}

//////////////////////////////////////////////////////////////////////////////

template <typename THost>
inline bool
atBegin(Segment<THost, InfixSegment> & segment)
{
    return (beginPosition(segment) == endPosition(segment));
}
template <typename THost>
inline bool
atBegin(Segment<THost, InfixSegment> const & segment)
{
    return (beginPosition(segment) == endPosition(segment));
}

//////////////////////////////////////////////////////////////////////////////

template <typename THost>
inline bool
atEnd(Segment<THost, InfixSegment> & segment)
{
    return (endPosition(segment) - beginPosition(segment)) > length(host(segment));
}
template <typename THost>
inline bool
atEnd(Segment<THost, InfixSegment> const & segment)
{
    return (endPosition(segment) - beginPosition(segment)) > length(host(segment));
}


//////////////////////////////////////////////////////////////////////////////

template <typename THost>
inline void
goBegin(Segment<THost, InfixSegment> & segment)
{
    setBeginPosition(segment, 0);
    setEndPosition(segment, 1);
}
template <typename THost, typename THost2>
inline void
goBegin(Segment<THost, InfixSegment> & segment,
        THost2 &)
{
    goBegin(segment);
}
template <typename THost, typename THost2>
inline void
goBegin(Segment<THost, InfixSegment> & segment,
        THost2 const &)
{
    goBegin(segment);
}

//////////////////////////////////////////////////////////////////////////////


template <typename THost>
inline void
goEnd(Segment<THost, InfixSegment> & segment)
{
    setBeginPosition(segment, 0);
    setEndPosition(segment, length(host(segment)));
}
template <typename THost, typename THost2>
inline void
goEnd(Segment<THost, InfixSegment> & segment,
      THost2 &)
{
    goEnd(segment);
}
template <typename THost, typename THost2>
inline void
goEnd(Segment<THost, InfixSegment> & segment,
      THost2 const &)
{
    goEnd(segment);
}

//////////////////////////////////////////////////////////////////////////////

template <typename THost>
inline Segment<THost, InfixSegment> &
operator ++(Segment<THost, InfixSegment> & segment)
{
    if (endPosition(segment) == length(host(segment)))
    {
        setEndPosition(segment, endPosition(segment) - beginPosition(segment) + 1);
        setBeginPosition(segment, 0);
    }
    else
    {
        setBeginPosition(segment, beginPosition(segment) + 1);
        setEndPosition(segment, endPosition(segment) + 1);
    }
    return segment;
}

//////////////////////////////////////////////////////////////////////////////

template <typename THost>
inline Segment<THost, InfixSegment> &
operator --(Segment<THost, InfixSegment> & segment)
{
    if (!beginPosition(segment))
    {
        typename Size<THost>::Type host_length = length(host(segment));

        setBeginPosition(segment, host_length - endPosition(segment) + beginPosition(segment) + 1);
        setEndPosition(segment, host_length);
    }
    else
    {
        setBeginPosition(segment, beginPosition(segment) - 1);
        setEndPosition(segment, endPosition(segment) - 1);
    }
    return segment;
}

//////////////////////////////////////////////////////////////////////////////

template <typename THost, typename TSpec, typename TPos>
inline typename Reference< Segment<THost, TSpec> >::Type
value(Segment<THost, TSpec> & me,
      TPos pos)
{
    SEQAN_ASSERT_LT_MSG(pos, static_cast<TPos>(length(me)), "Trying to acces an element behind the last one!");
    return *(begin(me, Standard()) + pos);
}

template <typename THost, typename TSpec, typename TPos>
inline typename Reference< Segment<THost, TSpec> const >::Type
value(Segment<THost, TSpec> const & me,
      TPos pos)
{
    SEQAN_ASSERT_LT_MSG(pos, static_cast<TPos>(length(me)), "Trying to acces an element behind the last one!");
    return *(begin(me, Standard()) + pos);
}

//////////////////////////////////////////////////////////////////////////////

template <typename T, typename TPosBegin, typename TPosEnd>
inline typename Infix<T>::Type
infix(T & t, TPosBegin pos_begin, TPosEnd pos_end)
{
    return typename Infix<T>::Type(t, pos_begin, pos_end);
}

template <typename T, typename TPosBegin, typename TPosEnd>
inline typename Infix<T *>::Type
infix(T * t, TPosBegin pos_begin, TPosEnd pos_end)
{
    return typename Infix<T *>::Type (t, pos_begin, pos_end);
}

template <typename T, typename TSpec, typename TPosBegin, typename TPosEnd>
inline typename Infix<Segment<T, TSpec> >::Type
infix(Segment<T, TSpec> & t, TPosBegin pos_begin, TPosEnd pos_end)
{
    return typename Infix<Segment<T, TSpec> >::Type (
        host(t),
        beginPosition(t) + pos_begin,
        beginPosition(t) + pos_end);
}

template <typename T, typename TSpec, typename TPosBegin, typename TPosEnd>
inline typename Infix<Segment<T, TSpec> const>::Type
infix(Segment<T, TSpec> const & t, TPosBegin pos_begin, TPosEnd pos_end)
{
    return typename Infix<Segment<T, TSpec> const>::Type (
        host(t),
        beginPosition(t) + pos_begin,
        beginPosition(t) + pos_end);
}

// infix() with iterators

template <typename T, typename TSpec, typename TIterSpec>
inline typename Infix<Segment<T, TSpec> >::Type
infix(Segment<T, TSpec> & t,
      Iter<Segment<T, TSpec>, TIterSpec> const & iterBegin,
      Iter<Segment<T, TSpec>, TIterSpec> const & iterEnd)
{
    return typename Infix<Segment<T, TSpec> >::Type (
        host(t),
        iterBegin,
        iterEnd);
}

template <typename T, typename TSpec, typename TIterSpec>
inline typename Infix<Segment<T, TSpec> const>::Type
infix(Segment<T, TSpec> const & t,
      Iter<Segment<T, TSpec> const, TIterSpec> const & iterBegin,
      Iter<Segment<T, TSpec> const, TIterSpec> const & iterEnd)
{
    return typename Infix<Segment<T, TSpec> >::Type (
        host(t),
        iterBegin,
        iterEnd);
}

//////////////////////////////////////////////////////////////////////////////

template <typename T, typename TPosBegin, typename TSize>
inline typename Infix<T>::Type
infixWithLength(T && t, TPosBegin pos_begin, TSize length)
{
    return infix(t, pos_begin, pos_begin + length);
}

template <typename T, typename TPosBegin, typename TSize>
inline typename Infix<T *>::Type
infixWithLength(T * t, TPosBegin pos_begin, TSize length)
{
    return infix(*t, pos_begin, pos_begin + length);
}

template <typename T, typename TSpec, typename TPosBegin, typename TSize>
inline typename Infix<Segment<T, TSpec> >::Type
infixWithLength(Segment<T, TSpec> & t, TPosBegin pos_begin, TSize length)
{
    return infix(host(t), beginPosition(t) + pos_begin, beginPosition(t) + pos_begin + length);
}

template <typename T, typename TSpec, typename TPosBegin, typename TSize>
inline typename Infix<Segment<T, TSpec> const>::Type
infixWithLength(Segment<T, TSpec> const & t, TPosBegin pos_begin, TSize length)
{
    return infix(host(t), beginPosition(t) + pos_begin, beginPosition(t) + pos_begin + length);
}

//////////////////////////////////////////////////////////////////////////////
//setBegin


template <typename TIterator>
inline void
setBegin(TIterator new_begin)
{
    setBegin(container(new_begin), hostIterator(new_begin));
}


//////////////////////////////////////////////////////////////////////////////
//setEnd

template <typename TIterator>
inline void
setEnd(TIterator new_end)
{
    setEnd(container(new_end), new_end);
}

//////////////////////////////////////////////////////////////////////////////

} //namespace seqan

#endif //#ifndef SEQAN_HEADER_...
