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

#ifndef INCLUDE_SEQAN_JOURNALED_STRING_TREE_DELTA_MAP_EXTENSION_ITERATOR_H_
#define INCLUDE_SEQAN_JOURNALED_STRING_TREE_DELTA_MAP_EXTENSION_ITERATOR_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// Bidirectional Iterator.
template <typename TContainer>
class Iter<TContainer, DeltaMapExtensionIterSpec>
{
public:
    typedef typename Host<TContainer>::Type                         THost;
    typedef typename Member<TContainer, ExtensionMapMember>::Type   TExtensionTable;
    typedef typename Value<TContainer>::Type                        TExtensionEntry;

    typedef typename Iterator<THost, Standard>::Type                THostIter;
    typedef typename Iterator<TExtensionTable, Standard>::Type      TExtensionIter;

    TContainer*     _contPtr;
    THostIter       _hostMapIter;
    TExtensionIter  _extTableIter;

    mutable TExtensionEntry _tmp;

    Iter()
    {}


    Iter(TContainer & cont) :
        _contPtr(&cont),
        _hostMapIter(begin(host(cont), Standard())),
        _extTableIter(begin(cont._extTable, Standard()))
    {}

    template <typename TOtherCont>
    Iter(Iter<TOtherCont, DeltaMapExtensionIterSpec> const & other,
         SEQAN_CTOR_ENABLE_IF(IsConstructible<TContainer, TOtherCont>)) :
         _contPtr(other._contPtr),
         _hostMapIter(other._hostMapIter),
         _extTableIter(other._extTableIter)
    {
        ignoreUnusedVariableWarning(dummy);
    }

    template <typename TOtherCont>
    SEQAN_FUNC_ENABLE_IF(IsConstructible<TContainer, TOtherCont>, Iter&)
    operator=(Iter<TOtherCont, DeltaMapExtensionIterSpec> const & other)
    {
        if (this != &other)
        {
            _contPtr(other._contPtr);
            _hostMapIter = other._hostMapIter;
            _extTableIter = other._extTableIter;
        }
        return *this;
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

template <typename TContainer>
struct Position<Iter<TContainer, DeltaMapExtensionIterSpec> >
{
    typedef typename Iter<TContainer, DeltaMapExtensionIterSpec>::THostIter THostIter_;
    typedef typename Position<THostIter_>::Type                             Type;
};

// ============================================================================
// Functions
// ============================================================================

// Special function to obtain the position.
//template <typename TContainer>
//inline typename Position<Iter<TContainer, DeltaMapExtensionIterSpec> >::Type
//getDeltaPosition(Iter<TContainer, DeltaMapExtensionIterSpec> const & iter)
//{
//    if (getDeltaPosition(*iter._hostMapIter) < iter._extTableIter->i1)
//        return getDeltaPosition(*iter._hostMapIter);
//    return iter._extTableIter->i1;
//}
//
//template <typename TContainer>
//inline bool
//isEndPoint(Iter<TContainer, DeltaMapExtensionIterSpec> const & iter)
//{
//    if (getDeltaPosition(*iter._hostMapIter) < iter._extTableIter->i1)
//        return false;
//    return true;
//}

namespace impl
{
template <typename TContainer>
inline typename HostIter_<typename Container<Iter<TContainer, DeltaMapExtensionIterSpec> >::Type>::Type
hostIter(Iter<TContainer, DeltaMapExtensionIterSpec> & iter)
{
    return (*iter).hostIter;
}

template <typename TContainer>
inline typename HostIter_<typename Container<Iter<TContainer, DeltaMapExtensionIterSpec> const>::Type>::Type
hostIter(Iter<TContainer, DeltaMapExtensionIterSpec> const & iter)
{
    return (*iter).hostIter;
}

template <typename TContainer>
inline bool
isEndPoint(Iter<TContainer, DeltaMapExtensionIterSpec> const & iter)
{
    return (*iter).info == ExtensionInfo::IS_END;
}
}

template <typename TContainer>
inline typename Container<Iter<TContainer, DeltaMapExtensionIterSpec> >::Type
container(Iter<TContainer, DeltaMapExtensionIterSpec> & iter)
{
    return *iter._contPtr;
}

template <typename TContainer>
inline typename Container<Iter<TContainer, DeltaMapExtensionIterSpec> const>::Type
container(Iter<TContainer, DeltaMapExtensionIterSpec> const & iter)
{
    return *iter._contPtr;
}

template <typename TContainer>
inline bool
atEnd(Iter<TContainer, DeltaMapExtensionIterSpec> const & iter)
{
    return iter._hostMapIter == end(host(*iter._contPtr), Standard()) &&
           iter._extTableIter == end(iter._contPtr->_extTable);
}

template <typename TContainer>
inline typename Reference<Iter<TContainer, DeltaMapExtensionIterSpec> >::Type
operator*(Iter<TContainer, DeltaMapExtensionIterSpec> & iter)
{
    typedef typename Value<TContainer>::Type TEntry;

    if (SEQAN_UNLIKELY(iter._extTableIter == end(iter._contPtr->_extTable, Standard())))
    {
        iter._tmp = TEntry(iter._hostMapIter, getDeltaPosition(*iter._hostMapIter), ExtensionInfo::IS_BEGIN);
        return iter._tmp;
    }
    else if (SEQAN_UNLIKELY(iter._hostMapIter == end(host(*iter._contPtr), Standard())))
    {
        return *iter._extTableIter;
    }

    if (getDeltaPosition(*iter._hostMapIter) < iter._extTableIter->deltaPos)
    {
        iter._tmp = TEntry(iter._hostMapIter, getDeltaPosition(*iter._hostMapIter), ExtensionInfo::IS_BEGIN);
        return iter._tmp;
    }
    return *iter._extTableIter;  // But we cannot find the position
}

template <typename TContainer>
inline typename Reference<Iter<TContainer, DeltaMapExtensionIterSpec> const>::Type
operator*(Iter<TContainer, DeltaMapExtensionIterSpec> const & iter)
{
    typedef typename Value<TContainer>::Type TEntry;

    if (getDeltaPosition(*iter._hostMapIter) < iter._extTableIter->deltaPos)
    {
        iter._tmp = TEntry(iter._hostMapIter, getDeltaPosition(*iter._hostMapIter), ExtensionInfo::IS_BEGIN);
        return iter._tmp;
    }
    return *iter._extTableIter;  // But we cannot find the position

}

template <typename TContainerLhs,
          typename TContainerRhs>
inline bool
operator==(Iter<TContainerLhs, DeltaMapExtensionIterSpec> const & iterLhs,
           Iter<TContainerRhs, DeltaMapExtensionIterSpec> const & iterRhs)
{
    return iterLhs._hostMapIter == iterRhs._hostMapIter && iterLhs._extTableIter == iterRhs._extTableIter;
}

template <typename TContainerLhs,
          typename TContainerRhs>
inline bool
operator!=(Iter<TContainerLhs, DeltaMapExtensionIterSpec> const & iterLhs,
           Iter<TContainerRhs, DeltaMapExtensionIterSpec> const & iterRhs)
{
    return !(iterLhs == iterRhs);
}

template <typename TContainer>
inline Iter<TContainer, DeltaMapExtensionIterSpec>&
operator++(Iter<TContainer, DeltaMapExtensionIterSpec> & iter)
{
    if (SEQAN_UNLIKELY(iter._hostMapIter == end(host(*iter._contPtr), Standard())))
    {
        ++iter._extTableIter;
        return iter;
    }
    if (SEQAN_UNLIKELY(iter._extTableIter == end(iter._contPtr->_extTable, Standard())))
    {
        ++iter._hostMapIter;
        return iter;
    }

    if (getDeltaPosition(*iter._hostMapIter) < iter._extTableIter->deltaPos)
        ++iter._hostMapIter;
    else
        ++iter._extTableIter;
    return iter;
}

template <typename TContainer>
inline Iter<TContainer, DeltaMapExtensionIterSpec>
operator++(Iter<TContainer, DeltaMapExtensionIterSpec> & iter, int)
{
    Iter<TContainer, DeltaMapExtensionIterSpec> tmp(iter);
    ++iter;
    return tmp;
}

template <typename TContainer,
          typename TIntegral>
inline Iter<TContainer, DeltaMapExtensionIterSpec> &
operator+=(Iter<TContainer, DeltaMapExtensionIterSpec> & iter, TIntegral steps)
{
    while (steps-- != 0)
        ++iter;
    return iter;
}

template <typename TContainer, typename TIntegral>
inline Iter<TContainer, DeltaMapExtensionIterSpec>
operator+(Iter<TContainer, DeltaMapExtensionIterSpec> iter,
          TIntegral step)
{
    return (iter += step);
}

template <typename TContainer>
inline Iter<TContainer, DeltaMapExtensionIterSpec>&
operator--(Iter<TContainer, DeltaMapExtensionIterSpec> & iter)
{
    if (SEQAN_UNLIKELY(iter._hostMapIter == begin(host(*iter._contPtr), Standard())))
    {
        --iter._extTableIter;
        return iter;
    }
    if (SEQAN_UNLIKELY(iter._extTableIter == begin(iter._contPtr->_extTable, Standard())))
    {
        --iter._hostMapIter;
        return iter;
    }

    if (SEQAN_LIKELY(getDeltaPosition(*(iter._hostMapIter - 1)) < (iter._extTableIter - 1)->deltaPos))
        --iter._extTableIter;
    else
        --iter._hostMapIter;
    return iter;
}

template <typename TContainer>
inline Iter<TContainer, DeltaMapExtensionIterSpec>
operator--(Iter<TContainer, DeltaMapExtensionIterSpec> & iter, int)
{
    Iter<TContainer, DeltaMapExtensionIterSpec> tmp(iter);
    --iter;
    return tmp;
}

template <typename TContainer,
          typename TIntegral>
inline Iter<TContainer, DeltaMapExtensionIterSpec> &
operator-=(Iter<TContainer, DeltaMapExtensionIterSpec> & iter, TIntegral steps)
{
    while (steps-- != 0)
        --iter;
    return iter;
}

template <typename TContainer, typename TIntegral>
inline Iter<TContainer, DeltaMapExtensionIterSpec>
operator-(Iter<TContainer, DeltaMapExtensionIterSpec> iter,
          TIntegral step)
{
    return (iter -= step);
}

template <typename TContainer>
inline typename Size<Iter<TContainer, DeltaMapExtensionIterSpec> >::Type
operator-(Iter<TContainer, DeltaMapExtensionIterSpec> const & lhs,
          Iter<TContainer, DeltaMapExtensionIterSpec> const & rhs)
{
    return (lhs._hostMapIter - rhs._hostMapIter) + (lhs._extTableIter - rhs._extTableIter);
}
}

#endif  // #ifndef INCLUDE_SEQAN_JOURNALED_STRING_TREE_DELTA_MAP_EXTENSION_ITERATOR_H_
