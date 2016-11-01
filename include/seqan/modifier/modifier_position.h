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
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================
// A Position ModifiedString represents a permutation of the host string.

#ifndef SEQAN_MODIFIER_MODIFIER_POSITION_H_
#define SEQAN_MODIFIER_MODIFIER_POSITION_H_

namespace seqan {

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Tag ModPos
// ----------------------------------------------------------------------------

template <typename TPositions>
struct ModPos {};

// ============================================================================
// Metafunctions
// ============================================================================

// --------------------------------------------------------------------------
// Metafunction Cargo
// --------------------------------------------------------------------------

template <typename THost, typename TPositions>
struct Cargo<ModifiedString<THost, ModPos<TPositions> > >
{
    typedef typename Pointer_<TPositions>::Type  Type;
};

template <typename THost, typename TPositions>
struct Cargo<ModifiedString<THost, ModPos<TPositions> > const>
{
    typedef typename Pointer_<TPositions>::Type const Type;
};

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

template <typename THost, typename TPositions>
struct Value<ModifiedString<THost, ModPos<TPositions> > > : Value<THost> {};

template <typename THost, typename TPositions>
struct Value<ModifiedString<THost, ModPos<TPositions> > const> : Value<THost> {};

// ----------------------------------------------------------------------------
// Metafunction GetValue
// ----------------------------------------------------------------------------

template <typename THost, typename TPositions>
struct GetValue<ModifiedString<THost, ModPos<TPositions> > > : GetValue<THost> {};

template <typename THost, typename TPositions>
struct GetValue<ModifiedString<THost, ModPos<TPositions> > const> : GetValue<THost> {};

// ----------------------------------------------------------------------------
// Metafunction Reference
// ----------------------------------------------------------------------------

template <typename THost, typename TPositions>
struct Reference<ModifiedString<THost, ModPos<TPositions> > > : Reference<THost> {};

template <typename THost, typename TPositions>
struct Reference<ModifiedString<THost, ModPos<TPositions> > const> : Reference<THost> {};

// ----------------------------------------------------------------------------
// Metafunction Difference
// ----------------------------------------------------------------------------

template <typename THost, typename TPositions>
struct Difference<ModifiedString<THost, ModPos<TPositions> > > : Difference<TPositions> {};

// ----------------------------------------------------------------------------
// Metafunction Size
// ----------------------------------------------------------------------------

template <typename THost, typename TPositions>
struct Size<ModifiedString<THost, ModPos<TPositions> > > : Size<TPositions> {};

// ----------------------------------------------------------------------------
// Metafunction Position
// ----------------------------------------------------------------------------

template <typename THost, typename TPositions>
struct Position<ModifiedString<THost, ModPos<TPositions> > > : Position<TPositions> {};

// ----------------------------------------------------------------------------
// Metafunction Iterator
// ----------------------------------------------------------------------------

template <typename THost, typename TPositions>
struct Iterator<ModifiedString<THost, ModPos<TPositions> >, Standard>
{
    typedef Iter<ModifiedString<THost, ModPos<TPositions> >, PositionIterator>   Type;
};

template <typename THost, typename TPositions>
struct Iterator<ModifiedString<THost, ModPos<TPositions> > const, Standard>
{
    typedef Iter<ModifiedString<THost, ModPos<TPositions> > const, PositionIterator>   Type;
};

template <typename THost, typename TPositions>
struct Iterator<ModifiedString<THost, ModPos<TPositions> >, Rooted>
{
    typedef Iter<ModifiedString<THost, ModPos<TPositions> >, PositionIterator>   Type;
};

template <typename THost, typename TPositions>
struct Iterator<ModifiedString<THost, ModPos<TPositions> > const, Rooted>
{
    typedef Iter<ModifiedString<THost, ModPos<TPositions> > const, PositionIterator>   Type;
};

// ----------------------------------------------------------------------------
// Metafunction Prefix
// ----------------------------------------------------------------------------
// NOTE(esiragusa): not implemented.

// ----------------------------------------------------------------------------
// Metafunction Suffix
// ----------------------------------------------------------------------------
// NOTE(esiragusa): not implemented.

// ----------------------------------------------------------------------------
// Metafunction Infix
// ----------------------------------------------------------------------------

template <typename THost, typename TPositions>
struct Infix<ModifiedString<THost, ModPos<TPositions> > >
{
    typedef ModifiedString<THost, ModPos<typename Infix<TPositions>::Type> >    Type;
};

template <typename THost, typename TPositions>
struct Infix<ModifiedString<THost, ModPos<TPositions> > const>
{
    typedef ModifiedString<THost, ModPos<typename Infix<TPositions>::Type> > const   Type;
};

// --------------------------------------------------------------------------
// Metafunction AllowsFastRandomAccess
// --------------------------------------------------------------------------

template <typename THost, typename TPositions>
struct AllowsFastRandomAccess<ModifiedString<THost, ModPos<TPositions> > > :
    And<AllowsFastRandomAccess<THost>, AllowsFastRandomAccess<TPositions> >
{};

// ============================================================================
// Functions
// ============================================================================

// --------------------------------------------------------------------------
// Function cargo()
// --------------------------------------------------------------------------

template <typename THost, typename TPositions>
inline typename Parameter_<TPositions>::Type
cargo(ModifiedString<THost, ModPos<TPositions> > & me)
{
    return _toParameter<TPositions>(me._cargo);
}

template <typename THost, typename TPositions>
inline typename Parameter_<TPositions>::Type
cargo(ModifiedString<THost, ModPos<TPositions> > const & me)
{
    return _toParameter<TPositions>(me._cargo);
}

// --------------------------------------------------------------------------
// Function setCargo()
// --------------------------------------------------------------------------

template <typename THost, typename TPositions>
inline void
setCargo(ModifiedString<THost, ModPos<TPositions> > & me, typename Parameter_<TPositions>::Type _cargo)
{
    me._cargo = _toPointer(_cargo);
}

// ----------------------------------------------------------------------------
// Function begin()
// ----------------------------------------------------------------------------

template <typename THost, typename TPositions, typename TTagSpec>
inline typename Iterator<ModifiedString<THost, ModPos<TPositions> >, Tag<TTagSpec> const>::Type
begin(ModifiedString<THost, ModPos<TPositions> > & me, Tag<TTagSpec> const /* tag */)
{
    return typename Iterator<ModifiedString<THost, ModPos<TPositions> >, Tag<TTagSpec> const>::Type(me);
}

template <typename THost, typename TPositions, typename TTagSpec>
inline typename Iterator<ModifiedString<THost, ModPos<TPositions> > const, Tag<TTagSpec> const>::Type
begin(ModifiedString<THost, ModPos<TPositions> > const & me, Tag<TTagSpec> const /* tag */)
{
    return typename Iterator<ModifiedString<THost, ModPos<TPositions> > const, Tag<TTagSpec> const>::Type(me);
}

// ----------------------------------------------------------------------------
// Function end()
// ----------------------------------------------------------------------------

template <typename THost, typename TPositions, typename TTagSpec>
inline typename Iterator<ModifiedString<THost, ModPos<TPositions> >, Tag<TTagSpec> const>::Type
end(ModifiedString<THost, ModPos<TPositions> > & me, Tag<TTagSpec> const /* tag */)
{
    return typename Iterator<ModifiedString<THost, ModPos<TPositions> >, Tag<TTagSpec> const>::Type(me, length(me));
}

template <typename THost, typename TPositions, typename TTagSpec>
inline typename Iterator<ModifiedString<THost, ModPos<TPositions> > const, Tag<TTagSpec> const>::Type
end(ModifiedString<THost, ModPos<TPositions> > const & me, Tag<TTagSpec> const /* tag */)
{
    return typename Iterator<ModifiedString<THost, ModPos<TPositions> > const, Tag<TTagSpec> const>::Type(me, length(me));
}

// ----------------------------------------------------------------------------
// Function value()
// ----------------------------------------------------------------------------

template <typename THost, typename TPositions, typename TPos>
inline typename Reference<ModifiedString<THost, ModPos<TPositions> > >::Type
value(ModifiedString<THost, ModPos<TPositions> > & me, TPos pos)
{
    return value(host(me), getValue(cargo(me), pos));
}

template <typename THost, typename TPositions, typename TPos>
inline typename Reference<ModifiedString<THost, ModPos<TPositions> > const>::Type
value(ModifiedString<THost, ModPos<TPositions> > const & me, TPos pos)
{
    return value(host(me), getValue(cargo(me), pos));
}

// ----------------------------------------------------------------------------
// Function getValue()
// ----------------------------------------------------------------------------

template <typename THost, typename TPositions, typename TPos>
inline typename GetValue<ModifiedString<THost, ModPos<TPositions> > >::Type
getValue(ModifiedString<THost, ModPos<TPositions> > & me, TPos pos)
{
    return getValue(host(me), getValue(cargo(me), pos));
}

template <typename THost, typename TPositions, typename TPos>
inline typename GetValue<ModifiedString<THost, ModPos<TPositions> > const>::Type
getValue(ModifiedString<THost, ModPos<TPositions> > const & me, TPos pos)
{
    return getValue(host(me), getValue(cargo(me), pos));
}

// ----------------------------------------------------------------------------
// Function length()
// ----------------------------------------------------------------------------

template <typename THost, typename TPositions>
inline typename Size<ModifiedString<THost, ModPos<TPositions> > >::Type
length(ModifiedString<THost, ModPos<TPositions> > const & me)
{
    return length(cargo(me));
}

// ----------------------------------------------------------------------------
// Function empty()
// ----------------------------------------------------------------------------

template <typename THost, typename TPositions>
inline bool empty(ModifiedString<THost, ModPos<TPositions> > const & me)
{
    return empty(cargo(me));
}

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

template <typename THost, typename TPositions>
inline void clear(ModifiedString<THost, ModPos<TPositions> > & me)
{
    clear(cargo(me));
}

// ----------------------------------------------------------------------------
// Function resize()
// ----------------------------------------------------------------------------
// NOTE(esiragusa): a dummy implementation like ContainerView could be necessary.

// ----------------------------------------------------------------------------
// Function prefix()
// ----------------------------------------------------------------------------
// NOTE(esiragusa): not implemented.

// ----------------------------------------------------------------------------
// Function suffix()
// ----------------------------------------------------------------------------
// NOTE(esiragusa): not implemented.

// ----------------------------------------------------------------------------
// Function infix()
// ----------------------------------------------------------------------------

template <typename THost, typename TPositions, typename TPosBegin, typename TPosEnd>
inline typename Infix<ModifiedString<THost, ModPos<TPositions> > const>::Type
infix(ModifiedString<THost, ModPos<TPositions> > const & me, TPosBegin pos_begin, TPosEnd pos_end)
{
    typename Infix<ModifiedString<THost, ModPos<TPositions> > const>::Type other;
    setHost(other, host(me));
    setCargo(other, infix(cargo(me), pos_begin, pos_end));
    return other;
}

template <typename THost, typename TPositions, typename TPosBegin, typename TPosEnd>
inline typename Infix<ModifiedString<THost, ModPos<TPositions> > >::Type
infix(ModifiedString<THost, ModPos<TPositions> > & me, TPosBegin pos_begin, TPosEnd pos_end)
{
    typename Infix<ModifiedString<THost, ModPos<TPositions> > >::Type other;
    setHost(other, host(me));
    setCargo(other, infix(cargo(me), pos_begin, pos_end));
    return other;
}

// ----------------------------------------------------------------------------
// Function position()
// ----------------------------------------------------------------------------

template <typename THost, typename TPositions, typename TPos>
inline typename Position<ModifiedString<THost, ModPos<TPositions> > const>::Type
position(ModifiedString<THost, ModPos<TPositions> > const & me, TPos i)
{
    return getValue(cargo(me), i);
}

// ----------------------------------------------------------------------------
// Function setPosition()
// ----------------------------------------------------------------------------

template <typename THost, typename TPositions, typename TPos>
inline void setPosition(ModifiedString<THost, ModPos<TPositions> > & me, TPos i, TPos j)
{
    assignValue(cargo(me), i, j);
}

template <typename THost, typename TPositions, typename TPos>
inline void setPosition(ModifiedString<THost, ModPos<TPositions> > const & me, TPos i, TPos j)
{
    assignValue(cargo(me), i, j);
}

// ----------------------------------------------------------------------------
// Functor PosLess_
// ----------------------------------------------------------------------------

template <typename THost, typename TPos = typename Position<THost>::Type, typename TPredicate = std::less<TPos> >
struct PosLess_ : public std::binary_function<TPos, TPos, bool>
{
    THost const & _host;
    TPredicate pred;

    PosLess_(THost const & _host) :
        _host(_host)
    {}

    PosLess_(THost const & _host, TPredicate const & pred) :
        _host(_host),
        pred(pred)
    {}

    bool operator()(TPos a, TPos b)
    {   
        return pred(getValue(_host, a), getValue(_host, b));
    }   
};

// ----------------------------------------------------------------------------
// Function sort()
// ----------------------------------------------------------------------------

template <typename THost, typename TPositions, typename TBinaryPredicate, typename TParallelTag>
inline void sort(ModifiedString<THost, ModPos<TPositions> > & me, TBinaryPredicate p, Tag<TParallelTag> const & tag)
{
    typedef typename Position<ModifiedString<THost, ModPos<TPositions> > >::Type TPos;

    sort(cargo(me), PosLess_<THost, TPos, TBinaryPredicate>(host(me), p), tag);
}

template <typename THost, typename TPositions, typename TParallelTag>
inline void sort(ModifiedString<THost, ModPos<TPositions> > & me, Tag<TParallelTag> const & tag)
{
    typedef typename Position<ModifiedString<THost, ModPos<TPositions> > >::Type TPos;

    sort(cargo(me), PosLess_<THost, TPos>(host(me)), tag);
}

template <typename THost, typename TPositions, typename TBinaryPredicate, typename TParallelTag>
inline void sort(ModifiedString<THost, ModPos<TPositions> > const & me, TBinaryPredicate p, Tag<TParallelTag> const & tag)
{
    typedef typename Position<ModifiedString<THost, ModPos<TPositions> > >::Type TPos;

    sort(cargo(me), PosLess_<THost, TPos, TBinaryPredicate>(host(me), p), tag);
}

template <typename THost, typename TPositions, typename TParallelTag>
inline void sort(ModifiedString<THost, ModPos<TPositions> > const & me, Tag<TParallelTag> const & tag)
{
    typedef typename Position<ModifiedString<THost, ModPos<TPositions> > >::Type TPos;

    sort(cargo(me), PosLess_<THost, TPos>(host(me)), tag);
}

}  // namespace seqan

#endif  // SEQAN_MODIFIER_MODIFIER_POSITION_H_
