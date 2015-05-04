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
    typedef TPositions    Type;
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

//template <typename THost, typename TPositions>
//struct Prefix<ModifiedString<THost, ModPos<TPositions> > >
//{
//    typedef ModifiedString<THost, ModPos<TPositions> >    Type;
//};
//
//template <typename THost, typename TPositions>
//struct Prefix<ModifiedString<THost, ModPos<TPositions> > const> :
//    Prefix<ModifiedString<THost, ModPos<TPositions> > > {};

// ----------------------------------------------------------------------------
// Metafunction Suffix
// ----------------------------------------------------------------------------

//template <typename THost, typename TPositions>
//struct Suffix<ModifiedString<THost, ModPos<TPositions> > >
//{
//    typedef ModifiedString<THost, ModPos<TPositions> >    Type;
//};
//
//template <typename THost, typename TPositions>
//struct Suffix<ModifiedString<THost, ModPos<TPositions> > const> :
//    Suffix<ModifiedString<THost, ModPos<TPositions> > > {};

// ----------------------------------------------------------------------------
// Metafunction Infix
// ----------------------------------------------------------------------------

//template <typename THost, typename TPositions>
//struct Infix<ModifiedString<THost, ModPos<TPositions> > >
//{
//    typedef ModifiedString<THost, ModPos<TPositions> >    Type;
//};
//
//template <typename THost, typename TPositions>
//struct Infix<ModifiedString<THost, ModPos<TPositions> > const> :
//    Infix<ModifiedString<THost, ModPos<TPositions> > > {};

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
// Function resize()
// ----------------------------------------------------------------------------
// this function doesn't do anything as we are not allowed to change the host (only its elements)
// it is, however, implemented for algorithms that get a sequence to work on
// and need to make sure that it has a certain length

//template <typename THost, typename TPositions, typename TSize, typename TValue, typename TExpand>
//inline typename Size< ModifiedString<THost, ModPos<TPositions> > >::Type
//resize(ModifiedString<THost, ModPos<TPositions> > & me, TSize new_length, TValue /* val */, Tag<TExpand>)
//{
//    ignoreUnusedVariableWarning(new_length);
//
//    SEQAN_ASSERT_EQ(new_length, (TSize)length(me));
//    return length(me);
//}
//
//template <typename THost, typename TPositions, typename TSize, typename TExpand>
//inline typename Size< ModifiedString<THost, ModPos<TPositions> > >::Type
//resize(ModifiedString<THost, ModPos<TPositions> > & me, TSize new_length, Tag<TExpand> tag)
//{
//    return resize(me, new_length, Nothing(), tag);
//}

// ----------------------------------------------------------------------------
// Function prefix()
// ----------------------------------------------------------------------------

//template <typename THost, typename TPositions, typename TPosEnd>
//inline typename Prefix<ModifiedString<THost, ModPos<TPositions> > const>::Type
//prefix(ModifiedString<THost, ModPos<TPositions> > const & me, TPosEnd pos_end)
//{
//    return typename Prefix<ModifiedString<THost, ModPos<TPositions> > const>::Type(me._begin, me._begin + pos_end);
//}
//
//template <typename THost, typename TPositions, typename TPosEnd>
//inline typename Prefix<ModifiedString<THost, ModPos<TPositions> > >::Type
//prefix(ModifiedString<THost, ModPos<TPositions> > & me, TPosEnd pos_end)
//{
//    return prefix(reinterpret_cast<ModifiedString<THost, ModPos<TPositions> > const &>(me), pos_end);
//}

// ----------------------------------------------------------------------------
// Function suffix()
// ----------------------------------------------------------------------------

//template <typename THost, typename TPositions, typename TPosBegin>
//inline typename Suffix<ModifiedString<THost, ModPos<TPositions> > const>::Type
//suffix(ModifiedString<THost, ModPos<TPositions> > const & me, TPosBegin pos_begin)
//{
//    return typename Suffix<ModifiedString<THost, ModPos<TPositions> > const>::Type(me._begin + pos_begin, me._end);
//}
//
//template <typename THost, typename TPositions, typename TPosBegin>
//inline typename Suffix<ModifiedString<THost, ModPos<TPositions> > >::Type
//suffix(ModifiedString<THost, ModPos<TPositions> > & me, TPosBegin pos_begin)
//{
//    return suffix(reinterpret_cast<ModifiedString<THost, ModPos<TPositions> > const &>(me), pos_begin);
//}

// ----------------------------------------------------------------------------
// Function infix()
// ----------------------------------------------------------------------------

//template <typename THost, typename TPositions, typename TPosBegin, typename TPosEnd>
//inline typename Infix<ModifiedString<THost, ModPos<TPositions> > const>::Type
//infix(ModifiedString<THost, ModPos<TPositions> > const & me, TPosBegin pos_begin, TPosEnd pos_end)
//{
//    return typename Infix<ModifiedString<THost, ModPos<TPositions> > >::Type(me._begin + pos_begin, me._begin + pos_end);
//}
//
//template <typename THost, typename TPositions, typename TPosBegin, typename TPosEnd>
//inline typename Infix<ModifiedString<THost, ModPos<TPositions> > >::Type
//infix(ModifiedString<THost, ModPos<TPositions> > & me, TPosBegin pos_begin, TPosEnd pos_end)
//{
//    return infix(reinterpret_cast<ModifiedString<THost, ModPos<TPositions> > const &>(me), pos_begin, pos_end);
//}

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

}  // namespace seqan

#endif  // SEQAN_MODIFIER_MODIFIER_POSITION_H_
