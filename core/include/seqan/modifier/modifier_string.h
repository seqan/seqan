// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_MODIFIER_MODIFIER_STRING_H_
#define SEQAN_MODIFIER_MODIFIER_STRING_H_

namespace SEQAN_NAMESPACE_MAIN
{

/**
.Class.ModifiedString:
..summary:Allows to modify arbitrary strings by specializing what differs from an origin.
..cat:Modifier
..signature:ModifiedString<THost[, TSpec]>
..param.THost:Original sequence type.
...type:Concept.ContainerConcept
..param.TSpec:The modifier type.
...metafunction:Metafunction.Spec
..implements:Concept.ContainerConcept
..remarks:$THost$ can also be a modified string, so you can create custom strings by combining predefined ones.
..include:seqan/modifier.h
*/

template < typename THost, typename TSpec = void >
class ModifiedString {
public:
    Holder<THost>							data_host;
    typename Cargo<ModifiedString>::Type	data_cargo;

    ModifiedString() {
        SEQAN_CHECKPOINT;
    }

    ModifiedString(ModifiedString &_origin):
            data_host(_origin.data_host),
            data_cargo(_origin.data_cargo) {
        SEQAN_CHECKPOINT;
    }

    ModifiedString(ModifiedString const &_origin):
            data_host(_origin.data_host),
            data_cargo(_origin.data_cargo) {
        SEQAN_CHECKPOINT;
    }
    
    ModifiedString(THost &_origin) {
        SEQAN_CHECKPOINT;
        setHost(*this, _origin);
    }

    template <typename T>
    ModifiedString(T & _origin) {
        SEQAN_CHECKPOINT;
        setValue(*this, _origin);
    }

    template <typename T>
    ModifiedString(T const & _origin) {
        SEQAN_CHECKPOINT;
        setValue(*this, _origin);
    }

    template <typename T>
    inline ModifiedString const &
    operator=(T & _origin) {
        SEQAN_CHECKPOINT;
        assign(*this, _origin);
        return *this;
    }
};

template < typename THost, typename TSpec >
struct Spec< ModifiedString<THost, TSpec> > {
    typedef TSpec Type;
};

template < typename THost, typename TSpec >
struct Spec< ModifiedString<THost, TSpec> const > {
    typedef TSpec Type;
};


// use Value, GetValue, Reference, Size, ... from corresponding iterator
template < typename THost, typename TSpec >
struct Value< ModifiedString<THost, TSpec> >:
    Value< typename Iterator< ModifiedString<THost, TSpec>, Rooted >::Type > {};

template < typename THost, typename TSpec >
struct Value< ModifiedString<THost, TSpec> const >:
    Value< typename Iterator< ModifiedString<THost, TSpec> const, Rooted >::Type > {};


template < typename THost, typename TSpec >
struct GetValue< ModifiedString<THost, TSpec> >:
    GetValue< typename Iterator< ModifiedString<THost, TSpec>, Rooted >::Type > {};

template < typename THost, typename TSpec >
struct GetValue< ModifiedString<THost, TSpec> const >:
    GetValue< typename Iterator< ModifiedString<THost, TSpec> const, Rooted >::Type > {};


template < typename THost, typename TSpec >
struct Reference< ModifiedString<THost, TSpec> >:
    Reference< typename Iterator< ModifiedString<THost, TSpec>, Rooted >::Type > {};

template < typename THost, typename TSpec >
struct Reference< ModifiedString<THost, TSpec> const >:
    Reference< typename Iterator< ModifiedString<THost, TSpec> const, Rooted >::Type > {};

///.Metafunction.Size.param.T.type:Class.ModifiedString
///.Metafunction.Size.class:Class.ModifiedString

template < typename THost, typename TSpec >
struct Size< ModifiedString<THost, TSpec> >:
    Size< typename Iterator< ModifiedString<THost, TSpec>, Rooted >::Type > {};

template < typename THost, typename TSpec >
struct Position< ModifiedString<THost, TSpec> >:
    Position< typename Iterator< ModifiedString<THost, TSpec>, Rooted >::Type > {};

template < typename THost, typename TSpec >
struct Difference< ModifiedString<THost, TSpec> >:
    Difference< typename Iterator< ModifiedString<THost, TSpec>, Rooted >::Type > {};


///.Metafunction.Iterator.param.T.type:Class.ModifiedString
///.Metafunction.Iterator.class:Class.ModifiedString

template <typename THost, typename TSpec>
struct Iterator< ModifiedString<THost, TSpec>, Standard > {
    typedef ModifiedIterator<typename Iterator<THost, Standard>::Type, TSpec> Type;
};

template <typename THost, typename TSpec >
struct Iterator< ModifiedString<THost, TSpec> const, Standard > {
    typedef ModifiedIterator<typename Iterator<THost const, Standard>::Type, TSpec> Type;
};

template <typename THost, typename TSpec>
struct Iterator< ModifiedString<THost, TSpec>, Rooted > {
    typedef ModifiedIterator<typename Iterator<THost, Rooted>::Type, TSpec> Type;
};

template <typename THost, typename TSpec >
struct Iterator< ModifiedString<THost, TSpec> const, Rooted > {
    typedef ModifiedIterator<typename Iterator<THost const, Rooted>::Type, TSpec> Type;
};

///.Metafunction.Host.param.T.type:Class.ModifiedString
///.Metafunction.Host.class:Class.ModifiedString

template < typename THost, typename TSpec >
struct Host< ModifiedString<THost, TSpec> > {
    typedef THost Type;
};

template < typename THost, typename TSpec >
struct Host< ModifiedString<THost, TSpec> const > {
    typedef THost const Type;
};



template < typename THost, typename TSpec >
struct Parameter_< ModifiedString<THost, TSpec> > {
    typedef ModifiedString<THost, TSpec> Type;
};

template < typename THost, typename TSpec >
struct Parameter_< ModifiedString<THost, TSpec> const > {
    typedef ModifiedString<THost, TSpec> const Type;
};


///.Metafunction.IsSequence.param.T.type:Class.ModifiedString

template < typename THost, typename TSpec >
struct IsSequence< ModifiedString<THost, TSpec> > {
    typedef True Type;
    enum { VALUE = true };
};

template < typename THost, typename TSpec >
struct AllowsFastRandomAccess< ModifiedString<THost, TSpec> >:
    AllowsFastRandomAccess<THost> {};

//////////////////////////////////////////////////////////////////////////////
// host interface
//////////////////////////////////////////////////////////////////////////////

template <typename THost, typename TSpec>
inline Holder<THost> &
_dataHost(ModifiedString<THost, TSpec> & me) 
{
    SEQAN_CHECKPOINT;
    return me.data_host;
}

template <typename THost, typename TSpec>
inline Holder<THost> const &
_dataHost(ModifiedString<THost, TSpec> const & me) 
{
    SEQAN_CHECKPOINT;
    return me.data_host;
}

template <typename THost, typename TSpec>
inline typename Reference< typename Cargo<ModifiedString<THost, TSpec> >::Type >::Type
cargo(ModifiedString<THost, TSpec> & me) 
{
    SEQAN_CHECKPOINT;
    return me.data_cargo;
}

template <typename THost, typename TSpec>
inline typename Reference< typename Cargo<ModifiedString<THost, TSpec> const>::Type >::Type
cargo(ModifiedString<THost, TSpec> const & me) 
{
    SEQAN_CHECKPOINT;
    return me.data_cargo;
}

//////////////////////////////////////////////////////////////////////////////
// assign
//////////////////////////////////////////////////////////////////////////////

template <typename TDest, typename TSource>
inline void _copyCargoImpl(TDest &, TSource &, False const) {}

template <typename TDest, typename TSource>
inline void _copyCargoImpl(TDest & me, TSource & _origin, True const) {
    SEQAN_CHECKPOINT;
    cargo(me) = cargo(_origin);
}

template <typename TDest, typename TSource>
inline void _copyCargo(TDest & me, TSource & _origin) {
    SEQAN_CHECKPOINT;
    _copyCargoImpl(me, _origin, typename IsSameType<
                   typename RemoveConst_<typename Cargo<TDest>::Type >::Type, 
                   typename RemoveConst_<typename Cargo<TSource>::Type>::Type >::Type());
}


template <typename THost, typename TSpec, typename THost2>
inline ModifiedString<THost, TSpec> const &
assign(ModifiedString<THost, TSpec> & me, ModifiedString<THost2, TSpec> & _origin) {
    SEQAN_CHECKPOINT;
   host(me) = host(_origin);
    _copyCargo(me, _origin);
    return me;
}

template <typename THost, typename TSpec, typename THost2>
inline ModifiedString<THost, TSpec> const &
assign(ModifiedString<THost, TSpec> & me, ModifiedString<THost2, TSpec> const & _origin) {
    SEQAN_CHECKPOINT;
    host(me) = host(_origin);
    _copyCargo(me, _origin);
    return me;
}

template <typename THost, typename TSpec, typename T>
inline ModifiedString<THost, TSpec> const &
assign(ModifiedString<THost, TSpec> & me, T & _origin) {
    SEQAN_CHECKPOINT;
    host(me) = _origin;
    return me;
}

template <typename THost, typename TSpec, typename T>
inline ModifiedString<THost, TSpec> const &
assign(ModifiedString<THost, TSpec> & me, T const & _origin) {
    SEQAN_CHECKPOINT;
    host(me) = _origin;
    return me;
}

//////////////////////////////////////////////////////////////////////////////
// value
//////////////////////////////////////////////////////////////////////////////

template <typename THost, typename TSpec, typename TPos>
inline typename Reference<ModifiedString<THost, TSpec> >::Type 
value(ModifiedString<THost, TSpec> & me, TPos pos)
{
    SEQAN_CHECKPOINT;
    return value(begin(me, Standard()) + pos);
}

template <typename THost, typename TSpec, typename TPos>
inline typename Reference<ModifiedString<THost, TSpec> const >::Type 
value(ModifiedString<THost, TSpec> const & me, TPos pos)
{
    SEQAN_CHECKPOINT;
    return value(begin(me, Standard()) + pos);
}

//////////////////////////////////////////////////////////////////////////////
// setValue
//////////////////////////////////////////////////////////////////////////////

template <typename THost, typename TSpec>
inline ModifiedString<THost, TSpec> const &
setValue(ModifiedString<THost, TSpec> & me, ModifiedString<THost, TSpec> const & _origin) {
    SEQAN_CHECKPOINT;
    setHost(me, host(_origin));
    _copyCargo(me, _origin);
    return me;
}

template <typename THost, typename TSpec>
inline ModifiedString<THost, TSpec> const &
setValue(ModifiedString<THost, TSpec> & me, ModifiedString<THost, TSpec> & _origin) {
    SEQAN_CHECKPOINT;
    setHost(me, host(_origin));
    _copyCargo(me, _origin);
    return me;
}

// pass _origin to parent modifier
template <typename THost, typename THostSpec, typename TSpec, typename THost2>
inline ModifiedString< ModifiedString<THost, THostSpec>, TSpec> const &
setValue(
    ModifiedString< ModifiedString<THost, THostSpec>, TSpec> & me, 
    THost2 const & _origin) 
{
    SEQAN_CHECKPOINT;
    setValue(host(me), _origin);
    return me;
}

template <typename THost, typename THostSpec, typename TSpec, typename THost2>
inline ModifiedString< ModifiedString<THost, THostSpec>, TSpec> const &
setValue(
    ModifiedString< ModifiedString<THost, THostSpec>, TSpec> & me, 
    THost2 & _origin) 
{
    SEQAN_CHECKPOINT;
    setValue(host(me), _origin);
    return me;
}

// set _origin at the innermost modifier
template <typename THost, typename TSpec>
inline ModifiedString<THost, TSpec> const &
setValue(ModifiedString<THost, TSpec> & me, THost const & _origin) {
    SEQAN_CHECKPOINT;
    assignHost(me, _origin);
    return me;
}

template <typename THost, typename TSpec>
inline ModifiedString<THost, TSpec> const &
setValue(ModifiedString<THost, TSpec> & me, THost & _origin) {
    SEQAN_CHECKPOINT;
    setHost(me, _origin);
    return me;
}

// allow _origin conversion at the innermost modifier
template <typename THost, typename TSpec, typename THost2>
inline ModifiedString<THost, TSpec> const &
setValue(ModifiedString<THost, TSpec> & me, THost2 & _origin) {
    SEQAN_CHECKPOINT;
    assignHost(me, _origin);
    return me;
}

template <typename THost, typename TSpec, typename THost2>
inline ModifiedString<THost, TSpec> const &
setValue(ModifiedString<THost, TSpec> & me, THost2 const & _origin) {
    SEQAN_CHECKPOINT;
    assignHost(me, _origin);
    return me;
}

//////////////////////////////////////////////////////////////////////////////
// length
//////////////////////////////////////////////////////////////////////////////

template < typename THost, typename TSpec >
inline typename Size< ModifiedString<THost, TSpec> >::Type 
length(ModifiedString<THost, TSpec> const & me) {
    SEQAN_CHECKPOINT;
    return length(host(me));
}

//////////////////////////////////////////////////////////////////////////////
// begin
//////////////////////////////////////////////////////////////////////////////

template < typename THost, typename TSpec >
inline typename Iterator< ModifiedString<THost, TSpec> const >::Type 
begin(ModifiedString<THost, TSpec> const & me) {
    SEQAN_CHECKPOINT;
    typename Iterator< ModifiedString<THost, TSpec> const >::Type temp_(begin(host(me)));
    _copyCargo(temp_, me);
    return temp_;
}

template < typename THost, typename TSpec >
inline typename Iterator< ModifiedString<THost, TSpec> >::Type 
begin(ModifiedString<THost, TSpec> & me) {
    SEQAN_CHECKPOINT;
    typename Iterator< ModifiedString<THost, TSpec> >::Type temp_(begin(host(me)));
    _copyCargo(temp_, me);
    return temp_;
}

template < typename THost, typename TSpec, typename TTagSpec >
inline typename Iterator< ModifiedString<THost, TSpec> const, Tag<TTagSpec> const >::Type 
begin(ModifiedString<THost, TSpec> const & me, Tag<TTagSpec> const tag_) {
    SEQAN_CHECKPOINT;
    typename Iterator< 
        ModifiedString<THost, TSpec> const, 
        Tag<TTagSpec> const 
    >::Type temp_(begin(host(me), tag_));
    _copyCargo(temp_, me);
    return temp_;
}

template < typename THost, typename TSpec, typename TTagSpec >
inline typename Iterator< ModifiedString<THost, TSpec>, Tag<TTagSpec> const >::Type 
begin(ModifiedString<THost, TSpec> & me, Tag<TTagSpec> const tag_) {
    SEQAN_CHECKPOINT;
    typename Iterator< 
        ModifiedString<THost, TSpec>, 
        Tag<TTagSpec> const 
    >::Type temp_(begin(host(me), tag_));
    _copyCargo(temp_, me);
    return temp_;
}

//////////////////////////////////////////////////////////////////////////////
// end
//////////////////////////////////////////////////////////////////////////////

template < typename THost, typename TSpec >
inline typename Iterator< ModifiedString<THost, TSpec> const >::Type 
end(ModifiedString<THost, TSpec> const & me) {
    SEQAN_CHECKPOINT;
    typename Iterator< ModifiedString<THost, TSpec> const >::Type temp_(end(host(me)));
    _copyCargo(temp_, me);
    return temp_;
}

template < typename THost, typename TSpec >
inline typename Iterator< ModifiedString<THost, TSpec> >::Type 
end(ModifiedString<THost, TSpec> & me) {
    SEQAN_CHECKPOINT;
    typename Iterator< ModifiedString<THost, TSpec> >::Type temp_(end(host(me)));
    _copyCargo(temp_, me);
    return temp_;
}

template < typename THost, typename TSpec, typename TTagSpec >
inline typename Iterator< ModifiedString<THost, TSpec> const, Tag<TTagSpec> const >::Type 
end(ModifiedString<THost, TSpec> const & me, Tag<TTagSpec> const tag_) {
    SEQAN_CHECKPOINT;
    typename Iterator< ModifiedString<THost, TSpec> const, Tag<TTagSpec> const >::Type temp_(end(host(me), tag_));
    _copyCargo(temp_, me);
    return temp_;
}

template < typename THost, typename TSpec, typename TTagSpec >
inline typename Iterator< ModifiedString<THost, TSpec>, Tag<TTagSpec> const >::Type 
end(ModifiedString<THost, TSpec> & me, Tag<TTagSpec> const tag_) {
    SEQAN_CHECKPOINT;
    typename Iterator< ModifiedString<THost, TSpec>, Tag<TTagSpec> const >::Type temp_(end(host(me), tag_));
    _copyCargo(temp_, me);
    return temp_;
}

//////////////////////////////////////////////////////////////////////////////

template <typename THost, typename TSpec, typename TRight >
inline bool
operator == (ModifiedString<THost, TSpec> const & left, 
			TRight const & right)
{
SEQAN_CHECKPOINT
	typename Comparator<ModifiedString<THost, TSpec> >::Type _lex(left, right);
    return isEqual(_lex);
}
template <typename TLeftValue, typename THost, typename TSpec >
inline bool
operator == (TLeftValue * left,
			 ModifiedString<THost, TSpec> const & right)
{
SEQAN_CHECKPOINT
	typename Comparator<ModifiedString<THost, TSpec> >::Type _lex(left, right);
    return isEqual(_lex);
}

//////////////////////////////////////////////////////////////////////////////

template <typename THost, typename TSpec, typename TRight >
inline bool
operator !=(ModifiedString<THost, TSpec> const & left, 
			TRight const & right)
{
SEQAN_CHECKPOINT
	typename Comparator<ModifiedString<THost, TSpec> >::Type _lex(left, right);
    return isNotEqual(_lex);
}
template <typename TLeftValue, typename THost, typename TSpec >
inline bool
operator != (TLeftValue * left,
			 ModifiedString<THost, TSpec> const & right)
{
SEQAN_CHECKPOINT
	typename Comparator<ModifiedString<THost, TSpec> >::Type _lex(left, right);
    return isNotEqual(_lex);
}

//////////////////////////////////////////////////////////////////////////////

template <typename THost, typename TSpec, typename TRight>
inline bool
operator < (ModifiedString<THost, TSpec> const & left, 
			TRight const & right)
{
SEQAN_CHECKPOINT
	return isLess(left, right, typename DefaultPrefixOrder<ModifiedString<THost, TSpec> >::Type());
}
template <typename TLeftValue, typename THost, typename TSpec >
inline bool
operator < (TLeftValue * left,
			 ModifiedString<THost, TSpec> const & right)
{
SEQAN_CHECKPOINT
	return isLess(left, right, typename DefaultPrefixOrder<TLeftValue *>::Type());
}

//////////////////////////////////////////////////////////////////////////////

template <typename THost, typename TSpec, typename TRight>
inline bool
operator <= (ModifiedString<THost, TSpec> const & left, 
			 TRight const & right)
{
SEQAN_CHECKPOINT
	return isLessOrEqual(left, right, typename DefaultPrefixOrder<ModifiedString<THost, TSpec> >::Type());
}
template <typename TLeftValue, typename THost, typename TSpec >
inline bool
operator <= (TLeftValue * left,
			 ModifiedString<THost, TSpec> const & right)
{
SEQAN_CHECKPOINT
	return isLessOrEqual(left, right, typename DefaultPrefixOrder<TLeftValue *>::Type());
}

//////////////////////////////////////////////////////////////////////////////

template <typename THost, typename TSpec, typename TRight>
inline bool
operator > (ModifiedString<THost, TSpec> const & left, 
		TRight const & right)
{
SEQAN_CHECKPOINT
	return isGreater(left, right, typename DefaultPrefixOrder<ModifiedString<THost, TSpec> >::Type());
}
template <typename TLeftValue, typename THost, typename TSpec >
inline bool
operator > (TLeftValue * left,
			 ModifiedString<THost, TSpec> const & right)
{
SEQAN_CHECKPOINT
	return isGreater(left, right, typename DefaultPrefixOrder<TLeftValue *>::Type());
}

//////////////////////////////////////////////////////////////////////////////

template <typename THost, typename TSpec, typename TRight>
inline bool
operator >= (ModifiedString<THost, TSpec> const & left, 
		TRight const & right)
{
SEQAN_CHECKPOINT
	return isGreaterOrEqual(left, right, typename DefaultPrefixOrder<ModifiedString<THost, TSpec> >::Type());
}
template <typename TLeftValue, typename THost, typename TSpec >
inline bool
operator >= (TLeftValue * left,
			 ModifiedString<THost, TSpec> const & right)
{
SEQAN_CHECKPOINT
	return isGreaterOrEqual(left, right, typename DefaultPrefixOrder<TLeftValue *>::Type());
}

//////////////////////////////////////////////////////////////////////////////
// stream operators
//////////////////////////////////////////////////////////////////////////////

template < typename TStream, typename THost, typename TSpec >
inline TStream &
operator << (TStream & target, ModifiedString<THost, TSpec> const & source)
{
    SEQAN_CHECKPOINT;
    write(target, source);
    return target;
}

//////////////////////////////////////////////////////////////////////////////

template < typename TStream, typename THost, typename TSpec >
inline TStream &
operator >> (TStream & source, ModifiedString<THost, TSpec> & target)
{
    SEQAN_CHECKPOINT;
    read(source, target);
    return source;
}

//////////////////////////////////////////////////////////////////////////////

template <typename THost, typename TSpec>
inline void const *
getObjectId(ModifiedString<THost, TSpec> & me) 
{
    SEQAN_CHECKPOINT;
    return getObjectId(host(me));
}
template <typename THost, typename TSpec>
inline void const *
getObjectId(ModifiedString<THost, TSpec> const & me) 
{
    SEQAN_CHECKPOINT;
    return getObjectId(host(me));
}

}  // namespace seqan

#endif  // SEQAN_MODIFIER_MODIFIER_STRING_H_
