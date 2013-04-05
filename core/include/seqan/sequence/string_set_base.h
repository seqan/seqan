// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
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
// Author: Andreas Gogol-Doering <andreas.doering@mdc-berlin.de>
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_SEQUENCE_STRING_SET_BASE_H_
#define SEQAN_SEQUENCE_STRING_SET_BASE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TSpec = Default>
struct Owner {};

/**
.Class.StringSet:
..cat:Sequences
..summary:A container class for a set of strings.
..signature:StringSet<TString, TSpec>
..param.TString:The string type.
...type:Class.String
..param.TSpec:The specializing type for the StringSet.
...metafunction:Metafunction.Spec
...default:$Owner<Generous>$.
..include:sequence.h
 */
template <typename TString, typename TSpec = Owner<> >
class StringSet;

// ============================================================================
// Metafunctions
// ============================================================================

// --------------------------------------------------------------------------
// Metafunction Concatenator
// --------------------------------------------------------------------------

/**
.Metafunction.Concatenator:
..class:Class.StringSet
..summary:Returns the type of the concatenation sequence of all sequences in a @Class.StringSet@.
..cat:Sequences
..signature:Concatenator<TStringSet>::Type
..param.TStringSet:The @Class.StringSet@ type.
...type:Class.StringSet
..returns:The type of a container that can be iterated like the concatenation string of all sequences in a @Class.StringSet@.
..include:seqan/sequence.h
*/

// TODO(holtgrew): Why is this specialized for all types?
template <typename TObject>
struct Concatenator
{
    typedef TObject Type;
};

template <typename TObject>
struct Concatenator<TObject const>
{
    typedef typename Concatenator<TObject>::Type const Type;
};

template <typename TString, typename TSpec >
struct Concatenator<StringSet<TString, TSpec> >
{
    typedef ConcatenatorManyToOne<StringSet<TString, TSpec> > Type;
};

// --------------------------------------------------------------------------
// Metafunction StringSetLimits
// --------------------------------------------------------------------------

// TODO(holtgrew): Document these metafunctions.
// TODO(holtgrew): Default specializations necessary?
template <typename TString>
struct StringSetLimits
{
    typedef Nothing Type;
};

template <typename TString>
struct StringSetLimits<TString const>
{
    typedef typename StringSetLimits<TString>::Type const Type;
};

template <typename TString, typename TSpec>
struct StringSetLimits<StringSet<TString, TSpec> >
{
    typedef typename Size<TString>::Type TSize_;
    typedef String<TSize_> Type;
};

// --------------------------------------------------------------------------
// Metafunction StringSetPosition
// --------------------------------------------------------------------------

// TODO(holtgrew): Default specializations necessary?
template <typename TString>
struct StringSetPosition
{
    typedef typename Size<TString>::Type Type;
};

template <typename TString, typename TSpec>
struct StringSetPosition<StringSet<TString, TSpec> >
{
    typedef typename Size<TString>::Type TSize_;
    typedef Pair<TSize_> Type;
};

// --------------------------------------------------------------------------
// Metafunction LengthSum
// --------------------------------------------------------------------------

template <typename TString>
struct LengthSum
{
    typedef typename Size<TString>::Type Type;
};

template <typename TString, typename TSpec>
struct LengthSum<StringSet<TString, TSpec> >
{
    typedef StringSet<TString, TSpec>                   TStringSet;
    typedef typename StringSetLimits<TStringSet>::Type  TLimits;
    typedef typename Value<TLimits>::Type               Type;
};

template <typename T>
struct LengthSum<T const> :
    public LengthSum<T> {};

// --------------------------------------------------------------------------
// Metafunction GetSequenceNo
// --------------------------------------------------------------------------

// TODO(holtgrew): Default specializations necessary?
template <typename TString>
struct GetSequenceByNo
{
    typedef TString & Type;
};

template <typename TString, typename TSpec>
struct GetSequenceByNo<StringSet<TString, TSpec> >
{
    typedef typename Reference< StringSet<TString, TSpec> >::Type Type;
};

template <typename TString, typename TSpec>
struct GetSequenceByNo<StringSet<TString, TSpec> const>
{
    typedef typename Reference< StringSet<TString, TSpec> const>::Type Type;
};

// --------------------------------------------------------------------------
// Metafunction Value
// --------------------------------------------------------------------------

template < typename TString, typename TSpec >
struct Value< StringSet< TString, TSpec > >
{
    typedef TString Type;
};

template < typename TString, typename TSpec >
struct Value< StringSet< TString, TSpec > const>
{
    typedef TString Type;
};

// --------------------------------------------------------------------------
// Metafunction Iterator
// --------------------------------------------------------------------------

template < typename TString, typename TSpec, typename TIteratorSpec>
struct Iterator< StringSet< TString, TSpec >, TIteratorSpec>
{
    typedef Iter< StringSet< TString, TSpec >, PositionIterator> Type;
};

template < typename TString, typename TSpec, typename TIteratorSpec >
struct Iterator< StringSet< TString, TSpec> const, TIteratorSpec>
{
    typedef Iter< StringSet< TString, TSpec > const, PositionIterator> Type;
};

// --------------------------------------------------------------------------
// Metafunction Size
// --------------------------------------------------------------------------

template <typename TString, typename TSpec>
struct Size< StringSet< TString, TSpec > >
    : Size<typename StringSetLimits< StringSet<TString, TSpec> >::Type > {};
// Default Size<T const> redirects to non-const.

// --------------------------------------------------------------------------
// Metafunction Prefix
// --------------------------------------------------------------------------
// TODO(holtgrew): Do Prefix, Suffix, Infix make sense if defined in this way for all StringSet classes?
// TODO(holtgrew): However, if this works nicely then it shows that implementing segments as Strings would not be advantageous since they now work for arbitrary sequential-access containers.

template <typename TString, typename TSpec>
struct Prefix< StringSet< TString, TSpec > >
    : Prefix<TString > {};

template <typename TString, typename TSpec>
struct Prefix<StringSet< TString, TSpec > const>
    : Prefix<TString const > {};

// --------------------------------------------------------------------------
// Metafunction Suffix
// --------------------------------------------------------------------------

template <typename TString, typename TSpec>
struct Suffix<StringSet< TString, TSpec> >
    : Suffix<TString> {};

template <typename TString, typename TSpec>
struct Suffix<StringSet< TString, TSpec> const>
    : Suffix<TString const> {};

// --------------------------------------------------------------------------
// Metafunction Infix
// --------------------------------------------------------------------------

template <typename TString, typename TSpec>
struct Infix<StringSet< TString, TSpec> >
    : Infix<TString> {};

template <typename TString, typename TSpec>
struct Infix<StringSet< TString, TSpec > const>
    : Infix< TString const > {};

// --------------------------------------------------------------------------
// Metafunction AllowsFastRandomAccess
// --------------------------------------------------------------------------

template <typename TString, typename TSpec>
struct AllowsFastRandomAccess<StringSet<TString, TSpec> >
    : AllowsFastRandomAccess<TString> {};
// Default AllowsFastRandomAccess<T const> redirects to non-const.

// --------------------------------------------------------------------------
// Metafunction DefaultOverflowImplicit
// --------------------------------------------------------------------------

template < typename TString, typename TSpec >
struct DefaultOverflowImplicit<StringSet< TString, TSpec> >
{
    typedef Generous Type;
};

template < typename TString, typename TSpec >
struct DefaultOverflowImplicit<StringSet< TString, TSpec> const>
{
    typedef Generous Type;
};

// ============================================================================
// Functions
// ============================================================================


// --------------------------------------------------------------------------
// Function swap()
// --------------------------------------------------------------------------

///.Function.swap.param.left.type:Class.StringSet
///.Function.swap.param.right.type:Class.StringSet
///.Function.swap.class:Class.StringSet

template <typename TString, typename TSpec>
inline void
swap(StringSet<TString, TSpec> & left,
     StringSet<TString, TSpec> & right)
{
    SEQAN_CHECKPOINT;
    typedef StringSet<TString, TSpec> TStringSet;

    TStringSet tmp(left, Move());
    move(left, right);
    move(right, tmp);
}

// --------------------------------------------------------------------------
// Function stringSetLimits()
// --------------------------------------------------------------------------

/**
.Function.stringSetLimits:
..cat:Sequences
..class:Class.String
..class:Class.StringSet
..summary:Retrieves a string of delimiter positions of a @Class.StringSet@ which is needed for local<->global position conversions.
..signature:stringSetLimits(me)
..param.me:A string or string set.
...type:Class.String
...type:Class.StringSet
..returns:A reference to a string.
...remarks:If $me$ is a @Class.StringSet@ then the returned string is of size $length(me)+1$ and contains the ascending (virtual) delimiter positions of the concatenation of all strings in the string set.
...remarks:If $me$ is a @Class.String@, @Tag.Nothing@ is returned.
..include:seqan/sequence.h
*/

// TODO(holtgrew): Default implementation necessary?!
template <typename TStringSet>
inline typename StringSetLimits<TStringSet>::Type
stringSetLimits(TStringSet &)
{
    return typename StringSetLimits<TStringSet>::Type();
}

template <typename TString, typename TSpec>
inline typename StringSetLimits< StringSet<TString, TSpec> >::Type &
stringSetLimits(StringSet<TString, TSpec> & stringSet)
{
    if (!_validStringSetLimits(stringSet))
        _refreshStringSetLimits(stringSet);
    return stringSet.limits;
}

template <typename TString, typename TSpec>
inline typename StringSetLimits< StringSet<TString, TSpec> const>::Type &
stringSetLimits(StringSet<TString, TSpec> const & stringSet)
{
    if (!_validStringSetLimits(stringSet))
        _refreshStringSetLimits(const_cast< StringSet<TString, TSpec>& >(stringSet));
    return stringSet.limits;
}

// --------------------------------------------------------------------------
// Function getSeqNo()
// --------------------------------------------------------------------------

/**
.Function.getSeqNo:
..cat:Sequences
..summary:Returns the sequence number of a position.
..signature:getSeqNo(pos[, limits])
..param.pos:A position.
...type:Class.Pair
..param.limits:The limits string returned by @Function.stringSetLimits@.
..returns:A single integer value that identifies the string within the stringset $pos$ points at.
...remarks:If $limits$ is omitted or @Tag.Nothing@ $getSeqNo$ returns 0.
...remarks:If $pos$ is a local position (of class @Class.Pair@) then $i1$ is returned.
...remarks:If $pos$ is a global position (integer type and $limits$ is a @Class.String@) then $pos$ is converted to a local position and $i1$ is returned.
..include:seqan/sequence.h
*/

// TODO(holtgrew): Auto-sequences should go away!
template <typename TPosition>
inline TPosition
getSeqNo(TPosition const &, Nothing const &)
{
    return 0;
}

// TODO(holtgrew): Auto-sequences should go away!
template <typename TPosition>
inline TPosition
getSeqNo(TPosition const &)
{
    return 0;
}

// n sequences (position type is Pair)
template <typename T1, typename T2, typename TPack, typename TLimitsString>
inline T1 getSeqNo(Pair<T1, T2, TPack> const & pos, TLimitsString const &)
{
    return getValueI1(pos);
}

// n sequences (position type is Pair)
template <typename T1, typename T2, typename TPack>
inline T1 getSeqNo(Pair<T1, T2, TPack> const & pos)
{
    return getValueI1(pos);
}

// n sequences (position type is an integral type)
template <typename TPos, typename TLimitsString>
inline TPos getSeqNo(TPos const & pos, TLimitsString const & limits)
{
    typedef typename Iterator<TLimitsString const, Standard>::Type TIter;
    typedef typename Value<TLimitsString>::Type TSize;
    TIter _begin = begin(limits, Standard());
    TIter _upper = ::std::upper_bound(_begin, end(limits, Standard()), (TSize)pos) - 1;
    return difference(_begin, _upper);
}

// --------------------------------------------------------------------------
// Function getSeqOffset()
// --------------------------------------------------------------------------

/**
.Function.getSeqOffset:
..cat:Sequences
..summary:Returns the local sequence offset of a position.
..signature:getSeqOffset(pos[, limits])
..param.pos:A position.
...type:Class.Pair
..param.limits:The limits string returned by @Function.stringSetLimits@.
..returns:A single integer value that identifies the position within the string $pos$ points at.
...remarks:If $limits$ is omitted or @Tag.Nothing@ $getSeqNo$ returns $pos$.
...remarks:If $pos$ is a local position (of class @Class.Pair@) then $i2$ is returned.
...remarks:If $pos$ is a global position (integer type and $limits$ is a @Class.String@) then $pos$ is converted to a local position and $i2$ is returned.
..include:seqan/sequence.h
*/

// TODO(holtgrew): Auto-sequences should go away!
template <typename TPosition>
inline TPosition
getSeqOffset(TPosition const & pos, Nothing const &)
{
    return pos;
}

// TODO(holtgrew): Auto-sequences should go away!
template <typename TPosition>
inline TPosition
getSeqOffset(TPosition const & pos)
{
    return pos;
}

// n sequences (position type is Pair)
template <typename T1, typename T2, typename TPack, typename TLimitsString>
inline T2 getSeqOffset(Pair<T1, T2, TPack> const & pos, TLimitsString const &) {
    return getValueI2(pos);
}

// n sequences (position type is Pair)
template <typename T1, typename T2, typename TPack>
inline T2 getSeqOffset(Pair<T1, T2, TPack> const & pos) {
    return getValueI2(pos);
}

// n sequences (position type is an integral type)
template <typename TPos, typename TLimitsString>
inline TPos getSeqOffset(TPos const & pos, TLimitsString const & limits) {
    typedef typename Iterator<TLimitsString const, Standard>::Type TIter;
    typedef typename Value<TLimitsString>::Type TSize;
    TIter _begin = begin(limits, Standard());
    TIter _upper = ::std::upper_bound(_begin, end(limits, Standard()), (TSize)pos) - 1;
    return pos - *_upper;
}

// --------------------------------------------------------------------------
// Function setSeqOffset()
// --------------------------------------------------------------------------

// TODO(esiragusa): Implement a spec for global positions.
template <typename TPosition, typename TSeqOffset>
inline void
setSeqOffset(TPosition & pos, TSeqOffset seqOffset)
{
    pos = seqOffset;
}

template <typename T1, typename T2, typename TPack, typename TSeqOffset>
inline void
setSeqOffset(Pair<T1, T2, TPack> & pos, TSeqOffset seqOffset)
{
    setValueI2(pos, seqOffset);
}


// --------------------------------------------------------------------------
// Function posGlobalize()
// --------------------------------------------------------------------------

/**
.Function.posGlobalize:
..cat:Sequences
..summary:Converts a local/global to a global position.
..signature:posGlobalize(pos, limits)
..param.pos:A local or global position (pair or integer value).
...type:Class.Pair
..param.limits:The limits string returned by @Function.stringSetLimits@.
..returns:The corresponding global position of $pos$.
...remarks:If $pos$ is an integral type $pos$ is returned.
...remarks:If not, $limits[getSeqNo(pos, limits)] + getSeqOffset(pos, limits)$ is returned.
..include:seqan/sequence.h
*/

// any_position and no limits_string -> any_position
template <typename TPosition>
inline TPosition posGlobalize(TPosition const & pos, Nothing const &)
{
    return pos;
}

// local_position (0,x) and no limits_string -> global_position x
template <typename T1, typename T2, typename TPack>
inline T2 posGlobalize(Pair<T1, T2, TPack> const & pos, Nothing const &)
{
    return getSeqOffset(pos);
}

// any_position and no limits_string -> any_position
template <typename TLimitsString, typename TPosition>
inline TPosition posGlobalize(TPosition const & pos, TLimitsString const &)
{
    return pos;
}

// local_position and limits_string -> global_position
template <typename TLimitsString, typename T1, typename T2, typename TPack>
inline typename Value<TLimitsString>::Type
posGlobalize(Pair<T1, T2, TPack> const & pos, TLimitsString const & limits)
{
    return limits[getSeqNo(pos, limits)] + getSeqOffset(pos, limits);
}

// --------------------------------------------------------------------------
// Function posLocalToX()
// --------------------------------------------------------------------------

/**
.Function.posLocalToX:
..cat:Sequences
..summary:Converts a local to a local/global position.
..signature:posGlobalize(dst, localPos, limits)
..param.dst:Destination value. A local or global position (pair or integer value).
...type:Class.Pair
..param.localPos:A local position (pair).
...type:Class.Pair
..param.limits:The limits string returned by @Function.stringSetLimits@.
..include:seqan/sequence.h
*/

template <typename TDest, typename TLimitsString, typename T1, typename T2, typename TPack>
inline void
posLocalToX(TDest & dst, Pair<T1, T2, TPack> const & localPos, TLimitsString const & limits)
{
    dst = posGlobalize(localPos, limits);
}

template <typename TD1, typename TD2, typename TDPack, typename TLimitsString, typename T1, typename T2, typename TPack>
inline void
posLocalToX(Pair<TD1, TD2, TDPack> & dst, Pair<T1, T2, TPack> const & localPos, TLimitsString const &)
{
    dst = localPos;
}

// --------------------------------------------------------------------------
// Function posLocalize()
// --------------------------------------------------------------------------

/**
.Function.posLocalize:
..cat:Sequences
..summary:Converts a local/global to a local position.
..signature:posLocalize(result, pos, limits)
..param.pos:A local or global position (pair or integer value).
...type:Class.Pair
..param.limits:The limits string returned by @Function.stringSetLimits@.
..param.result:Reference to the resulting corresponding local position of $pos$.
...remarks:If $pos$ is an integral type and $limits$ is omitted or @Tag.Nothing@, $pos$ is returned.
...remarks:If $pos$ is a local position (of class @Class.Pair@) then $pos$ is returned.
...remarks:If $pos$ is a global position (integer type and $limits$ is a @Class.String@) then $pos$ is converted to a local position.
..include:seqan/sequence.h
*/

// any_position and no limits_string -> any_position
template <typename TResult, typename TPosition>
inline void posLocalize(TResult & result, TPosition const & pos, Nothing const &) {
    result = pos;
}

template <typename T1, typename T2, typename TPack, typename TPosition>
inline void posLocalize(Pair<T1, T2, TPack> & result, TPosition const & pos, Nothing const &) {
    result.i1 = 0;
    result.i2 = pos;
}

// global_position and limits_string -> local_position
template <typename TResult, typename TSize, typename TSpec, typename TPosition>
inline void posLocalize(TResult & result, TPosition const & pos, String<TSize, TSpec> const & limits) {
    typedef typename Iterator<String<TSize, TSpec> const, Standard>::Type TIter;
    TIter _begin = begin(limits, Standard());
    TIter _upper = ::std::upper_bound(_begin, end(limits, Standard()), (TSize)pos) - 1;
    result.i1 = difference(_begin, _upper);
    result.i2 = pos - *_upper;
}

// local_position -> local_position
template <typename TResult, typename TSize, typename TSpec, typename T1, typename T2, typename TPack>
inline void posLocalize(TResult & result, Pair<T1, T2, TPack> const & pos, String<TSize, TSpec> const &/*limits*/) {
    result = pos;
}

// --------------------------------------------------------------------------
// Function prefix()
// --------------------------------------------------------------------------

///.Function.prefix.param.host.type:Class.StringSet
///.Function.prefix.class:Class.StringSet

template < typename TString, typename TSpec, typename TPosition >
inline typename Prefix<TString>::Type
prefix(StringSet< TString, TSpec > & me, TPosition const & pos)
{
    typedef StringSet<TString, TSpec>               TStringSet;
    typedef typename Size<TStringSet>::Type         TSetSize;
    typedef typename Size<TString>::Type            TStringSize;
    typedef Pair<TSetSize, TStringSize, Pack> TPair;

    TPair lPos;
    posLocalize(lPos, pos, stringSetLimits(me));
    return prefix(me[getSeqNo(lPos)], getSeqOffset(lPos));
}

template < typename TString, typename TSpec, typename TPosition >
inline typename Prefix<TString const>::Type
prefix(StringSet< TString, TSpec > const & me, TPosition const & pos)
{
    typedef StringSet<TString, TSpec>               TStringSet;
    typedef typename Size<TStringSet>::Type         TSetSize;
    typedef typename Size<TString>::Type            TStringSize;
    typedef Pair<TSetSize, TStringSize, Pack> TPair;

    TPair lPos;
    posLocalize(lPos, pos, stringSetLimits(me));
    return prefix(me[getSeqNo(lPos)], getSeqOffset(lPos));
}

// --------------------------------------------------------------------------
// Function suffix()
// --------------------------------------------------------------------------

///.Function.suffix.param.host.type:Class.StringSet
///.Function.suffix.class:Class.StringSet

template < typename TString, typename TSpec, typename TPosition >
inline typename Suffix<TString>::Type
suffix(StringSet< TString, TSpec > & me, TPosition const & pos)
{
    typedef StringSet<TString, TSpec>               TStringSet;
    typedef typename Size<TStringSet>::Type         TSetSize;
    typedef typename Size<TString>::Type            TStringSize;
    typedef Pair<TSetSize, TStringSize, Pack> TPair;

    TPair lPos;
    posLocalize(lPos, pos, stringSetLimits(me));
    return suffix(me[getSeqNo(lPos)], getSeqOffset(lPos));
}

template < typename TString, typename TSpec, typename TPosition >
inline typename Suffix<TString const>::Type
suffix(StringSet< TString, TSpec > const & me, TPosition const & pos)
{
    typedef StringSet<TString, TSpec>               TStringSet;
    typedef typename Size<TStringSet>::Type         TSetSize;
    typedef typename Size<TString>::Type            TStringSize;
    typedef Pair<TSetSize, TStringSize, Pack> TPair;

    TPair lPos;
    posLocalize(lPos, pos, stringSetLimits(me));
    return suffix(me[getSeqNo(lPos)], getSeqOffset(lPos));
}

// --------------------------------------------------------------------------
// Function infixWithLength()
// --------------------------------------------------------------------------

///.Function.infixWithLength.param.host.type:Class.StringSet
///.Function.infixWithLength.class:Class.StringSet

template < typename TString, typename TSpec, typename TPosition, typename TSize >
inline typename Infix<TString>::Type
infixWithLength(StringSet< TString, TSpec > & me, TPosition const & pos, TSize length)
{
    typedef StringSet<TString, TSpec>               TStringSet;
    typedef typename Size<TStringSet>::Type         TSetSize;
    typedef typename Size<TString>::Type            TStringSize;
    typedef Pair<TSetSize, TStringSize, Pack> TPair;

    TPair lPos;
    posLocalize(lPos, pos, stringSetLimits(me));
    return infixWithLength(me[getSeqNo(lPos)], getSeqOffset(lPos), length);
}

template < typename TString, typename TSpec, typename TPosition, typename TSize >
inline typename Infix<TString const>::Type
infixWithLength(StringSet< TString, TSpec > const & me, TPosition const & pos, TSize length)
{
    typedef StringSet<TString, TSpec>               TStringSet;
    typedef typename Size<TStringSet>::Type         TSetSize;
    typedef typename Size<TString>::Type            TStringSize;
    typedef Pair<TSetSize, TStringSize, Pack> TPair;

    TPair lPos;
    posLocalize(lPos, pos, stringSetLimits(me));
    return infixWithLength(me[getSeqNo(lPos)], getSeqOffset(lPos), length);
}

// --------------------------------------------------------------------------
// Function infix()
// --------------------------------------------------------------------------

///.Function.infix.param.host.type:Class.StringSet
///.Function.infix.class:Class.StringSet

template < typename TString, typename TSpec, typename TPosBegin, typename TPosEnd >
inline typename Infix<TString>::Type
infix(StringSet< TString, TSpec > & me, TPosBegin const & posBegin, TPosEnd const & posEnd)
{
    typedef StringSet<TString, TSpec>               TStringSet;
    typedef typename Size<TStringSet>::Type         TSetSize;
    typedef typename Size<TString>::Type            TStringSize;
    typedef Pair<TSetSize, TStringSize, Pack> TPair;

    TPair localPosBegin, localPosEnd;
    posLocalize(localPosBegin, posBegin, stringSetLimits(me));
    posLocalize(localPosEnd, posEnd, stringSetLimits(me));
    return infix(me[getSeqNo(localPosBegin)], getSeqOffset(localPosBegin), getSeqOffset(localPosEnd));
}

template < typename TString, typename TSpec, typename TPosBegin, typename TPosEnd >
inline typename Infix<TString const>::Type
infix(StringSet< TString, TSpec > const & me, TPosBegin const & posBegin, TPosEnd const & posEnd)
{
    typedef StringSet<TString, TSpec>               TStringSet;
    typedef typename Size<TStringSet>::Type         TSetSize;
    typedef typename Size<TString>::Type            TStringSize;
    typedef Pair<TSetSize, TStringSize, Pack> TPair;

    TPair localPosBegin, localPosEnd;
    posLocalize(localPosBegin, posBegin, stringSetLimits(me));
    posLocalize(localPosEnd, posEnd, stringSetLimits(me));
    return infix(me[getSeqNo(localPosBegin)], getSeqOffset(localPosBegin), getSeqOffset(localPosEnd));
}

// --------------------------------------------------------------------------
// Function posAtFirstLocal()
// --------------------------------------------------------------------------

template <typename TPos, typename TLimitsString>
inline bool posAtFirstLocal(TPos const & pos, TLimitsString const & limits) {
    return getSeqOffset(pos, limits) == 0;
}
template <typename TPos>
inline bool posAtFirstLocal(TPos const & pos) {
    return getSeqOffset(pos) == 0;
}

// --------------------------------------------------------------------------
// Function posAtEnd()
// --------------------------------------------------------------------------

template <typename T1, typename T2, typename TPack, typename TSequence, typename TSpec>
inline bool posAtEnd(Pair<T1, T2, TPack> const & pos, StringSet<TSequence, TSpec> const & stringSet) {
    return pos.i2 == sequenceLength(pos.i1, stringSet);
}
template <typename TPos, typename TSequence, typename TSpec>
inline bool posAtEnd(TPos pos, StringSet<TSequence, TSpec> const & stringSet) {
    return getSeqOffset(pos, stringSetLimits(stringSet)) == 0;
}
template <typename TPos, typename TSequence>
inline bool posAtEnd(TPos pos, TSequence const & seq) {
    return pos == length(seq);
}

// --------------------------------------------------------------------------
// Function posPrev()
// --------------------------------------------------------------------------

template <typename TPos>
inline TPos posPrev(TPos pos) {
    return pos - 1;
}

template <typename T1, typename T2, typename TPack>
inline Pair<T1, T2, TPack> posPrev(Pair<T1, T2, TPack> const & pos) {
    return Pair<T1, T2, TPack>(getValueI1(pos), getValueI2(pos) - 1);
}

// --------------------------------------------------------------------------
// Function posNext()
// --------------------------------------------------------------------------

template <typename TPos>
inline TPos posNext(TPos pos) {
    return pos + 1;
}

template <typename T1, typename T2, typename TPack>
inline Pair<T1, T2, TPack>
posNext(Pair<T1, T2, TPack> const & pos) {
    return Pair<T1, T2, TPack>(getValueI1(pos), getValueI2(pos) + 1);
}

// --------------------------------------------------------------------------
// Function posAdd()
// --------------------------------------------------------------------------

template <typename TPos, typename TDelta>
inline TPos posAdd(TPos pos, TDelta delta) {
    return pos + delta;
}

template <typename T1, typename T2, typename TPack, typename TDelta>
inline Pair<T1, T2, TPack>
posAdd(Pair<T1, T2, TPack> const & pos, TDelta delta) {
    return Pair<T1, T2, TPack>(getValueI1(pos), getValueI2(pos) + delta);
}


// --------------------------------------------------------------------------
// Function posAddAndCheck()
// --------------------------------------------------------------------------

template <typename TPos, typename TDelta, typename TSequence>
inline TPos posAddAndCheck(TPos & pos, TDelta delta, TSequence const & sequence) {
    return (pos += delta) < length(sequence);
}

template <typename TPos, typename TDelta, typename TSequence, typename TSpec>
inline TPos posAddAndCheck(TPos & pos, TDelta delta, StringSet<TSequence, TSpec> const & stringSet)
{
    typedef StringSet<TSequence, TSpec> TStringSet;
    typedef typename StringSetLimits<TStringSet const>::Type TLimits;
    typedef typename Iterator<TLimits, Standard>::Type TIter;
    typedef typename Value<TLimits>::Type TSize;

    TLimits & limits = stringSetLimits(stringSet);
    TIter _end = end(limits, Standard());
    TIter _endMark = ::std::upper_bound(begin(limits, Standard()), _end, (TSize)pos);
    pos += delta;
    if (_endMark < _end)
        return pos < *_endMark;
    else
        return false;
}

template <typename T1, typename T2, typename TPack, typename TDelta, typename TSequence, typename TSpec>
inline bool
posAddAndCheck(Pair<T1, T2, TPack> & pos, TDelta delta, StringSet<TSequence, TSpec> const & stringSet) {
    return (pos.i2 += delta) < length(stringSet[pos.i1]);
}

// --------------------------------------------------------------------------
// Function posSub()
// --------------------------------------------------------------------------

template <typename TA, typename TB>
inline TA posSub(TA a, TB b) {
    return a - b;
}

template <
    typename TA1, typename TA2, typename TAPack,
    typename TB1, typename TB2, typename TBPack
>
inline TA2
posSub(Pair<TA1, TA2, TAPack> const & a, Pair<TB1, TB2, TBPack> const & b) {
    return getValueI2(a) - getValueI2(b);
}

// --------------------------------------------------------------------------
// Function posLess()
// --------------------------------------------------------------------------

template <typename TPos>
inline bool posLess(TPos const & a, TPos const & b) {
    return a < b;
}

template <typename T1, typename T2, typename TPack>
inline bool posLess(Pair<T1, T2, TPack> const & a, Pair<T1, T2, TPack> const & b) {
    return
         (getValueI1(a) <  getValueI1(b)) ||
        ((getValueI1(a) == getValueI1(b)) && (getValueI2(a) < getValueI2(b)));
}

// --------------------------------------------------------------------------
// Function posCompare()
// --------------------------------------------------------------------------

template <typename TPos>
inline int posCompare(TPos const & a, TPos const & b) {
    if (a < b) return -1;
    if (a > b) return 1;
    return 0;
}

template <typename T1, typename T2, typename TPack>
inline int posCompare(Pair<T1, T2, TPack> const & a, Pair<T1, T2, TPack> const & b) {
    if (getValueI1(a) < getValueI1(b)) return -1;
    if (getValueI1(a) > getValueI1(b)) return 1;
    return posCompare(getValueI2(a), getValueI2(b));
}

// --------------------------------------------------------------------------
// Function suffixLength()
// --------------------------------------------------------------------------

template <typename TPos, typename TString>
inline typename Size<TString>::Type
suffixLength(TPos pos, TString const & string) {
    return length(string) - pos;
}

template <typename TPos, typename TString, typename TSpec>
inline typename Size<TString>::Type
suffixLength(TPos pos, StringSet<TString, TSpec> const & stringSet) {
    return length(stringSet[getSeqNo(pos, stringSetLimits(stringSet))]) - getSeqOffset(pos, stringSetLimits(stringSet));
}

// --------------------------------------------------------------------------
// Function countSequences()
// --------------------------------------------------------------------------

template <typename TString>
inline unsigned
countSequences(TString const &) {
    return 1;
}

template <typename TString, typename TSpec>
inline typename Size<StringSet<TString, TSpec> >::Type
countSequences(StringSet<TString, TSpec> const & stringSet) {
    return length(stringSet);
}

// --------------------------------------------------------------------------
// Function getSequenceByNo()
// --------------------------------------------------------------------------

template <typename TSeqNo, typename TString>
inline typename GetSequenceByNo<TString>::Type
getSequenceByNo(TSeqNo /*seqNo*/, TString & string)
{
    return string;
}

template <typename TSeqNo, typename TString, typename TSpec>
inline typename GetSequenceByNo< StringSet<TString, TSpec> >::Type
getSequenceByNo(TSeqNo seqNo, StringSet<TString, TSpec> & stringSet)
{
    return stringSet[seqNo];
}

template <typename TSeqNo, typename TString, typename TSpec>
inline typename GetSequenceByNo< StringSet<TString, TSpec> const>::Type
getSequenceByNo(TSeqNo seqNo, StringSet<TString, TSpec> const & stringSet)
{
    return stringSet[seqNo];
}

// --------------------------------------------------------------------------
// Function sequenceLength()
// --------------------------------------------------------------------------

template <typename TSeqNo, typename TText>
inline typename Size< typename GetSequenceByNo<TText const>::Type>::Type
sequenceLength(TSeqNo seqNo, TText const & text)
{
    return length(getSequenceByNo(seqNo, text));
}

// --------------------------------------------------------------------------
// Function _validStringSetLimits
// --------------------------------------------------------------------------

// TODO(holtgrew): Anti auto-stringset
template < typename T >
inline bool _validStringSetLimits(T const &) {
    return true;
}

template < typename TString, typename TSpec >
inline bool _validStringSetLimits(StringSet< TString, TSpec > const & me) {
    return me.limitsValid;
}

// --------------------------------------------------------------------------
// Function _refreshStringSetLimits
// --------------------------------------------------------------------------

template < typename T >
inline void _refreshStringSetLimits(T &) {}

template < typename TString, typename TSpec >
void _refreshStringSetLimits(StringSet< TString, TSpec > & me)
{
    typedef StringSet< TString, TSpec >                 TStringSet;
    typedef typename StringSetLimits<TStringSet>::Type  TLimits;

    typename Value<TLimits>::Type   sum = 0;
    typename Size<TStringSet>::Type len = length(me);
    typename Size<TStringSet>::Type i = 0;

//      SEQAN_ASSERT_EQ(length(me.limits), len + 1);
    resize(me.limits, len + 1, Generous());
    for(; i < len; ++i)
    {
        me.limits[i] = sum;
        sum += length(me[i]);
        SEQAN_ASSERT_LEQ(me.limits[i], sum);
    }
    me.limits[i] = sum;
    me.limitsValid = true;
}

// --------------------------------------------------------------------------
// Function _findIthNonZeroValue()
// --------------------------------------------------------------------------

// find the i-th non-zero value of a string me
template < typename TValue, typename TSpec, typename TPos >
inline typename Size< String<TValue, TSpec> >::Type
_findIthNonZeroValue(String<TValue, TSpec> const & me, TPos i)
{
    typename Iterator< String<TValue, TSpec> const, Standard >::Type it = begin(me, Standard());
    typename Iterator< String<TValue, TSpec> const, Standard >::Type itEnd = end(me, Standard());

    for(; it != itEnd; ++it)
        if (*it)
        {
            if (i)
                --i;
            else
                return position(it, me);
        }
    return length(me);
}

// --------------------------------------------------------------------------
// Function _countNonZeroValues()
// --------------------------------------------------------------------------

// count non-zero values before position i
template < typename TValue, typename TSpec, typename TPos >
inline typename Size< String<TValue, TSpec> >::Type
_countNonZeroValues(String<TValue, TSpec> const & me, TPos i)
{
    typename Iterator< String<TValue, TSpec> const, Standard >::Type it = begin(me, Standard());
    typename Iterator< String<TValue, TSpec> const, Standard >::Type itEnd = begin(me, Standard()) + i;
    typename Size< String<TValue, TSpec> >::Type counter = 0;

    for(; it != itEnd; ++it)
        if (*it) ++counter;
    return counter;
}

// --------------------------------------------------------------------------
// Function lengthSum()
// --------------------------------------------------------------------------

template <typename TString>
inline typename LengthSum<TString>::Type 
lengthSum(TString const & me)
{
    return length(me);
}

template <typename TString, typename TSpec>
inline typename LengthSum<StringSet<TString, TSpec> >::Type
lengthSum(StringSet<TString, TSpec> & me)
{
    if (!_validStringSetLimits(me))
        _refreshStringSetLimits(me);
    return back(stringSetLimits(me));
}

template < typename TString, typename TSpec >
inline typename LengthSum<StringSet<TString, TSpec> >::Type
lengthSum(StringSet<TString, TSpec> const & me)
{
    if (!_validStringSetLimits(me))
        _refreshStringSetLimits(me);
    return back(stringSetLimits(me));
}


// --------------------------------------------------------------------------
// Function length()
// --------------------------------------------------------------------------

///.Function.appendValue.param.target.type:Class.StringSet
///.Function.appendValue.class:Class.StringSet
///.Function.clear.param.object.type:Class.StringSet
///.Function.clear.class:Class.StringSet
///.Function.resize.param.object.type:Class.StringSet
///.Function.resize.class:Class.StringSet
///.Function.length.param.object.type:Class.StringSet
///.Function.length.class:Class.StringSet

template <typename TString, typename TSpec >
inline typename Size<StringSet<TString, TSpec > >::Type
length(StringSet<TString, TSpec > const & me)
{
    return length(me.strings);
}

// --------------------------------------------------------------------------
// Function resize()
// --------------------------------------------------------------------------

// TODO(rmaerker): This belongs to string_set_base.h. Move it!
template <typename TString, typename TSpec, typename TSize, typename TExpand >
inline typename Size<StringSet<TString, TSpec > >::Type
resize(StringSet<TString, TSpec > & me, TSize new_size, Tag<TExpand> tag)
{
    resize(me.limits, new_size + 1, tag);
    me.limitsValid = (new_size == 0);

    //     the following would not work as changing the size of
    //     a single string cannot be recognized by the stringset
    //
    //        if (_validStringSetLimits(me))
    //            resize(me.limits, new_size + 1, back(me.limits), tag);
    
    return resize(me.strings, new_size, tag);
}

// --------------------------------------------------------------------------
// Function reserve()
// --------------------------------------------------------------------------

template <typename TString, typename TSpec, typename TSize, typename TExpand>
inline typename Size<StringSet<TString, TSpec > >::Type
reserve(StringSet<TString, TSpec > & me,
        TSize const & new_capacity,
        Tag<TExpand> tag)
{
    reserve(me.limits, new_capacity + 1, tag);
    return reserve(me.strings, new_capacity, tag);
}

// --------------------------------------------------------------------------
// Function iter()
// --------------------------------------------------------------------------

///.Function.iter.param.object.type:Class.StringSet
template <typename TString, typename TSpec, typename TPos, typename TTag>
inline typename Iterator<StringSet<TString, TSpec >, Tag<TTag> const>::Type
iter(StringSet<TString, TSpec > & me,
     TPos pos,
     Tag<TTag> const &)
{
    typedef StringSet<TString, TSpec > TStringSet;
    typedef typename Iterator<TStringSet, Tag<TTag> const>::Type TIterator;
    typedef typename Position<TStringSet>::Type TPosition;
    return TIterator(me, (TPosition) pos);
}

template <typename TString, typename TSpec, typename TPos, typename TTag>
inline typename Iterator<StringSet<TString, TSpec > const, Tag<TTag> const>::Type
iter(StringSet<TString, TSpec > const & me,
     TPos pos,
     Tag<TTag> const &)
{
    typedef StringSet<TString, TSpec > const TStringSet;
    typedef typename Iterator<TStringSet, Tag<TTag> const>::Type TIterator;
    typedef typename Position<TStringSet>::Type TPosition;
    return TIterator(me, (TPosition) pos);
}

// --------------------------------------------------------------------------
// Function begin()
// --------------------------------------------------------------------------

///.Function.begin.param.object.type:Class.StringSet
///.Function.begin.class:Class.StringSet

template <typename TString, typename TSpec, typename TTag>
inline typename Iterator<StringSet<TString, TSpec >, Tag<TTag> const>::Type
begin(StringSet<TString, TSpec > & me,
      Tag<TTag> const & tag)
{
    return iter(me, 0, tag);
}

template <typename TString, typename TSpec, typename TTag>
inline typename Iterator<StringSet<TString, TSpec > const, Tag<TTag> const>::Type
begin(StringSet<TString, TSpec > const & me,
      Tag<TTag> const & tag)
{
    return iter(me, 0, tag);
}

// --------------------------------------------------------------------------
// Function end()
// --------------------------------------------------------------------------

///.Function.end.param.object.type:Class.StringSet
///.Function.end.class:Class.StringSet

template <typename TString, typename TSpec, typename TTag>
inline typename Iterator<StringSet<TString, TSpec >, Tag<TTag> const>::Type
end(StringSet<TString, TSpec > & me,
Tag<TTag> const & tag)
{
return iter(me, length(me), tag);
}

template <typename TString, typename TSpec, typename TTag>
inline typename Iterator<StringSet<TString, TSpec > const, Tag<TTag> const>::Type
end(StringSet<TString, TSpec > const & me,
Tag<TTag> const & tag)
{
return iter(me, length(me), tag);
}

// --------------------------------------------------------------------------
// Function value()
// --------------------------------------------------------------------------

///.Function.value.param.object.type:Class.StringSet
///.Function.value.class:Class.StringSet

// --------------------------------------------------------------------------
// Function getValueById()
// --------------------------------------------------------------------------

/**
.Function.getValueById:
..cat:Sequences
..class:Class.StringSet
..summary:Retrieves a string from the StringSet given an id.
..signature:getValueById(me, id)
..param.me:A StringSet.
...type:Class.StringSet
..param.id:An id.
...type:Metafunction.Id
..returns:A reference to a string.
..see:Function.assignValueById
..see:Function.valueById
..include:seqan/sequence.h
*/

// TODO(holtgrew): Why is there no generic implementation for StringSets??

// --------------------------------------------------------------------------
// Function valueById()
// --------------------------------------------------------------------------

/**
.Function.valueById:
..cat:Sequences
..class:Class.StringSet
..summary:Retrieves a string from the StringSet given an id.
..signature:valueById(me, id)
..param.me:A StringSet.
...type:Class.StringSet
..param.id:An id.
...type:Metafunction.Id
..returns:A reference to a string.
..see:Function.assignValueById
..see:Function.getValueById
..include:seqan/sequence.h
*/

template<typename TString, typename TSpec, typename TId>
inline typename Reference<StringSet<TString, TSpec> >::Type
valueById(StringSet<TString, TSpec> & me,
        TId const id)
{
    SEQAN_CHECKPOINT;
    return getValueById(me, id);
}

// --------------------------------------------------------------------------
// Function assignValueById()
// --------------------------------------------------------------------------

/**
.Function.assignValueById:
..cat:Sequences
..class:Class.StringSet
..summary:Adds a new string to the StringSet and returns an id.
..signature:assignValueById(dest, str, [id])
..signature:assignValueById(dest, source, id)
..param.dest:A StringSet.
...type:Class.StringSet
..param.source:A StringSet.
...type:Class.StringSet
..param.str:A new string.
...type:Metafunction.Value
..param.id:An associated id.
...type:Metafunction.Id
..returns:A new id
...type:Metafunction.Id
..see:Function.getValueById
..see:Function.valueById
..include:seqan/sequence.h
*/

template<typename TString, typename TSpec, typename TString2>
inline typename Id<StringSet<TString, TSpec> >::Type
assignValueById(StringSet<TString, TSpec>& me,
                TString2& obj)
{
    SEQAN_CHECKPOINT;
    appendValue(me, obj);
    SEQAN_ASSERT_EQ(length(me.limits), length(me) + 1);
    return length(me.strings) - 1;
}

template<typename TString, typename TSpec1, typename TSpec2, typename TId>
inline typename Id<StringSet<TString, TSpec1> >::Type
assignValueById(StringSet<TString, TSpec1>& dest,
                StringSet<TString, TSpec2>& source,
                TId id)
{
    SEQAN_CHECKPOINT;
    return assignValueById(dest, getValueById(source, id), id);
}

// --------------------------------------------------------------------------
// Function removeValueById()
// --------------------------------------------------------------------------

/**
.Function.removeValueById:
..cat:Sequences
..class:Class.StringSet
..summary:Removes a string from the StringSet given an id.
..signature:removeValueById(me, id)
..param.me:A StringSet.
...type:Class.StringSet
..param.id:An id.
...type:Metafunction.Id
..returns:void
..see:Function.assignValueById
..include:seqan/sequence.h
*/

// TODO(holtgrew): Why is there no generic implementation for StringSets??

// --------------------------------------------------------------------------
// Function positionToId()
// --------------------------------------------------------------------------

/**
.Function.positionToId:
..cat:Sequences
..class:Class.StringSet
..summary:Retrieves the id of a string in the StringSet given a position.
..signature:positionToId(string_set, pos)
..param.string_set:A StringSet.
...type:Class.StringSet
..param.pos:A position that is transfored into an id.
..returns:An id that corresponds to $pos$ within $string_set$
..see:Function.assignValueById
..see:Function.valueById
..include:seqan/sequence.h
*/

// TODO(holtgrew): Why is there no generic implementation for StringSets??

// --------------------------------------------------------------------------
// Function concat()
// --------------------------------------------------------------------------

/**
.Function.concat:
..summary:Returns the concatenation sequence of all sequences in a @Class.StringSet@.
..cat:Sequences
..class:Class.StringSet
..signature:concat(stringSet)
..param.stringSet:A @Class.StringSet@ object.
...type:Class.StringSet
..returns:A container that can be iterated like the concatenation string of all sequences in a @Class.StringSet@.
..remarks:If $stringSet$ is a @Spec.ConcatDirect@ StringSet a reference to $stringSet.concat$ is returned.
For all other StringSets a @Class.ConcatenatorManyToOne@ object is returned.
...type:Metafunction.Concatenator
..include:seqan/sequence.h
*/

// TODO(holtgrew): Why default concat() for any class?
template <typename TString>
inline typename Concatenator<TString>::Type &
concat(TString & string)
{
    return string;
}

// TODO(holtgrew): Why default concat() for any class?
template <typename TString>
inline typename Concatenator<TString const>::Type &
concat(TString const & string)
{
    return string;
}

template <typename TString, typename TSpec>
inline typename Concatenator<StringSet<TString, TSpec> >::Type &
concat(StringSet<TString, TSpec> & me)
{
    me.concat.set = &me;
    return me.concat;
}

template <typename TString, typename TSpec>
inline typename Concatenator<StringSet<TString, TSpec> const>::Type &
concat(StringSet<TString, TSpec> const & constMe)
{
    StringSet<TString, TSpec> &me = const_cast<StringSet<TString, TSpec> &>(constMe);
    me.concat.set = &me;
    return me.concat;
}

// --------------------------------------------------------------------------
// Function strSplit()
// --------------------------------------------------------------------------

/**
.Function.stringSplit:
..summary:Append a list of the words in the string, using sep as the delimiter string @Class.StringSet@.
..cat:Sequences
..class:Class.StringSet
..signature:strSplit(stringSet, sequence)
..signature:strSplit(stringSet, sequence, sep)
..signature:strSplit(stringSet, sequence, sep, allowEmptyStrings)
..signature:strSplit(stringSet, sequence, sep, allowEmptyStrings, maxSplit)
..param.stringSet:The @Class.StringSet@ object the words are appended to.
...type:Class.StringSet
..param.sequence:A sequence of words.
..param.sep:Word separator (default: ' ').
..param.allowEmptyStrings:Boolean to specify whether empty words should be considered (default: true, iff sep is given).
..param.maxSplit:If maxsplit is given, at most maxsplit splits are done.
..include:seqan/sequence.h
*/

template <typename TString, typename TSpec, typename TSequence, typename TSeparator, typename TSize>
inline void
strSplit(StringSet<TString, TSpec> & result, TSequence const &sequence, TSeparator sep, bool allowEmptyStrings, TSize maxSplit)
{
    typedef typename Iterator<TSequence const, Standard>::Type TIter;
    
    TIter itBeg = begin(sequence, Standard());
    TIter itEnd = end(sequence, Standard());
    TIter itFrom = itBeg;
    
    if (maxSplit == 0)
    {
        appendValue(result, sequence);
        return;
    }
    
    for (TIter it = itBeg; it != itEnd; ++it)
        if (*it == sep)
        {
            if (allowEmptyStrings || itFrom != it)
            {
                appendValue(result, infix(sequence, itFrom - itBeg, it - itBeg));
                if (--maxSplit == 0)
                {
                    if (!allowEmptyStrings)
                    {
                        while (it != itEnd && *it == sep)
                            ++it;
                    }
                    else
                        ++it;
                    
                    if (it != itEnd)
                        appendValue(result, infix(sequence, it - itBeg, itEnd - itBeg));
                    
                    return;
                }
            }
            itFrom = it + 1;
        }
    
    if (allowEmptyStrings || itFrom != itEnd)
        appendValue(result, infix(sequence, itFrom - itBeg, itEnd - itBeg));
}

template <typename TString, typename TSpec, typename TSequence, typename TSeparator>
inline void
strSplit(StringSet<TString, TSpec> & result, TSequence const &sequence, TSeparator sep, bool allowEmptyStrings)
{
    strSplit(result, sequence, sep, allowEmptyStrings, maxValue<typename Size<TSequence>::Type>());
}

template <typename TString, typename TSpec, typename TSequence, typename TSeparator>
inline void
strSplit(StringSet<TString, TSpec> & result, TSequence const &sequence, TSeparator sep)
{
    strSplit(result, sequence, sep, true);
}

template <typename TString, typename TSpec, typename TSequence>
inline void
strSplit(StringSet<TString, TSpec> & result, TSequence const &sequence)
{
    strSplit(result, sequence, ' ', false);
}

// --------------------------------------------------------------------------
// Function idToPosition()
// --------------------------------------------------------------------------

/**
.Function.idToPosition:
..cat:Sequences
..class:Class.StringSet
..summary:Retrieves the position of a string in the StringSet given an id.
..signature:idToPosition(me, id)
..param.me:A StringSet.
...type:Class.StringSet
..param.id:An id.
...type:Metafunction.Id
..returns:A reference to a string.
..see:Function.assignValueById
..see:Function.valueById
..include:seqan/sequence.h
*/

// TODO(holtgrew): Why is there no generic implementation for StringSets??


// TODO(holtgrew): Should the following code be thrown away?

//template <typename TString, typename TSpec, typename TDestSpec, typename TIds, typename TLength>
//inline void
//subset(StringSet<TString, Owner<TSpec> >& source,
//    StringSet<TString, TDestSpec>& dest,
//    TIds ids,
//    TLength len)
//{
//SEQAN_CHECKPOINT
//}

//template <typename TString, typename TIds, typename TLength>
//inline void
//subset(StringSet<TString, Dependent<Generous> >& source,
//    StringSet<TString, Dependent<Generous> >& dest,
//    TIds ids,
//    TLength len)
//{
//SEQAN_CHECKPOINT
//    typedef StringSet<TString, Dependent<Generous> > TStringSet;
//    typedef typename Id<TStringSet>::Type TId;
//    typedef typename Size<TStringSet>::Type TSize;

//    clear(dest);
//    resize(dest.limits, len + 1);
//    dest.limitsValid = (len == 0);
//    resize(dest.strings, length(source.strings), (TString*) 0);
//    for(TSize i = 0; i < len; ++i)
//        dest.strings[ids[i]] = source.strings[ids[i]];
//}

//template <typename TString, typename TIds, typename TLength>
//inline void
//subset(StringSet<TString, Dependent<Tight> >& source,
//    StringSet<TString, Dependent<Tight> >& dest,
//    TIds ids,
//    TLength len)
//{
//SEQAN_CHECKPOINT
//    typedef StringSet<TString, Dependent<Tight> > TStringSet;
//    typedef typename Id<TStringSet>::Type TId;
//    typedef typename Size<TStringSet>::Type TSize;

//    clear(dest);
//    resize(dest.limits, len + 1);
//    dest.limitsValid = (len == 0);
//    TLength upperBound = length(source.ids);
//    for(TSize i=0;i<len;++i) {
//        TId id = ids[i];
//        if ((upperBound > id) &&
//            (source.ids[id] == id)) {
//                appendValue(dest.strings, source.strings[id]);
//                appendValue(dest.ids, id);
//        } else {
//            typedef String<TId> TIdString;
//            typedef typename Iterator<TIdString, Rooted>::Type TIter;
//            TIter it = begin(source.ids);
//            for(;!atEnd(it);goNext(it)) {
//                if (*it == id) {
//                    appendValue(dest.strings, source.strings[position(it)]);
//                    appendValue(dest.ids, id);
//                }
//            }
//        }
//    }
//}

//template <typename TString, typename TSpec, typename TIds>
//inline void
//subset(StringSet<TString, TSpec>& source,
//    StringSet<TString, TSpec>& dest,
//    TIds ids)
//{
//SEQAN_CHECKPOINT
//    subset(source, dest, ids, length(ids));
//}

}  // namespace seqan

#endif  // #ifndef SEQAN_SEQUENCE_STRING_SET_BASE_H_
