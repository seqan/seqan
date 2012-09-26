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

#ifndef SEQAN_MODIFIER_MODIFIER_VIEW_H_
#define SEQAN_MODIFIER_MODIFIER_VIEW_H_

namespace seqan {

/**
.Spec.ModView:
..summary:Transforms the characters of the $THost$ string/iterator using a custom function.
..cat:Modifier
..general:Class.ModifiedIterator
..general:Class.ModifiedString
..signature:ModifiedIterator<THost, ModView<TFunctor> >
..signature:ModifiedString<THost, ModView<TFunctor> >
..param.THost:Original string/iterator.
...type:Concept.RandomAccessIteratorConcept
..param.TFunctor:A unary function (see STL's $unary_function$).
...remarks:The argument type of $TFunctor$ must be $VALUE<THost>::Type$.
..remarks:The @Metafunction.Value@ type of this modifier is the result type of $TFunctor$.
..include:seqan/modifier.h
*/

template <typename TFunctor>
struct ModView {};

template <typename TFunctor>
struct ModViewCargo {
    TFunctor	func;
};


//////////////////////////////////////////////////////////////////////////////
// view iterator
//////////////////////////////////////////////////////////////////////////////


template <typename THost, typename TFunctor>
struct Cargo< ModifiedIterator<THost, ModView<TFunctor> > > {
    typedef ModViewCargo<TFunctor>	Type;
};


template <typename THost, typename TFunctor>
class ModifiedIterator<THost, ModView<TFunctor> > {
public:
    Holder<THost, Simple>					data_host;
    typename Cargo<ModifiedIterator>::Type	data_cargo;

    mutable typename Value<ModifiedIterator>::Type	tmp_value;

    ModifiedIterator() {}

    explicit ModifiedIterator(TFunctor &_func) {
        SEQAN_CHECKPOINT;
        assignModViewFunctor(*this, _func);
    }

    explicit ModifiedIterator(TFunctor const &_func) {
        SEQAN_CHECKPOINT;
        assignModViewFunctor(*this, _func);
    }

    ModifiedIterator(ModifiedIterator &_origin):
        data_host(_origin.data_host),
        data_cargo(_origin.data_cargo) {
        SEQAN_CHECKPOINT;
    }

    ModifiedIterator(ModifiedIterator const &_origin):
        data_host(_origin.data_host),
        data_cargo(_origin.data_cargo) {
        SEQAN_CHECKPOINT;
    }

    template <typename T>
    ModifiedIterator(T & _origin) {
        SEQAN_CHECKPOINT;
        assign(*this, _origin);
    }

    template <typename T>
    ModifiedIterator(T const & _origin) {
        SEQAN_CHECKPOINT;
        assign(*this, _origin);
    }
//____________________________________________________________________________

    template <typename T>
    inline ModifiedIterator const &
    operator = (T & _origin) {
        SEQAN_CHECKPOINT;
        assign(*this, _origin);
        return *this;
    }

    template <typename T>
    inline ModifiedIterator const &
    operator = (T const & _origin) {
        SEQAN_CHECKPOINT;
        assign(*this, _origin);
        return *this;
    }
};

template <typename THost, typename TFunctor>
struct Value< ModifiedIterator<THost, ModView<TFunctor> > > {
    typedef typename TFunctor::result_type			TResult;
    typedef typename RemoveConst_<TResult>::Type	Type;
};

template <typename THost, typename TFunctor>
struct GetValue< ModifiedIterator<THost, ModView<TFunctor> > >:
    Value< ModifiedIterator<THost, ModView<TFunctor> > > {};

template <typename THost, typename TFunctor>
struct Reference< ModifiedIterator<THost, ModView<TFunctor> > > {
    typedef typename Value< ModifiedIterator<THost, ModView<TFunctor> > >::Type & Type;
};


//////////////////////////////////////////////////////////////////////////////
// value
//////////////////////////////////////////////////////////////////////////////

template <typename THost, typename TFunctor>
inline typename Reference<ModifiedIterator<THost, ModView<TFunctor> > >::Type 
value(ModifiedIterator<THost, ModView<TFunctor> > & me)
{
    SEQAN_CHECKPOINT;
    me.tmp_value = cargo(me).func(getValue(host(me)));
    return me.tmp_value;
}

template <typename THost, typename TFunctor>
inline typename Reference<ModifiedIterator<THost, ModView<TFunctor> > const>::Type 
value(ModifiedIterator<THost, ModView<TFunctor> > const & me)
{
    SEQAN_CHECKPOINT;
    me.tmp_value = cargo(me).func(getValue(host(me)));
    return me.tmp_value;
}


//////////////////////////////////////////////////////////////////////////////
// getValue
//////////////////////////////////////////////////////////////////////////////

template <typename THost, typename TFunctor>
inline typename GetValue<ModifiedIterator<THost, ModView<TFunctor> > >::Type 
getValue(ModifiedIterator<THost, ModView<TFunctor> > & me)
{
    SEQAN_CHECKPOINT;
    return cargo(me).func(getValue(host(me)));
}

template <typename THost, typename TFunctor>
inline typename GetValue<ModifiedIterator<THost, ModView<TFunctor> > const>::Type 
getValue(ModifiedIterator<THost, ModView<TFunctor> > const & me)
{
    SEQAN_CHECKPOINT;
    return cargo(me).func(getValue(host(me)));
}


//////////////////////////////////////////////////////////////////////////////
// assignModViewFunctor
//////////////////////////////////////////////////////////////////////////////

template <typename THost, typename TFunctor>
inline void
assignModViewFunctor(ModifiedIterator<THost, ModView<TFunctor> > & me, TFunctor const & _func) 
{
    SEQAN_CHECKPOINT;
    cargo(me).func = _func;
}


//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
// view string
//////////////////////////////////////////////////////////////////////////////


template <typename THost, typename TFunctor>
struct Cargo< ModifiedString<THost, ModView<TFunctor> > > {
    typedef ModViewCargo<TFunctor>	Type;
};

template <typename THost, typename TFunctor>
class ModifiedString<THost, ModView<TFunctor> > {
public:
    Holder<THost>							data_host;
    typename Cargo<ModifiedString>::Type	data_cargo;

    mutable typename Value<ModifiedString>::Type	tmp_value;

    ModifiedString() {}

    explicit ModifiedString(TFunctor &_func) {
        SEQAN_CHECKPOINT;
        cargo(*this).func = _func;
    }

    explicit ModifiedString(TFunctor const &_func) {
        SEQAN_CHECKPOINT;
        cargo(*this).func = _func;
    }

    explicit ModifiedString(ModifiedString const &_origin, TFunctor const &_func):
        data_host(_origin.data_host)
    {
        SEQAN_CHECKPOINT;
        cargo(*this).func = _func;
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
    operator = (T & _origin) {
        SEQAN_CHECKPOINT;
        assign(*this, _origin);
        return *this;
    }

    template <typename TPos>
    inline typename Reference<ModifiedString>::Type 
    operator [] (TPos pos)
    {
        SEQAN_CHECKPOINT;
        return value(*this, pos);
    }

    template <typename TPos>
    inline typename Reference<ModifiedString const>::Type 
    operator [] (TPos pos) const
    {
        SEQAN_CHECKPOINT;
        return value(*this, pos);
    }
};


//////////////////////////////////////////////////////////////////////////////
// value
//////////////////////////////////////////////////////////////////////////////

template <typename THost, typename TFunctor, typename TPos>
inline typename Reference<ModifiedString<THost, ModView<TFunctor> > >::Type 
value(ModifiedString<THost, ModView<TFunctor> > & me, TPos pos)
{
    SEQAN_CHECKPOINT;
    me.tmp_value = cargo(me).func(getValue(host(me), pos));
    return me.tmp_value;
}

template <typename THost, typename TFunctor, typename TPos>
inline typename Reference<ModifiedString<THost, ModView<TFunctor> > const>::Type 
value(ModifiedString<THost, ModView<TFunctor> > const & me, TPos pos)
{
    SEQAN_CHECKPOINT;
    me.tmp_value = cargo(me).func(getValue(host(me), pos));
    return me.tmp_value;
}


//////////////////////////////////////////////////////////////////////////////
// getValue
//////////////////////////////////////////////////////////////////////////////

template <typename THost, typename TFunctor, typename TPos>
inline typename GetValue<ModifiedString<THost, ModView<TFunctor> > >::Type 
getValue(ModifiedString<THost, ModView<TFunctor> > & me, TPos pos)
{
    SEQAN_CHECKPOINT;
    return cargo(me).func(getValue(host(me), pos));
}

template <typename THost, typename TFunctor, typename TPos>
inline typename GetValue<ModifiedString<THost, ModView<TFunctor> > const>::Type 
getValue(ModifiedString<THost, ModView<TFunctor> > const & me, TPos pos)
{
    SEQAN_CHECKPOINT;
    return cargo(me).func(getValue(host(me), pos));
}


//////////////////////////////////////////////////////////////////////////////
// assignModViewFunctor
//////////////////////////////////////////////////////////////////////////////

template <typename THost, typename TFunctor>
inline void
assignModViewFunctor(ModifiedString<THost, ModView<TFunctor> > & me, TFunctor const & _func)
{
    SEQAN_CHECKPOINT;
    cargo(me).func = _func;
}


//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
// convert
//////////////////////////////////////////////////////////////////////////////

template < typename TSequence, typename TFunctor >
inline void
convert(TSequence & sequence, TFunctor const &F)
{
	SEQAN_CHECKPOINT;
#if defined (_OPENMP) && defined (SEQAN_PARALLEL)
	// OpenMP does not support for loop with iterators. Therefore use index variables.
	typedef typename Position<TSequence>::Type	TPos;
	typedef typename MakeSigned_<TPos>::Type	TSignedPos;

	#pragma omp parallel for if(length(sequence) > 1000000)
	for(TSignedPos p = 0; p < (TSignedPos)length(sequence); ++p)
		sequence[p] = F(sequence[p]);
	
#else
	typedef typename Iterator<TSequence, Standard>::Type	TIter;

	TIter it = begin(sequence, Standard());
	TIter itEnd = end(sequence, Standard());
	for(; it != itEnd; ++it)
		*it = F(*it);
#endif
}

template < typename TSequence, typename TFunctor >
inline void
convert(TSequence const & sequence, TFunctor const &F)
{
	SEQAN_CHECKPOINT;
#if defined (_OPENMP) && defined (SEQAN_PARALLEL)
	// OpenMP does not support for loop with iterators. Therefore use index variables.
	typedef typename Position<TSequence>::Type	TPos;
	typedef typename MakeSigned_<TPos>::Type	TSignedPos;

	#pragma omp parallel for if(length(sequence) > 1000000)
	for(TSignedPos p = 0; p < (TSignedPos)length(sequence); ++p)
		sequence[p] = F(sequence[p]);
	
#else
	typedef typename Iterator<TSequence, Standard>::Type	TIter;

	TIter it = begin(sequence, Standard());
	TIter itEnd = end(sequence, Standard());
	for(; it != itEnd; ++it)
		*it = F(*it);
#endif
}

}  // namespace seqan

#endif  // SEQAN_MODIFIER_MODIFIER_VIEW_H_
