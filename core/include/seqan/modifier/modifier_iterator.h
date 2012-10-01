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

#ifndef SEQAN_HEADER_MODIFIER_ITERATOR_H
#define SEQAN_HEADER_MODIFIER_ITERATOR_H

namespace SEQAN_NAMESPACE_MAIN
{


//////////////////////////////////////////////////////////////////////////////

/**
.Class.ModifiedIterator:
..summary:Allows to modify arbitrary iterators by specializing what differs from an origin.
..cat:Modifier
..signature:ModifiedIterator<THost[, TSpec]>
..param.THost:Original iterator.
...type:Concept.RandomAccessIteratorConcept
..param.TSpec:The modifier type.
...metafunction:Metafunction.Spec
..implements:Concept.RandomAccessIteratorConcept
..remarks:$THost$ can also be a modified iterator, so you can create custom iterators by combining predefined ones.
..include:seqan/modifier.h
*/

	template < typename THost, typename TSpec = void >
	class ModifiedIterator {
	public:
		Holder<THost, Simple>					data_host;
		typename Cargo<ModifiedIterator>::Type	data_cargo;

		ModifiedIterator() {
            SEQAN_CHECKPOINT;
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

	template < typename THost, typename TSpec >
	struct Spec< ModifiedIterator<THost, TSpec> > {
		typedef TSpec Type;
	};

	template < typename THost, typename TSpec >
	struct Spec< ModifiedIterator<THost, TSpec> const > {
		typedef TSpec Type;
	};


	// an iterator is not the owner of the values pointing at
	// it can be constant while
	// - pointing to an alterable object
	// - returning an non-constant value
	// - being an iterator of an alterable container

	template < typename THost, typename TSpec >
	struct Value< ModifiedIterator<THost, TSpec> >:
		Value<THost> {};

	template < typename THost, typename TSpec >
	struct Value< ModifiedIterator<THost, TSpec> const >:
		Value< ModifiedIterator<THost, TSpec> > {};


	template < typename THost, typename TSpec >
	struct GetValue< ModifiedIterator<THost, TSpec> >:
		GetValue<THost> {};

	template < typename THost, typename TSpec >
	struct GetValue< ModifiedIterator<THost, TSpec> const >:
		GetValue< ModifiedIterator<THost, TSpec> > {};


	template < typename THost, typename TSpec >
	struct Reference< ModifiedIterator<THost, TSpec> >:
		Reference<THost> {};

	template < typename THost, typename TSpec >
	struct Reference< ModifiedIterator<THost, TSpec> const >:
		Reference< ModifiedIterator<THost, TSpec> > {};

	template < typename THost, typename TSpec >
	struct Size< ModifiedIterator<THost, TSpec> >:
		Size<THost> {};

	template < typename THost, typename TSpec >
	struct Position< ModifiedIterator<THost, TSpec> >:
		Position<THost> {};

	template < typename THost, typename TSpec >
	struct Difference< ModifiedIterator<THost, TSpec> >:
		Difference<THost> {};


	template < typename THost, typename TSpec >
	struct Host< ModifiedIterator<THost, TSpec> > {
		typedef THost Type;
	};

	template < typename THost, typename TSpec >
	struct Host< ModifiedIterator<THost, TSpec> const > {
		typedef THost const Type;
	};


	//template < typename THost, typename TSpec >
	//struct Container< ModifiedIterator<THost, TSpec> >:
	//	Container<THost> {};

	//template < typename THost, typename TSpec >
	//struct Container< ModifiedIterator<THost, TSpec> const >:
	//	Container< ModifiedIterator<THost, TSpec> > {};

	template <typename THost, typename TSpec>
	class ModifiedString;

	template <typename THost, typename TSpec >
	struct Container< ModifiedIterator<THost, TSpec> >
	{
		typedef typename Container<THost>::Type THostContainer;
		typedef ModifiedString<THostContainer, TSpec> Type;
	};
	template <typename THost, typename TSpec >
	struct Container< ModifiedIterator<THost, TSpec> const>
	{
		typedef typename Container<THost>::Type THostContainer;
		typedef ModifiedString<THostContainer, TSpec> Type;
	};

	//////////////////////////////////////////////////////////////////////////////
	// host interface
	//////////////////////////////////////////////////////////////////////////////

	template <typename THost, typename TSpec>
	inline Holder<THost, Simple> &
	_dataHost(ModifiedIterator<THost, TSpec> & me) 
	{
        SEQAN_CHECKPOINT;
		return me.data_host;
	}
	
	template <typename THost, typename TSpec>
	inline Holder<THost, Simple> const &
	_dataHost(ModifiedIterator<THost, TSpec> const & me) 
	{
        SEQAN_CHECKPOINT;
		return me.data_host;
	}

	template <typename THost, typename TSpec>
	inline typename Reference< typename Cargo<ModifiedIterator<THost, TSpec> >::Type >::Type
	cargo(ModifiedIterator<THost, TSpec> & me) 
	{
        SEQAN_CHECKPOINT;
		return me.data_cargo;
	}

	template <typename THost, typename TSpec>
	inline typename Reference< typename Cargo<ModifiedIterator<THost, TSpec> const>::Type >::Type
	cargo(ModifiedIterator<THost, TSpec> const & me) 
	{
        SEQAN_CHECKPOINT;
		return me.data_cargo;
	}

	//////////////////////////////////////////////////////////////////////////////
	// container/setContainer interface
	//////////////////////////////////////////////////////////////////////////////

	template <typename THost, typename TSpec>
	inline typename Container<ModifiedIterator<THost, TSpec> >::Type //no reference
	container(ModifiedIterator<THost, TSpec> & me) 
	{
        SEQAN_CHECKPOINT;
		typedef typename Container<ModifiedIterator<THost, TSpec> >::Type TContainer;
		TContainer temp_(container(host(me)));
		_copyCargo(temp_, me);
		return temp_;
	}

	template <typename THost, typename TSpec>
	inline typename Container<ModifiedIterator<THost, TSpec> const>::Type //no reference
	container(ModifiedIterator<THost, TSpec> const & me) 
	{
        SEQAN_CHECKPOINT;
		typedef typename Container<ModifiedIterator<THost, TSpec> const>::Type TContainer;
		TContainer temp_(container(host(me)));
		_copyCargo(temp_, me);
		return temp_;
	}

	//////////////////////////////////////////////////////////////////////////////

	template <typename TIteratorHost, typename TSpec, typename TStringHost>
	inline void
	setContainer(
		ModifiedIterator<TIteratorHost, TSpec> & me, 
		ModifiedString<TStringHost, TSpec> & cont) 
	{
        SEQAN_CHECKPOINT;
		setContainer(host(me), host(cont));
		_copyCargo(me, cont);
	}
	template <typename TIteratorHost, typename TSpec, typename TStringHost>
	inline void
	setContainer(
		ModifiedIterator<TIteratorHost, TSpec> & me, 
		ModifiedString<TStringHost, TSpec> const & cont) 
	{
        SEQAN_CHECKPOINT;
		setContainer(host(me), host(const_cast<ModifiedString<TStringHost, TSpec> &>(cont)));
		_copyCargo(me, cont);
	}
	template <typename THost, typename TSpec, typename TContainer>
	inline void
	setContainer(ModifiedIterator<THost, TSpec> & me, TContainer & cont) 
	{
        SEQAN_CHECKPOINT;
		setContainer(host(me), cont);
	}
/*	template <typename THost, typename TSpec, typename TContainer>
	inline void
	setContainer(ModifiedIterator<THost, TSpec> & me, TContainer const & cont) 
	{
	SEQAN_CHECKPOINT
		THost &_host = host(me);
		setContainer(_host, host(cont));
		_copyCargo(me, cont);
	}
*/
	//////////////////////////////////////////////////////////////////////////////
	// assign
	//////////////////////////////////////////////////////////////////////////////
    
    template <typename TTarget, typename TSource>
    inline void 
    _assignModifiedIterator(TTarget &me, TSource &_origin, True)
    {
		host(me) = _origin;
    }

    template <typename TTarget, typename TSource>
    inline void 
    _assignModifiedIterator(TTarget &me, TSource &_origin, False)
    {
		host(me) = host(_origin);
		cargo(me) = cargo(_origin);
    }

	template <typename THost, typename TSpec, typename THost2>
	inline ModifiedIterator<THost, TSpec> const &
	assign(ModifiedIterator<THost, TSpec> & me, ModifiedIterator<THost2, TSpec> & _origin) {
        SEQAN_CHECKPOINT;
        _assignModifiedIterator(me, _origin, typename IsSameType<THost, ModifiedIterator<THost2, TSpec> >::Type());
		return me;
	}

	template <typename THost, typename TSpec, typename THost2>
	inline ModifiedIterator<THost, TSpec> const &
	assign(ModifiedIterator<THost, TSpec> & me, ModifiedIterator<THost2, TSpec> const & _origin) {
        SEQAN_CHECKPOINT;
        _assignModifiedIterator(me, _origin, typename IsSameType<THost, ModifiedIterator<THost2, TSpec> >::Type());
		return me;
	}

	template <typename THost, typename TSpec, typename T>
	inline ModifiedIterator<THost, TSpec> const &
	assign(ModifiedIterator<THost, TSpec> & me, T & _origin) {
        SEQAN_CHECKPOINT;
		host(me) = _origin;
		return me;
	}

	template <typename THost, typename TSpec, typename T>
	inline ModifiedIterator<THost, TSpec> const &
	assign(ModifiedIterator<THost, TSpec> & me, T const & _origin) {
        SEQAN_CHECKPOINT;
		host(me) = _origin;
		return me;
	}

	//////////////////////////////////////////////////////////////////////////////
	// value
	//////////////////////////////////////////////////////////////////////////////

	// redefinition candidate
	template <typename THost, typename TSpec>
	inline typename Reference<ModifiedIterator<THost, TSpec> >::Type 
	value(ModifiedIterator<THost, TSpec> & me)
	{
        SEQAN_CHECKPOINT;
		return value(host(me));
	}

	// redefinition candidate
	template <typename THost, typename TSpec>
	inline typename Reference<ModifiedIterator<THost, TSpec> const>::Type 
	value(ModifiedIterator<THost, TSpec> const & me)
	{
        SEQAN_CHECKPOINT;
		return value(host(me));
	}

	template <typename THost, typename TSpec>
	inline typename Reference<ModifiedIterator<THost, TSpec> >::Type 
	operator * (ModifiedIterator<THost, TSpec> & me)
	{
        SEQAN_CHECKPOINT;
		return value(me);
	}

	template <typename THost, typename TSpec>
	inline typename Reference<ModifiedIterator<THost, TSpec> const>::Type 
	operator * (ModifiedIterator<THost, TSpec> const & me)
	{
        SEQAN_CHECKPOINT;
		return value(me);
	}


	//////////////////////////////////////////////////////////////////////////////
	// operator ++
	//////////////////////////////////////////////////////////////////////////////

	// redefinition candidate
	template <typename THost, typename TSpec>
	inline void
	goNext(ModifiedIterator<THost, TSpec> & me)
	{
	SEQAN_CHECKPOINT
		goNext(host(me));
	}

	template <typename THost, typename TSpec>
	inline ModifiedIterator<THost, TSpec> const &
	operator ++ (ModifiedIterator<THost, TSpec> & me)
	{
	SEQAN_CHECKPOINT
		goNext(me);
		return me;
	}

	template <typename THost, typename TSpec>
	inline ModifiedIterator<THost, TSpec>
	operator ++ (ModifiedIterator<THost, TSpec> & me, int)
	{
	SEQAN_CHECKPOINT
		ModifiedIterator<THost, TSpec> temp_(me);
		goNext(me);
		return temp_;
	}

	//////////////////////////////////////////////////////////////////////////////
	// operator --
	//////////////////////////////////////////////////////////////////////////////

	// redefinition candidate
	template <typename THost, typename TSpec>
	inline void
	goPrevious(ModifiedIterator<THost, TSpec> & me)
	{
        SEQAN_CHECKPOINT;
		goPrevious(host(me));
	}

	template <typename THost, typename TSpec>
	inline ModifiedIterator<THost, TSpec> const &
	operator -- (ModifiedIterator<THost, TSpec> & me)
	{
        SEQAN_CHECKPOINT;
		goPrevious(me);
		return me;
	}

	template <typename THost, typename TSpec>
	inline ModifiedIterator<THost, TSpec>
	operator -- (ModifiedIterator<THost, TSpec> & me, int)
	{
        SEQAN_CHECKPOINT;
		ModifiedIterator<THost, TSpec> temp_(me);
		goPrevious(me);
		return temp_;
	}

	//////////////////////////////////////////////////////////////////////////////
	// operator +
	//////////////////////////////////////////////////////////////////////////////

	// redefinition candidate
	template <typename THost, typename TSpec, typename TDelta>
	inline ModifiedIterator<THost, TSpec> &
	operator += (ModifiedIterator<THost, TSpec> & me, TDelta delta) {
        SEQAN_CHECKPOINT;
		host(me) += delta;
		return me;
	}

	template <typename THost, typename TSpec, typename TDelta>
	inline ModifiedIterator<THost, TSpec>
	operator + (ModifiedIterator<THost, TSpec> const & me, TDelta delta) {
        SEQAN_CHECKPOINT;
		ModifiedIterator<THost, TSpec> temp_(me);
		temp_ += delta;
		return temp_;
	}

	//////////////////////////////////////////////////////////////////////////////
	// operator -
	//////////////////////////////////////////////////////////////////////////////

	// redefinition candidate
	template <typename THost, typename TSpec, typename TDelta>
	inline ModifiedIterator<THost, TSpec> &
	operator -= (ModifiedIterator<THost, TSpec> & me, TDelta delta) {
        SEQAN_CHECKPOINT;
		host(me) -= delta;
		return me;
	}

	template <typename THost, typename TSpec, typename TDelta>
	inline ModifiedIterator<THost, TSpec>
	operator - (ModifiedIterator<THost, TSpec> const & me, TDelta delta) {
        SEQAN_CHECKPOINT;
		ModifiedIterator<THost, TSpec> temp_(me);
		temp_ -= delta;
		return temp_;
	}

	// redefinition candidate
	template <typename THost, typename TSpec>
	inline typename Difference< ModifiedIterator<THost, TSpec> >::Type
	operator - (ModifiedIterator<THost, TSpec> const & a, ModifiedIterator<THost, TSpec> const & b) {
        SEQAN_CHECKPOINT;
		return host(a) - host(b);
	}

	//////////////////////////////////////////////////////////////////////////////
	// goBegin
	//////////////////////////////////////////////////////////////////////////////

	template <typename THost, typename TSpec, typename TContainer>
	inline void
	goBegin(ModifiedIterator<THost, TSpec> & me,
			TContainer const & container)
	{
        SEQAN_CHECKPOINT;
		host(me) = begin(container);
	}

	template <typename THost, typename TSpec>
	inline void
	goBegin(ModifiedIterator<THost, TSpec> & me)
	{
        SEQAN_CHECKPOINT;
		goBegin(me, container(me));
	}

	//////////////////////////////////////////////////////////////////////////////
	// goEnd
	//////////////////////////////////////////////////////////////////////////////

	template <typename THost, typename TSpec, typename TContainer>
	inline void
	goEnd(ModifiedIterator<THost, TSpec> & me,
			TContainer const & container)
	{
        SEQAN_CHECKPOINT;
		host(me) = end(container);
	}

	template <typename THost, typename TSpec>
	inline void
	goEnd(ModifiedIterator<THost, TSpec> & me)
	{
        SEQAN_CHECKPOINT;
		goEnd(me, container(me));
	}

	//////////////////////////////////////////////////////////////////////////////
	// position
	//////////////////////////////////////////////////////////////////////////////

	// redefinition candidate
	template <typename THost, typename TSpec>
	inline typename Position<ModifiedIterator<THost, TSpec> const>::Type 
	position(ModifiedIterator<THost, TSpec> const & me)
	{
        SEQAN_CHECKPOINT;
		return position(host(me));
	}

	// redefinition candidate
	template <typename THost, typename TSpec, typename TContainer>
	inline typename Position<ModifiedIterator<THost, TSpec> const>::Type 
	position(ModifiedIterator<THost, TSpec> const & me, TContainer const &cont)
	{
        SEQAN_CHECKPOINT;
		return position(host(me), cont);
	}

	//////////////////////////////////////////////////////////////////////////////
	// operator ==
	//////////////////////////////////////////////////////////////////////////////

	// redefinition candidate
	template <typename THost, typename TSpec>
	inline bool
	operator == (ModifiedIterator<THost, TSpec> const & a, ModifiedIterator<THost, TSpec> const & b) {
		return host(a) == host(b);
	}

	template <typename THost, typename TSpec>
	inline bool
	operator != (ModifiedIterator<THost, TSpec> const & a, ModifiedIterator<THost, TSpec> const & b) {
		return !(a == b);
	}

	//////////////////////////////////////////////////////////////////////////////
	// operator <
	//////////////////////////////////////////////////////////////////////////////

	// redefinition candidate
	template <typename THost, typename TSpec>
	inline bool
	operator < (ModifiedIterator<THost, TSpec> const & a, ModifiedIterator<THost, TSpec> const & b) {
        SEQAN_CHECKPOINT;
		return host(a) < host(b);
	}

	template <typename THost, typename TSpec>
	inline bool
	operator > (ModifiedIterator<THost, TSpec> const & a, ModifiedIterator<THost, TSpec> const & b) {
        SEQAN_CHECKPOINT;
		return b < a;
	}

	//////////////////////////////////////////////////////////////////////////////
	// atBegin
	//////////////////////////////////////////////////////////////////////////////

	template <typename THost, typename TSpec, typename TContainer>
	inline bool
	atBegin(ModifiedIterator<THost, TSpec> & me,
			TContainer const & container)
	{
        SEQAN_CHECKPOINT;
		return atBegin(const_cast<ModifiedIterator<THost, TSpec> const &>(me), container);
	}

	// redefinition candidate
	template <typename THost, typename TSpec, typename TContainer>
	inline bool
	atBegin(ModifiedIterator<THost, TSpec> const & me,
			TContainer const & container)
	{
        SEQAN_CHECKPOINT;
		return atBegin(host(me), container);
	}

	template <typename THost, typename TSpec>
	inline bool
	atBegin(ModifiedIterator<THost, TSpec> & me)
	{
        SEQAN_CHECKPOINT;
		return atBegin(const_cast<ModifiedIterator<THost, TSpec> const &>(me));
	}

	// redefinition candidate
	template <typename THost, typename TSpec>
	inline bool
	atBegin(ModifiedIterator<THost, TSpec> const & me)
	{
        SEQAN_CHECKPOINT;
		return atBegin(host(me));
	}

	//////////////////////////////////////////////////////////////////////////////
	// atEnd
	//////////////////////////////////////////////////////////////////////////////

	template <typename THost, typename TSpec, typename TContainer>
	inline bool
	atEnd(ModifiedIterator<THost, TSpec> & me,
			TContainer const & container)
	{
        SEQAN_CHECKPOINT;
		return atEnd(const_cast<ModifiedIterator<THost, TSpec> const &>(me), container);
	}

	// redefinition candidate
	template <typename THost, typename TSpec, typename TContainer>
	inline bool
	atEnd(ModifiedIterator<THost, TSpec> const & me,
			TContainer const & container)
	{
        SEQAN_CHECKPOINT;
		return atEnd(host(me), container);
	}

	template <typename THost, typename TSpec>
	inline bool
	atEnd(ModifiedIterator<THost, TSpec> & me)
	{
        SEQAN_CHECKPOINT;
		return atEnd(const_cast<ModifiedIterator<THost, TSpec> const &>(me));
	}

	// redefinition candidate
	template <typename THost, typename TSpec>
	inline bool
	atEnd(ModifiedIterator<THost, TSpec> const & me)
	{
        SEQAN_CHECKPOINT;
		return atEnd(host(me));
	}

}

// Adapt SeqAn modified to std.
namespace std
{
	template<typename THost, typename TSpec>
	struct iterator_traits<seqan::ModifiedIterator<THost, TSpec> >
	{
		typedef ::seqan::ModifiedIterator<THost, TSpec> TIter;

		typedef random_access_iterator_tag iterator_category;
		typedef typename ::seqan::Value<TIter>::Type value_type;
		typedef typename ::seqan::Difference<TIter>::Type difference_type;
		typedef typename ::seqan::Value<TIter>::Type * pointer;
		typedef typename ::seqan::Reference<TIter>::Type reference;
	};
}

#endif
