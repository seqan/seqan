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

#ifndef SEQAN_HEADER_MODIFIER_REVERSE_H
#define SEQAN_HEADER_MODIFIER_REVERSE_H

#ifdef _OPENMP
#include <omp.h>
#endif

namespace SEQAN_NAMESPACE_MAIN
{


//////////////////////////////////////////////////////////////////////////////

/**
.Spec.ModReverse:
..summary:Mirrors the characters from begin to end.
..cat:Modifier
..general:Class.ModifiedIterator
..general:Class.ModifiedString
..signature:ModifiedIterator<THost, ModReverse>
..signature:ModifiedString<THost, ModReverse>
..param.THost:Original string/iterator.
...type:Concept.RandomAccessIteratorConcept
..include:seqan/modifier.h
*/

	struct ModReverse {};

//////////////////////////////////////////////////////////////////////////////


	//////////////////////////////////////////////////////////////////////////////
	// reverse iterator
	//////////////////////////////////////////////////////////////////////////////


	template <typename THost>
	struct Cargo< ModifiedIterator<THost, ModReverse> > {
		typedef Cargo Type;		// to reduce namespace pollution
		bool _atEnd;
		Cargo(): _atEnd(false) {}
	};

	template <typename THost>
	class ModifiedIterator<THost, ModReverse> {
	public:
		Holder<THost, Simple>					data_host;
		typename Cargo<ModifiedIterator>::Type	data_cargo;

		ModifiedIterator() {}
		ModifiedIterator(ModifiedIterator &_origin):
			data_host(_origin.data_host),
			data_cargo(_origin.data_cargo) {}

		ModifiedIterator(ModifiedIterator const &_origin):
			data_host(_origin.data_host),
			data_cargo(_origin.data_cargo) {}

		template <typename T>
		ModifiedIterator(T & _origin) {
			assign(*this, _origin);
		}

		template <typename T>
		ModifiedIterator(T const & _origin) {
			assign(*this, _origin);
		}
//____________________________________________________________________________

		template <typename T>
		inline ModifiedIterator const &
		operator = (T & _origin) {
			assign(*this, _origin);
			return *this;
		}

		template <typename T>
		inline ModifiedIterator const &
		operator = (T const & _origin) {
			assign(*this, _origin);
			return *this;
		}
	};

	//////////////////////////////////////////////////////////////////////////////
	// operator ++
	//////////////////////////////////////////////////////////////////////////////

	template <typename THost>
	inline void
	goNext(ModifiedIterator<THost, ModReverse> & me)
	{
	SEQAN_CHECKPOINT
		if (atBegin(host(me)))
			cargo(me)._atEnd = true;
		else
			goPrevious(host(me));
	}

	//////////////////////////////////////////////////////////////////////////////
	// operator --
	//////////////////////////////////////////////////////////////////////////////

	template <typename THost>
	inline void
	goPrevious(ModifiedIterator<THost, ModReverse> & me)
	{
	SEQAN_CHECKPOINT
		if (cargo(me)._atEnd)
			cargo(me)._atEnd = false;
		else
			goNext(host(me));
	}

	//////////////////////////////////////////////////////////////////////////////
	// goEnd
	//////////////////////////////////////////////////////////////////////////////

	template <typename THost>
	inline void
	goEnd(ModifiedIterator<THost, ModReverse> & me)
	{
	SEQAN_CHECKPOINT
		goBegin(host(me));
		cargo(me)._atEnd = true;
	}

	//////////////////////////////////////////////////////////////////////////////
	// goBegin
	//////////////////////////////////////////////////////////////////////////////

	template <typename THost>
	inline void
	goBegin(ModifiedIterator<THost, ModReverse> & me)
	{
	SEQAN_CHECKPOINT
		goEnd(host(me));
		if (atBegin(host(me)))
			cargo(me)._atEnd = true;
		else
		{
			cargo(me)._atEnd = false;
			goPrevious(host(me));
		}
	}

	//////////////////////////////////////////////////////////////////////////////
	// operator +
	//////////////////////////////////////////////////////////////////////////////

	template <typename THost, typename TDelta>
	inline ModifiedIterator<THost, ModReverse> &
	operator += (ModifiedIterator<THost, ModReverse> & me, TDelta delta_) 
	{
		typedef ModifiedIterator<THost, ModReverse> TIterator;
		typedef typename Position<TIterator>::Type TPosition;
		TPosition delta = delta_;

		if (delta == 0)
		{
			return me;
		}
		if (delta > 0)
		{
			if (position(host(me)) < delta) {
				cargo(me)._atEnd = true;
				--delta;
			}
			host(me) -= delta;
		}
		else
		{
			if (cargo(me)._atEnd) {
				cargo(me)._atEnd = false;
				++delta;
			}
			host(me) -= delta;
		} 
		return me;
	}

	//////////////////////////////////////////////////////////////////////////////
	// operator -
	//////////////////////////////////////////////////////////////////////////////

	template <typename THost, typename TDelta>
	inline ModifiedIterator<THost, ModReverse> &
	operator -= (ModifiedIterator<THost, ModReverse> & me, TDelta delta) {
		if (delta > 0) {
			if (cargo(me)._atEnd) {
				cargo(me)._atEnd = false;
				--delta;
			}
			host(me) += delta;
		} else {
			if (position(host(me)) < -delta) {
				cargo(me)._atEnd = true;
				++delta;
			}
			host(me) -= -delta;
		}
		return me;
	}

	template <typename THost>
	inline typename Difference< ModifiedIterator<THost, ModReverse> >::Type
	operator - (ModifiedIterator<THost, ModReverse> const & a, ModifiedIterator<THost, ModReverse> const & b) {
		typename Difference< ModifiedIterator<THost, ModReverse> >::Type diff = host(b) - host(a);
		if (cargo(a)._atEnd) ++diff;
		if (cargo(b)._atEnd) --diff;
		return diff;
	}

	//////////////////////////////////////////////////////////////////////////////
	// position
	//////////////////////////////////////////////////////////////////////////////

	template <typename THost>
	inline typename Position<ModifiedIterator<THost, ModReverse> const>::Type 
	position(ModifiedIterator<THost, ModReverse> const & me)
	{
	SEQAN_CHECKPOINT
		if (cargo(me)._atEnd)
			return length(container(host(me)));
		else
			return length(container(host(me))) - 1 - position(host(me));
	}

	template <typename THost, typename TContainer>
	inline typename Position<ModifiedIterator<THost, ModReverse> const>::Type 
	position(ModifiedIterator<THost, ModReverse> const & me, TContainer const &cont)
	{
	SEQAN_CHECKPOINT
		if (cargo(me)._atEnd)
			return length(cont);
		else
			return length(cont) - 1 - position(host(me), cont);
	}

	//////////////////////////////////////////////////////////////////////////////
	// setPosition
	//////////////////////////////////////////////////////////////////////////////

	template <typename THost, typename TPosition>
	inline void
	setPosition(ModifiedIterator<THost, ModReverse> const & me, TPosition pos)
	{
	SEQAN_CHECKPOINT
		setPosition(host(me), length(container(host(me))) - 1 - pos);
	}

	//////////////////////////////////////////////////////////////////////////////
	// operator ==
	//////////////////////////////////////////////////////////////////////////////

	template <typename THost>
	inline bool
	operator == (ModifiedIterator<THost, ModReverse> const & a, ModifiedIterator<THost, ModReverse> const & b) {
		return cargo(a)._atEnd == cargo(b)._atEnd && host(a) == host(b);
	}

	//////////////////////////////////////////////////////////////////////////////
	// operator <
	//////////////////////////////////////////////////////////////////////////////

	// redefinition candidate
	template <typename THost>
	inline bool
	operator < (ModifiedIterator<THost, ModReverse> const & a, ModifiedIterator<THost, ModReverse> const & b) {
		return (!cargo(a)._atEnd && cargo(b)._atEnd) ||
			   (!cargo(a)._atEnd && !cargo(b)._atEnd && host(a) > host(b));
	}

	//////////////////////////////////////////////////////////////////////////////
	// atBegin
	//////////////////////////////////////////////////////////////////////////////

	template <typename THost, typename TContainer>
	inline bool
	atBegin(ModifiedIterator<THost, ModReverse> const & me,
			TContainer const & container)
	{
	SEQAN_CHECKPOINT
		return position(me, container) == 0;
	}

	template <typename THost>
	inline bool
	atBegin(ModifiedIterator<THost, ModReverse> const & me)
	{
	SEQAN_CHECKPOINT
		return position(me) == 0;
	}

	//////////////////////////////////////////////////////////////////////////////
	// atEnd
	//////////////////////////////////////////////////////////////////////////////

	template <typename THost, typename TContainer>
	inline bool
	atEnd(ModifiedIterator<THost, ModReverse> const & me,
			TContainer const & /*container*/)
	{
	SEQAN_CHECKPOINT
		return cargo(me)._atEnd;
	}

	template <typename THost>
	inline bool
	atEnd(ModifiedIterator<THost, ModReverse> const & me)
	{
	SEQAN_CHECKPOINT
		return cargo(me)._atEnd;
	}


//////////////////////////////////////////////////////////////////////////////


	//////////////////////////////////////////////////////////////////////////////
	// reverse string
	//////////////////////////////////////////////////////////////////////////////


	template <typename THost>
	class ModifiedString<THost, ModReverse> {
	public:
		Holder<THost>							data_host;
		typename Cargo<ModifiedString>::Type	data_cargo;

		ModifiedString() {}

		ModifiedString(ModifiedString &_origin):
			data_host(_origin.data_host),
			data_cargo(_origin.data_cargo) {}

		ModifiedString(ModifiedString const &_origin):
			data_host(_origin.data_host),
			data_cargo(_origin.data_cargo) {}

		template <typename THostHost, typename THostSpec>
		ModifiedString(ModifiedString<THostHost, THostSpec> &_origin):
			data_host(_origin.data_host) {}

		ModifiedString(THost &_origin) {
			setHost(*this, _origin);
		}

		template <typename T>
		ModifiedString(T & _origin) {
			setValue(*this, _origin);
		}

		template <typename T>
		ModifiedString(T const & _origin) {
			setValue(*this, _origin);
		}

		template <typename T>
		inline ModifiedString const &
		operator = (T & _origin) {
			assign(*this, _origin);
			return *this;
		}

		template <typename TPos>
		inline typename Reference<ModifiedString>::Type 
		operator [] (TPos pos)
		{
		SEQAN_CHECKPOINT
			return value(*this, pos);
		}

		template <typename TPos>
		inline typename Reference<ModifiedString const>::Type 
		operator [] (TPos pos) const
		{
		SEQAN_CHECKPOINT
			return value(*this, pos);
		}
	};


	template <typename THost>
	struct Iterator< ModifiedString<THost, ModReverse>, Standard > {
		typedef ModifiedIterator<typename Iterator<THost, Rooted>::Type, ModReverse> Type;
	};

	template <typename THost>
	struct Iterator< ModifiedString<THost, ModReverse> const, Standard > {
		typedef ModifiedIterator<typename Iterator<THost const, Rooted>::Type, ModReverse> Type;
	};

	template <typename THost>
	struct DefaultIteratorSpec< ModifiedString<THost, ModReverse> >
	{
		typedef Rooted Type;
	};




	//////////////////////////////////////////////////////////////////////////////
	// value
	//////////////////////////////////////////////////////////////////////////////

	template <typename THost, typename TPos>
	inline typename Reference<ModifiedString<THost, ModReverse> >::Type 
	value(ModifiedString<THost, ModReverse> & me, TPos pos)
	{
	SEQAN_CHECKPOINT
		return value(host(me), (length(host(me)) - 1) - pos);
	}

	template <typename THost, typename TPos>
	inline typename Reference<ModifiedString<THost, ModReverse> const>::Type 
	value(ModifiedString<THost, ModReverse> const & me, TPos pos)
	{
	SEQAN_CHECKPOINT
		return value(host(me), (length(host(me)) - 1) - pos);
	}

	//////////////////////////////////////////////////////////////////////////////
	// begin
	//////////////////////////////////////////////////////////////////////////////

	template < typename THost, typename TTag >
	inline typename Iterator< ModifiedString<THost, ModReverse> const >::Type 
	begin(ModifiedString<THost, ModReverse> const & me) {
		typename Iterator< ModifiedString<THost, ModReverse> const >::Type temp_(end(host(me), Rooted()));
		_copyCargo(temp_, me);
		goNext(temp_);
		return temp_;
	}

	template < typename THost >
	inline typename Iterator< ModifiedString<THost, ModReverse> >::Type 
	begin(ModifiedString<THost, ModReverse> & me) {
		typename Iterator< ModifiedString<THost, ModReverse> >::Type temp_(end(host(me), Rooted()));
		_copyCargo(temp_, me);
		goNext(temp_);
		return temp_;
	}

	template < typename THost, typename TTagSpec >
	inline typename Iterator< ModifiedString<THost, ModReverse> const, Tag<TTagSpec> const >::Type 
	begin(ModifiedString<THost, ModReverse> const & me, Tag<TTagSpec> const) {
		typename Iterator< ModifiedString<THost, ModReverse> const, Tag<TTagSpec> const >::Type temp_(end(host(me), Rooted()));
		_copyCargo(temp_, me);
		goNext(temp_);
		return temp_;
	}

	template < typename THost, typename TTagSpec >
	inline typename Iterator< ModifiedString<THost, ModReverse>, Tag<TTagSpec> const >::Type 
	begin(ModifiedString<THost, ModReverse> & me, Tag<TTagSpec> const) {
		typename Iterator< ModifiedString<THost, ModReverse>, Tag<TTagSpec> const >::Type temp_(end(host(me), Rooted()));
		_copyCargo(temp_, me);
		goNext(temp_);
		return temp_;
	}

	//////////////////////////////////////////////////////////////////////////////
	// end
	//////////////////////////////////////////////////////////////////////////////

	template < typename THost >
	inline typename Iterator< ModifiedString<THost, ModReverse> const >::Type 
	end(ModifiedString<THost, ModReverse> const & me) {
		typename Iterator< ModifiedString<THost, ModReverse> const >::Type temp_(begin(host(me), Rooted()));
		_copyCargo(temp_, me);
		goNext(temp_);
		return temp_;
	}

	template < typename THost >
	inline typename Iterator< ModifiedString<THost, ModReverse> >::Type 
	end(ModifiedString<THost, ModReverse> & me) {
		typename Iterator< ModifiedString<THost, ModReverse> >::Type temp_(begin(host(me), Rooted()));
		_copyCargo(temp_, me);
		goNext(temp_);
		return temp_;
	}

	template < typename THost, typename TTagSpec >
	inline typename Iterator< ModifiedString<THost, ModReverse> const, Tag<TTagSpec> const >::Type 
	end(ModifiedString<THost, ModReverse> const & me, Tag<TTagSpec> const) {
		typename Iterator< ModifiedString<THost, ModReverse> const, Tag<TTagSpec> const >::Type temp_(begin(host(me), Rooted()));
		_copyCargo(temp_, me);
		goNext(temp_);
		return temp_;
	}

	template < typename THost, typename TTagSpec >
	inline typename Iterator< ModifiedString<THost, ModReverse>, Tag<TTagSpec> const >::Type 
	end(ModifiedString<THost, ModReverse> & me, Tag<TTagSpec> const) {
		typename Iterator< ModifiedString<THost, ModReverse>, Tag<TTagSpec> const >::Type temp_(begin(host(me), Rooted()));
		_copyCargo(temp_, me);
		goNext(temp_);
		return temp_;
	}


//////////////////////////////////////////////////////////////////////////////


	//////////////////////////////////////////////////////////////////////////////
	// reverse
	//////////////////////////////////////////////////////////////////////////////

/**
.Function.reverse
..summary:Reverse an object/container in-place.
..cat:Modifier
..signature:reverse(object)
..param.object:The object/container whose elements to reverse.
...type:Concept.ContainerConcept
...type:Adaption.std::list
..include:seqan/modifier.h
*/
	template < typename TSequence >
	inline void
	reverse(TSequence & sequence) 
	{
		typedef typename Value<TSequence>::Type					TValue;

#if defined (_OPENMP) && defined (SEQAN_PARALLEL)
		// OpenMP does not support for loop with iterators. Therefore use index variables.
		typedef typename Position<TSequence>::Type				TPos;
		typedef typename MakeSigned_<TPos>::Type				TSignedPos;

		TSignedPos pMid = length(sequence) / 2;

		#pragma omp parallel for if(length(sequence) > 1000000)
		for(TSignedPos p1 = 0; p1 < pMid; ++p1) {
			TPos p2 = length(sequence) - 1 - p1;
			TValue tmp = sequence[p1];
			sequence[p1] = sequence[p2];
			sequence[p2] = tmp;
		}
#else
		typedef typename Iterator<TSequence, Standard>::Type	TIter;
		TIter it1 = begin(sequence, Standard());
		TIter it2 = it1 + (length(sequence) - 1);
		TIter itMid = it1 + length(sequence) / 2;

		for(; it1 != itMid; ++it1, --it2) {
			TValue tmp = *it1;
			*it1 = *it2;
			*it2 = tmp;
		}
#endif
	}

	template < typename TSequence >
	inline void
	reverse(TSequence const & sequence) 
	{
	    reverse(const_cast<TSequence &>(sequence));
    }

	template < typename TSequence, typename TSpec >
	inline void
	reverse(StringSet<TSequence, TSpec> & stringSet) 
	{
		unsigned seqCount = length(stringSet);
		for(unsigned seqNo = 0; seqNo < seqCount; ++seqNo)
			reverse(stringSet[seqNo]);
	}

	template < typename TSequence, typename TSpec >
	inline void
	reverse(StringSet<TSequence, TSpec> const & stringSet) 
	{
		unsigned seqCount = length(stringSet);
		for(unsigned seqNo = 0; seqNo < seqCount; ++seqNo)
			reverse(stringSet[seqNo]);
	}

    template <typename TValue>
    inline void
    reverse(std::list<TValue> & list)
    {
        SEQAN_CHECKPOINT;
        list.reverse();
    }

//////////////////////////////////////////////////////////////////////////////
// shortcut

template <typename THost>
inline ModifiedString<THost, ModReverse>
reverseString(THost const & host)
{
	return ModifiedString<THost, ModReverse>(host);
}

//////////////////////////////////////////////////////////////////////////////
}

#endif
