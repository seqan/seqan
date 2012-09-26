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


#ifndef SEQAN_HEADER_SKIPLIST_TYPE_H
#define SEQAN_HEADER_SKIPLIST_TYPE_H

namespace SEQAN_NAMESPACE_MAIN
{

//SEQAN_NO_DDDOC: do not generate documentation for this file


/**
.Metafunction.Key:
..summary:Type of the theKeyattribute of an object. 
..signature:Key<T>::Type
..param.T:Type for which the key type is determined.
..returns.param.Type:Key type of $T$.
..remarks.text:The theKeytype of an object is used for sorting and searching.
..include:seqan/chaining.h
*/

/* moved to basic_h, see #6
	template< typename T >
	struct Key
	{
		typedef T Type;
	};
*/
		// default for Pair
	template < typename T1_, typename T2_, typename Pack >
    struct Pair;

	// auxiliary functions for Objects that are plugged into the 
		// SkipList
		// Specialized for std::pair
	template< typename TKey, typename TVal > inline
	TKey & key( Pair<TKey, TVal> & p )
	{
		return p.i1;
	}

	template< typename TKey, typename TVal > inline
	TVal & getValue( Pair<TKey, TVal> & p )
	{
		return p.i2;
	}

	template< typename TKey2, typename TVal >
	void setKey( Pair<TKey2, TVal> & p, TKey2 theKey)
	{
		p.i1 = theKey;
	}

	template< typename TKey, typename TVal >
	struct Value< Pair< TKey, TVal> >
	{
		typedef TVal Type;
	};

	template< typename TKey, typename TVal >
	struct Key< Pair< TKey, TVal > >
	{
		typedef TKey Type;
	};

		// specialization for std::pair
	template< typename TKey, typename TVal > inline
	TKey key( std::pair< TKey, TVal > & p )
	{
		return p.first;
	}

	template< typename TKey, typename TVal >
	void setKey( std::pair<TKey, TVal> & p, TKey theKey )
	{
		p.first = theKey;
	}

	template< typename TKey, typename TVal >
	struct Value< std::pair< TKey, TVal> >
	{
		typedef TVal Type;
	};

	template< typename TKey, typename TVal >
	struct Key< std::pair< TKey, TVal > >
	{
		typedef TKey Type;
	};

	
	/**
	.Metafunction.Cargo:
	..summary:An additional cargo member of an object. 
	..signature:Cargo<T>::Type
	..param.T:Type for which the cargo type is determined.
	..returns.param.Type:Cargo type of $T$.
	..remarks.text:The cargo type is used for additional cargo information.
	*/

	struct EmptyCargo_
	{};
/*
	template< typename T >
	struct Cargo
	{
		typedef EmptyCargo_ Type;
	};
*/

//////////////////////////////////////////////////////////////////////////////
// Tags
//////////////////////////////////////////////////////////////////////////////

struct SkipListDynamic;

struct SkipListStatic;

struct Complete;

struct Deferred;

//////////////////////////////////////////////////////////////////////////////
// Forward declarations
//////////////////////////////////////////////////////////////////////////////


template< typename TObject, typename TModus = SkipListDynamic, typename TSpec = Default, typename TStructuring = Complete >
struct SkipElement;

template< typename TObject, typename TModus = SkipListDynamic, typename TSpec = Default, typename TStructuring  = Complete >
struct SkipBaseElement;

template< typename TObject, typename TModus = SkipListDynamic, typename TSpec = Default, typename TStructuring  = Complete >
struct SkipList;


//////////////////////////////////////////////////////////////////////////////
//
// METAFUNCTIONS
//
//////////////////////////////////////////////////////////////////////////////


template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct Size< SkipList< TObject, TModus, TSpec, TStructuring > >
{
	typedef size_t Type;
};

template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct Size< SkipList< TObject, TModus, TSpec, TStructuring > const >
{
	typedef typename Size< SkipList< TObject, TModus, TSpec, TStructuring > >::Type const Type;
};

template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct Size< SkipElement< TObject, TModus, TSpec, TStructuring > >
{
	typedef typename Size< SkipList< TObject, TModus, TSpec, TStructuring > >::Type Type;
};

template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct Size< SkipElement< TObject, TModus, TSpec, TStructuring > const >
{
	typedef typename Size< SkipList< TObject, TModus, TSpec, TStructuring > >::Type const Type;
};

template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct Size< SkipBaseElement< TObject, TModus, TSpec, TStructuring > >
{
	typedef typename Size< SkipList< TObject, TModus, TSpec, TStructuring > >::Type Type;
};

template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct Size< SkipBaseElement< TObject, TModus, TSpec, TStructuring > const >
{
	typedef typename Size< SkipList< TObject, TModus, TSpec, TStructuring > >::Type const Type;
};

///////////////////////////////////////////////////////////////////////////////////////////////
//
//		Position Type
//
///////////////////////////////////////////////////////////////////////////////////////////////


template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct Position< SkipList< TObject, TModus, TSpec, TStructuring > >
{
	typedef SkipBaseElement< TObject, TModus, TSpec, TStructuring > * Type;
};


///////////////////////////////////////////////////////////////////////////////////////////////
//
//		GetValue Type
//
///////////////////////////////////////////////////////////////////////////////////////////////


template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct GetValue< SkipList< TObject, TModus, TSpec, TStructuring > const >
{
	typedef typename GetValue< SkipList< TObject, TModus, TSpec, TStructuring > >::Type const Type;
};

template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct GetValue< SkipElement< TObject, TModus, TSpec, TStructuring > >
{
	typedef typename GetValue< SkipList< TObject, TModus, TSpec, TStructuring > >::Type Type;
};
template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct GetValue< SkipElement< TObject, TModus, TSpec, TStructuring > const >
{
	typedef typename GetValue< SkipList< TObject, TModus, TSpec, TStructuring > >::Type const Type;
};

template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct GetValue< SkipBaseElement< TObject, TModus, TSpec, TStructuring > >
{
	typedef typename GetValue< SkipList< TObject, TModus, TSpec, TStructuring > >::Type Type;
};

template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct GetValue< SkipBaseElement< TObject, TModus, TSpec, TStructuring > const >
{
	typedef typename GetValue< SkipList< TObject, TModus, TSpec, TStructuring > >::Type const Type;
};


template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct GetValue< SkipBaseElement< TObject, TModus, TSpec, TStructuring > * >
{
	typedef SkipBaseElement< TObject, TModus, TSpec, TStructuring > Type;
};


///////////////////////////////////////////////////////////////////////////////////////////////
//
//		Value Type
//
///////////////////////////////////////////////////////////////////////////////////////////////


template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct Value< SkipList< TObject, TModus, TSpec, TStructuring > >
{
	typedef TObject Type;
};

template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct Value< SkipList< TObject, TModus, TSpec, TStructuring > const >
{
	typedef typename Value< SkipList< TObject, TModus, TSpec, TStructuring > >::Type const Type;
};

template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct Value< SkipElement< TObject, TModus, TSpec, TStructuring > >
{
	typedef typename Value< SkipList< TObject, TModus, TSpec, TStructuring > >::Type Type;
};
template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct Value< SkipElement< TObject, TModus, TSpec, TStructuring > const >
{
	typedef typename Value< SkipList< TObject, TModus, TSpec, TStructuring > >::Type const Type;
};

template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct Value< SkipBaseElement< TObject, TModus, TSpec, TStructuring > >
{
	typedef typename Value< SkipList< TObject, TModus, TSpec, TStructuring > >::Type Type;
};
template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct Value< SkipBaseElement< TObject, TModus, TSpec, TStructuring > const >
{
	typedef typename Value< SkipList< TObject, TModus, TSpec, TStructuring > >::Type const Type;
};


///////////////////////////////////////////////////////////////////////////////////////////////
//
//		Key Type
//
///////////////////////////////////////////////////////////////////////////////////////////////


template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct Key< SkipList< TObject, TModus, TSpec, TStructuring > >
{
	typedef typename Key< TObject >::Type Type;
};

template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct Key< SkipList< TObject, TModus, TSpec, TStructuring > const >
{
	typedef typename Key< SkipList< TObject, TModus, TSpec, TStructuring > >::Type const Type;
};

template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct Key< SkipElement< TObject, TModus, TSpec, TStructuring > >
{
	typedef typename Key< SkipList< TObject, TModus, TSpec, TStructuring > >::Type Type;
};

template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct Key< SkipElement< TObject, TModus, TSpec, TStructuring > const >
{
	typedef typename Key< SkipList< TObject, TModus, TSpec, TStructuring > >::Type const Type;
};

template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct Key< SkipBaseElement< TObject, TModus, TSpec, TStructuring > >
{
	typedef typename Key< SkipList< TObject, TModus, TSpec, TStructuring > >::Type Type;
};

template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct Key< SkipBaseElement< TObject, TModus, TSpec, TStructuring > const >
{
	typedef typename Key< SkipList< TObject, TModus, TSpec, TStructuring > >::Type const Type;
};


///////////////////////////////////////////////////////////////////////////////////////////////
//
//		Cargo Type
//
///////////////////////////////////////////////////////////////////////////////////////////////


template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct Cargo< SkipList< TObject, TModus, TSpec, TStructuring > >
{
	typedef EmptyCargo_ Type;	// default: no cargo
};

template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct Cargo< SkipElement< TObject, TModus, TSpec, TStructuring > >
{
	typedef typename Cargo< SkipList< TObject, TModus, TSpec, TStructuring > >::Type Type;
};

template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct Cargo< SkipBaseElement< TObject, TModus, TSpec, TStructuring > >
{
	typedef typename Cargo< SkipList< TObject, TModus, TSpec, TStructuring > >::Type Type;
};


template< typename TTag, typename TCargo > inline
void
_initCargo( TTag * /*tag*/, 
		   TCargo & /*_cargo*/ )
{}


}

#endif // SEQAN_HEADER_SKIPLIST_TYPE_H
