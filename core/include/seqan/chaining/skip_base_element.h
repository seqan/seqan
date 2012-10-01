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


/*	2006 Hendrik Woehrle
*
*	Deferred Skip List Datastructure
*
*	Elements in the base layer of the Skip List
*
*	Specializations:
*
*	TModus:
*		* SkipListDynamic: contains predecessor/successor pointers to have the properties of an double linked list in the base layer
*		* SkipListStatic: comes without these pointers to reduce size
*
*	TSpec:
*		* Complete: Not using Deferred Data Structuring
*		* Deferred: For use with Deferred Data Structuring
*
*/

//SEQAN_NO_DDDOC: do not generate documentation for this file

#ifndef SEQAN_HEADER_SKIP_BASE_ELEMENT_H
#define SEQAN_HEADER_SKIP_BASE_ELEMENT_H

namespace seqan
{


//////////////////////////////////////////////////////////////////////////////
//
//		SkipBaseElement
//
//////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////
//
//	SkipListDynamic <-> SkipListStatic and Complete <-> Deferred member wrapper structs
//
//	depending on this specializations, the elements contain different members
//
/////////////////////////////////////////////////////////////////

// special members of static/dynamic elements

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
	struct
	DynamicStruct_
	{};

	template< typename TObject, typename TSpec, typename TStructuring >
	struct
	DynamicStruct_< TObject, SkipListDynamic, TSpec, TStructuring >
	{
			// predecessing element in the skip list
		SkipBaseElement< TObject, SkipListDynamic, TSpec, TStructuring > * _pred;
			// succeeding element
		SkipBaseElement< TObject, SkipListDynamic, TSpec, TStructuring > * _succ;

		DynamicStruct_()
			: _pred( NULL )
			, _succ( NULL )
		{}

		DynamicStruct_( SkipBaseElement< TObject, SkipListDynamic, TSpec, TStructuring > * goPrevious, 
						SkipBaseElement< TObject, SkipListDynamic, TSpec, TStructuring > * succ )
			: _pred( goPrevious )
			, _succ( succ )
		{}
	};

// special members of complete/deferred elements

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
	struct
	DeferredStruct_
	{
		union{
				// number of unsorted elements to the right of this
			typename Size< SkipBaseElement< TObject, TModus, TSpec, TStructuring > >::Type _count;
				// height of related tower
			typename Size< SkipBaseElement< TObject, TModus, TSpec, TStructuring > >::Type _height;
		};

		DeferredStruct_()
			: _height( 0 )
		{}

		DeferredStruct_( typename Size< SkipBaseElement< TObject, TModus, TSpec, TStructuring > >::Type h )
			: _height( h )
		{}
	};
	
	template< typename TObject, typename TSpec >
	struct
	DeferredStruct_< TObject, SkipListDynamic, TSpec, Deferred >
	{
			// next sorted element on the left side
		SkipBaseElement< TObject, SkipListDynamic, TSpec, Deferred > * _left;
			// next sorted element on the right side
		SkipBaseElement< TObject, SkipListDynamic, TSpec, Deferred > * _right;
			// number of unsorted elements to the right of this
		typename Size< SkipBaseElement< TObject, SkipListDynamic, TSpec, Deferred > >::Type _count;
			// height of related tower
		typename Size< SkipBaseElement< TObject, SkipListDynamic, TSpec, Deferred > >::Type _height;

		DeferredStruct_()
			: _left( NULL )
			, _right( NULL )
			, _count( 0 )
			, _height( 0 )
		{}

		DeferredStruct_( SkipBaseElement< TObject, SkipListDynamic, TSpec, Deferred > * left, 
							SkipBaseElement< TObject, SkipListDynamic, TSpec, Deferred > * right )
			: _left( left )
			, _right( right )
			, _count( 0 )
			, _height( 0 )
		{}

		template< typename TSize >
		DeferredStruct_( SkipBaseElement< TObject, SkipListDynamic, TSpec, Deferred > * left, 
							SkipBaseElement< TObject, SkipListDynamic, TSpec, Deferred > * right,
							TSize _count,
							TSize _height )
			: _left( left )
			, _right( right )
			, _count( _count )
			, _height( _height )
		{}
	};


	template< typename TObject, typename TSpec >
	struct
	DeferredStruct_< TObject, SkipListStatic, TSpec, Deferred >
	{
		//	// next sorted element on the left side
		//SkipBaseElement< TObject, SkipListStatic, TSpec, Deferred > * _left;
			// number of unsorted elements to the right of this
		typename Size< SkipBaseElement< TObject, SkipListStatic, TSpec, Deferred > >::Type _count;
			// height of related tower
		typename Size< SkipBaseElement< TObject, SkipListStatic, TSpec, Deferred > >::Type _height;

		DeferredStruct_()
			: _count( 0 )
			, _height( 0 )
		{}

	};

		// overloading for iterators

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
	inline typename Reference< SkipBaseElement< TObject, TModus, TSpec, TStructuring > >::Type
	value( SkipBaseElement< TObject, TModus, TSpec, TStructuring > & me)
	{
	SEQAN_CHECKPOINT
		return *me._obj;
	}

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
	inline typename Reference< SkipBaseElement< TObject, TModus, TSpec, TStructuring > >::Type
	value( SkipBaseElement< TObject, TModus, TSpec, TStructuring > const & me)
	{
	SEQAN_CHECKPOINT
		return *me._obj;
	}

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
	inline typename Reference< SkipBaseElement< TObject, TModus, TSpec, TStructuring > >::Type
	value( SkipBaseElement< TObject, TModus, TSpec, TStructuring > * me)
	{
	SEQAN_CHECKPOINT
		return *me->_obj;
	}
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
	inline typename Reference< SkipBaseElement< TObject, TModus, TSpec, TStructuring > >::Type
	value( SkipBaseElement< TObject, TModus, TSpec, TStructuring > const * me)
	{
	SEQAN_CHECKPOINT
		return *me->_obj;
	}

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
	inline typename GetValue< SkipBaseElement< TObject, TModus, TSpec, TStructuring > >::Type
	getValue( SkipBaseElement< TObject, TModus, TSpec, TStructuring > & me )
	{
	SEQAN_CHECKPOINT
		return value(me);
	}

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
	inline typename GetValue< SkipBaseElement< TObject, TModus, TSpec, TStructuring > const >::Type
	getValue( SkipBaseElement< TObject, TModus, TSpec, TStructuring > const & me )
	{
	SEQAN_CHECKPOINT
		return value(me);
	}

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
	inline typename GetValue< SkipBaseElement< TObject, TModus, TSpec, TStructuring > * >::Type
	getValue( SkipBaseElement< TObject, TModus, TSpec, TStructuring > * me )
	{
	SEQAN_CHECKPOINT
		return *me;
	}


	template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
	inline void
	valueDestruct( SkipBaseElement< TObject, TModus, TSpec, TStructuring > * it )
	{
	SEQAN_CHECKPOINT
		it->~SkipBaseElement< TObject, TModus, TSpec, TStructuring >();
	}

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
	inline void
	valueConstruct( SkipBaseElement< TObject, TModus, TSpec, TStructuring > * it )
	{
	SEQAN_CHECKPOINT
		new( it ) SkipBaseElement< TObject, TModus, TSpec, TStructuring >;
	}

	template < typename TObject, typename TModus, typename TSpec, typename TStructuring >
	inline void
	valueConstruct( SkipBaseElement< TObject, TModus, TSpec, TStructuring > * it,
				   TObject * obj )
	{
	SEQAN_CHECKPOINT
		new( it ) SkipBaseElement< TObject, TModus, TSpec, TStructuring >( obj );
	}

	template < typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TKey >
	inline void
	valueConstruct( SkipBaseElement< TObject, TModus, TSpec, TStructuring > * it,
				   TObject * obj,
				   TKey key )
	{
	SEQAN_CHECKPOINT
		new( it ) SkipBaseElement< TObject, TModus, TSpec, TStructuring >( obj, key );
	}


/*
.Internal.getObject:
..summary:Get the saved object which stores key-value information.
..cat:SkipList
..signature:getObject(element)
..param.element:The desired element.
...type:Class.SkipBaseElement
*/


	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
	TObject * 
	getObject( SkipBaseElement< TObject, TModus, TSpec, TStructuring > * me )
	{
		SEQAN_CHECKPOINT
		return me->_obj;
	}


/*
.Internal._setObject:
..summary:Set the saved object of a SkipBaseElement which stores key-value information.
..remarks: For Internal use only. To insert an object into a dynamic skip list, use insert.
..cat:SkipList
..signature:getObject(element, object)
..param.element:The desired element.
...type:Class.SkipBaseElement
..param.element:The desired element.
...type:Class.SkipBaseElement
*/

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
	void 
	_setObject( SkipBaseElement< TObject, TModus, TSpec, TStructuring > & me,
				TObject * obj )
	{
		SEQAN_CHECKPOINT
		me._obj = obj;
	}


	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
	void 
	_setObject( SkipBaseElement< TObject, TModus, TSpec, TStructuring > & me,
				SkipBaseElement< TObject, TModus, TSpec, TStructuring > * data )
	{
		SEQAN_CHECKPOINT
		me._obj = getObject( data );
	}


/*
.Internal._setUp:
..summary:Set the pointer to the lowest SkipElement of the related tower.
..cat:SkipList
..signature:_setUp(element, up)
..param.element:The desired object.
...type:Class.SkipBaseElement
..param.up:The SkipElement of the tower.
...type:Class.SkipElement
*/

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline 
	void 
	_setUp( SkipBaseElement< TObject, TModus, TSpec, TStructuring > & me, 
			SkipElement< TObject, TModus, TSpec, TStructuring > & up )
	{
		SEQAN_CHECKPOINT
		me._up = &up;
	}


/*
.Internal._getUp:
..summary:_get a pointer to the lowest SkipElement of the related tower.
..cat:SkipList
..signature:_getUp(element)
..param.element:The desired object.
...type:Class.SkipBaseElement
..returns:A @Class.SkipElement.$SkipElement*$@ pointing to the lowest SkipElement in the related tower. NULL if no tower is related to $element$.
*/

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline 
	SkipElement< TObject, TModus, TSpec, TStructuring > & 
	_getUp( SkipBaseElement< TObject, TModus, TSpec, TStructuring > & me )
	{
		SEQAN_CHECKPOINT
		return *me._up;
	}


/*
.Internal._getSucc:
..summary:Get a pointer to the succeeding SkipBaseElement.
..cat:SkipList
..signature:_getSucc(element)
..param.element:The desired object.
...type:Class.SkipBaseElement
..returns:A @Class.SkipBaseElement.$SkipBaseElement*$@, pointing to $element1$ succeeding $element$. NULL if no such $element1$ exists.
..remarks:If the containing SkipList is Deferred and not completely sorted, it is unlikely that the ordering of the keys is already established. 
..So $key( _getSucc( element ) ) >= key( element )$ does not hold for sure.
*/

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
	SkipBaseElement< TObject, TModus, TSpec, TStructuring > *
	_getSucc( SkipBaseElement< TObject, TModus, TSpec, TStructuring > & me )
	{
		SEQAN_CHECKPOINT
		return ( &me + 1 );
	}

	template< typename TObject, typename TSpec, typename TStructuring > inline
	SkipBaseElement< TObject, SkipListDynamic, TSpec, TStructuring > *
	_getSucc( SkipBaseElement< TObject, SkipListDynamic, TSpec, TStructuring > & me )
	{
		SEQAN_CHECKPOINT
		return me._dynStruct._succ;
	}

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
	void
	goNext( SkipBaseElement< TObject, TModus, TSpec, TStructuring > *& me )
	{
		SEQAN_CHECKPOINT
		++me;
	}

	template< typename TObject, typename TSpec, typename TStructuring > inline
	void
	goNext( SkipBaseElement< TObject, SkipListDynamic, TSpec, TStructuring > *& me )
	{
		SEQAN_CHECKPOINT
		me = _getSucc( *me );
	}


/*
.Internal._setSucc:
..summary:Set the pointer to the succeeding SkipBaseElement.
..cat:SkipList
..signature:_setSucc(element, succ)
..param.element:The element.
...type:Class.SkipBaseElement
..param.succ:The succeeding SkipBaseElement.
...type:Class.SkipBaseElement
*/

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
	void 
	_setSucc( SkipBaseElement< TObject, TModus, TSpec, TStructuring > & /*me*/, 
				SkipBaseElement< TObject, TModus, TSpec, TStructuring > * /*succ*/ )
	{
		SEQAN_ASSERT_FAIL("No dynamic related members in default mode");
	}

	template< typename TObject, typename TSpec, typename TStructuring > inline
	void 
	_setSucc( SkipBaseElement< TObject, SkipListDynamic, TSpec, TStructuring > & me, 
				SkipBaseElement< TObject, SkipListDynamic, TSpec, TStructuring > * succ )
	{
		SEQAN_CHECKPOINT
		me._dynStruct._succ = succ;
	}

	
/*
.Function._getPred:
..summary:Get a pointer to the preceding SkipBaseElement.
..cat:SkipList
..signature:_getPred(element)
..param.element:The desired object.
...type:Class.SkipBaseElement
..returns:A @Class.SkipBaseElement.$SkipBaseElement*$@, pointing to $element1$ preceding $element$. NULL if no such $element1$ exists.
..remarks:If the containing SkipList is Deferred and not completely sorted, it is unlikely that the ordering of the keys is already established. 
..So $key( _getPred( element ) ) >= key( element )$ does not hold for sure.
..include:seqan/chaining.h
*/
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
	SkipBaseElement< TObject, TModus, TSpec, TStructuring > * 
	_getPred( SkipBaseElement< TObject, TModus, TSpec, TStructuring > & me );

	template< typename TObject, typename TSpec, typename TStructuring > inline
	SkipBaseElement< TObject, SkipListStatic, TSpec, TStructuring > * 
	_getPred( SkipBaseElement< TObject, SkipListStatic, TSpec, TStructuring > & me )
	{
		SEQAN_CHECKPOINT
		return ( &me - 1 );
	}

	template< typename TObject, typename TSpec, typename TStructuring > inline
	SkipBaseElement< TObject, SkipListDynamic, TSpec, TStructuring > * 
	_getPred( SkipBaseElement< TObject, SkipListDynamic, TSpec, TStructuring > & me )
	{
		SEQAN_CHECKPOINT
		return me._dynStruct._pred;
	}

	template< typename TObject, typename TSpec, typename TStructuring > inline
	void
	goPrevious( SkipBaseElement< TObject, SkipListStatic, TSpec, TStructuring > *& me )
	{
		SEQAN_CHECKPOINT
		--me;
	}

	template< typename TObject, typename TSpec, typename TStructuring > inline
	void
	goPrevious( SkipBaseElement< TObject, SkipListDynamic, TSpec, TStructuring > *& me )
	{
		SEQAN_CHECKPOINT
		me = _getPred( *me );
	}


/*
.Internal._setPred:
..summary:Set the pointer to the preceding SkipBaseElement.
..cat:SkipList
..signature:_setPred(element, goPrevious)
..param.element:The element.
...type:Class.SkipBaseElement
..param.goPrevious:The succeeding SkipBaseElement.
...type:Class.SkipBaseElement
*/
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
	void 
	_setPred(	SkipBaseElement< TObject, TModus, TSpec, TStructuring > & /*me*/, 
				SkipBaseElement< TObject, TModus, TSpec, TStructuring > * /*goPrevious*/ )
	{
		SEQAN_ASSERT_FAIL("No dynamic related members in default mode");
	}

	template< typename TObject, typename TSpec, typename TStructuring > inline
	void 
	_setPred(	SkipBaseElement< TObject, SkipListDynamic, TSpec, TStructuring > & me, 
				SkipBaseElement< TObject, SkipListDynamic, TSpec, TStructuring > * goPrevious )
	{
		SEQAN_CHECKPOINT
		me._dynStruct._pred = goPrevious;
	}

/*
.Internal._getLeft:
..summary:Get a pointer to the next sorted SkipBaseElement on the left side.
..cat:SkipList
..signature:_getLeft(element)
..param.element:The desired object.
...type:Class.SkipBaseElement
..returns:A @Class.SkipBaseElement.$SkipBaseElement*$@, pointing to the next sorted element on the left side.
*/

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline 
	SkipBaseElement< TObject, TModus, TSpec, TStructuring > * 
	_getLeft( SkipBaseElement< TObject, TModus, TSpec, TStructuring > & me )
	{
		SEQAN_CHECKPOINT
		return _getPred( me );
	}

	template< typename TObject, typename TModus, typename TSpec > inline 
	SkipBaseElement< TObject, TModus, TSpec, Deferred > * 
	_getLeft( SkipBaseElement< TObject, TModus, TSpec, Deferred > & me )
	{
		SEQAN_CHECKPOINT
		return me._defStruct._left;
	}


/*
.Internal._setLeft:
..summary:Get a pointer to the next sorted SkipBaseElement on the left side.
..cat:SkipList
..signature:_getLeft(element)
..param.element:The desired object.
...type:Class.SkipBaseElement
..returns:A @Class.SkipBaseElement.$SkipBaseElement*$@, pointing to the next sorted element on the left side.
*/

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline 
	void 
	_setLeft( SkipBaseElement< TObject, TModus, TSpec, TStructuring > & /*me*/,
							SkipBaseElement< TObject, TModus, TSpec, TStructuring > * /*left*/ )
	{
		SEQAN_ASSERT_FAIL("No deferred related members in default mode");
	}

	template< typename TObject, typename TModus, typename TSpec > inline 
	void 
	_setLeft( SkipBaseElement< TObject, TModus, TSpec, Deferred > & me,
							SkipBaseElement< TObject, TModus, TSpec, Deferred > * left )
	{
		SEQAN_CHECKPOINT
		me._defStruct._left = left;
	}


/*
.Internal._getRight:
..summary:Get a pointer to the next sorted SkipBaseElement on the right side.
..cat:SkipList
..signature:_getRight(element)
..param.element:The desired object.
...type:Class.SkipBaseElement
..returns:A @Class.SkipBaseElement.$SkipBaseElement*$@, pointing to the next sorted element on the right side.
*/


	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline 
	SkipBaseElement< TObject, TModus, TSpec, TStructuring > * 
	_getRight( SkipBaseElement< TObject, TModus, TSpec, TStructuring > & me )
	{
		SEQAN_CHECKPOINT
		return _getSucc( me );
	}


	template< typename TObject, typename TSpec > inline 
	SkipBaseElement< TObject, SkipListDynamic, TSpec, Deferred > * 
	_getRight( SkipBaseElement< TObject, SkipListDynamic, TSpec, Deferred > & me )
	{
		SEQAN_CHECKPOINT
		return me._defStruct._right;
	}

	template< typename TObject, typename TSpec > inline 
	SkipBaseElement< TObject, SkipListStatic, TSpec, Deferred > * 
	_getRight( SkipBaseElement< TObject, SkipListStatic, TSpec, Deferred > & me )
	{
		SEQAN_CHECKPOINT
		return &me +( _getCount( me ) + 1 );
	}


	template< typename TObject, typename TSpec, typename TStructuring > inline
	void
	right( SkipBaseElement< TObject, SkipListStatic, TSpec, TStructuring > *& me )
	{
		SEQAN_CHECKPOINT
		me += ( _getCount( *me ) + 1 );
	}

	template< typename TObject, typename TSpec, typename TStructuring > inline
	void
	right( SkipBaseElement< TObject, SkipListDynamic, TSpec, TStructuring > *& me )
	{
		SEQAN_CHECKPOINT
		me = _getRight( *me );
	}


/*
.Internal._setRight:
..summary:Get a pointer to the next sorted SkipBaseElement on the right side.
..cat:SkipList
..signature:_getRight(element)
..param.element:The desired object.
...type:Class.SkipBaseElement
..returns:A @Class.SkipBaseElement.$SkipBaseElement*$@, pointing to the next sorted element on the right side.
*/

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline 
	void 
	_setRight( SkipBaseElement< TObject, TModus, TSpec, TStructuring > & /*me*/,
				SkipBaseElement< TObject, TModus, TSpec, TStructuring > * /*right*/ )
	{
	}

	template< typename TObject, typename TSpec > inline 
	void 
	_setRight( SkipBaseElement< TObject, SkipListStatic, TSpec, Deferred > & /*me*/,
				SkipBaseElement< TObject, SkipListStatic, TSpec, Deferred > * /*right*/ )
	{
		// do nothing
	}

	template< typename TObject, typename TSpec > inline 
	void 
	_setRight( SkipBaseElement< TObject, SkipListDynamic, TSpec, Deferred > & me,
				SkipBaseElement< TObject, SkipListDynamic, TSpec, Deferred > * right )
	{
		SEQAN_CHECKPOINT
		me._defStruct._right = right;
	}


	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline 
	SkipBaseElement< TObject, TModus, TSpec, TStructuring > * 
	_getNext( SkipBaseElement< TObject, TModus, TSpec, TStructuring > & me )
	{
		SEQAN_CHECKPOINT
		return me._next;
	}


	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline 
	void 
	_setNext(	SkipBaseElement< TObject, TModus, TSpec, TStructuring > & me, 
				SkipBaseElement< TObject, TModus, TSpec, TStructuring > * next )
	{
		SEQAN_CHECKPOINT
		me._next = next;
	}


/*
.Internal._getHeight:
..summary:Get the height of the related tower.
..cat:SkipList
..signature:_getHeight(element)
..param.element:The element.
...type:Class.SkipBaseElement
..returns:Height of the associated tower.
...type:
*/

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
	typename Size< SkipList< TObject, TModus, TSpec, TStructuring > >::Type
	_getHeight( SkipBaseElement< TObject, TModus, TSpec, TStructuring > & me )
	{
		SEQAN_CHECKPOINT
		return me._defStruct._height;	
	}

/*
.Internal._setHeight:
..summary:Set the height of the related tower.
..cat:SkipList
..signature:_setHeight(element, height)
..param.element:The desired object.
...type:Class.SkipBaseElement
..param.height:The height.
...type:@Metafunction.Height.$Size< element-type >::Type$@
*/

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TSize > inline
	void 
	_setHeight(	SkipBaseElement< TObject, TModus, TSpec, TStructuring > & me,
				TSize height )
	{
		SEQAN_CHECKPOINT
		me._defStruct._height = height;
	}

// key

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TKey > inline 
	void 
	setKey( SkipBaseElement< TObject, TModus, TSpec, TStructuring > & me,
			TKey theKey)
	{
		SEQAN_CHECKPOINT
		me._key= theKey;
	}

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TKey > inline 
	void 
	setKey( SkipBaseElement< TObject, TModus, TSpec, TStructuring > * me,
			TKey theKey)
	{
		SEQAN_CHECKPOINT
		me->_key= theKey;
	}

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline 
	typename Key< TObject >::Type 
	key( SkipBaseElement< TObject, TModus, TSpec, TStructuring > & me )
	{
		SEQAN_CHECKPOINT
		return me._key;
	}

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TParam > inline 
	typename Key< TObject >::Type 
	key( SkipBaseElement< TObject, TModus, TSpec, TStructuring > & me,
				TParam & /*p*/ )
	{
		SEQAN_CHECKPOINT
		return me._key;
	}


/*
.Internal._getCount:
..summary:_get the count-value of the SkipBaseElement.
..cat:SkipList
..signature:_getCount(element)
..param.element:The desired object.
...type:Class.SkipBaseElement
..returns:The count-value of element. Type is $TKey$.
..remarks:The count value is defined as the number of unsorted elements between $element$ and $_getRight(element)$
*/
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline 
	typename Size< SkipList< TObject, TModus, TSpec, TStructuring > >::Type
	_getCount( SkipBaseElement< TObject, TModus, TSpec, TStructuring > & me )
	{
		SEQAN_CHECKPOINT
		return me._defStruct._count;
	}
	

/*
.Internal._setCount:
..summary:Set the count value.
..cat:SkipList
..signature:_setCount(element, count)
..param.element:The desired object.
...type:Class.SkipBaseElement
..param.count:The new count-value of element.
*/

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TSize > inline
	void 
	_setCount(	SkipBaseElement< TObject, TModus, TSpec, TStructuring > & me,
				TSize count )
	{
		SEQAN_CHECKPOINT
		me._defStruct._count = count;
	}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// struct SkipBaseElement
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////


/*
.Internal.SkipBaseElement:
..summary:Type of elements in the base layer of a SkipList
..cat:SkipList
..signature:SkipBaseElement< TObject, TModus, TSpec, TStructuring >
.signature:SkipElement< TObject, TModus, TSpec, TStructuring >
..param.TObject:Type of the stored object.
..param.TModus:The TModus parameter of the SkipList.
..param.TSpec:The TSpec parameter of the SkipList.
..param.TStructuring:The TStructuring parameter of the SkipList.
*/

template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct SkipBaseElement
{
	
	DynamicStruct_< TObject, TModus, TSpec, TStructuring > _dynStruct;

	DeferredStruct_< TObject, TModus, TSpec, TStructuring > _defStruct;
	
		// pointer to tower in the upper layers
	union{
		SkipElement< TObject, TModus, TSpec, TStructuring > * _up;
		SkipBaseElement< TObject, TModus, TSpec, TStructuring > * _next;	// pointer needed for memory pool allocator, points to next free memory block 
	};
		// the actual key
	typename Key< TObject >::Type _key;
		// saved key/value pair
	TObject * _obj;
		

	typename Cargo< SkipBaseElement< TObject, TModus, TSpec, TStructuring > >::Type _cargo;

	
	friend inline
	typename Cargo< SkipBaseElement< TObject, TModus, TSpec, TStructuring > >::Type * 
	cargo( SkipBaseElement< TObject, TModus, TSpec, TStructuring > & me )
	{
		SEQAN_CHECKPOINT
		return &me._cargo;
	}

	template< typename TCargo > friend inline
	void
	_setCargo( SkipBaseElement< TObject, TModus, TSpec, TStructuring > & me,
					TCargo & cargo )
	{
		SEQAN_CHECKPOINT
		me._cargo = cargo;
	}
	
	
public:

			// standard constructor
	SkipBaseElement(void)
		: _up(NULL)
		, _obj(NULL)
	{
		_initCargo( this, _cargo );
	}

	template< typename TKey >
	SkipBaseElement( TObject * obj, TKey key )
		: _up(NULL)
		, _key( key )
		, _obj(obj)
	{
		_initCargo( this, _cargo );
	}
	
	~SkipBaseElement(void)
	{
	}
  
};	// struct SkipBaseElement

	
} // namespace seqan
#endif //SKIP_BASE_ELEMENT_H_


