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
*	Elements in the higher layers of the Skip List
*
*/

//SEQAN_NO_ DDDOC: do not generate documentation for this file

#ifndef SEQAN_HEADER_SKIP_ELEMENT_H
#define SEQAN_HEADER_SKIP_ELEMENT_H

namespace seqan
{

//////////////////////////////////////////////////////////////////////////////
//		SkipElement
//////////////////////////////////////////////////////////////////////////////

/*
.Internal.SkipElement:
..summary:Basic data structure for elements in higher levels of a SkipList
..cat:SkipList
..signature:SkipElement< TObject, TModus, TSpec, TStructuring >
..param.TObject:Type of the stored object.
..param.TModus:The TModus parameter of the SkipList.
..param.TSpec:The TSpec parameter of the SkipList.
..param.TStructuring:The TStructuring parameter of the SkipList.
*/

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
	struct SkipElement;


/*
.Internal.getObject:
..summary:Get a reference to the related object in the SkipList
..cat:SkipList
..signature:getObject( element )
..param.element:The object that has the desired value.
...type:Class.SkipElement
...type:Class.SkipBaseElement
..returns:A Reference to the related object. Type depends on the TObject parameter of the SkipList.
*/

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > __inline 
	TObject * 
	getObject( SkipElement< TObject, TModus, TSpec, TStructuring > * me )
	{
		SEQAN_CHECKPOINT
		return getObject( _getDown( me ));
	}

/*
.Internal._getRight:
..summary:Get the SkipElement to the right of object.
..cat:SkipList
..signature:_getRight(element)
..param.element:The object that references the desired SkipElement to the right.
...type:Class.SkipElement
...type:Class.SkipBaseElement
..returns:The pointer to the right element. Types are @Class.SkipBaseElement.$SkipBaseElement*$@ or @Class.SkipElement.$SkipElement*$@ respectively.
*/
	
		//default case
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > __inline 
	SkipElement< TObject, TModus, TSpec, TStructuring > * 
	_getRight(SkipElement< TObject, TModus, TSpec, TStructuring > & me)
	{
		SEQAN_CHECKPOINT
		return me._right;
	}

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > __inline
	void
	right( SkipElement< TObject, TModus, TSpec, TStructuring > *& me )
	{
		SEQAN_CHECKPOINT
		me = me->_right;
	}

/*
.Internal._getDown:
..summary:Get the underlying SkipBaseElement of a SkipElement.
..cat:SkipList
..signature:_getDown(element)
..param.element:The element.
...type:Class.SkipElement
..returns:A pointer to the down element. Type is @Class.SkipBaseElement.$SkipBaseElement*$@
*/

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > __inline 
	SkipBaseElement< TObject, TModus, TSpec, TStructuring > * 
	_getDown( SkipElement< TObject, TModus, TSpec, TStructuring > & me )
	{
		SEQAN_CHECKPOINT
		return me._down;
	}


/*
.Internal._setRight:
..summary:Set the $right$ member variable of a SkipElement to the SkipElement (or of a SkipBaseElement to the SkipBaseElement) that should be on the right.
..cat:SkipList
..signature:_setRight(element1, element2)
..param.element1:The desired object.
...type:Class.SkipElement
...type:Class.SkipBaseElement
..param.element2:Pointer to the element, that should become the element on the right.
...type:Class.SkipElement
...type:Class.SkipBaseElement
..remarks.text:A SkipElement can only have another SkipElement on the right, a SkipBaseElement only another SkipBaseElement.
*/

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > __inline 
	void 
	_setRight( SkipElement< TObject, TModus, TSpec, TStructuring > & me, 
			   SkipElement< TObject, TModus, TSpec, TStructuring > * right )
	{
		SEQAN_CHECKPOINT
		me._right = right;
	}

/*
.Internal._setDown:
..summary:Set the $down$ member variable to the related SkipBaseElement in the base layer.
..cat:SkipList
..signature:_setDown(element1, element2)
..param.element1:The element.
...type:Class.SkipElement
..param.element2:A pointer to the realted SkipBaseElement.
...type:Class.SkipBaseElement
*/

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > __inline 
	void 
	_setDown( SkipElement< TObject, TModus, TSpec, TStructuring > & me,
			  SkipBaseElement< TObject, TModus, TSpec, TStructuring > * down )
	{
		SEQAN_CHECKPOINT
		me._down = down;
	}

/*
.Internal._getNext:
..summary:Get the next-pointer of the element. This is only used in @Class.ClassPool@ to administrate memory blocks.
..cat:SkipList
..signature:_getNext(element)
..param.element:The desired object.
...type:Class.SkipElement
...type:Class.SkipBaseElement
...type:Class.SkipList
..returns:Pointer to the next free block in ClassPool. Type is $SkipElement*$ or $SkipBaseElement*$ respectively.
*/

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > __inline 
	SkipElement< TObject, TModus, TSpec, TStructuring > * 
	_getNext( SkipElement< TObject, TModus, TSpec, TStructuring > & me )
	{
		SEQAN_CHECKPOINT
		return me._next;
	}

/*
.Internal._setNext:
..summary:Set the next-pointer of the element. This is only used in @Class.ClassPool@ to govern memory blocks.
..cat:SkipList
..signature:setNext(element, next)
..param.element:The desired object.
...type:Class.SkipBaseElement
...type:Class.SkipElement
...type:Class.SkipList
..param.next:Pointer to next-object.
...type:Class.SkipBaseElement
...type:Class.SkipElement
...type:Class.SkipList
*/

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > __inline 
	void 
	_setNext( SkipElement< TObject, TModus, TSpec, TStructuring > & me,
			  SkipElement< TObject, TModus, TSpec, TStructuring > * next )
	{
		SEQAN_CHECKPOINT
		me._next = next;
	}


/*
.Internal._getHeight:
..summary:_get the height of the tower of which object is a part of.
..cat:SkipList
..signature:_getHeight(element)
..param.element:The desired object.
...type:Class.SkipElement
...type:Class.SkipBaseElement
..returns:The height of the tower related which contains $element$. The type is @Metafunction.Size.$Size< "element-type" >::Type$@
*/

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > __inline 
	typename Size< SkipElement< TObject, TModus, TSpec, TStructuring > >::Type 
	_getHeight( SkipElement< TObject, TModus, TSpec, TStructuring > & me )
	{
		SEQAN_CHECKPOINT
		return _getHeight( * _getDown( me ) );
	}

/**
.Function.key:
..summary:Get the the key of the element.
..cat:Map
..signature:key(element)
..param.element:The desired object.
..returns:The the key of the element.  Type is @Metafunction.Key.$Key< "element-type" >::Type$@.
..include:seqan/chaining.h
*/

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > __inline 
	typename Key< TObject >::Type
	key( SkipElement< TObject, TModus, TSpec, TStructuring > & me )
	{
		return me._key;
	}

	
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TKey > __inline
	void 
	setKey( SkipElement< TObject, TModus, TSpec, TStructuring > & me,
			TKey pos )
	{
		me._key= pos;
	}

/*
.Function.setKey:
..summary:Set the theKeyattribute of the element.
..cat:SkipList
..signature:setKey(element, theKey)
..param.element:The desired object.
...type:Class.SkipElement
..param.key:The desired object.
...type:Metafunction.Key.$Key< "element-type" >::Type
..include:seqan/chaining.h
*/

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TKey >
	void 
	setKey( SkipElement< TObject, TModus, TSpec, TStructuring > & me,
			TKey pos );


/*
.Internal._performDestructorAction:
..summary:Action that should be performed while destruction of the element.
..cat:SkipList
..signature:_performDestructorAction(element)
..param.element:The element.
...type:Class.SkipElement
...type:Class.SkipBaseElement
..Remarks:In the default case, no special operations are performed. In the $RangeTree$ spec, the destructor of the associated structure is called.
Can be used to perform operation s
*/


	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > __inline
	void
	_performDestructorAction( SkipElement< TObject, TModus, TSpec, TStructuring > & /*me*/ )
	{}

/*
.Function.cargo:
..summary:Get the cargo of the element.
..cat:SkipList
..signature:cargo(element)
..param.element:The desired element.
...type:Class.SkipElement
...type:Class.SkipBaseElement
..returns:A reference to the cargo of the element. Type is $TKey$.
..include:seqan/chaining.h
*/

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > __inline
	typename Cargo< SkipElement< TObject, TModus, TSpec, TStructuring > >::Type &
	cargo( SkipElement< TObject, TModus, TSpec, TStructuring > & me );

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TCargo > __inline
	void
	setCargo( SkipElement< TObject, TModus, TSpec, TStructuring > & me,
				TCargo cargo );

	
//___________________________________________________________________	
	
	
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > __inline 
	void 
	dump( SkipElement< TObject, TModus, TSpec, TStructuring > & me )
	{
		std::cout.width(5);
		if( key( me ) == minValue< typename Key< TObject >::Type >( ) )
			std::cout << std::left << "L";
		else
			std::cout<< key( me );
	}


	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > __inline 
	void 
	printElementValue( SkipElement< TObject, TModus, TSpec, TStructuring > & me )
	{
		SkipBaseElement< TObject, TModus, TSpec, TStructuring > * buffer =  _getDown( me );
		if( buffer != NULL )
			printElementValue( buffer );
		else
			std::cout << "_ ";
	}


	template < typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TKey >
	__inline void
	valueConstruct( SkipElement< TObject, TModus, TSpec, TStructuring > * it,
				    SkipElement< TObject, TModus, TSpec, TStructuring > * right,
					SkipBaseElement< TObject, TModus, TSpec, TStructuring > * down,
					TKey theKey )
	{
	SEQAN_CHECKPOINT
		new( it ) SkipElement< TObject, TModus, TSpec, TStructuring >( right, down, theKey );
	}


template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct SkipElement
{
public:

	union{
			// pointer to the SkipElement on the level of the tower on the right of this tower
		SkipElement< TObject, TModus, TSpec, TStructuring > * _right;
			// pointer needed for memory pool allocator, points to next free memory block
		SkipElement< TObject, TModus, TSpec, TStructuring > * _next;
	};
		
		// pointer to the underlying SkipBaseElement
	SkipBaseElement< TObject, TModus, TSpec, TStructuring > * _down;
		// the key
	typename Key< TObject >::Type _key;
		// related SkipList for entries of lower dimension
		// only constructed if needed by a range query, 
		// and number of elements is not under threshold
	typename Cargo< SkipElement< TObject, TModus, TSpec, TStructuring > >::Type _cargo;


		// don't allow copy assignment operator and copy constructor
		// (without the SkipList environment, SkipElements are senseless. Forbidding avoids unintended corruption of the
		// SkipList containing the original SkipElement by deleting the copy)	
	SkipElement( SkipElement & elem );
	void operator=( const SkipElement & elem );

	friend __inline
	typename Cargo< SkipElement< TObject, TModus, TSpec, TStructuring > >::Type *
	cargo( SkipElement< TObject, TModus, TSpec, TStructuring > & me )
	{
		return &me._cargo;
	}

	template< typename TCargo > friend __inline
	void
	setCargo( SkipElement< TObject, TModus, TSpec, TStructuring > & me ,
				TCargo cargo )
	{
		me._cargo = cargo;
	}
/*
	friend 
	void 
	dump<>( SkipElement< TObject, TModus, TSpec, TStructuring > & me );

	friend 
	void 
	printElementValue<>( SkipElement< TObject, TModus, TSpec, TStructuring > & me );

*/
public:

	SkipElement( void )
	: _right(NULL)
	, _down(NULL)
	{
		_initCargo( this, _cargo );
		#ifdef MEMTEST
			++created;
			++current;
			if( current > max )
				max = current;
		#endif
	}

	template< typename TKey >
	SkipElement( SkipElement< TObject, TModus, TSpec, TStructuring > * right,
				 SkipBaseElement< TObject, TModus, TSpec, TStructuring > * down,
				 TKey theKey )
	: _right(right)
	, _down(down)
	, _key(theKey)
	{
		_initCargo( this, _cargo );
	}


	~SkipElement( void )
	{	
		_performDestructorAction( *this );
		#ifdef MEMTEST
			--current;
		#endif
	}

	
}; // struct SkipElement

} // namespace seqan

#endif // SKIP_ELEMENT_H
