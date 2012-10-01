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
*	Contains struct SkipList< TObject, TModus, TSpec, TStructuring > implementation
*
*/

//SEQAN_NO_DDDOC: do not generate documentation for this file


#ifndef SEQAN_HEADER_SKIP_LIST_STATIC_H
#define SEQAN_HEADER_SKIP_LIST_STATIC_H

namespace seqan
{

/**
.Class.SkipList:
..cat:SkipList
..summary:A SkipList is a randomized alternative for rooted trees. Offers fast searching, insertion and deletion operations
objects. Saved objects are sorted with respect to their "key"-attribute.
..signature:SkipList< TObject, [ TModus, TSpec, TStructuring] >
..param.TObject:Type of stored objects.
..param.TModus:Modus of operation of a SkipList. A SkipList can either be dynamic or static. 
SkipListDynamic Skip Lists admit insertion and deletion operations of elements, but the construction time is higher 
compared to a static SkipList. Default is SkipListDynamic.
..param.TSpec:Specialization of the SkipList.
..param.TStructuring:Parameter to specify whether the SkipList uses Deferred Data Structuring or not.
..remarks:The object given to the SkipList should offer the following functions:
..remarks:$key( obj )$: returns the key of the object.
..remarks:$setKey( obj, k )$: set the key of the object to k.
..remarks:In contrast to STL-like containers, the objects are not cloned by the SkipList. It only supports searching operations on a set of objects. This set must be handled by the user.
..include:seqan/chaining.h
*/
	
		// no append operations on skip lists
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TParam > inline 
	void
	append( SkipList< TObject, TModus, TSpec, TStructuring > & /*me*/ ,
			TParam & /*param*/)
	{
		// do nothing
	}

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TParam1, typename TParam2 > inline 
	void
	append( SkipList< TObject, TModus, TSpec, TStructuring > & /*me*/ ,
			TParam1 & /*param1*/,
			TParam2 & /*param2*/ )
	{
		// do nothing
	}

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TParam > inline 
	void
	appendValue( SkipList< TObject, TModus, TSpec, TStructuring > & /*me*/ ,
					TParam & /*param*/)
	{
		// do nothing
	}

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TParam1, typename TParam2 > inline 
	void
	appendValue( SkipList< TObject, TModus, TSpec, TStructuring > & /*me*/ ,
				TParam1 & /*param1*/,
				TParam2 & /*param2*/ )
	{
		// do nothing
	}

		// no assign operations on skip lists
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TParam > inline 
	void
	assign( SkipList< TObject, TModus, TSpec, TStructuring > & /*me*/,
			TParam & /*param*/)
	{
		// do nothing
	}

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TParam1, typename TParam2 > inline 
	void
	assign( SkipList< TObject, TModus, TSpec, TStructuring > & /*me*/,
			TParam1 & /*param1*/,
			TParam2 & /*param2*/ )
	{
		// do nothing
	}

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TParam > inline 
	void
	assignValue( SkipList< TObject, TModus, TSpec, TStructuring > & /*me*/ ,
					TParam & /*param*/)
	{
		// do nothing
	}

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TParam1, typename TParam2 > inline 
	void
	assignValue( SkipList< TObject, TModus, TSpec, TStructuring > & /*me*/ ,
				TParam1 & /*param1*/,
				TParam2 & /*param2*/ )
	{
		// do nothing
	}


	// get the beginning of the skiplist
		// i.e. the left bording element with score - infinity
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
	inline typename Iterator< SkipList< TObject, TModus, TSpec, TStructuring > , Standard >::Type 
	_beginDefault( SkipList< TObject, TModus, TSpec, TStructuring > & me,
				   Standard)
	{
	SEQAN_CHECKPOINT
		return _getSucc( *me._baseStore );
	}

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
	inline typename Iterator< SkipList< TObject, TModus, TSpec, TStructuring > , Standard >::Type 
	_beginDefault( SkipList< TObject, TModus, TSpec, TStructuring > const & me,
				   Standard)
	{
	SEQAN_CHECKPOINT
		return _getSucc( *me._baseStore );
	}


		// capacity
		// static case: the size of the list
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline 
	typename Size< SkipList< TObject, TModus, TSpec, TStructuring > >::Type
	capacity( SkipList< TObject, TModus, TSpec, TStructuring > & me)
	{
	SEQAN_CHECKPOINT
		return length(me);
	}

		// dynamic case: the list can hold an infinite number of objevts
	template< typename TObject, typename TSpec, typename TStructuring > inline 
	typename Size< SkipList< TObject, SkipListDynamic, TSpec, TStructuring > >::Type
	capacity( SkipList< TObject, SkipListDynamic, TSpec, TStructuring > & /*me*/)
	{
	SEQAN_CHECKPOINT
		return maxValue< typename Size< SkipList< TObject, SkipListDynamic, TSpec, TStructuring > >::Type >();
	}


	template < typename TObject, typename TModus, typename TSpec, typename TStructuring >
	inline typename Iterator< SkipList< TObject, TModus, TSpec, TStructuring >, Standard>::Type 
	_endDefault( SkipList< TObject, TModus, TSpec, TStructuring > & me,
					Standard)
	{
	SEQAN_CHECKPOINT
		return _getRightBorder( me );
	}
	template < typename TObject, typename TModus, typename TSpec, typename TStructuring >
	inline typename Iterator< SkipList< TObject, TModus, TSpec, TStructuring > const, Standard>::Type 
	_endDefault( SkipList< TObject, TModus, TSpec, TStructuring > const & me,
					Standard)
	{
	SEQAN_CHECKPOINT
		return _getRightBorder( me );
	}

		// the length == size
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
	typename Size< SkipList< TObject, TModus, TSpec, TStructuring > >::Type
	length( SkipList< TObject, TModus, TSpec, TStructuring > & me )
	{
	SEQAN_CHECKPOINT
		return me._numOfElements;
	}

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TSize > inline
	void 
	_setLength( SkipList< TObject, TModus, TSpec, TStructuring > & me,
				TSize num_elems )
	{
		SEQAN_CHECKPOINT
		me._numOfElements = num_elems;
	}

		// no moving objects in skip lists
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TPos >
	void
	moveValue( SkipList< TObject, TModus, TSpec, TStructuring > & /*me*/, 
				TPos & /*pos*/)
	{
		// do nothing
	}

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TPos >
	void
	moveValue( SkipList< TObject, TModus, TSpec, TStructuring > const & /*me*/, 
				TPos & /*pos*/)
	{
		// do nothing
	}

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TPos >
	void
	moveValue( SkipList< TObject, TModus, TSpec, TStructuring > & /*me*/, 
				TPos const & /*pos*/)
	{
		// do nothing
	}

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TPos >
	void
	moveValue( SkipList< TObject, TModus, TSpec, TStructuring > const & /*me*/, 
				TPos const & /*pos*/)
	{
		// do nothing
	}


			// no replacements in skip lists
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TPos >
	void
	replace( SkipList< TObject, TModus, TSpec, TStructuring > & /*me*/, 
				TPos & /*pos*/)
	{
	SEQAN_CHECKPOINT
		// do nothing
	}

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TPos >
	void
	replace( SkipList< TObject, TModus, TSpec, TStructuring >  const & /*me*/, 
				TPos & /*pos*/)
	{
	SEQAN_CHECKPOINT
		// do nothing
	}

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TPos >
	void
	replace( SkipList< TObject, TModus, TSpec, TStructuring > & /*me*/, 
				TPos const & /*pos*/)
	{
	SEQAN_CHECKPOINT
		// do nothing
	}


		// the value at a given position
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TPos>
	inline typename Reference< SkipList< TObject, TModus, TSpec, TStructuring > >::Type
	value( SkipList< TObject, TModus, TSpec, TStructuring > & me, 
			TPos pos)
	{
	SEQAN_CHECKPOINT
		SkipBaseElement< TObject, TModus, TSpec, TStructuring > * base = me._baseStore;
		goNext( base );
		while( pos > 0 )
		{
			goNext( base );
			--pos;
		}
		return value( *base );
	}


	template< typename TObject, typename TSpec, typename TStructuring, typename TPos>
	inline typename Reference< SkipList< TObject, SkipListStatic, TSpec, TStructuring > >::Type
	value( SkipList< TObject, SkipListStatic, TSpec, TStructuring > & me, 
			TPos pos)
	{
	SEQAN_CHECKPOINT
		return value( me._baseStore + pos + 1 );
	} 

///////////////////////////////////////////////////////
//
//	member accessors (internal)
//
///////////////////////////////////////////////////////
	

		// get root node
		// i.e. the skip element in left bording tower in the current layer
		// current layer == 0 at the beginning
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline 
	SkipElement< TObject, TModus, TSpec, TStructuring > * 
	_getRoot( SkipList< TObject, TModus, TSpec, TStructuring > & me )
	{
		SEQAN_CHECKPOINT
		return _getLeftSideStore( me )+ _getCurrentLayer( me ) - 1;
	}

		// get the element allocator of a skiplist
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline 
	Allocator< ClassPool< SkipElement< TObject, TModus, TSpec, TStructuring >, Limited, SimpleAllocator > > & 
	_getElementAlloc( SkipList< TObject, TModus, TSpec, TStructuring > & me )
	{
		SEQAN_CHECKPOINT
		SEQAN_ASSERT_MSG(&me._elementAlloc != NULL, "Allocator for layer elements not initialized");
		return me._elementAlloc;
	}

		// get the base element allocator of skiplist
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline 
	Allocator< ClassPool< SkipBaseElement< TObject, TModus, TSpec, TStructuring >, Unlimited, SimpleAllocator> > & 
	_getBaseAlloc( SkipList< TObject, TModus, TSpec, TStructuring > & me )
	{
		SEQAN_CHECKPOINT
		SEQAN_ASSERT_MSG(&me._baseAlloc != NULL, "Allocator for base elements not initialized");
		return me._baseAlloc;
	}

		
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline 
	void 
	_setNext(	SkipList< TObject, TModus, TSpec, TStructuring > & me, 
				SkipList< TObject, TModus, TSpec, TStructuring > * next )
	{
	SEQAN_CHECKPOINT
		me._next = next;
	}


		//get the right bording element
		//i.e. the lement with theKey== + infinity
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline 
	SkipBaseElement< TObject, TModus, TSpec, TStructuring > * 
	_getRightBorder( SkipList< TObject, TModus, TSpec, TStructuring > & me )
	{
		SEQAN_CHECKPOINT
		SEQAN_ASSERT_MSG(me._rightBorder != NULL, "Right border not created");
		return me._rightBorder;
	}

		// get the left side store of the skip list
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
	SkipElement< TObject, TModus, TSpec, TStructuring > * 
	_getLeftSideStore( SkipList< TObject, TModus, TSpec, TStructuring > & me )
	{
		SEQAN_CHECKPOINT
		SEQAN_ASSERT_MSG(me._leftSideStore != NULL, "_leftSideStore not initialized");
		return me._leftSideStore;
	}

		// get the base layer of the skip list
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
	SkipBaseElement< TObject, TModus, TSpec, TStructuring > *
	_getBaseStore( SkipList< TObject, TModus, TSpec, TStructuring > & me )
	{
		SEQAN_CHECKPOINT
		SEQAN_ASSERT_MSG(me._baseStore != NULL, "_baseStore not initialized");
		return me._baseStore;
	}

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
	SkipElement< TObject, TModus, TSpec, TStructuring > ** 
	_getSearchPath( SkipList< TObject, TModus, TSpec, TStructuring > & me )
	{
		SEQAN_CHECKPOINT
		return _getSearchPath( me._sp );
	}


	
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
	typename Size< SkipList< TObject, TModus, TSpec, TStructuring > >::Type 
	_getCurrentLayer( SkipList< TObject, TModus, TSpec, TStructuring > & me )
	{
		SEQAN_CHECKPOINT
		return _getHeight( *me._baseStore );
	}

		// renew current layer
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TSize > inline 
	void 
	_setCurrentLayer(	SkipList< TObject, TModus, TSpec, TStructuring > & me,
						TSize layer )
	{
		SEQAN_CHECKPOINT
		SEQAN_ASSERT_GT_MSG(layer, _getCurrentLayer( me ), "Setting smaller value for _currentLayer");
		_setHeight( *me._baseStore, layer );
	}


/////////////////////////////////////////////////////////////////////////////////
//
//								 struct SkipList
//
/////////////////////////////////////////////////////////////////////////////////


template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct SkipList
{	
	// private members:
		// container for elements in the lowest layer
	union{
		SkipBaseElement< TObject, TModus, TSpec, TStructuring > * _baseStore;
		SkipList< TObject, TModus, TSpec, TStructuring > * _next;
	};
		// container for elements on the left side
		// search operations start from this bording elements
	SkipElement< TObject, TModus, TSpec, TStructuring > * _leftSideStore;

		// number of elements in the list
	typename Size< SkipList< TObject, TModus, TSpec, TStructuring > >::Type _numOfElements;

		// pool allocators
	Allocator< ClassPool< SkipBaseElement< TObject, TModus, TSpec, TStructuring >, Unlimited > > _baseAlloc;
	Allocator< ClassPool< SkipElement< TObject, TModus, TSpec, TStructuring >, Limited > > _elementAlloc;	
	
		// search path
	SearchPath_< TObject, TModus, TSpec, TStructuring > _sp;
		
		// border element for the right side
	SkipBaseElement< TObject, TModus, TSpec, TStructuring > * _rightBorder;

		// border objects
	TObject l_border_obj;
	TObject r_border_obj;

		// 
	bool _initialState;
	
//*************************************** internal help functions: ************************************************


private:

	SkipList & operator=( const SkipList & /*old*/ )
	{}

	SkipList( const SkipList & old )
		: _numOfElements( length( old ) ) // bording elements are not included in Entries array => numOfElements + 2 elements in base layer	
		, _baseAlloc( _numOfElements + 2 )	// guessing needed space: a basic layer must be possible
		, _elementAlloc( _numOfElements )
		, _sp( _getMaximalSLTowerHeight( _numOfElements ) )
		, _initialState( true )
	{
			// construct bording elements
		setKey( l_border_obj, minValue< typename Key< TObject >::Type >() );
		setKey( r_border_obj, maxValue< typename Key< TObject >::Type >() );

		SkipBaseElement< TObject, TModus, TSpec, TStructuring > * base_right;
		
		_initBases( *this, _baseStore, base_right, begin( old ), end( old ), _numOfElements ); 

		_initSL( *this, _baseStore, base_right, _numOfElements );
		_rightBorder = base_right;
	}

public:
				
	SkipList(void)
		: _numOfElements( 0 ) // bording elements are not included in Entries array => numOfElements + 2 elements in base layer	
		, _baseAlloc( 100  )	// guessing needed space: a basic layer must be possible
		, _elementAlloc( 100 )
		, _sp( 10 )
		, _initialState( true )
	{
		mtRandInit();

			// construct bording elements
		setKey( l_border_obj, minValue< typename Key< TObject >::Type >() );
		setKey( r_border_obj, maxValue< typename Key< TObject >::Type >() );
			
	}

	template< typename TContainer >
	SkipList( TContainer & cont )
		: _numOfElements( length( cont ) ) // bording elements are not included in Entries array => numOfElements + 2 elements in base layer	
		, _baseAlloc( _numOfElements + 2 )	// guessing needed space: a basic layer must be possible
		, _elementAlloc( _numOfElements )
		, _sp( _getMaximalSLTowerHeight( _numOfElements ) )
		, _initialState( true )
	{
		mtRandInit();

			// construct bording elements
		setKey( l_border_obj, minValue< typename Key< TObject >::Type >() );
		setKey( r_border_obj, maxValue< typename Key< TObject >::Type >() );

		SkipBaseElement< TObject, TModus, TSpec, TStructuring > * base_right;
		
		_initBases( *this, _baseStore, base_right, begin( cont ), end( cont ), _numOfElements ); 

		_initSL( *this, _baseStore, base_right, _numOfElements );
		_rightBorder = base_right;
	}

		// Constructor
		// Needs a range, defined by to iterators

		// TODO: für allgmeinen iterator verfügbar machen
	template< typename TIterator >
	SkipList( TIterator beg, TIterator end )
		: _numOfElements( end - beg ) // bording elements are not included in Entries array => numOfElements + 2 elements in base layer	
		, _baseAlloc( _numOfElements + 2 )	// guessing needed space: a basic layer must be possible
		, _elementAlloc( _numOfElements )
		, _sp( _getMaximalSLTowerHeight( _numOfElements ) )
		, _initialState( true )
	{			
		mtRandInit();

		setKey( l_border_obj, minValue< typename Key< TObject >::Type >() );
		setKey( r_border_obj, maxValue< typename Key< TObject >::Type >() );

			// construct bording elements
		SkipBaseElement< TObject, TModus, TSpec, TStructuring > * base_right;

		_initBases( *this, _baseStore, base_right, beg, end, _numOfElements ); 

		_initSL( *this, _baseStore, base_right, _numOfElements );
		_rightBorder = base_right;
	}
	

	~SkipList(void)
	{
		_clearSearchPath( this->_sp, _getMaximalSLTowerHeight( _numOfElements ) );
	}

	// Debug print methods
private:

	template< typename TSize1, typename TSize2 > friend
	void 
	printLayer(	SkipList< TObject, TModus, TSpec, TStructuring > & me,
				TSize1 layer,
				TSize2 column )
	{
		if( layer == 0 )
		{
			for( typename seqan::Size< SkipList< TObject, TModus, TSpec, TStructuring > >::Type j = 0; j < 11; ++j )
			{
				std::cout<< "______";
			}
			std::cout<<std::endl;
			typename Iterator< SkipList< TObject, TModus, TSpec, TStructuring > >::Type temp = begin( me );
			goPrevious( temp );
			goFurther( temp, column );
			typename Iterator< SkipList< TObject, TModus, TSpec, TStructuring > >::Type tempEnd = temp;
			goFurther( tempEnd, 11 );
			while( temp != end( me ) && temp != tempEnd )
			{
				std::cout.width(7);
				if( key( temp ) == minValue< typename Key< TObject >::Type >( ) )
					std::cout << std::left << "L";
				else
					std::cout << std::left << key( temp );
				goNext( temp );
			}
			//std::cout<<std::endl;
			//printCounts( me );
			//std::cout<<std::endl;
			//printSorting( me );
			//std::cout<<std::endl;
			//printHeights( me );
			//std::cout<<std::endl;
			//printValues( me );
		}
		else if( layer <= _getCurrentLayer( me ) && !_getInitialState( me ) )
		{
			typename Iterator< SkipList< TObject, TModus, TSpec, TStructuring > >::Type temp = begin( me );
			goPrevious( temp );
			goFurther( temp, column );
			typename Iterator< SkipList< TObject, TModus, TSpec, TStructuring > >::Type tempEnd = temp;
			goFurther( tempEnd, 11 );
			while( temp != end( me ) && temp != tempEnd )
			{
				if( _getHeight( *hostIterator( temp ) ) >= layer )
				{
					if( _getRight( *( &_getUp( *hostIterator( temp ) ) + layer - 1) ) )
					{
						std::stringstream s;
						if( key( temp ) == minValue< typename Key< TObject >::Type >( ) )
							s << std::left << "L";
						else
							s << key( temp );
						s << ">";
						std::cout.width(7);
						std::cout << std::left << s.str();
					}
					else
						dump( *( &_getUp( *hostIterator( temp ) ) + layer - 1 ) );
				}
				else 
				{
					std::cout.width(7);
					std::cout << " ";
				}
				goNext( temp );
			} 
		}
		std::cout<<std::endl;
	}

	friend
	void 
	dump( SkipList< TObject, TModus, TSpec, TStructuring > & me )
	{
		int column = 0;

		while( column < length( me ) )
		{
			typename seqan::Size< SkipList< TObject, TModus, TSpec, TStructuring > >::Type j = _getMaximalSLTowerHeight( me._numOfElements );
			while( j > 0 ){
				printLayer( me, --j, column );
			}
			std::cout<< std::endl;
			column += 14;
		}
	}

	friend
	void 
	printCounts( SkipList< TObject, TModus, TSpec, TStructuring > & me )
	{
		SkipBaseElement< TObject, TModus, TSpec, TStructuring > * temp = &me._baseStore[0];
		for( typename seqan::Size< SkipList< TObject, TModus, TSpec, TStructuring > >::Type i = 0; i < me._numOfElements + 1; ++i ){
			std::cout.width(5);
			std::cout << std::left << _getCount( *temp );
			goNext( temp );
		}
	}

	friend
	void 
	printSorting( SkipList< TObject, TModus, TSpec, TStructuring > & me )
	{
		SkipBaseElement< TObject, TModus, TSpec, TStructuring >* temp = &me._baseStore[0];
		for( typename seqan::Size< SkipList< TObject, TModus, TSpec, TStructuring > >::Type i = 0; i < me._numOfElements + 1; ++i ){
			std::cout.width(5);
			if( _getRight( *temp ) != NULL )
				std::cout << std::left << "x";
			else
				std::cout << std::left << "o";
			goNext( temp );
		}
	}

	friend
	void
	printHeights( SkipList< TObject, TModus, TSpec, TStructuring > & me )
	{
		SkipBaseElement< TObject, TModus, TSpec, TStructuring >* temp = &me._baseStore[0];
		for( typename seqan::Size< SkipList< TObject, TModus, TSpec, TStructuring > >::Type i = 0; i < me._numOfElements + 1; ++i ){
			std::cout.width(5);
			std::cout << std::left << _getHeight( *temp );
			goNext( temp );
		}
	}

	friend
	void
	printValues( SkipList< TObject, TModus, TSpec, TStructuring > & /*me*/ )
	{
		//SkipBaseElement< TObject, TModus, TSpec, TStructuring >* temp = &me._baseStore[0];
		//for( typename seqan::Size< SkipList< TObject, TModus, TSpec, TStructuring > >::Type i = 0; i < me._numOfElements + 1; ++i ){
		//	std::cout.width(5);
		//	std::cout << std::left << getValue( getObject( *temp ) );
		//	temp = _getSucc( *temp );
		//}
	}
	

}; // struct SkipList

} // namespace seqan
#endif //SKIP_LIST_STATIC_H_

