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
*	Skip List Datastructure
*
*	Algorithms, helper functions, accessor functions of the SkipList
*
*/

//SEQAN_NO_DDDOC: do not generate documentation for this file


#ifndef SEQAN_HEADER_SKIP_LIST_BASE_H
#define SEQAN_HEADER_SKIP_LIST_BASE_H

namespace seqan
{

//////////////////////////////////////////////////////////////////////////////////////////
//
//	container adaptor methods - specs for the seqan default container functions
//
//////////////////////////////////////////////////////////////////////////////////////////

		
	template< typename TObject, typename TSpec, typename TStructuring >
	struct IsContiguous< SkipList< TObject, SkipListStatic, TSpec, TStructuring > >
	{
		enum { VALUE = true };
	};

	template< typename TObject, typename TSpec, typename TStructuring >
	struct IsContiguous< SkipList< TObject, SkipListDynamic, TSpec, TStructuring > >
	{
		enum { VALUE = false };
	};


////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//		dependend member classes
//
//
////////////////////////////////////////////////////////////////////////////////////////////////
	

////////////////////////////////////////////////////////////////////////////////////////////////
//
//		search path: saves elements which where traversed by a search operation
//		(needed by several other functions, e.g. insertion, deletion )
//
////////////////////////////////////////////////////////////////////////////////////////////////

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
	struct
	SearchPath_
	{
			// buffer for elements on searchng path
		SkipElement< TObject, TModus, TSpec, TStructuring > ** _searchPath;

		SearchPath_( )
		{
		}

		template< typename TSize >
		SearchPath_( TSize size )
		{
		SEQAN_CHECKPOINT
			allocate( this, _searchPath, size );
		}

		~SearchPath_( void )
		{
		SEQAN_CHECKPOINT
			SEQAN_ASSERT(_searchPath == NULL);
		}


	};

		// get the search path of a skip list
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
	SkipElement< TObject, TModus, TSpec, TStructuring > **
	_getSearchPath( SearchPath_< TObject, TModus, TSpec, TStructuring > & sp )
	{
	SEQAN_CHECKPOINT
		return sp._searchPath;
	}

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TSize > inline
	void
	_clearSearchPath( SearchPath_< TObject, TModus, TSpec, TStructuring > & sp,
						TSize size )
	{
	SEQAN_CHECKPOINT
		deallocate( sp, sp._searchPath, size );
		sp._searchPath = NULL;
	}

		// change the maximal size of the skip list
		// ( e.g. when the maximal tower height has changed after insertion operations )
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TSize > inline
	void
	_renewSearchPath( SkipList< TObject, TModus, TSpec, TStructuring > & list,
						TSize old_size,
						TSize new_size )
	{
	SEQAN_CHECKPOINT
		deallocate( list, list._sp._searchPath, old_size );
		allocate( list, list._sp._searchPath, new_size );
	}



///////////////////////////////////////////////////////
//
//	state of skip list
//
///////////////////////////////////////////////////////


	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
	bool
	_getInitialState( SkipList< TObject, TModus, TSpec, TStructuring > & list )
	{
	SEQAN_CHECKPOINT
		return list._initialState;
	}

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
	void
	_noInitialState( SkipList< TObject, TModus, TSpec, TStructuring > & list )
	{
	SEQAN_CHECKPOINT
		list._initialState = false;
	}

		// get the dimension
		// in default case, the dimension is 0
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
	typename Size< SkipList< TObject, TModus, TSpec, TStructuring > >::Type 
	dimension( SkipList< TObject, TModus, TSpec, TStructuring > & me );


///////////////////////////////////////////////////////////////////////////////////////////
//
//		helper functions
//
//	several usefull basic helper functions for different purposes 
//
///////////////////////////////////////////////////////////////////////////////////////////

		// get the maximal tower height, that depends on the number of elemnts in the skip list
	template< typename TNumber > inline
	TNumber
	_getMaximalSLTowerHeight( TNumber elements )
	{
		SEQAN_CHECKPOINT
		return log2( elements ) + 6;
	}

		// get the maximal tower height, that depends on the size type of skip list
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
	typename Size< SkipList< TObject, TModus, TSpec, TStructuring > >::Type
	_getMaximalSLTowerHeight( SkipList< TObject, TModus, TSpec, TStructuring > & me )
	{
		SEQAN_CHECKPOINT
		return _getMaximalSLTowerHeight( length( me ) );
	}

	// _throwCoin
		// generates height with geometric distribution of a fair coin
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename THeight > inline 
	typename seqan::Size< SkipList< TObject, TModus, TSpec, TStructuring > >::Type 
	_throwCoin( SkipList< TObject, TModus, TSpec, TStructuring > & /*list*/,
				THeight max_height)
	{
		SEQAN_CHECKPOINT
		typename Size< SkipList< TObject, TModus, TSpec> >::Type height = _geomRand< typename Size< SkipList< TObject, TModus, TSpec, TStructuring > >::Type >();
		//if( height < 4 )
		//	height = 4;
		//height -= 4;
		if( height > max_height )
			height = max_height;
		return height;
	}

///////////////////////////////////////////////////////////////////////////////////////////
//
//		helper functions
//
//	several usefull basic helper functions for different purposes 
//
///////////////////////////////////////////////////////////////////////////////////////////

		// swapBases
		// swaps the related objects of two skip base elements
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline 
	void 
	_swapBases(	SkipBaseElement< TObject, TModus, TSpec, TStructuring > & first,
				SkipBaseElement< TObject, TModus, TSpec, TStructuring > & second )
	{
		SEQAN_CHECKPOINT
		SEQAN_ASSERT_MSG( &first != NULL && &second != NULL, "An SkipBaseElements's object is corrupted");
		typename Key< TObject >::Type buffer_key= key( first );
		setKey( first, key( second ) );
		setKey( second, buffer_key );
		TObject * bufferObject = getObject( &first );
		_setObject( first, getObject( &second ) );
		_setObject( second, bufferObject );
	}


		//	_calcPivot
		//	calculates the pivot element from a given element on left side
		//	depends on the TModus spec of the SkipList

		// dynamic: base layer is a list
	template< typename TObject, typename TSpec, typename TStructuring > inline
	SkipBaseElement< TObject, SkipListDynamic, TSpec, TStructuring > *
	_calcPivot(	SkipBaseElement< TObject, SkipListDynamic, TSpec, TStructuring > & base )
	{
		SEQAN_CHECKPOINT
		SEQAN_ASSERT( &base != NULL );
		typename Size< SkipList< TObject, SkipListDynamic, TSpec, TStructuring > >::Type offset = mtRand() % _getCount( base ) + 1;
		SkipBaseElement< TObject, SkipListDynamic, TSpec, TStructuring > * pivot = &base;
		while( offset > 0 )
		{
			SEQAN_ASSERT( _getSucc( *pivot ) != NULL );
			goNext( pivot );
			--offset;
		}
		return pivot;
	}

		// static: base layer is an array
	template< typename TObject, typename TSpec, typename TStructuring > inline
	SkipBaseElement< TObject, SkipListStatic, TSpec, TStructuring > *
	_calcPivot(	SkipBaseElement< TObject, SkipListStatic, TSpec, TStructuring > & base )
	{
		SEQAN_CHECKPOINT
		SEQAN_ASSERT( &base != NULL );
		SEQAN_ASSERT_MSG( _getCount( base ) != 0, "Tried to sort element that is already sorted" );
		return  ( &base + ( mtRand() % _getCount( base ) + 1 ) );
	}


///////////////////////////////////////////////////////
//
//		sorting
//
///////////////////////////////////////////////////////

		//	_sort
		//	sort elements between elem and right with a step of quicksort
		//	depends on the TStructuring spec of the SkipList

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
	SkipBaseElement< TObject, TModus, TSpec, TStructuring > * 
	_sort(	SkipList< TObject, TModus, TSpec, TStructuring > & list,
			SkipBaseElement< TObject, TModus, TSpec, TStructuring > * elem,
			SkipBaseElement< TObject, TModus, TSpec, TStructuring > * right )
	{
		SEQAN_CHECKPOINT
		return _sort( list, elem, right, list );
	}


		// complete (=default) case
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TExtraParam > 
	SkipBaseElement< TObject, TModus, TSpec, TStructuring > * 
	_sort(	SkipList< TObject, TModus, TSpec, TStructuring > & /*list*/,
			SkipBaseElement< TObject, TModus, TSpec, TStructuring > * elem,
			SkipBaseElement< TObject, TModus, TSpec, TStructuring > * right_border,
			TExtraParam & param )
	{
			// security checking
		SEQAN_CHECKPOINT
		SEQAN_ASSERT_MSG( elem != NULL, "Tried to sort element that is NULL" );
			// helper variables
		SkipBaseElement< TObject, TModus, TSpec, TStructuring > * pivot = _calcPivot( *elem );
		typename Key< TObject >::Type pivot_key= key( *pivot, param );
			
			// buffer elements
		SkipBaseElement< TObject, TModus, TSpec, TStructuring > * index = elem;
		goNext( index );
		SkipBaseElement< TObject, TModus, TSpec, TStructuring > * pointer = index;
		goPrevious( right_border );
		_swapBases( *pivot, *right_border );
		
		typename seqan::Size< SkipList< TObject, TModus, TSpec, TStructuring > >::Type i_count = 0;
		while( pointer != right_border )
		{
			if( key( *pointer, param ) < pivot_key){
				_swapBases( *index, *pointer );
				goNext( index );
				i_count++;
			}
			goNext( pointer );
		}
		_swapBases( *index, *right_border );
		pivot = index;
		_renewConnects( *elem, *pivot, _getCount( *elem ), i_count );
		return pivot;
	}


	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
	void
	_sortEquals(	SkipList< TObject, TModus, TSpec, TStructuring > & list,
					SkipBaseElement< TObject, TModus, TSpec, TStructuring > * elem )
	{
		SEQAN_CHECKPOINT
		_sortEquals( list, elem, key( *elem ) );
	}

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TKey > inline
	void
	_sortEquals(	SkipList< TObject, TModus, TSpec, TStructuring > & /*list*/,
					SkipBaseElement< TObject, TModus, TSpec, TStructuring > * /*elem*/,
					TKey /*theKey*/)
	{	
		// do nothing
	}

	template< typename TObject, typename TModus, typename TSpec, typename TKey > 
	void
	_sortEquals(	SkipList< TObject, TModus, TSpec, Deferred > & /*list*/,
					SkipBaseElement< TObject, TModus, TSpec, Deferred > * elem,
					TKey theKey)
	{
			// security checking
		SEQAN_CHECKPOINT
		SEQAN_ASSERT_MSG( elem != NULL, "Tried to sort element that is NULL" );
			
		while( key( *_getSucc( *elem ) ) == theKey&& _getCount( *elem ) == 0 )
			goNext( elem );

		typename Size< SkipList< TObject, TModus, TSpec, Deferred > >::Type i_count = _getCount( *elem );
		SkipBaseElement< TObject, TModus, TSpec, Deferred > * right_border = _getRight( *elem );
		SkipBaseElement< TObject, TModus, TSpec, Deferred > * index = elem;
		goNext( index );
		SkipBaseElement< TObject, TModus, TSpec, Deferred > * pointer = index;
		
		while( pointer != right_border )
		{
			if( key( *pointer ) == theKey){
				_swapBases( *index, *pointer );
				_setRight( *elem, index );
				_setCount( *elem, 0 );
				goNext( elem );
				goNext( index );
				--i_count;
			}
			goNext( pointer );
		}
		if( index != elem ){
			goPrevious( index );
			_setCount( *index, i_count );
			_setRight( *index, right_border );
		}
	}

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TSize > inline
	void
	_renewDynConnects( SkipBaseElement< TObject, TModus, TSpec, TStructuring > & /*elem*/,
						SkipBaseElement< TObject, TModus, TSpec, TStructuring > & /*pivot*/,
						TSize /*elem_count*/,
						TSize /*i_count*/ )
	{
	}

	template< typename TObject, typename TModus, typename TSpec, typename TSize > inline
	void
	_renewDynConnects( SkipBaseElement< TObject, TModus, TSpec, Deferred > & elem,
						SkipBaseElement< TObject, TModus, TSpec, Deferred > & pivot,
						TSize elem_count,
						TSize i_count )
	{
		SEQAN_CHECKPOINT
		_setCount( pivot, elem_count - i_count - 1 );
		_setCount( elem, i_count );
	}


	template< typename TObject, typename TSpec, typename TSize > inline
	void
	_renewDynConnects( SkipBaseElement< TObject, SkipListDynamic, TSpec, Deferred > & elem,
						SkipBaseElement< TObject, SkipListDynamic, TSpec, Deferred > & pivot,
						TSize elem_count,
						TSize i_count )
	{
		SEQAN_CHECKPOINT
		_setLeft( *_getRight( elem ), &pivot );
		_setLeft( pivot, &elem );
		_setRight( pivot, _getRight( elem ) );
		_setRight( elem, &pivot );
		_setCount( pivot, elem_count - i_count - 1 );
		_setCount( elem, i_count );
	}


	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TSize > inline
	void
	_renewConnects( SkipBaseElement< TObject, TModus, TSpec, TStructuring > & elem,
					SkipBaseElement< TObject, TModus, TSpec, TStructuring > & pivot,
					TSize elem_count,
					TSize i_count )
	{
		SEQAN_CHECKPOINT
		_setCount( pivot, elem_count - i_count - 1 );
		_setCount( elem, i_count );
	}


	template< typename TObject, typename TSpec, typename TSize > inline
	void
	_renewConnects( SkipBaseElement< TObject, SkipListDynamic, TSpec, Deferred > & elem,
					SkipBaseElement< TObject, SkipListDynamic, TSpec, Deferred > & pivot,
					TSize elem_count,
					TSize i_count )
	{
		SEQAN_CHECKPOINT
		_setLeft( *_getRight( elem ), &pivot );
		_setLeft( pivot, &elem );
		_setRight( pivot, _getRight( elem ) );
		_setRight( elem, &pivot );
		_setCount( pivot, elem_count - i_count - 1 );
		_setCount( elem, i_count );
	}

		// quick sort recursive step
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
	void
	_sortRecursive( SkipList< TObject, TModus, TSpec, TStructuring > & list,
					SkipBaseElement< TObject, TModus, TSpec, TStructuring > * elem,
					SkipBaseElement< TObject, TModus, TSpec, TStructuring > * right )
	{
		SEQAN_CHECKPOINT
		if( _getCount( *elem ) != 0 )
		{
			SkipBaseElement< TObject, TModus, TSpec, TStructuring > * pivot = _sort( list, elem, right );
			_sortRecursive( list, elem, pivot );
			_sortRecursive( list, pivot, right );
		}	
	}


		// construct the towers for the complete skip list
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
	void
	_buildTowers( SkipList< TObject, TModus, TSpec, TStructuring > & list )
	{
		SEQAN_CHECKPOINT
		SkipBaseElement< TObject, TModus, TSpec, TStructuring > * buffer = _getBaseStore( list );
		SkipBaseElement< TObject, TModus, TSpec, TStructuring > * last = list._rightBorder;
		SkipElement< TObject, TModus, TSpec, TStructuring > ** search_path = new SkipElement< TObject, TModus, TSpec, TStructuring >*[ _getMaximalSLTowerHeight( length( list ) ) ];
		search_path[0] = &_getUp( *buffer ); 
		typename Size< SkipList< TObject, TModus, TSpec, TStructuring > >::Type height, maxHeight = _getMaximalSLTowerHeight( list );
		typename Key< TObject >::Type buffer_key = key( *buffer );
		goNext( buffer );

		while( buffer != last )
		{
			if( key( *buffer ) != buffer_key )
			{
				height = _throwCoin( list, maxHeight );
				if( height > 0 ){
					_add( list, buffer, height, search_path );
					_connectUpdate( list, buffer, height, search_path, list );
				}
			}
			buffer_key= key( *buffer );
			goNext( buffer  );
		}
		
		for( typename Size< SkipList< TObject, TModus, TSpec, TStructuring > >::Type i = _getCurrentLayer( list ) + 1; i > 0; --i )
			search_path[ i-1 ] = NULL;
	}


		// complete sorting and building of towers
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
	void
	_completeBuild( SkipList< TObject, TModus, TSpec, TStructuring > & list )
	{
		SEQAN_CHECKPOINT
		_sortRecursive( list, _getBaseStore( list ), list._rightBorder );
		_setHeight( *_getBaseStore( list ), 1 );
		_buildTowers( list );
	}

	template< typename TObject, typename TModus, typename TSpec > inline
	void
	_completeBuild( SkipList< TObject, TModus, TSpec, Deferred > & /*list*/ )
	{}

////////////////////////////////////////////////////////////////////////////////////////////////
//
//			searching
//
////////////////////////////////////////////////////////////////////////////////////////////////


		// searching from a given element in a tower
		// deferred skip list
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TParam, typename TKey >
	SkipBaseElement< TObject, TModus, TSpec, TStructuring > *
	_searchFrom(	SkipList< TObject, TModus, TSpec, TStructuring > & list,
					SkipElement< TObject, TModus, TSpec, TStructuring > * layer_element, 
					TKey theKey,
					SkipElement< TObject, TModus, TSpec, TStructuring > ** search_path,
					TParam & param )
	{
		SEQAN_CHECKPOINT
		SEQAN_ASSERT_MSG( theKey < maxValue< typename Key< TObject >::Type >( ), "search theKeyexceeds supremum" );
		SEQAN_ASSERT_MSG( theKey > minValue< typename Key< TObject >::Type >( ), "search theKeyexceeds infimum" );
	
		typename Size< SkipList< TObject, TModus, TSpec, TStructuring > >::Type height = layer_element - &_getUp( *_getDown( *layer_element ) ) + 1;
		SkipElement< TObject, TModus, TSpec, TStructuring > * temp_right = _getRight( *layer_element );
		SkipElement< TObject, TModus, TSpec, TStructuring > ** sp = &search_path[ height - 1];
			// search in higher layers		
		while( height > 0 ){
			while( key( *temp_right ) <= theKey){ 
				layer_element = temp_right;
				right( temp_right );
			}
			--height;
			*sp = layer_element;
			--sp;
			--layer_element;
			temp_right = _getRight( *layer_element );
		}
		++layer_element;
			// in the lowest layer, searching to the right
		SkipBaseElement< TObject, TModus, TSpec, TStructuring > * base_element = _getDown( *layer_element );
		SkipBaseElement< TObject, TModus, TSpec, TStructuring > * right_base = _getRight( *base_element );
		while( key( *right_base, param ) < theKey ){
			base_element = right_base;
			right( right_base );
		}
		if( key( *right_base, param ) == theKey && key( *base_element, param ) != theKey )
		{
			base_element = right_base;
		}

			// if the element found is in correct place and no unsorted elements follow,
			// it is the correct element. Otherwise search following intervall
		if( key( *base_element ) < theKey&& _getCount( *base_element ) > 0 )
		{
			base_element = _splitAction( list, base_element, theKey, search_path, param );
		}
			// handling for multiple elements
		right_base = _getRight( *base_element );
		if( key( *right_base, param ) <= theKey&& key( *base_element, param ) < theKey )
			return right_base;
		return  base_element;
	}

		// complete skip list
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TParam, typename TKey >
	SkipBaseElement< TObject, TModus, TSpec, TStructuring > *
	_searchFrom(	SkipList< TObject, TModus, TSpec, TStructuring > & /*list*/,
						SkipElement< TObject, TModus, TSpec, TStructuring > * layer_element, 
						TKey theKey,
						TParam & param )
	{
		SEQAN_CHECKPOINT
		SEQAN_ASSERT_MSG( theKey!= maxValue< typename Key< TObject >::Type >( ), "search theKeyexceeds supremum" );
		SEQAN_ASSERT_MSG( theKey!= minValue< typename Key< TObject >::Type >( ), "search theKeyexceeds infimum" );
	
		typename Size< SkipList< TObject, TModus, TSpec, TStructuring > >::Type height = layer_element - &_getUp( *_getDown( *layer_element ) ) + 1;

		SkipElement< TObject, TModus, TSpec, TStructuring > * temp_right = _getRight( *layer_element );

			// search in higher layers		
		while( height > 0 ){
			while( key( *temp_right ) <= theKey){ 
				layer_element = temp_right;
				right( temp_right );
			}
			--height;
			--layer_element;
			temp_right = _getRight( *layer_element );
		}
		++layer_element;
			// in the lowest layer, searching to the right
		SkipBaseElement< TObject, TModus, TSpec, TStructuring > * base_element = _getDown( *layer_element );
		SkipBaseElement< TObject, TModus, TSpec, TStructuring > * succ_element = _getSucc( *base_element );
		while( key( *succ_element, param ) < theKey){
			base_element = succ_element;
			goNext( succ_element );
		}
		if( key( *succ_element, param ) == theKey&& key( *base_element, param ) != theKey)
		{
			return succ_element;
		}
		return  base_element;
	}


		// search wrapper functions
	template< typename TObject, typename TModus, typename TSpec, typename TKey > inline
	SkipBaseElement< TObject, TModus, TSpec, Deferred > *
	_specialSearch( SkipList< TObject, TModus, TSpec, Deferred > & list,
					SkipElement< TObject, TModus, TSpec, Deferred > * layer_element, 
					TKey theKey)
	{ 
		SEQAN_CHECKPOINT
		SkipElement< TObject, TModus, TSpec, Deferred > ** search_path = _getSearchPath( list );
		return _searchFrom( list, layer_element, theKey, search_path, list );
	}

	
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TKey > inline
	SkipBaseElement< TObject, TModus, TSpec, TStructuring > *
	_specialSearch( SkipList< TObject, TModus, TSpec, TStructuring > & list,
					SkipElement< TObject, TModus, TSpec, TStructuring > * layer_element, 
					TKey theKey)
	{
		SEQAN_CHECKPOINT
		SkipElement< TObject, TModus, TSpec, TStructuring > ** search_path = _getSearchPath( list );
		return _searchFrom( list, layer_element, theKey, search_path, list );
	}

/**
.Function.searchElement:
..class:Class.SkipList
..summary:Get the leftmost SkipBaseElement element in a SkipList list, where $key(element) <= search_key$
..cat:SkipList
..signature:searchElement(list, search_key)
..remarks:If no Element exists which fullfills the conditions, the returned element is the right(maximal) bounding element.
..param.list:The list to be searched.
...type:Class.SkipList
..param.search_key:The key.
..include:seqan/chaining.h
*/

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TKey > inline
	typename Iterator< SkipList< TObject, TModus, TSpec, TStructuring > >::Type 
	searchElement(	SkipList< TObject, TModus, TSpec, TStructuring > & list, 
					TKey theKey )
	{
		SEQAN_ASSERT( theKey < maxValue< typename Key< TObject >::Type >() );
		SEQAN_CHECKPOINT
		if( _getInitialState( list ) )
		{
			_completeBuild( list );
			_noInitialState( list );
		}
			// starting in the root, search to the right and downwards
		return typename Iterator< SkipList< TObject, TModus, TSpec, TStructuring > >::Type( list, _specialSearch( list, _getRoot( list ), theKey ) );
	}

////////////////////////////////////////////////////////////////////////////////////////////////
//
//		structure building
//
////////////////////////////////////////////////////////////////////////////////////////////////

/*
.internal._split:
..summary:Manages sorting status of elements in the base layer 
..cat:SkipList
..signature:_split( list, elem, theKey, search_path )
..param.list:The affected skip list.
...type:Class.SkipList
..param.elem:Element in the base layer.
...type:Class.SkipBaseElement
..param.key:The search key
...type:Key type of TObject
..param.search_path:Buffer of traversed elements
...type:Ayrray of Pointers to SkipElement
*/

		// wrapper functions for split (split should only be performed on deferred skip lists)
	template< typename TObject, typename TModus, typename TSpec, typename TParam, typename TKey > inline
	SkipBaseElement< TObject, TModus, TSpec, Deferred > *
	_splitAction( SkipList< TObject, TModus, TSpec, Deferred > & list,
					SkipBaseElement< TObject, TModus, TSpec, Deferred > * base,
					TKey theKey, 
					SkipElement< TObject, TModus, TSpec, Deferred > ** search_path,
					TParam & param )
	{
		return _split( list, base, theKey, search_path, param );
	}

		// general case (not deferred)
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TParam, typename TKey > inline
	SkipBaseElement< TObject, TModus, TSpec, TStructuring > *
	_splitAction( SkipList< TObject, TModus, TSpec, TStructuring > & /*list*/,
					SkipBaseElement< TObject, TModus, TSpec, TStructuring > * base,
					TKey /*theKey*/, 
					SkipElement< TObject, TModus, TSpec, TStructuring > ** /*search_path*/,
					TParam & /*param*/ )
	{
		return base;
	}

		// spec for deferred sl		
	template< typename TObject, typename TModus, typename TSpec, typename TParam, typename TKey > 
	SkipBaseElement< TObject, TModus, TSpec, Deferred > * 
	_split( SkipList< TObject, TModus, TSpec, Deferred > & list,
			SkipBaseElement< TObject, TModus, TSpec, Deferred > * base,
			TKey theKey, 
			SkipElement< TObject, TModus, TSpec, Deferred > ** search_path,
			TParam & param )
	{	 
		SEQAN_CHECKPOINT
		SkipBaseElement< TObject, TModus, TSpec, Deferred > * elem = base;
		SEQAN_ASSERT_MSG( _getCount( *elem ) > 0, "count is < 0" );
		SkipBaseElement< TObject, TModus, TSpec, Deferred > * pivot ;
		typename Key< SkipList< TObject, TModus, TSpec, Deferred > >::Type pivot_key;
		typename Size< SkipList< TObject, TModus, TSpec, Deferred > >::Type height = 0;
		typename Size< SkipList< TObject, TModus, TSpec, Deferred > >::Type maxHeight = _getMaximalSLTowerHeight( list );
		while( _getCount( *elem ) > 0 ){

				// sorting
			pivot = _sort( list, elem, _getRight( *elem ), param );
			
				// helper variables
			pivot_key= key( *pivot, param );

			if( key( ( *elem ), param ) != pivot_key)
			{
				height = _throwCoin( list, maxHeight );
				if( pivot_key<= theKey){
					if( height > 0){
						_add( list, pivot, height, search_path );
						_connectUpdate( list, pivot, height, search_path, param );
					}
					elem = pivot;
				}
				else{
					if( height > 0){
						_add( list, pivot, height, search_path );
						_connect( list, pivot, height, search_path, param );
					}
				}
			}
			else
				elem = pivot;
			if( key( *elem, param ) == theKey)
				break;
		}
		if( key( *elem, param ) != key( *base, param ) )
			return elem;
		else
			return base;
	}

	
		// reconnect the pointers after building a tower
		// 
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TParam, typename TSize > inline
	void 
	_connectUpdate(	SkipList< TObject, TModus, TSpec, TStructuring > & /*list*/,
						SkipBaseElement< TObject, TModus, TSpec, TStructuring > * base,
						TSize height,
						SkipElement< TObject, TModus, TSpec, TStructuring > ** search_path,
						TParam & param )
	{
		SEQAN_CHECKPOINT
		SkipElement< TObject, TModus, TSpec, TStructuring > * buffer = &_getUp( *base );
		SkipElement< TObject, TModus, TSpec, TStructuring > * tower_top = buffer + height;
		typename Key< TObject >::Type theKey= key( *base, param );
		SkipElement< TObject, TModus, TSpec, TStructuring > ** sp = search_path;
		while( buffer != tower_top ){
			new( buffer ) SkipElement< TObject, TModus, TSpec, TStructuring >( _getRight( **sp ), base, theKey);
			_setRight( **sp, buffer );
			*sp = buffer;
			++buffer;
			++sp;
		}
	}

		// reconnect the pointers after building a tower
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TParam, typename TSize > inline
	void 
	_connect(	SkipList< TObject, TModus, TSpec, TStructuring > & /*list*/,
				SkipBaseElement< TObject, TModus, TSpec, TStructuring > * base,
				TSize height,
				SkipElement< TObject, TModus, TSpec, TStructuring > ** search_path,
				TParam & param )
	{
		SEQAN_CHECKPOINT
		SkipElement< TObject, TModus, TSpec, TStructuring > * buffer = &_getUp( *base );
		SkipElement< TObject, TModus, TSpec, TStructuring > * tower_top = buffer + height;
		typename Key< TObject >::Type theKey= key( *base, param );
		SkipElement< TObject, TModus, TSpec, TStructuring > ** sp = search_path;
		while( buffer != tower_top ){
			new( buffer ) SkipElement< TObject, TModus, TSpec, TStructuring >( _getRight( **sp ), base, theKey);
			_setRight( **sp, buffer );
			++buffer;
			++sp;
		}
	}


		// add a new tower to a base element
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename THeight > inline
	void 
	_add(	SkipList< TObject, TModus, TSpec, TStructuring >  & list,
			SkipBaseElement< TObject, TModus, TSpec, TStructuring > * base,
			THeight height,
			SkipElement< TObject, TModus, TSpec, TStructuring > ** search_path )
	{			
		SEQAN_CHECKPOINT
		if( height > _getCurrentLayer( list ) ){
			for( THeight i = _getCurrentLayer( list ); i < height; ++i )
				search_path[i] = &( _getLeftSideStore( list )[i] );
			_setCurrentLayer( list, height );
		}
		SkipElement< TObject, TModus, TSpec, TStructuring > * tower;
		allocate( _getElementAlloc( list ), tower, height );
		_setUp( *base, *tower );
		_setHeight( *base, height );
	}

///////////////////////////////////////////////////////
//
//	initialization
//
///////////////////////////////////////////////////////

	
		// setting pointers between elements in base layer
		
			// static case
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
	void
	_setDynConnects( SkipBaseElement< TObject, TModus, TSpec, TStructuring > * /*pred*/,
						SkipBaseElement< TObject, TModus, TSpec, TStructuring > * /*succ*/ )
	{
		// do nothing
	}
			
			// dynamic case
	template< typename TObject, typename TSpec, typename TStructuring > inline
	void
	_setDynConnects( SkipBaseElement< TObject, SkipListDynamic, TSpec, TStructuring > * pred,
						SkipBaseElement< TObject, SkipListDynamic, TSpec, TStructuring > * succ )
	{
		SEQAN_CHECKPOINT
		_setSucc( *pred, succ );
		_setPred( *succ, pred );
	}


	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
	void
	_setDefConnects( SkipBaseElement< TObject, TModus, TSpec, TStructuring > * /*left*/,
						SkipBaseElement< TObject, TModus, TSpec, TStructuring > * /*right*/ )
	{
		// do nothing
	}


	template< typename TObject, typename TSpec > inline
	void
	_setDefConnects( SkipBaseElement< TObject, SkipListDynamic, TSpec, Deferred > * left,
						SkipBaseElement< TObject, SkipListDynamic, TSpec, Deferred > * right )
	{
		SEQAN_CHECKPOINT
		_setRight( *left, right );
		_setLeft( *right, left );
	}


	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TIter, typename TSize > inline
	void
	_initBases( SkipList< TObject, TModus, TSpec, TStructuring > & list, 
				SkipBaseElement< TObject, TModus, TSpec, TStructuring > & firstBase,
				SkipBaseElement< TObject, TModus, TSpec, TStructuring > *& lastBase,
				TIter & firstData,
				TIter & lastData,
				TSize numEntries );

	template< typename TObject, typename TSpec, typename TStructuring, typename TIter, typename TSize > inline
	void
	_initBases( SkipList< TObject, SkipListStatic, TSpec, TStructuring > & list, 
				SkipBaseElement< TObject, SkipListStatic, TSpec, TStructuring > * firstBase,
				SkipBaseElement< TObject, SkipListStatic, TSpec, TStructuring > *& lastBase,
				TIter & firstData,
				TIter & lastData,
				TSize numEntries )
	{
		SEQAN_CHECKPOINT
		allocate( _getBaseAlloc( list ), firstBase, numEntries + 2 );
		list._baseStore = firstBase;

		valueConstruct( firstBase, &list.l_border_obj, minValue< typename Key< TObject >::Type >() );

		++firstBase;

		while( firstData != lastData )
		{
			valueConstruct( firstBase, &value( firstData ), key( *firstData ) );
			++firstBase;
			++firstData;
		}
		lastBase = firstBase;
		valueConstruct( lastBase, &list.r_border_obj, maxValue< typename Key< TObject >::Type >() );
		firstBase = list._baseStore;
		_setDefConnects( lastBase, lastBase );
		_setDefConnects( firstBase, lastBase );
	}

	template< typename TObject, typename TSpec, typename TStructuring, typename TIter, typename TSize > inline
	void
	_initBases( SkipList< TObject, SkipListDynamic, TSpec, TStructuring > & list, 
				SkipBaseElement< TObject, SkipListDynamic, TSpec, TStructuring > * firstBase,
				SkipBaseElement< TObject, SkipListDynamic, TSpec, TStructuring > *& lastBase,
				TIter & firstData,
				TIter & lastData,
				TSize numEntries )
	{
		SEQAN_CHECKPOINT
		allocate( _getBaseAlloc( list ), firstBase, numEntries + 2 );
		list._baseStore = firstBase;

		valueConstruct( firstBase, &list.l_border_obj, minValue< typename Key< TObject >::Type >() );
		SkipBaseElement< TObject, SkipListDynamic, TSpec, TStructuring > * previous = firstBase;
		++firstBase;

		while( firstData != lastData )
		{
			valueConstruct( firstBase, &value( firstData ), key( *firstData ) );
			_setDynConnects( previous, firstBase );
			++firstBase;
			++previous;
			++firstData;
		}
		lastBase = firstBase;
		valueConstruct( lastBase, &list.r_border_obj, maxValue< typename Key< TObject >::Type >() );
		_setDynConnects( previous, lastBase );
		_setSucc( *lastBase, lastBase );
		firstBase = list._baseStore;
		_setDefConnects( lastBase, lastBase );
		_setDefConnects( firstBase, lastBase );
	}


	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TSize > inline
	void
	_initSL(	SkipList< TObject, TModus, TSpec, TStructuring > & list,
				SkipBaseElement< TObject, TModus, TSpec, TStructuring > * first_base,
				SkipBaseElement< TObject, TModus, TSpec, TStructuring > * /*last_base*/,
				TSize numEntries )
	{
			// allocate space for base elements and base layer
		SEQAN_CHECKPOINT
		SkipElement< TObject, TModus, TSpec, TStructuring  > * _rightBorder;
		allocate( _getElementAlloc( list ), _rightBorder, 2 );
		arrayConstruct( _rightBorder, _rightBorder + 2 );
		allocate( _getElementAlloc( list ), list._leftSideStore, _getMaximalSLTowerHeight( numEntries ) );
				// construct bording elements
				// right side
		_setDown( *_rightBorder, &list._baseStore[ numEntries + 1 ] );
		_setHeight( list._baseStore[ numEntries + 1 ], 1 );
		_setRight( * _rightBorder, _rightBorder );
		setKey( *_rightBorder, maxValue< typename Key< TObject >::Type >() );
		_setUp( list._baseStore[ numEntries + 1 ], *_rightBorder );
		
				// ... left side
		_setDown( *list._leftSideStore, first_base );
		_setHeight( *first_base, 1 );
		_setUp( *first_base, *list._leftSideStore );

		typename Key< TObject >::Type left_border_key= key( list.l_border_obj );
		SkipElement< TObject, TModus, TSpec, TStructuring > * buffer = list._leftSideStore;
		for( typename Size< SkipList< TObject, TModus, TSpec, TStructuring > >::Type i = 0; i < _getMaximalSLTowerHeight( numEntries ); ++i )
		{
			new( buffer ) SkipElement< TObject, TModus, TSpec, TStructuring >( _rightBorder, list._baseStore, left_border_key);
			++buffer;
		}
		
		_setCount( *first_base, numEntries );

	}
	

} // namespace seqan
#endif //SKIP_LIST_BASE_H_
