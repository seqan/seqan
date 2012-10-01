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


//SEQAN_NO_DDDOC: do not generate documentation for this file

#ifndef SEQAN_RT_SKIP_LIST_BASE_H
#define SEQAN_RT_SKIP_LIST_BASE_H

#include "skip_list.h"
#include "rt_skip_element.h"
#include "rt_skip_base_element.h"

namespace seqan{


/////////////////////////////////////////////////////////////////////////////////////////
//
//	basic accessor functions
//
/////////////////////////////////////////////////////////////////////////////////////////

		// get thr right border element
		// e.g. the skip base element with supremum key
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
	SkipElement< TObject, TModus, RT< TSpec >, TStructuring > *
	_getRightBorder( SkipList< TObject, TModus, RT< TSpec >, TStructuring > & me );
		
	
		// specialization for static case
	template< typename TObject, typename TSpec, typename TStructuring > inline
	SkipElement< TObject, SkipListStatic, RT< TSpec >, TStructuring > * 
	_getRightBorder( SkipList< TObject, SkipListStatic, RT< TSpec >, TStructuring > & me )
	{
		SEQAN_CHECKPOINT
		return &_getUp( me._baseStore[ me._numOfElements + 1 ] );
	}

		// specialization for dynamic case
	//template< typename TObject, typename TSpec, typename TStructuring > inline
	//SkipElement< TObject, SkipListDynamic, RT< TSpec >, TStructuring > * 
	//_getRightBorder( SkipList< TObject, SkipListDynamic, RT< TSpec >, TStructuring > & me )
	//{
	//	SEQAN_CHECKPOINT
	//	SkipElement< TObject, SkipListDynamic, RT< TSpec >, TStructuring > * border = _getRoot( me );
	//	while( key( *border ) < maxValue< typename Key< TObject >::Type >() )
	//		border = _getRight( *border );
	//	return border;
	//}

		// get the element allocator of a skip list
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
	Allocator< ClassPool< SkipElement< TObject, TModus, RT< TSpec >, TStructuring >, Limited, SimpleAllocator > > &
	_getElementAlloc( SkipList< TObject, TModus, RT< TSpec >, TStructuring > & me )
	{
		SEQAN_CHECKPOINT
		return _getElementAlloc( *_getMainTree( me ) );
	}

		// get the list allocator of a skip list
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
	Allocator< ClassPool< SkipList< TObject, TModus, RT< TSpec >, TStructuring >, Unlimited, SimpleAllocator > > &
	_getListAlloc( SkipList< TObject, TModus, RT< TSpec >, TStructuring > & me )
	{
		SEQAN_CHECKPOINT
		return _getListAlloc( *_getMainTree( me ) );
	}

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
	Allocator< SimpleAlloc<> > &
	_getBaseAlloc( SkipList< TObject, TModus, RT< TSpec >, TStructuring > & me )
	{
		SEQAN_CHECKPOINT
		return _getBaseAlloc( *_getMainTree( me ) );
	}

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
	RangeTree< TObject, TModus, RT< TSpec >, TStructuring > *
	_getMainTree( SkipList< TObject, TModus, RT< TSpec >, TStructuring > & me );

	
///////////////////////////////////////////////////////
//
//	initialization
//
///////////////////////////////////////////////////////

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TContainer, typename TSize >
	void
	_create( SkipList< TObject, TModus, RT< TSpec >, TStructuring > & me,
				TContainer & data,
				RangeTree< TObject,	TModus, RT< TSpec >, TStructuring > & Tree,
				TSize dim )
	{
		SEQAN_CHECKPOINT
		me._mainTree = &Tree;
		me._dim = dim;
		me._numOfElements = length( data );

		typename Iterator< TContainer >::Type first = begin( data );
		typename Iterator< TContainer >::Type last = end( data );

		SkipBaseElement< TObject, TModus, RT< TSpec >, TStructuring > * base_right;

		_initBases( me, me._baseStore, base_right, first, last, me._numOfElements, dim ); 

		_initSL( me, me._baseStore, base_right, me._numOfElements, me._dim );

		_completeBuild( me, dim );
	}

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TSize >
	void
	_create( SkipList< TObject, TModus, RT< TSpec >, TStructuring > & me,
				SkipBaseElement< TObject, TModus, RT< TSpec >, TStructuring > * first,
				SkipBaseElement< TObject, TModus, RT< TSpec >, TStructuring > * last,
				RangeTree< TObject,	TModus, RT< TSpec >, TStructuring > & Tree,
				TSize dim )
	{
		SEQAN_CHECKPOINT
		me._mainTree = &Tree;
		me._dim = dim;
		me._numOfElements = last - first;

		SkipBaseElement< TObject, TModus, RT< TSpec >, TStructuring > * base_right;

		_initBases( me, me._baseStore, base_right, first, last, me._numOfElements, dim ); 

		_initSL( me, me._baseStore, base_right, me._numOfElements, dim );

		_completeBuild( me, dim );
	}


	template< typename TObject, typename TSpec, typename TStructuring, typename TIter, typename TSize > inline
	void
	_initBases( SkipList< TObject, SkipListStatic, RT< TSpec >, TStructuring > & list, 
				SkipBaseElement< TObject, SkipListStatic, RT< TSpec >, TStructuring > * firstBase,
				SkipBaseElement< TObject, SkipListStatic, RT< TSpec >, TStructuring > *& lastBase,
				TIter & firstData,
				TIter & lastData,
				TSize numEntries,
				TSize dim )
	{
		SEQAN_CHECKPOINT
		allocate( _getBaseAlloc( list ), firstBase, numEntries + 2 );
		list._baseStore = firstBase;

		valueConstruct( firstBase, _getLBorderObj( *_getMainTree( list ) ), minValue< typename Key< TObject >::Type >() );
		
		++firstBase;

		while( firstData != lastData )
		{
			valueConstruct( firstBase, &value( firstData ), key( value( firstData ), dim ) );
			++firstBase;
			++firstData;
		}
		lastBase = firstBase;
		firstBase = list._baseStore;
		valueConstruct( lastBase, _getRBorderObj( *_getMainTree( list ) ), maxValue< typename Key< TObject >::Type >() );
	}


	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TSize > inline
	void
	_initSL(	SkipList< TObject, TModus, RT< TSpec >, TStructuring > & list,
				SkipBaseElement< TObject, TModus, RT< TSpec >, TStructuring > * first_base,
				SkipBaseElement< TObject, TModus, RT< TSpec >, TStructuring > * /*last_base*/,
				TSize numEntries,
				TSize /*dim*/ )
	{
		SEQAN_CHECKPOINT
			// allocate space for bording elements
		SkipElement< TObject, TModus, RT< TSpec >, TStructuring  > * _rightBorder;
		allocate( _getElementAlloc( list ), _rightBorder, 2 );
		arrayConstruct( _rightBorder, _rightBorder + 2 );
		allocate( _getElementAlloc( list ), list._leftSideStore, _getMaximalSLTowerHeight( numEntries ) );
				
			// set values
			// left side ...
		_setHeight( *first_base, 1 );
		_setUp( *first_base, *list._leftSideStore );

		typename Key< TObject >::Type left_border_key = minValue< typename Key< TObject >::Type >();
		SkipElement< TObject, TModus, RT< TSpec >, TStructuring > * buffer = list._leftSideStore;
		for( typename Size< SkipList< TObject, TModus, RT< TSpec >, TStructuring > >::Type i = 0; i < _getMaximalSLTowerHeight( numEntries ); ++i )
		{
			valueConstruct( buffer, _rightBorder, first_base, left_border_key );
			++buffer;
		}
			// ... right side
		SkipBaseElement< TObject, TModus, RT< TSpec >, TStructuring > * r_base = &list._baseStore[ numEntries + 1 ];
		_setDown( *_rightBorder, r_base);
		_setHeight( *r_base, 1 );
		_setRight( * _rightBorder, _rightBorder );
		setKey( *_rightBorder, maxValue< typename Key< TObject >::Type >() );
		_setUp( *r_base, *_rightBorder );

		_setCount( *first_base, numEntries );
	}

/////////////////////////////////////////////////////////////////////////////////////////
//
//	structure building algorithms
//
/////////////////////////////////////////////////////////////////////////////////////////

	template< typename TObject, typename TModus, typename TSpec, typename TSize > inline
	void 
	_connectUpdate(	SkipList< TObject, TModus, RT< TSpec >, Deferred > & /*list*/,
						SkipBaseElement< TObject, TModus, RT< TSpec >, Deferred > * base,
						TSize height,
						SkipElement< TObject, TModus, RT< TSpec >, Deferred > ** search_path,
						TSize dim )
	{
		SEQAN_CHECKPOINT
		SkipElement< TObject, TModus, RT< TSpec >, Deferred > * buffer = &_getUp( *base );
		SkipElement< TObject, TModus, RT< TSpec >, Deferred > * top = buffer + height;
		typename Key< TObject >::Type theKey = key( *base, dim );
		while( buffer != top )
		{
			new( buffer ) SkipElement< TObject, TModus, RT< TSpec >, Deferred >( _getRight( **search_path ), base, theKey );
			_setRight( **search_path, buffer );
			_deleteAssocStruct( *search_path );
			*search_path = buffer;
			++buffer;
			++search_path;
		}
	}

	// reconnect the pointers after building a tower
	template< typename TObject, typename TModus, typename TSpec, typename TSize > inline
	void 
	_connect(	SkipList< TObject, TModus, RT< TSpec >, Deferred > & /*list*/,
				SkipBaseElement< TObject, TModus, RT< TSpec >, Deferred > * base,
				TSize height,
				SkipElement< TObject, TModus, RT< TSpec >, Deferred > ** search_path,
				TSize dim )
	{
		SEQAN_CHECKPOINT
		SkipElement< TObject, TModus, RT< TSpec >, Deferred > * buffer = &_getUp( *base );
		SkipElement< TObject, TModus, RT< TSpec >, Deferred > * top = buffer + height;
		typename Key< TObject >::Type theKey = key( *base, dim );
		while( buffer != top )
		{
			new( buffer ) SkipElement< TObject, TModus, RT< TSpec >, Deferred >( _getRight( **search_path ), base, theKey );
			_setRight( **search_path, buffer );
			_deleteAssocStruct( *search_path );
			++buffer;
			++search_path;
		}
	}

	template< typename TObject, typename TModus, typename TSpec, typename TSize > inline
	void 
	_connectUpdate(	SkipList< TObject, TModus, RT< TSpec >, SemiDeferred > & /*list*/,
						SkipBaseElement< TObject, TModus, RT< TSpec >, SemiDeferred > * base,
						TSize height,
						SkipElement< TObject, TModus, RT< TSpec >, SemiDeferred > ** search_path,
						TSize dim )
	{
		SEQAN_CHECKPOINT
		SkipElement< TObject, TModus, RT< TSpec >, SemiDeferred > * buffer = &_getUp( *base );
		SkipElement< TObject, TModus, RT< TSpec >, SemiDeferred > * top = buffer + height;
		typename Key< TObject >::Type theKey = key( *base, dim );
		while( buffer != top )
		{
			new( buffer ) SkipElement< TObject, TModus, RT< TSpec >, SemiDeferred >( _getRight( **search_path ), base, theKey );
			_setRight( **search_path, buffer );
			*search_path = buffer;
			++buffer;
			++search_path;
		}
	}


		// reconnect the pointers after building special for
		// the complete case
	template< typename TObject, typename TModus, typename TSpec, typename TSize > inline
	void 
	_connectUpdate(	SkipList< TObject, TModus, RT< TSpec >, Complete > & list,
						SkipBaseElement< TObject, TModus, RT< TSpec >, Complete > * base,
						TSize height,
						SkipElement< TObject, TModus, RT< TSpec >, Complete > ** search_path,
						TSize dim )
	{
		SEQAN_CHECKPOINT
		SkipElement< TObject, TModus, RT< TSpec >, Complete > * buffer = &_getUp( *base );
		typename Key< TObject >::Type theKey = key( *base, dim );
		typename Size< SkipList< TObject, TModus, RT< TSpec >, Complete > >::Type current_height = 0;
		while( current_height < height - 1 ){
			new( buffer ) SkipElement< TObject, TModus, RT< TSpec >, Complete >( _getRight( **search_path ), base, theKey );
			_setRight( **search_path, buffer );
			if( _getHeight( **search_path ) == current_height + 1)
				_createAssocStruct( *search_path, &list, dim );
			*search_path = buffer;
			++current_height;
			++buffer;
			++search_path;
		}
		new( buffer ) SkipElement< TObject, TModus, RT< TSpec >, Complete >( _getRight( **search_path ), base, theKey );
		_setRight( **search_path, buffer );	
		_createAssocStruct( *search_path, &list, dim );
		*search_path = buffer;
	}


	// linker rand muss getestet werden -> darf nicht mit eigefügt werden
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TSize > inline
	void
	_buildAssocStruct_left(  SkipList< TObject, TModus, RT< TSpec >, TStructuring >  * me,
							 SkipElement< TObject, TModus, RT< TSpec >, TStructuring > * left,
							 SkipBaseElement< TObject, TModus, RT< TSpec >, TStructuring > * right,
							 TSize dim )
	{
		SEQAN_CHECKPOINT
		if( _getDown( *left ) != _getBaseStore( *me ) )
			_createAssocStruct( left, me, _getDown( *left ), right, dim );
		else
			_createAssocStruct( left, me, _getSucc( *_getDown( *left ) ), _getDown( *_getRight( *_getRight( *left ) ) ), dim );
	}


	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TSize > inline
	void
	_buildAssocStruct_right(  SkipList< TObject, TModus, RT< TSpec >, TStructuring >  * me,
							  SkipElement< TObject, TModus, RT< TSpec >, TStructuring > * left,
							  SkipBaseElement< TObject, TModus, RT< TSpec >, TStructuring > * right,
							  TSize dim )
	{
		SEQAN_CHECKPOINT
		_createAssocStruct( left, me, _getDown( *left ), right, dim );
	}


}

#endif
