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

#ifndef SEQAN_HEADER_RT_SL_COMPL_ALGOS_H
#define SEQAN_HEADER_RT_SL_COMPL_ALGOS_H

namespace seqan{

/////////////////////////////////////////////////////////////////////////////////////////
//
//	structure building algorithms
//
/////////////////////////////////////////////////////////////////////////////////////////
	
		// quick sort recursive step
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TSize > inline
	void
	_sortRecursive( SkipList< TObject, TModus, RT< TSpec >, TStructuring > & list,
					SkipBaseElement< TObject, TModus, RT< TSpec >, TStructuring > * elem,
					SkipBaseElement< TObject, TModus, RT< TSpec >, TStructuring > * right,
					TSize dim )
	{
		SEQAN_CHECKPOINT
		if( _getCount( *elem ) != 0 )
		{
			SkipBaseElement< TObject, TModus, RT< TSpec >, TStructuring > * pivot = _sort( list, elem, right, dim );
			_sortRecursive( list, elem, pivot, dim );
			_sortRecursive( list, pivot, right, dim );
		}	
	}
		// construct the towers for the complete skip list
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TSize >
	void
	_buildTowers( SkipList< TObject, TModus, RT< TSpec >, TStructuring > & list,
					TSize dim )
	{
		SEQAN_CHECKPOINT
		SkipBaseElement< TObject, TModus, RT< TSpec >, TStructuring > * buffer = _getBaseStore( list );
		SkipBaseElement< TObject, TModus, RT< TSpec >, TStructuring > * end = buffer + length( list ) + 1;
		SkipElement< TObject, TModus, RT< TSpec >, TStructuring > ** search_path = new SkipElement< TObject, TModus, RT< TSpec >, TStructuring >*[ _getMaximalSLTowerHeight( length( list ) ) ];
		*search_path = &_getUp( *buffer ); 
		typename Size< SkipList< TObject, TModus, RT< TSpec >, TStructuring  > >::Type height;
		typename Size< SkipList< TObject, TModus, RT< TSpec >, TStructuring  > >::Type max_height = _getMaximalSLTowerHeight( list );
		typename Key< TObject >::Type buffer_theKey = key( *_getBaseStore( list ), dim );
		goNext( buffer );

		while( buffer != end )
		{
			SEQAN_ASSERT_LEQ(buffer_theKey, key(*buffer, dim));
			if( key( *buffer, dim ) != buffer_theKey )
			{
				height = _throwCoin< TObject, TModus, RT< TSpec >, TStructuring  >( list, max_height );
				if( height > 0 ){
					_add( list, buffer, height, search_path );
					_connectUpdate( list, buffer, height, search_path, dim );
				}
				buffer_theKey = key( *buffer, dim );
			}
			goNext( buffer  );
		}
		
		delete[] search_path;
	}

	

/////////////////////////////////////////////////////////////////////////////////////////
//
//	search algorithms
//
/////////////////////////////////////////////////////////////////////////////////////////

	
	
	
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TResultSet, typename TSize >
	void
	_fingerSearch(	SkipList< TObject, TModus, RT< TSpec >, TStructuring > * list,
					TObject * left_border,
					TObject * right_border,
					TSize dim,
					TResultSet & results )
	{
		SEQAN_CHECKPOINT
		typename Key< TObject >::Type left_theKey = key( *left_border, dim );
		typename Key< TObject >::Type right_theKey = key( *right_border, dim );
			
			// search for the left base element
		SkipBaseElement< TObject, TModus, RT< TSpec >, TStructuring > * base = _searchFrom( *list, _getRoot( *list ), left_theKey, dim );

		base = _checkBaseElementsLeft( base, left_border, right_border, dim, left_theKey, right_theKey, results );
		if( key( *base, dim ) > right_theKey )
			return;
	
				//	1 ) searching for highest layer,
				//		on-line search in associated structures
		SkipElement< TObject, TModus, RT< TSpec >, TStructuring > * tower_buffer = _findTowerTop( base, list, left_border, right_border, dim, right_theKey, results );
				//	2 ) on-line search in associated structures of the higher layers
		base = _collectAssocStructs( tower_buffer, left_border, right_border, dim, right_theKey, list, results );
				//	3 ) check the remaining base elements
		_checkBaseElementsRight( base, left_border, right_border, dim, right_theKey, results );
	}

	
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TResultSet >
	void
	_bottomSearch(	SkipList< TObject, TModus, RT< TSpec >, TStructuring > * list,
					TObject * left_border,
					TObject * right_border,
					TResultSet & results )
	{
		SEQAN_CHECKPOINT
		typename Size< SkipList< TObject, TModus, RT< TSpec >, TStructuring > >::Type dim = 0;
		typename Key< TObject >::Type left_theKey = key( *left_border, dim );
		typename Key< TObject >::Type right_theKey = key( *right_border, dim );
			// search for the left base element
		SkipBaseElement< TObject, TModus, RT< TSpec >, TStructuring > * base = _searchFrom( *list, _getRoot( *list ), left_theKey, dim );
		
		base = _checkBaseElementsLeftBottom( base, left_border, right_border, left_theKey, right_theKey, results );
		if( key( *base, dim ) > right_theKey )
			return;

		while( key( *base ) <= right_theKey )
		{
			_pushBack( results, getObject( base ) );
			goNext( base );
		}
	}

}

#endif // SEQAN_HEADER_RT_SL_COMPL_ALGOS_H
