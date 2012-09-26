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

#ifndef SEQAN_HEADER_RT_SL_COMMON_ALGOS_H
#define SEQAN_HEADER_RT_SL_COMMON_ALGOS_H

namespace seqan{

/////////////////////////////////////////////////////////////////////////////////////////
//
//	structure building algorithms
//
/////////////////////////////////////////////////////////////////////////////////////////


	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TSize > inline
	void
	_completeBuild( SkipList< TObject, TModus, RT< TSpec >, TStructuring > & list,
					TSize dim )
	{
		_sortRecursive( list, _getBaseStore( list ), _getBaseStore( list ) + length( list ) + 1, dim );
		_setHeight( *_getBaseStore( list ), 1 );
		_buildTowers( list, dim );
	}

	template< typename TObject, typename TModus, typename TSpec, typename TSize > inline
	void
	_completeBuild( SkipList< TObject, TModus, RT< TSpec >, Deferred > & /*list*/,
					TSize /*dim*/ )
	{
	} 

/////////////////////////////////////////////////////////////////////////////////////////
//
//	searching algorithms
//
/////////////////////////////////////////////////////////////////////////////////////////

	// _findTowerTop searches for the highest layer of a tower, 
		// whose right-pointer points to an element beyond the given border
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TResultSet, typename TSize, typename TKey > inline
	SkipElement< TObject, TModus, RT< TSpec >, TStructuring > *
	_findTowerTop(	SkipBaseElement< TObject, TModus, RT< TSpec >, TStructuring > * base_buffer,
					SkipList< TObject, TModus, RT< TSpec >, TStructuring > * list,
					TObject * left_border,
					TObject * right_border,
					TSize dim,
					TKey theKey,
					TResultSet & results )
	{
		SkipElement< TObject, TModus, RT< TSpec >, TStructuring > * tower_buffer = &_getUp( *base_buffer ) + _getHeight( * base_buffer ) - 1;
		SkipElement< TObject, TModus, RT< TSpec >, TStructuring > * right_buffer = _getRight( *tower_buffer );
		while( key( *right_buffer ) <= theKey )
		{
			_scanAssocStruct( tower_buffer, list, left_border, right_border, dim, results );
			base_buffer = _getDown( *right_buffer );
			tower_buffer = &_getUp( *base_buffer ) + _getHeight( * base_buffer ) - 1;
			right_buffer = _getRight( *tower_buffer );
		}
		return tower_buffer;	
	}


	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TSize, typename TKey, typename TResultSet > inline
	SkipBaseElement< TObject, TModus, RT< TSpec >, TStructuring > *
	_collectAssocStructs(	SkipElement< TObject, TModus, RT< TSpec >, TStructuring > * tower_buffer,
							TObject * left_border,
							TObject * right_border,
							TSize dim,
							TKey search_theKey,
							SkipList< TObject, TModus, RT< TSpec >, TStructuring > * list,
							TResultSet & results )
	{
		typename Size< SkipList< TObject, TModus, RT< TSpec >, TStructuring > >::Type height = tower_buffer - &_getUp( * _getDown( * tower_buffer ) ) + 1;
		SkipElement< TObject, TModus, RT< TSpec >, TStructuring > * right_buffer = _getRight( *tower_buffer );
		while( height > 0 )
		{
			while( key( *right_buffer ) <= search_theKey )
			{
				_scanAssocStruct( tower_buffer, list, left_border, right_border, dim, results );
				tower_buffer = right_buffer;
				right( right_buffer );
			}
			--height;
			--tower_buffer;
			right_buffer = _getRight( *tower_buffer );
		}
		++tower_buffer;
		return _getDown( *tower_buffer );
	}

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TSize, typename TKey, typename TResultSet > inline
	SkipBaseElement< TObject, TModus, RT< TSpec >, TStructuring > *
	_checkBaseElementsLeft( SkipBaseElement< TObject, TModus, RT< TSpec >, TStructuring > * base,
							TObject * left_border,
							TObject * right_border,
							TSize dim,
							TKey left_theKey,
							TKey right_theKey,
							TResultSet & results )
	{
			// check elemnts in base layer, until
			//		1 ) element with tower is reached
			//		2 ) sorted element is reached, which has a larger key
			//			( all elements on right side of a sorted elemnt have a larger theKey )
		
		while( key( *base ) < left_theKey )
			goNext( base );

		typename Size< TObject >::Type l_dim = dim - 1;
		while( _getHeight( *base ) == 0 )
		{
			if( key( *base, dim ) <= right_theKey ){
				if( _testRange( *getObject( base ), *left_border, *right_border, l_dim ) )
					_pushBack( results, getObject( base ) );
			}
			else if( _getCount( *base ) != 0 ){
				break;
			}
			goNext( base );
		}
		return base;
	}


	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TKey, typename TResultSet > inline
	SkipBaseElement< TObject, TModus, RT< TSpec >, TStructuring > *
	_checkBaseElementsLeftBottom( SkipBaseElement< TObject, TModus, RT< TSpec >, TStructuring > * base,
									TObject * /*left_border*/,
									TObject * /*right_border*/,
									TKey left_theKey,
									TKey right_theKey,
									TResultSet & results )
	{
			// check elemnts in base layer, until
			//		1 ) element with tower is reached
			//		2 ) sorted element is reached, which has a larger key
			//			( all elements on right side of a sorted elemnt have a larger theKey )
		
		while( key( *base ) < left_theKey )
			goNext( base );

		while( _getHeight( *base ) == 0 )
		{
			if( key( *base ) <= right_theKey ){
				_pushBack( results, getObject( base ) );
			}
			else if( _getCount( *base ) != 0 ){
				break;
			}
			goNext( base );
		}
		return base;
	}

	

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TSize, typename TKey, typename TResultSet > inline
	void
	_checkBaseElementsRight(	SkipBaseElement< TObject, TModus, RT< TSpec >, TStructuring > * base,
								TObject * left_border,
								TObject * right_border,
								TSize dim,
								TKey right_theKey,
								TResultSet & results )
	{
		TSize l_dim = dim-1;
		while( key( *base ) <= right_theKey )
		{
			if( _testRange( *getObject( base ), *left_border, *right_border, l_dim ) )
			{
				_pushBack( results, getObject( base ) );
			}
			goNext( base );
		}
	}

}

#endif

