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

#ifndef SEQAN_HEADER_RT_SL_DEF_ALGOS_H
#define SEQAN_HEADER_RT_SL_DEF_ALGOS_H

namespace seqan{


	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TParam > inline
	SkipElement< TObject, TModus, RT< TSpec >, TStructuring > ** 
	_getSearchPath( SkipList< TObject, TModus, RT< TSpec >, TStructuring > & me,
					TParam dim )
	{
		SEQAN_CHECKPOINT
		return _getSearchPath( *_getMainTree( me ), dim );
	}



/////////////////////////////////////////////////////////////////////////////////////////
//
//	search algorithms
//
/////////////////////////////////////////////////////////////////////////////////////////


	template< typename TObject, typename TModus, typename TSpec, typename TResultSet, typename TSize >
	void
	_fingerSearch(	SkipList< TObject, TModus, RT< TSpec >, Deferred > * list,
					TObject * left_border,
					TObject * right_border,
					TSize dim,
					TResultSet & results )
	{
		SEQAN_CHECKPOINT
		typename Key< TObject >::Type left_theKey = key( *left_border, dim );
		typename Key< TObject >::Type right_theKey = key( *right_border, dim );
			// search for the left base element
		SkipElement< TObject, TModus, RT< TSpec >, Deferred > ** search_path = _getSearchPath( *list, dim );
		SkipBaseElement< TObject, TModus, RT< TSpec >, Deferred > * base = _searchFrom( *list, _getRoot( *list ), left_theKey, search_path, dim );
		
		base = _checkBaseElementsLeft( base, left_border, right_border, dim, left_theKey, right_theKey, results );
		if( key( *base, dim ) > right_theKey )
			return;

			// search for the right border
			//	the cooresponding element might be sorted or unsorted
			//	=> search for the corresponding element

				// I. search for the smallest element, which is outside the search intervall
				//	1 ) search for the element, of which the right pointer of the highest tower
				//		points to an element outside of the search intervall
		SkipElement< TObject, TModus, RT< TSpec >, Deferred > * tower_buffer = _findTowerTop( base, list, left_border, right_border, dim, right_theKey, results );					
		
		base = _getDown( *tower_buffer );

				//  2 ) search downwards from highest element, sort all equal keys behind found object
		_sortEquals( *list, _searchFrom( *list, tower_buffer, right_theKey, search_path, dim ), right_theKey );
		
				// II. the range is now limited by a sorted element
				//	1 ) searching for highest layer again( may have changed ),
				//		on-line search in associated structures
		tower_buffer = _findTowerTop( base, list, left_border, right_border, dim, right_theKey, results );

				//	2 ) on-line search in associated structures of the higher layers
		base = _collectAssocStructs( tower_buffer, left_border, right_border, dim, right_theKey, list, results );

				//	3 ) check the remaining base elements
		_checkBaseElementsRight( base, left_border, right_border, dim, right_theKey, results );
	}

	template< typename TObject, typename TModus, typename TSpec, typename TResultSet >
	void
	_bottomSearch(	SkipList< TObject, TModus, RT< TSpec >, Deferred > * list,
					TObject * left_border,
					TObject * right_border,
					TResultSet & results )
	{
		SEQAN_CHECKPOINT
		typename Size< TObject >::Type dim = 0;
		typename Key< TObject >::Type left_theKey = key( *left_border, dim );
		typename Key< TObject >::Type right_theKey = key( *right_border, dim );
			// search for the left base element
		SkipElement< TObject, TModus, RT< TSpec >, Deferred > ** search_path = _getSearchPath( *list, 0 );
		SkipBaseElement< TObject, TModus, RT< TSpec >, Deferred > * base = _searchFrom( *list, _getRoot( *list ), left_theKey, search_path, dim );

		while( key( *base ) < left_theKey )
			goNext( base );
		if( key( *base ) > right_theKey )
			return;

		SkipBaseElement< TObject, TModus, RT< TSpec >, Deferred > * right_base  = _searchFrom( *list, _getRoot( *list ), right_theKey, search_path, dim );
		_sortEquals( *list, right_base, right_theKey );

			//  2 ) search downwards from highest element
		while( base != right_base )
		{
			_pushBack( results, getObject( base ) );
			goNext( base );
		}
		
		while( key( *base ) <= right_theKey )
		{
			_pushBack( results, getObject( base ) );
			goNext( base );
		}
	}

}

#endif // SEQAN_HEADER_RT_SL_DEF_ALGOS_H
