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

#ifndef SEQAN_HEADER_RMT_SL_DEF_ALGOS_H
#define SEQAN_HEADER_RMT_SL_DEF_ALGOS_H

namespace seqan{


/////////////////////////////////////////////////////////////////////////////////////////
//
//	structure building algorithms
//
/////////////////////////////////////////////////////////////////////////////////////////

		// construct the towers for the complete skip list
		// -> the deferred RMT comletely builds the list of the lowest layer
	template< typename TObject, typename TSpec > inline
	void
	_completeBuild( SkipList< TObject, SkipListStatic, RT< MaxTree< TSpec > >, Deferred > & list,
					typename Size< TObject >::Type dim )
	{
		if( dim == 0 )
		{
			_sortRecursive( list, _getBaseStore( list ), _getBaseStore( list ) + length( list ) + 1, dim );
			_setHeight( *_getBaseStore( list ), 1 );	
			_buildMaxTowers( list );
		}
	}
		
/////////////////////////////////////////////////////////////////////////////////////////
//
//	range maximum tree algorithms
//
/////////////////////////////////////////////////////////////////////////////////////////

		// _performRMQ adaption for the deferred RMT
		// The defferd RMT needs 2 searches per layer:
		// The first one assures that the element with the appropriate key is in an ordered state,
		// the second one is equivalent to _performRMQ of the common case
	template< typename TObject, typename TSpec, typename TBorder, typename TSize, typename TKey >
	void
	_performRMQ(	SkipList< TObject, SkipListStatic, RT< MaxTree< TSpec > >, Deferred > * list,
					TBorder & borderObj,
					TSize dim,
					TKey searchKey,
					TObject *& maxObject )
	{
		SkipElement< TObject, SkipListStatic, RT< MaxTree< TSpec > >, Deferred > * layer_element = _getRoot( *list );
		typename Size< SkipList< TObject, SkipListStatic, RT< MaxTree< TSpec > >, Deferred > >::Type height = _getCurrentLayer( *list );
		SkipElement< TObject, SkipListStatic, RT< MaxTree< TSpec > >, Deferred > ** search_path = _getSearchPath( *_getMainTree( *list ), 0 ) + height - 1;
				// first search operation
		_searchFrom( *list, layer_element, searchKey, search_path, dim );
				// search again
		height = _getCurrentLayer( *list );
		layer_element = _getRoot( *list );
		while( height > 0 )
		{
			while( key( *_getRight( *layer_element ) ) < searchKey )
			{ 
				if( !_hasAssocStruct( layer_element ) )
					_createAssocStruct( layer_element, list, dim );
				if( _hasSmallAssocStruct( layer_element ) )
					_performSmallRMQ( layer_element, dim - 1, borderObj, maxObject );
				else
					_processRMQ( _getAssocStruct( layer_element ), dim - 1, borderObj, maxObject );
				layer_element = _getRight( *layer_element );
			}
			--layer_element;
			--height;
		}
		++layer_element;
		SkipBaseElement< TObject, SkipListStatic, RT< MaxTree< TSpec > >, Deferred > * base_element = _getDown( *layer_element );
		SkipBaseElement< TObject, SkipListStatic, RT< MaxTree< TSpec > >, Deferred > * base_buffer = base_element;
		while( key( *base_buffer, dim ) < searchKey )
		{
			base_buffer = _getRight( *base_buffer );
		}
			// in the lowest layer, searching to the right
		while( base_element != base_buffer )
		{
			_testRangeMax( getObject( base_element), &borderObj, maxObject, dim - 1 );
			goNext( base_element );
		}
	}

	
}

#endif // SEQAN_HEADER_RMT_SL_DEF_ALGOS_H
