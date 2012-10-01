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

#ifndef SEQAN_HEADER_RMT_SL_COMMON_ALGOS_H
#define SEQAN_HEADER_RT_SL_COMMON_ALGOS_H

/*
 *  rmt_common_algos.h
 *  rmt
 *
 *  Created by Hendrik Woehrle
 *
 *	Contains functions for all specs of the RMT
 *
 */

namespace seqan{

/////////////////////////////////////////////////////////////////////////////////////////
//
//	tree algorithms
//
/////////////////////////////////////////////////////////////////////////////////////////


/*DISABLED
.Function.rangeMaxQuery:
..summary:Get the object with maximum priority in the RMT in a given intervall
..cat:Range Tree
..signature:rangeMaxQuery(tree, border)
..param.tree:A Range Tree with spec MaxTree.
...type:RangeMaximumTree
..param.border:The object that stores the borders for all dimensions.
..returns:A pointer to the object witch maximal priority in the given intervall.
..include:seqan/chaining.h
*/

	template< typename TObject, typename TSpec, typename TStructuring, typename TBorder >
	TObject *
	rangeMaxQuery( RangeTree< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > & tree,
					TBorder & border_obj )
	{
		SEQAN_ASSERT_GEQ(dimension(border_obj), dimension(tree));
		TObject * maxObject = _getLBorderObj( tree );
		_processRMQ( _getList( tree ), dimension( tree ) - 1, border_obj, maxObject );
		return maxObject;
	}

		// wrappper function to distinguis between higher dimensions or lower dimensions
		// in a RMQ
	template< typename TObject, typename TSpec, typename TStructuring, typename TBorder > inline
	void
	_processRMQ(	SkipList< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > * list,
					typename Size< TObject >::Type dim,
					TBorder & borderObj,
					TObject *& maxObject )
	{
		typename Key< TObject >::Type searchKey = key( borderObj, dim );
		if( dim > 0 )
			_performRMQ( list, borderObj, dim, searchKey, maxObject );
		else
			_performRMQ( list, borderObj, searchKey, maxObject );
	}

/*DISABLED
.Function.activate:
..summary:Update the internal pointer structure of the RMT. 
..cat:Range Tree
..signature:activate( tree, obj )
..signature:activate( tree, obj, prio )
..param.tree:The tree.
...type:Class.RangeMaximumTree
..param.obj:The object.
..param.prio:The new priority (optional). $prio > priority( obj )$ must hold.
...type:$Metafunction.Weight< TObject >::Type$.
..include:seqan/chaining.h
*/
	template< typename TObject, typename TSpec, typename TStructuring > inline
	void
	activate( RangeTree< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > & tree,
			   TObject & obj )
	{
		_activate( _getList( tree ), &obj, dimension( tree ) - 1 );
	}

	template< typename TObject, typename TSpec, typename TStructuring, typename TWeight > inline
	void
	activate( RangeTree< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > & tree,
			   TObject & obj,
			   TWeight prio )
	{
		SEQAN_ASSERT_GEQ(prio, priority( obj ) );
		setPriority( obj, prio );
		_activate( _getList( tree ), &obj, dimension( tree ) - 1 );
	}


		// activate for higher layers
		// -> searches lower layers that contain the object whose priority should be increased
	template< typename TObject, typename TSpec, typename TStructuring, typename TSize >
	void
	_activateHigherLayer(	SkipList< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > * list,
							TObject * obj,
							TSize dim )
	{
		typename Size< SkipList< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > >::Type height = _getCurrentLayer( *list );
		SkipElement< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > * layer_element = _getRoot( *list );
		SkipElement< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > * right_buffer = _getRight( *layer_element );

		typename Key< TObject >::Type searchKey = key( *obj, dim );
		while( height > 0 )
		{
			while( key( *right_buffer ) <= searchKey )
			{	
				layer_element = right_buffer;
				right_buffer = _getRight( *right_buffer );
			}
			if( _hasBigAssocStruct( layer_element ) )
				_activate( _getAssocStruct( layer_element ), obj, dim-1 );
			--height;
			--layer_element;
			right_buffer = _getRight( *layer_element );
		}

	}


		// the activation of the object in the lowest layer,
		// updates the max-pointers
	template< typename TObject, typename TSpec, typename TStructuring, typename TSize >
	void
	_activate(  SkipList< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > * list,
				TObject * obj,
				TSize dim )
	{
		if( dim > 0 ){
			_activateHigherLayer( list, obj, dim );
		}
		else {
			typename Size< SkipList< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > >::Type height = _getCurrentLayer( *list );
			SkipElement< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > ** sp = _getSearchPath( *_getMainTree( *list ), dim ) + height - 1;
			SkipElement< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > * layer_element = _getRoot( *list );
			typename Key< TObject >::Type search_key = key( *obj, dim );
			while( height > 0 )
			{
				while( key( *_getRight( *layer_element ) ) <= search_key )
				{	
					layer_element = _getRight( *layer_element );		
				}
				--height;
				*sp = layer_element;
				--layer_element;
				--sp;
			}
			++sp;
			++layer_element;
			++height;
			typename Weight< TObject >::Type new_score = priority( *obj );
			typename Size< SkipList< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > >::Type max_height = _getCurrentLayer( *list ) + 1;
			while( height < max_height ){
				if( priority( *sp ) <= new_score )
					_setMaxObject( *sp, obj );
				else
					break;
				++height;
				++sp;
			}
		}
	}

/////////////////////////////////////////////////////////////////////////////////////////
//
//	structure building algorithms
//
/////////////////////////////////////////////////////////////////////////////////////////
	

/*DISABLED
.Internal._connectUpdateMax:
..summary:As _connectUpdate, additionaly updates the max-pointers in the lowest layer of the RMT.
..cat:RangeMaximumTree
..signature:_connectUpdateMax( list, base, height, key, search_path, max_obj )
..param.list:The Skip List of the RMT representing this layer.
...type:Class.SkipList
..param.base:The base element.
...type:Class.SkipBaseElement*
..param.height:The height of the tower.
..param.key:The key of base.
..param.searchPath:The search path.
...type:SkipElement**
..param.max_obj:The object with maximal priority.
*/
	
	template< typename TObject, typename TSpec, typename TStructuring, typename TSize, typename TKey >
	void
	_connectUpdateMax(	SkipList< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > & list,
							SkipBaseElement< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > * base, 
							TSize height,
							TKey searchKey,
							SkipElement< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > ** search_path,
							TObject * max_obj )
	{
		SkipElement< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > * buffer = &_getUp( *base );
		SkipElement< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > * tower_top = buffer + height;
		
		while( buffer != tower_top ){
			new( buffer ) SkipElement< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring >( _getRight( **search_path ), base, searchKey );
			_setRight( **search_path, buffer );
			_setMaxObject( buffer, max_obj );
			*search_path = buffer;
			++buffer;
			++search_path;
		}
		typename Weight< TObject >::Type score = priority( *max_obj );
		typename Size< TObject >::Type max_height = _getMaximalSLTowerHeight( list );
		while( height < max_height )
		{
			SEQAN_ASSERT(max_obj != NULL);
			if( priority( *search_path ) < score )
				_setMaxObject( *search_path, max_obj );
			else break;
			++search_path;
			++height;
		}
	}
	
			// build towers with maximum pointers in the lowest layer of the rmt
	template< typename TObject, typename TSpec, typename TStructuring >
	void
	_buildMaxTowers( SkipList< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > & list )
	{
		SkipBaseElement< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > * buffer = _getBaseStore( list );
		SkipBaseElement< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > * end = buffer + length( list ) + 1;
		typename Size< SkipList< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > >::Type max_height = _getMaximalSLTowerHeight( list );
		SkipElement< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > ** search_path = new SkipElement< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring >*[ max_height ];
		SkipElement< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > * elem_buffer = &_getUp( *buffer );
		TObject * border_obj = getObject( buffer );
		
		typename Size< SkipList< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > >::Type height = 0;
		
		typename Key< TObject >::Type buffer_key = key( *buffer );

		for( typename Size< TObject >::Type i = 0; i < max_height; ++i ){
			*search_path = elem_buffer;
			_setMaxObject( elem_buffer, border_obj );
			++search_path;
			++elem_buffer;
		}
		search_path -= max_height;
		goNext( buffer );

		typename Key< TObject >::Type act_key = key( *buffer );

		while( buffer != end )
		{
			act_key = key( *buffer );
			if( act_key != buffer_key )
			{
				height = _throwCoin< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring >( list, max_height );
				if( height > 0 ){
					_addMax( list, buffer, height, search_path );
					_connectUpdateMax( list, buffer, height, act_key, search_path, getObject( buffer ) );
				}
				else
					_activateScoreBuild( buffer, &list, search_path );
				buffer_key = act_key;
			}
			else
				_activateScoreBuild( buffer, &list, search_path );
			goNext( buffer  );
		}
		delete[] search_path;
	}


		// _add adaption for the lowest layer in a rmt
		// adjusts the max pointers
	template< typename TObject, typename TSpec, typename TStructuring, typename THeight > inline
	void 
	_addMax(	SkipList< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring >  & list,
				SkipBaseElement< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > * base,
				THeight height,
				SkipElement< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > ** /*search_path*/ )
	{			
			// adding additional layers, if necessary 
		SEQAN_ASSERT_MSG( &_getUp( *base ) == NULL, "tried to build tower on bas element with tower" );
		if( height > _getCurrentLayer( list ) ){
			_setCurrentLayer( list, height );
		}
		SkipElement< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > * tower;
		allocate( _getElementAlloc( list ), tower, height );
		_setUp( *base, *tower );
		_setHeight( *base, height );
	}


}

#endif // SEQAN_HEADER_RMT_SL_COMMON_ALGOS_H

