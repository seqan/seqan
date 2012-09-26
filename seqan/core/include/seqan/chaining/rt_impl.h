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

#ifndef SEQAN_HEADER_COMPLETE_RANGE_TREE_H
#define SEQAN_HEADER_COMPLETE_RANGE_TREE_H


namespace seqan
{


/*DISABLED
.Class.RangeTree:
..cat:Range Tree
..summary:The RangeTree is a data structure to solve the orthogonal range searching problem.
..signature:RangeTree< TObject, [ TModus, TSpec, TStructuring] >
..param.TObject:Type of stored objects.
..param.TModus:Modus of operation of a RangeTree. A RangeTree is static.
..param.TSpec:Specialization of the RangeTree.
..param.TStructuring:Parameter to specify whether the RangeTree uses Deferred Data Structuring or not.
..remarks:The object given to the RangeTree should offer the following functions:
..remarks:$key( obj, dim )$: returns the key of the object for dimension $dim$.
..remarks:$setKey( obj, dim, k )$: set the key of the object to $k$ for dimension $dim$.
..remarks:In contrast to STL-like containers, the objects are not cloned by the RangeTree. It only supports searching operations on a set of objects. This set must be handled by the user.
..remarks:The $MaxTree$ specialization offers the abbility to perform Range Maximum Queries.
..include:seqan/chaining.h
*/

		
template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
class RangeTree< TObject, TModus, RT< TSpec >, TStructuring >
{
public:
	SkipList< TObject, TModus, RT< TSpec >, TStructuring > * _list;

	TObject RBorderObj_;
	TObject _LBorderObj;

	typename Size< TObject >::Type _dim;
	typename Size< SkipList< TObject, TModus, RT< TSpec >, TStructuring > >::Type _numOfElems;

	SearchPath_< TObject, TModus, RT< TSpec >, TStructuring > _sp;
	RangeTreeAllocators< TObject, TModus, RT< TSpec >, TStructuring > _allocs;
	
	friend inline
	typename Size< SkipList< TObject, TModus, RT< TSpec >, TStructuring > >::Type 
	length( RangeTree< TObject, TModus, RT< TSpec >, TStructuring > & me )
	{
		SEQAN_CHECKPOINT
		return me._numOfElems;
	}
/*
	friend
	Allocator< ClassPool< SkipElement< TObject, TModus, RT< TSpec >, TStructuring >, Limited > > &
	_getElementAlloc<>( RangeTree< TObject, TModus, RT< TSpec >, TStructuring > & me );

	friend
	Allocator< ClassPool< SkipList< TObject, TModus, RT< TSpec >, TStructuring >, Unlimited > > & 
	_getListAlloc<>( RangeTree< TObject, TModus, RT< TSpec >, TStructuring > & me );

	friend
	Allocator< SimpleAlloc<> > &
	_getBaseAlloc<>( RangeTree< TObject, TModus, RT< TSpec >, TStructuring > & me );

	friend
	TObject * 
	_getRBorderObj<>( RangeTree< TObject, TModus, RT< TSpec >, TStructuring > & me );

	friend
	TObject * 
	_getLBorderObj<>( RangeTree< TObject, TModus, RT< TSpec >, TStructuring > & me );

	friend
	SkipList< TObject, TModus, RT< TSpec >, TStructuring > *
	_getList<>( RangeTree< TObject, TModus, RT< TSpec >, TStructuring > & me );
*/
	template< typename TSize > friend inline
	SkipElement< TObject, TModus, RT< TSpec >, TStructuring > **
	_getSearchPath( RangeTree< TObject, TModus, RT< TSpec >, TStructuring > & me,
					TSize dim )
	{
		SEQAN_CHECKPOINT
		return me._sp._searchPath + _getMaximalSLTowerHeight( me ) * dim;
	}
	
	friend inline
	typename Size< TObject >::Type
	dimension( RangeTree< TObject, TModus, RT< TSpec >, TStructuring > & me )
	{
		SEQAN_CHECKPOINT
		return me._dim;
	}

	RangeTree()
	{}

public:
	
	template< typename TContainer, typename TSize >
	RangeTree(	TContainer & data,
				TSize dim )
	: _dim( dim )
	, _numOfElems( length( data ) )
	, _sp( _getMaximalSLTowerHeight( length( data ) ) * _dim )
	, _allocs( 2 * length( data ) * dim * dim * dim, length( data ) * dim * dim )
	{
		SEQAN_CHECKPOINT
		_setMinInfty( _LBorderObj, _dim );
		_setMaxInfty( RBorderObj_, _dim );

		allocate( _getListAlloc( *this ), _list, 1 );
		valueConstruct( _list );
		_create( *_list, data, *this, _dim - 1 );
	}


	~RangeTree()
	{	
		SEQAN_CHECKPOINT
		_clearSearchPath( this->_sp, _getMaximalSLTowerHeight( length( *this ) ) * _dim );
	}

}; // struct RangeTree

}

#endif // SEQAN_HEADER_COMPLETE_RANGE_TREE_H




