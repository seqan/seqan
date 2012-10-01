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


#ifndef SEQAN_RT_BASE_H
#define SEQAN_RT_BASE_H



namespace seqan{

/////////////////////////////////////////////////////////////////////////////////////////
//
//	declarations
//
/////////////////////////////////////////////////////////////////////////////////////////
		
		// standard tag struct for range tree
	template< typename TSpec = Default >
	struct RT
	{};

		// tag structs for grade of deferredness of the range tree
	struct SemiDeferred
	{};

	struct Complete
	{};

	struct Deferred
	{};

		// main classes
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
	struct SkipList;
	
	template< typename TObject, typename TModus, typename TSpec = Default, typename TStructuring = Complete >
	class RangeTree;


	template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
	struct Cargo< SkipList< TObject, TModus, RT< TSpec >, TStructuring > >
	{
		typedef EmptyCargo_ Type;
	};

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
	struct Cargo< SkipBaseElement< TObject, TModus, RT< TSpec >, TStructuring > >
	{
		typedef typename Cargo< SkipList< TObject, TModus, RT< TSpec >, TStructuring > >::Type Type;
	};

/////////////////////////////////////////////////////////////////////////////////////////
//
//	dependent members
//
/////////////////////////////////////////////////////////////////////////////////////////

		// the memory allocators of the range tree
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
	struct
	RangeTreeAllocators;

	template< typename TObject, typename TSpec, typename TStructuring >
	struct
	RangeTreeAllocators< TObject, SkipListStatic, TSpec, TStructuring >
	{
		
		Allocator< ClassPool< SkipElement< TObject, SkipListStatic, TSpec, TStructuring >, Limited > > _elementAlloc;
		Allocator< ClassPool< SkipList< TObject, SkipListStatic, TSpec, TStructuring >, Unlimited > > _listAlloc;
		Allocator< SimpleAlloc<> > _baseAlloc;

		RangeTreeAllocators()
			: _elementAlloc( NULL )
		{}

		template< typename TSize >
		RangeTreeAllocators( TSize size1, TSize size2 )
			:_elementAlloc( size1 )
			,_listAlloc( size2 )
		{
			SEQAN_CHECKPOINT
		}

		~RangeTreeAllocators( )
		{
			SEQAN_CHECKPOINT
		}
	};



/////////////////////////////////////////////////////////////////////////////////////////
//
//	utilities
//
/////////////////////////////////////////////////////////////////////////////////////////
		

	const size_t RANGE_TREE_THRESH_ = 16;

		
	template< typename TObject, typename TSpec, typename TStructuring > inline
	bool
	_checkAssocThresh( SkipBaseElement< TObject, SkipListStatic, TSpec, TStructuring > * first,
						SkipBaseElement< TObject, SkipListStatic, TSpec, TStructuring > * second )
	{
		SEQAN_CHECKPOINT
		return ( ( second - first ) > RANGE_TREE_THRESH_ );
	}

	template< typename TTarget, typename TSource > inline
	void
	_pushBack( TTarget & target, TSource const & source )
	{
		SEQAN_CHECKPOINT
#ifdef RTTIMETEST
		volatile TSource temp = source;
#else
		appendValue( target, source );
#endif
	}


			// setting the element to - infinity
	template< typename TObject, typename TSize > inline 
	void 
	_setMinInfty(	TObject & me,  
					TSize dim )
	{
		SEQAN_CHECKPOINT
		me = TObject( dim );
		typename Key< TObject>::Type infValue = minValue< typename Key< TObject>::Type >();
		for( typename Size< TObject >::Type i = 0; i < dimension( me ); ++i )
		{
			setKey( me, i, infValue );
		}
	}

		// setting the element to + infinity
	template< typename TObject, typename TSize > inline 
	void 
	_setMaxInfty(	TObject & me,  
					TSize dim )
	{
		
		SEQAN_CHECKPOINT
		me = TObject( dim );
		typename Key< TObject>::Type supValue = maxValue< typename Key< TObject>::Type >();
		for( typename Size< TObject >::Type i = 0; i < dimension( me ); ++i )
		{
			setKey( me, i, supValue );
		}
	}

/////////////////////////////////////////////////////////////////////////////////////////
//
//	basic accessor functions
//
/////////////////////////////////////////////////////////////////////////////////////////
	

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline 
	Allocator< ClassPool< SkipElement< TObject, TModus, RT< TSpec >, TStructuring >, Limited, SimpleAllocator > > & 
	_getElementAlloc( RangeTree< TObject, TModus, RT< TSpec >, TStructuring > & me )
	{
		SEQAN_CHECKPOINT
		return me._allocs._elementAlloc;
	}

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline 
	Allocator< ClassPool< SkipList< TObject, TModus, RT< TSpec >, TStructuring >, Unlimited, SimpleAllocator > > & 
	_getListAlloc( RangeTree< TObject, TModus, RT< TSpec >, TStructuring > & me )
	{
		SEQAN_CHECKPOINT
		return me._allocs._listAlloc;
	}

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline 
	Allocator< SimpleAlloc<> > &
	_getBaseAlloc( RangeTree< TObject, TModus, RT< TSpec >, TStructuring > & me )
	{
		SEQAN_CHECKPOINT
		return me._allocs._baseAlloc;
	}

	
		// accessor für grenzobjekte
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline 
	TObject * 
	_getLBorderObj( RangeTree< TObject, TModus, TSpec, TStructuring > & me )
	{
		SEQAN_CHECKPOINT
		return &me._LBorderObj;
	}


	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
	TObject * 
	_getRBorderObj( RangeTree< TObject, TModus, TSpec, TStructuring > & me )
	{
		SEQAN_CHECKPOINT
		return &me.RBorderObj_;
	}

		// skip list of the main tree
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
	SkipList< TObject, TModus, RT< TSpec >, TStructuring > *
	_getList( RangeTree< TObject, TModus, RT< TSpec >, TStructuring > & me )
	{
		SEQAN_CHECKPOINT
		return me._list;
	}



	
/////////////////////////////////////////////////////////////////////////////////////////
//
//	algorithms
//
/////////////////////////////////////////////////////////////////////////////////////////
	
	
/*DISABLED
.Function.rangeQuery:
..summary:Get the object with maximum priority in the RMT in a given intervall
..cat:Range Tree
..signature:rangeMaxQuery(tree, lower_border, upper_border, dest)
..param.tree:A Range Tree.
...type:RangeTree
..param.lower_border:The object that stores the lower borders for all dimensions, i.e. $key( lower_border ) <= key( point in range )$
..param.lower_border:The object that stores the upper borders for all dimensions, i.e. $key( point in range ) <= key( upper_border )$
..param.dest:A container to save the objects.
..remarks:The size of $dest$ should be sufficient.
..include:seqan/chaining.h
*/

		// perform a range query
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TResultSet >
	void
	rangeQuery( RangeTree< TObject, TModus, TSpec, TStructuring > & me, 
				TObject & first, 
				TObject & second,
				TResultSet & results )
	{
		SEQAN_CHECKPOINT
		if( dimension( me ) > 1 )
			_fingerSearch( _getList( me ), &first, &second, dimension( me ) - 1, results );
		else
			_bottomSearch( _getList( me ), &first, &second, results );
	}

	template< typename TObject, typename TSize > inline
	bool 
	_testBruteForce( TObject & elem,
					 TObject & first,
					 TObject & second,
					 TSize dim )
	{
		SEQAN_CHECKPOINT
		bool in_range = true;
		typename Size< TObject >::Type _dim = 0;
		while( in_range && _dim <= dim )
		{
			in_range = ( ( key( first, _dim ) <= key( elem, _dim ) ) && ( key( elem, _dim ) <= key( second, _dim ) ) );
			++_dim;
		}
		return in_range;
	}

		// test if an element is in range
		// from dim to dim - 1 to 0
	template< typename TObject, typename TSize > inline
	bool 
	_testRange(	TObject & elem,
				TObject & first,
				TObject & second,
				TSize dim )
	{
		SEQAN_CHECKPOINT
		bool in_range = true;
		typename Key< TObject >::Type theKey;
		while ( in_range && dim > 0 )
		{
			theKey = key( elem, dim );
			in_range = ( ( key( first, dim ) <= theKey ) && ( theKey <= key( second, dim ) ) );
			--dim;
		}
		theKey = key( elem, dim );
		in_range &= ( ( key( first, dim ) <= theKey ) && ( theKey <= key( second, dim ) ) );
		return in_range;
	}

	
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
	typename Size< SkipList< TObject, TModus, RT< TSpec >, TStructuring > >::Type
	_getMaximalSLTowerHeight( RangeTree< TObject, TModus, RT< TSpec >, TStructuring > & rt )
	{
		SEQAN_CHECKPOINT
		return log2( length( rt ) ) + 2;
	}


	template< typename TObject, typename TSpec, typename TStructuring, typename TSize >
	void 
	printLayerScores(	SkipList< TObject, SkipListStatic, RT< TSpec >, TStructuring > * /*list*/,
						TSize /*layer*/,
						TSize /*_dim*/ )
	{}

}

#endif // SEQAN_RT_BASE
