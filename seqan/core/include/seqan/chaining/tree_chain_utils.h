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


#ifndef SEQAN_HEADER_TREECHAIN_UTILS_H
#define SEQAN_HEADER_TREECHAIN_UTILS_H



namespace seqan{

//////////////////////////////////////////////////////////////////////////////////////////
//	helper functions
//////////////////////////////////////////////////////////////////////////////////////////

		// calculate faculty of a number
	template< typename TSimpleType > inline
	TSimpleType
	fac( TSimpleType data )
	{
		TSimpleType value = 1;
		for( TSimpleType i = 2; i <= data; ++i )
		{
			value *= i;
		}
		return value;
	}

		// initailize a permutation to 0, 1, ..., k
	template< typename TPerm, typename TSize > inline
	void
	_initPerm( TPerm & perm, 
				TSize length )
	{
		for( TSize i = 0; i < length; ++i )
		{
			appendValue( perm, i );
		}
	}
		
		// reset a permutation to 0, 1, ..., k
	template< typename TPerm > inline
	void
	_resetPerm( TPerm & perm )
	{
		typename Iterator< TPerm, Rooted >::Type permIt = begin( perm );
		int i = 0;
		while( permIt != end( perm ) )
		{
			assignValue( permIt, i );
			++i;
			goNext( permIt );
		}
	}

		// get the permutation e.g. order of to sequences of numbers
	template< typename TData, typename TItPerm, typename TSize > inline
	void
	_getPerm( TData * values,
				TItPerm perm,
				TSize length )
	{
		TData max = minValue< TData>();
		TSize maxIndex = 0;
		for( TSize i = 0; i < length; ++i )
		{
			for( TSize j = 0; j < length; ++j )
			{
				if( values[ j ] > max )
				{
					max = values[ j ];
					maxIndex = j;
					values[ j ] = minValue< TData>();
				}
			}
			assignValue( perm, maxIndex );
			maxIndex = 0;
			max = minValue< TData>();
			goNext( perm );
		}		
	}

		// get the differences of coordinates of two fragments
	template< typename TData, typename TSize, typename FragType >
	void
	_getPermDifference( TData * values,
						TSize dim,
						FragType & upper,
						FragType & lower )
	{
		for( TSize i = 0; i < dim; ++i )
		{
			values[ i ] = leftPosition( upper, i ) - rightPosition( lower, i );
		}		
	}

//////////////////////////////////////////////////////////////////////////////////////////
//	transformations
//////////////////////////////////////////////////////////////////////////////////////////

		// transform the coordinates of a chain point
	template< typename FragType, typename SpecType, typename TPerm > inline
	void
	_chainTransformCoords( ChainPoint_< FragType, SpecType > & point_src,
							ChainPoint_< FragType, SpecType > & point_dst,
							TPerm & perm )
	{
		typename Iterator< TPerm >::Type permIt = begin( perm );
		typename Iterator< TPerm >::Type permEnd = end( perm );
		goPrevious( permEnd );
		while( permIt != permEnd )
		{
			setKey( point_dst, *permIt, key( point_src, *permIt ) - key( point_src, *( permIt + 1 ) ) );
			++permIt;
		}
		setKey( point_dst, *permIt, key( point_src, *permIt ) );
	}

			// transform the coordinates of a chain point for searching
	template< typename FragType, typename SpecType, typename TPerm > inline
	void
	_chainTransformCoordsSearch( ChainPoint_< FragType, SpecType > & point_src,
									ChainPoint_< FragType, SpecType > & point_dst,
									TPerm & perm )
	{
		typename Iterator< TPerm >::Type permIt = begin( perm );
		typename Iterator< TPerm >::Type permEnd = end( perm );
		goPrevious( permEnd );
		while( permIt != permEnd )
		{
			setKey( point_dst, *permIt, key( point_src, *permIt ) - key( point_src, *( permIt + 1 ) ) + 1 );
			++permIt;
		}
		setKey( point_dst, *permIt, key( point_src, *permIt ) );
	}


//////////////////////////////////////////////////////////////////////////////////////////
// cost functions
//////////////////////////////////////////////////////////////////////////////////////////
	
		// calulate the score
	template< typename FragType, typename SpecType, typename TCostModell, typename TScore, typename TSize > inline
	typename Weight< FragType >::Type
	_maxPriority( MetaFragment_< FragType > & last_meta,
						MetaFragment_< FragType > & current_meta,
						ChainPoint_< FragType, SpecType > & point,
						TCostModell cost, 
						TScore const & score_,
						TSize dim );

	template< typename FragType, typename SpecType, typename TScore, typename TSize > inline
	typename Weight< FragType >::Type
	_maxPriority( MetaFragment_< FragType > & ,
					MetaFragment_< FragType > & current_meta,
					ChainPoint_< FragType, SpecType > & point,
					GZeroCost, 
					TScore const &,
					TSize )
	{
		return weight( current_meta ) + priority( point );
	}


		// calculate the maximum priority
	template< typename FragType, typename SpecType, typename TScore, typename TSize > inline
	typename Weight< FragType >::Type
	_maxPriority( MetaFragment_< FragType > &,
					MetaFragment_< FragType > & current_meta,
					ChainPoint_< FragType, SpecType > & point,
					GOneCost, 
					TScore const & score_,
					TSize dim )
	{
		typename Weight< FragType >::Type prio = weight( current_meta );
		prio += score( _meta( point ) );
		prio -= _costGOne( current_meta, _meta( point ), score_, dim );
		return prio;
	}

	
		// calculate the priority for activation
	template< typename FragType, typename SpecType, typename TCostModell, typename TScore, typename TSize > inline
	typename Weight< FragType >::Type
	_activatePriority( MetaFragment_< FragType > & last_meta,
						ChainPoint_< FragType, SpecType > & point,
						TCostModell cost, 
						TScore const & score_,
						TSize dim );


	template< typename FragType, typename SpecType, typename TScore, typename TSize > inline
	typename Weight< FragType >::Type
	_activatePriority( MetaFragment_< FragType > &,
						ChainPoint_< FragType, SpecType > &,
						GZeroCost, 
						TScore const &,
						TSize )
	{
		return 0;
	}


	template< typename FragType, typename SpecType, typename TScore, typename TSize > inline
	typename Weight< FragType >::Type
	_activatePriority( MetaFragment_< FragType > & last_meta,
						ChainPoint_< FragType, SpecType > & point,
						GOneCost, 
						TScore const & score_,
						TSize dim )
	{
		return _costGOne( last_meta, _meta( point ), score_, dim );
	}

	
		// the cost function for manhattan metric
	template< typename FragType, typename TScore, typename TSize >
	typename Weight< FragType >::Type
	_costGOne( MetaFragment_< FragType > & upper,
				MetaFragment_< FragType > & lower,
				TScore const & score,
				TSize dim )
	{
		typename Weight< FragType >::Type weight = 0;
		for( typename Size< FragType >::Type i = 0; i < dim; ++i )
		{
//!!!Change dist semantics
//			weight += ( scoreGapExtend( score ) * ( leftPosition( _getFrag( upper ), i ) - rightPosition( _getFrag( lower ), i ) ) );
			weight += ( scoreGap( score ) * ( leftPosition( _getFrag( upper ), i ) - rightPosition( _getFrag( lower ), i ) - 1) );
		}
		return weight;
	}

		// the cost function for SoP metric
	template< typename FragType, typename TItPerm, typename TScoreValue, typename TScoreType, typename TSize >
	typename Weight< FragType >::Type
	_costGSoP( MetaFragment_< FragType > & upper,
				MetaFragment_< FragType > & lower,
				Score< TScoreValue, TScoreType > const & score,
				TItPerm permBeg,
				TItPerm,
				TSize dim )
	{
		typename Weight< FragType >::Type weight = 0;
		typename Weight< FragType >::Type weight_buffer = 0;
		typename Weight< FragType >::Type delta;
		typename Size< FragType >::Type dim_factor = dim - 1;
		for( typename Size< FragType >::Type i = 0; i < dim; ++i )
		{
//!!!Change dist semantics
//			delta = static_cast< typename Weight< FragType >::Type >( leftPosition( _getFrag( upper ), value( permBeg ) ) - rightPosition( _getFrag( lower ), value( permBeg ) ) );
			delta = static_cast< typename Weight< FragType >::Type >( leftPosition( _getFrag( upper ), value( permBeg ) ) - rightPosition( _getFrag( lower ), value( permBeg ) ) -1);
			weight_buffer = ( scoreGap( score ) * delta );
			weight_buffer *= static_cast< typename Weight< FragType >::Type >( ( dim_factor - i ) );
			weight += weight_buffer;
			weight_buffer = ( delta * ( scoreMismatch( score ) - scoreGap( score ) ) );
			weight_buffer *= static_cast< typename Weight< FragType >::Type >( i );
			weight += weight_buffer;
			goNext( permBeg );
		}
		return weight;
	}

//////////////////////////////////////////////////////////////////////////////////////////
//	sorting
//////////////////////////////////////////////////////////////////////////////////////////


		// struct for std::sort of wrapper points
	template< typename T >
	struct
	ChainSorter_
	{
		inline bool
		operator()( T & first, T & second  )
		{
			if ( key( first ) < key( second ) )
				return true;
			else if ( key( first ) == key( second ) && _isBegin( first ) && _isEnd( second ) )
				return true;
			return false;
		}
		inline bool 
		operator()( const T & first, const T & second  )
		{
			if ( key( first ) < key( second ) )
				return true;
			else if ( key( first ) == key( second ) && _isBegin( first ) && _isEnd( second ) )
				return true;
			return false;
		}

		ChainSorter_()
		{}	

	};
	

//////////////////////////////////////////////////////////////////////////////////////////
//	dynamic programming helper functions
//////////////////////////////////////////////////////////////////////////////////////////


		// backtracking
	template< typename TDest, typename TMetas >
	typename Weight< typename Value< TDest >::Type >::Type
	_chainTrace( TDest & dest,
					TMetas & metas )
	{
		typedef typename Value< TDest >::Type FragType;
		typename Iterator< TMetas >::Type meta = end( metas );
		goPrevious( meta );
		MetaFragment_< FragType > * pMeta = & value( meta );
		typename Weight< FragType >::Type chain_score = score( *pMeta );
		//pMeta = &_getPred( *pMeta );
		while( pMeta != &value( begin( metas ) ) )
		{
			SEQAN_ASSERT(&_getFrag( *pMeta ) != 0);
			appendValue( dest, _getFrag( *pMeta ) );
			pMeta = &_getPred( *pMeta );
		}
		appendValue( dest, _getFrag( *pMeta ) );
		std::reverse( begin( dest ), end( dest ) );
//		typename Iterator< TDest >::Type destIt = begin( dest );
		return chain_score;
	}


//////////////////////////////////////////////////////////////////////////////////////////
//	initialization helper functions
//////////////////////////////////////////////////////////////////////////////////////////

		// init of starting frag (the origin)
	template< typename FragType, typename TSize > inline
	void
	_initStartingFrag( FragType & frag,
							TSize dim )
	{
		TSize dim_counter = 0;
		while( dim_counter != dim )
		{
			_setLeftPosition( frag, dim_counter, 0 );
			_setRightPosition( frag, dim_counter, 0 );
			++dim_counter;
		}
		setWeight( frag, 0 );
	}

		// init of needed variables
		// * the wrapper points
		// * the chain points
		// * the metainformation structures
		// spec for G0 and G1 metric
	template< typename FragType, typename TSource, typename TMetas, typename TWPoints, typename TCPoints, typename TSpec > inline
	void
	_buildChainEnvironment( TSource & source,  
								TMetas & metas, 
								TWPoints & wPoints,
								TCPoints & cPoints,
								FragType & startingFrag,
								FragType & endFrag,
								TSpec )
	{
		typedef typename Key< FragType >::Type KeyType;
		typedef typename Size< FragType >::Type SizeType;
		typedef typename TSpec::Type SpecType;
		SizeType dim = dimension( value( begin( source ) ) );
		SizeType lower_dim  = dim - 1;
				
			// initialize starting frag
		_initStartingFrag( startingFrag, dim );
		appendValue( metas, MetaFragment_< FragType >( startingFrag ) );
		wPoints.push_back( WrapperPoint_< FragType >( value( begin( metas ) ), true ) );
		appendValue( cPoints, ChainPoint_< FragType, SpecType >( value( begin( metas ) ) ) );

			// buffers to find the maximal coordinates
		KeyType * maxCoords;
		allocate( maxCoords, maxCoords, dim );
		for( SizeType i = 0; i < dim; ++i )
		{
			maxCoords[ i ] = 0;
		}

			// traverse the set of fragments
		typename Iterator< TSource, Rooted >::Type sourceIt = begin( source );
		typename Iterator< TMetas, Rooted >::Type metaIt = begin( metas );
		goNext( metaIt );
		while( sourceIt != end( source ) )
		{
			appendValue( metas, MetaFragment_< FragType >( value( sourceIt ) ) );
			wPoints.push_back( WrapperPoint_< FragType >( value( metaIt ), leftPosition( value( sourceIt ), lower_dim ), false ) );
			wPoints.push_back( WrapperPoint_< FragType >( value( metaIt ), rightPosition( value( sourceIt ), lower_dim ), true ) );
			appendValue( cPoints, ChainPoint_< FragType, SpecType >( value( metaIt ) ) );

			for( SizeType i = 0; i < dim; ++i )
			{
				if( maxCoords[ i ] < rightPosition( value( sourceIt ), i ) )
					maxCoords[ i ] = rightPosition( value( sourceIt ), i );
			}
			
			goNext( metaIt );
			goNext( sourceIt );
		}

			// set the coordinates of the terminus
		for( SizeType dim_counter = 0; dim_counter < dim; ++dim_counter )
		{
			_setLeftPosition( endFrag, dim_counter, maxCoords[ dim_counter ] + 1 );
			_setRightPosition( endFrag, dim_counter, maxCoords[ dim_counter ] + 1 );
		}
		appendValue( metas, MetaFragment_< FragType >( endFrag ) );
		wPoints.push_back( WrapperPoint_< FragType >( *metaIt, false ) );
		deallocate( maxCoords, maxCoords, dim );
		appendValue( cPoints, ChainPoint_< FragType, SpecType >( *metaIt ) );
	}


		// spec for G_SoP metric
	template< typename FragType, typename TSource, typename TMetas, typename TWPoints, typename TPoints, typename TPerm, typename TSpec > inline
	void
	_buildChainEnvironment( TSource & source, 
								TMetas & metas, 
								TWPoints & wPoints,
								TPoints  & tPoints,
								FragType & startingFrag,
								FragType & endFrag,
								TPerm & perm,
								typename Size< FragType >::Type fac,
								TSpec &)
	{
		typedef typename Key< FragType >::Type KeyType;
		typedef typename Size< FragType >::Type SizeType;
		typedef typename TSpec::Type SpecType;
		SizeType dim = dimension( value( begin( source ) ) );
		SizeType lower_dim = dim - 1;
		SizeType dim_counter = 0;
		
		_initStartingFrag( startingFrag, dim );

		appendValue( metas, MetaFragment_< FragType >( startingFrag ) );
		wPoints.push_back( WrapperPoint_< FragType >( value( begin( metas ) ), true ) );
		
		KeyType * maxCoords;
		allocate( maxCoords, maxCoords, dim );
		for( SizeType i = 0; i < dim; ++ i )
			maxCoords[ i ] = 0;

		typename Iterator< TSource, Rooted >::Type sourceIt = begin( source );
		typename Iterator< TMetas, Rooted >::Type metaIt = begin( metas );
		goNext( metaIt );

		while( sourceIt != end( source ) )
		{
			appendValue( metas, MetaFragment_< FragType >( value( sourceIt ) ) );
			wPoints.push_back( WrapperPoint_< FragType >( value( metaIt ), leftPosition( value( sourceIt ), lower_dim ), false ) );
			wPoints.push_back( WrapperPoint_< FragType >( value( metaIt ), rightPosition( value( sourceIt ), lower_dim ), true ) );

			for( dim_counter = 0; dim_counter < dim; ++dim_counter )
			{
				if( maxCoords[ dim_counter ] < rightPosition( value( sourceIt ), dim_counter ) )
					 maxCoords[ dim_counter ] = rightPosition( value( sourceIt ), dim_counter );
			}

			goNext( metaIt );
			goNext( sourceIt );
		}
		dim_counter = 0;
		while( dim_counter != dim )
		{
			_setLeftPosition( endFrag, dim_counter, maxCoords[ dim_counter ] + 1);
			_setRightPosition( endFrag, dim_counter, maxCoords[ dim_counter ] + 1 );
			++dim_counter;
		}
		appendValue( metas, MetaFragment_< FragType >( endFrag ) );
		deallocate( maxCoords, maxCoords, dim );
		wPoints.push_back( WrapperPoint_< FragType >( *metaIt, false ) );
		metaIt = begin( metas );
		
			// transformate the point coordinates
		typename Iterator< TPoints >::Type transIt =  begin( tPoints );
			// traverse all dim! permutations
		for( SizeType i = 0; i < fac; ++i )
		{
			assignValue( transIt, typename Value< TPoints >::Type() );
			reserve( value( transIt ), length( source ) + 2 );
			
				// transform points
			typename Iterator< typename Value< TPoints >::Type >::Type cPointIt = begin( value( transIt ) );
			for( SizeType pointCount = 0; pointCount < ( length( source ) + 2 ); ++pointCount )
			{
				ChainPoint_< FragType, SpecType > buffer( value( metaIt ), dim );
				appendValue( value( transIt ), buffer );
				_chainTransformCoords( buffer, value( cPointIt ), perm );
				goNext( metaIt );
				goNext( cPointIt );
			}
			std::next_permutation( begin( perm ), end( perm ) );
			metaIt = begin( metas );
			goNext( transIt );
		}
	}

		// construct dim! range trees for gSoP metric
	template< typename TTrees, typename TTPoints, typename TSize >
	void
	_buildChainTrees( TTrees & trees, 
						TTPoints & tPoints, 
						TSize dim,
						TSize facValue )
	{
		typename Iterator< TTPoints >::Type tPointIt = begin( tPoints );
		typename Iterator< TTrees >::Type treeIt = begin( trees );
		for( TSize i = 0; i < facValue; ++i )
		{
			allocate( value( treeIt ), value( treeIt ), 1 );
			new( value( treeIt ) ) typename Value< typename Value< TTrees >::Type >::Type ( value( tPointIt ), dim );
			goNext( tPointIt );
			goNext( treeIt );
		}
	}

		// delete all RMT's
	template< typename TTrees, typename TSize >
	void
	_deleteChainTrees( TTrees & trees, 
							TSize facValue )
	{
		typename Iterator< TTrees >::Type treeIt = begin( trees );
		for( TSize i = 0; i < facValue; ++i )
		{
			delete ( value( treeIt ) );
			goNext( treeIt );
		}
	}

}

#endif
