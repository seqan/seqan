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


#ifndef SEQAN_HEADER_TREECHAIN_H
#define SEQAN_HEADER_TREECHAIN_H

#include <algorithm>
#include <vector>

namespace seqan{

		// compute chain spec for G0 and G1 cost metric
	template< typename TSource, typename TDest, typename TScoreValue, typename TScoreType, typename TStructuring, typename TCostModell, typename TSpec >
	TScoreValue
	_computeChain( TSource & source, 
					TDest & dest, 
					TCostModell cost, 
					Score< TScoreValue, TScoreType > const & score_,
					TStructuring,
					TSpec spec )
	{
		SEQAN_ASSERT_NOT(empty(source));
			// define some basic types
		typedef typename Value< TSource >::Type FragType;
		typedef typename Weight< FragType >::Type WeightType;
		typedef typename Key< FragType >::Type PositionType;
		typedef typename Size< FragType >::Type SizeType;
		typedef typename TSpec::Type SpecType;

		SizeType dim = dimension( value( begin( source ) ) );
		
			// construct containers for classes
		String< MetaFragment_< FragType > > metas;
		reserve( metas, length( source ) + 2 );

		std::vector< WrapperPoint_< FragType > > points;
		points.reserve( 2 * ( length( source ) + 2 ) );
		
		//String< WrapperPoint_< FragType > > points;
		//reserve( points, 2 * ( length( source ) + 2 ) );

		String< ChainPoint_< FragType, SpecType > > end_points;
		reserve( end_points, length( source ) + 2 );
		
			// define origin and terminus fragment
		FragType startingFrag( dim );		
		FragType endFrag( dim );

			// build the environment (construct wrapper points, chain point, get coordinates of the terminus)
		_buildChainEnvironment( source, metas, points, end_points, startingFrag, endFrag, spec );

		typename Iterator< String< MetaFragment_< FragType > > >::Type lastMeta = end( metas );
		goPrevious( lastMeta );

			// set the score of the origin from -infinity to 0
		setScore( value( begin( metas ) ), 0 );

			// sort the wrapper points to apply the line sweep paradigma
		std::sort( points.begin(), points.end(), ChainSorter_< WrapperPoint_< FragType > >( ) );
	
			// build the RMT
		RangeTree< ChainPoint_< FragType, SpecType >, SkipListStatic, RT< MaxTree< > >, TStructuring > tree( end_points, dim-1 );

			// algorithm main loop
			// traverse wrapper points
		typename std::vector< WrapperPoint_< FragType > >::iterator pointIt = points.begin();
		while( pointIt != points.end() )
		{
				// actual point is the beginning of a frag
				// => search for preceding fragment
			MetaFragment_< FragType > & meta = _meta(*pointIt );
			if( !_isEnd( *pointIt ) )
			{
				ChainPoint_< FragType, SpecType > buffer( meta, dim - 1, true );
				ChainPoint_< FragType, SpecType > * result = rangeMaxQuery( tree, buffer );

				SEQAN_ASSERT(result != NULL);

				_setPred( meta, _meta( *result ) );
				setScore( meta, _maxPriority( value( lastMeta ), meta, *result, cost, score_, dim ) );
			}
			else{
					// point is the end of a frag
					// => activate it
				size_t offset = &meta - &_meta( *points.begin() );
				ChainPoint_< FragType, SpecType > & point = end_points[ offset ];

				setPriority( point, score( meta ) - _activatePriority( value( lastMeta ), point, cost, score_, dim ) );
				activate( tree, point );
			}
			goNext( pointIt );
		}
			// perform backtracking
		return _chainTrace( dest, metas );
	}
	
}

#endif

