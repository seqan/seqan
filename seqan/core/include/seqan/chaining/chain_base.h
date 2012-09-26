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

#ifndef SEQAN_HEADER_CHAIN_BASE_H
#define SEQAN_HEADER_CHAIN_BASE_H

/*
 *  chain_base.h
 *  chaining
 *
 *  Created by Hendrik Woehrle
 *
 *	Basis declarations & definitions for chaining algorithms
 *
 */

namespace seqan
{
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//			class forward declarations
//
////////////////////////////////////////////////////////////////////////////////////////////////////////
		
		// basic data structure for fragments (see seed_base.h)
	//template< typename TBorder, typename TSpec = Default >
	//struct Seed;

		// wrapper data structure for point in the RMT
	template< typename T, typename TSpec >
	struct ChainPoint_;

		// wrapper for points for the line sweep paradigma
	template< typename T >
	struct WrapperPoint_;

		// structure for metainformation about fragments in chains
	template< typename T >
	struct MetaFragment_;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//			tag classes
//
////////////////////////////////////////////////////////////////////////////////////////////////////////
		
	struct GenericChaining_;
	typedef Tag<GenericChaining_> const GenericChaining;

	struct RangeTreeChaining_;
	typedef Tag<RangeTreeChaining_> const RangetreeChaining;


	struct GZeroCost
	{};
	//typedef Tag<ZeroCost_> const GZeroCost;

	struct GOneCost
	{};
	//typedef Tag<G_1_Cost_> const GOneCost;

	struct GInftyCost
	{};
	//typedef Tag<G_Inf_Cost_> const GInftyCost;

	struct GSumOfPairCost
	{};
	//typedef Tag<F_SoP_Cost_> const GInftyCost;

	template< typename Spec >
	struct ChainSpecType_
	{
		typedef Spec Type;
	};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//			special metatype declarations
//
///////////////////////////////////////////////////////////////////////////////////////////////////////
// sizes

	//template< typename TBorder, typename TSpec >
	//struct Size< Seed< TBorder, TSpec > >
	//{
	//	typedef size_t Type;
	//};

	template< typename T, typename TSpec >
	struct Size< ChainPoint_< T, TSpec > >
	{
		typedef size_t Type;
	};

	template< typename T >
	struct Size< WrapperPoint_< T > >
	{
		typedef size_t Type;
	};

	template< typename T >
	struct Size< MetaFragment_< T > >
	{
		typedef size_t Type;
	};


//////////////////////////////////////////////////////////////////////////
// key types

	//template< typename TBorder, typename TSpec >
	//struct Key< Seed< TBorder, TSpec > >
	//{
	//	typedef TBorder Type;
	//};

	template< typename T, typename TSpec >
	struct Key< ChainPoint_< T, TSpec > >
	{
		typedef typename Key< T >::Type Type;
	};

	template< typename T >
	struct Key< WrapperPoint_< T > >
	{
		typedef typename Key< T >::Type Type;
	};

	template< typename T >
	struct Key< MetaFragment_< T > >
	{
		typedef typename Key< T >::Type Type;
	};



///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//			functions
//
///////////////////////////////////////////////////////////////////////////////////////////////////////


/*
	template< typename TSource, typename TDest, typename TScoreValue, typename TScoreType, typename TStructuring, typename TAlgorithm >
	TScoreValue
	computeChain( TSource & source, TDest & dest, Score< TScoreValue, TScoreType > const & score_, TAlgorithm tag, TStructuring structuring );

	template< typename TSource, typename TDest, typename TCostModel, typename TScoreValue, typename TScoreType, typename TStructuring, typename TAlgorithm, typename TSpec >
	TScoreValue
	_computeChain( TSource & source, TDest & dest, TCostModel type, Score< TScoreValue, TScoreType > const & score_, TAlgorithm tag, TStructuring structuring, TSpec spec );

*/
	template< typename TSource, typename TDest, typename TScoreValue, typename TScoreType, typename TStructuring >
	TScoreValue
	computeChain( TSource & source, TDest & dest, Score< TScoreValue, TScoreType > const & score_, TStructuring structuring )
	{

//		SEQAN_ASSERT_MSG( scoreGapOpen( score_ ) == scoreGapExtend( score_ ), "Chaining only supports linear gap costs" );
//		SEQAN_ASSERT_MSG( scoreGapOpen( score_ ) >= 0 && scoreGapExtend( score_ ) >= 0 && scoreMismatch( score_ ) >= 0, "Scores should be positive" );
		switch( dimension( value( begin( source ) ) ) )
		{
			case 1: SEQAN_ASSERT_FAIL("One dimensional chaining not supported");
					return 0;
			case 2: if( scoreMismatch( score_ ) == 0 && scoreGap( score_ ) == 0 )
					{
						return _computeChain( source, dest, GZeroCost(), score_, structuring, ChainSpecType_< Array< 1 > >() );
					}
					else if( scoreMismatch( score_ ) > 0 && scoreGap( score_ ) > 0 )
					{
						if( scoreMismatch( score_ ) == scoreGap( score_ ) )
							return _computeChain( source, dest, GOneCost(), score_, structuring, ChainSpecType_< Array< 1 > >() );
						else //if( scoreMismatch( score_ ) > 2 * scoreGapExtend( score_ ) )
							return _computeChain( source, dest, GSumOfPairCost(), score_, structuring,ChainSpecType_< Array< 2 > >() );
					}
					SEQAN_ASSERT_FAIL("Gap/mismatch model not supported");
					break;
			case 3: if( scoreMismatch( score_ ) == 0 && scoreGap( score_ ) == 0 )
					return _computeChain( source, dest, GZeroCost(), score_, structuring, ChainSpecType_< Array< 2 > >() );
					else if( scoreMismatch( score_ ) > 0 && scoreGap( score_ ) > 0 )
					{
						if( scoreMismatch( score_ ) == scoreGap( score_ ) )
							return _computeChain( source, dest, GOneCost(), score_, structuring, ChainSpecType_< Array< 2 > >() );
						else //if( scoreMismatch( score_ ) > 2 * scoreGapExtend( score_ ) )
							return _computeChain( source, dest, GSumOfPairCost(), score_, structuring,ChainSpecType_< Array< 3 > >() );
					}
					SEQAN_ASSERT_FAIL("Gap/mismatch model not supported");
					break;
			case 4:	if( scoreMismatch( score_ ) == 0 && scoreGap( score_ ) == 0 )
					return _computeChain( source, dest, GZeroCost(), score_, structuring, ChainSpecType_< Array< 3 > >() );
					else if( scoreMismatch( score_ ) > 0 && scoreGap( score_ ) > 0 )
					{
						if( scoreMismatch( score_ ) == scoreGap( score_ ) )
							return _computeChain( source, dest, GOneCost(), score_, structuring, ChainSpecType_< Array< 3 > >() );
						else //if( scoreMismatch( score_ ) > 2 * scoreGapExtend( score_ ) )
							return _computeChain( source, dest, GSumOfPairCost(), score_, structuring,ChainSpecType_< Array< 4 > >() );
					}
					SEQAN_ASSERT_FAIL("Gap/mismatch model not supported");
					break;
			case 5:	if( scoreMismatch( score_ ) == 0 && scoreGap( score_ ) == 0 )
					return _computeChain( source, dest, GZeroCost(), score_, structuring, ChainSpecType_< Array< 4 > >() );
					else if( scoreMismatch( score_ ) > 0 && scoreGap( score_ ) > 0 )
					{
						if( scoreMismatch( score_ ) == scoreGap( score_ ) )
							return _computeChain( source, dest, GOneCost(), score_, structuring, ChainSpecType_< Array< 4 > >() );
						else //if( scoreMismatch( score_ ) > 2 * scoreGapExtend( score_ ) )
							return _computeChain( source, dest, GSumOfPairCost(), score_, structuring,ChainSpecType_< Array< 5 > >() );
					}
					SEQAN_ASSERT_FAIL("Gap/mismatch model not supported");
					break;
			default:if( scoreMismatch( score_ ) == 0 && scoreGap( score_ ) == 0 )
					return _computeChain( source, dest, GZeroCost(), score_, structuring, ChainSpecType_< Default >() );
					else if( scoreMismatch( score_ ) > 0 && scoreGap( score_ ) > 0 )
					{
						if( scoreMismatch( score_ ) == scoreGap( score_ ) )
							return _computeChain( source, dest, GOneCost(), score_, structuring, ChainSpecType_< Default >() );
						else //if( scoreMismatch( score_ ) > 2 * scoreGapExtend( score_ ) )
							return _computeChain( source, dest, GSumOfPairCost(), score_, structuring, ChainSpecType_< Default >() );
					}
					SEQAN_ASSERT_FAIL("Gap/mismatch model not supported");

		}
		return minValue< TScoreValue >();
	}

///////////////////////////////////////////////////////////////////////////////////////////////////////

template< typename TSource, typename TDest, typename TScoring>
inline typename Value<TScoring>::Type
globalChaining(TSource & source,
	  TDest & dest,
	  TScoring const & scoring,
	  RangetreeChaining)
{
	typedef typename Iterator<TSource, Standard>::Type TSourceIterator;
	typedef typename Iterator<TDest, Standard>::Type TDestIterator;

	SEQAN_ASSERT_GT(length(source), 0u);

	//transform coordinates to old style ("end is last item")
	unsigned int dim = dimension(source[0]);
	for (TSourceIterator it = begin(source, Standard()); it < end(source, Standard()); ++it)
	{
		for (unsigned int i = 0; i < dim; ++i)
		{
			_setLeftPosition(*it, i, leftPosition(*it, i) + 1);
		}
	}

	//compute chain
	typename Value<TScoring>::Type ret_value = computeChain(source, dest, scoring, SemiDeferred()); //Deferred(), Complete()

	//retransform coordinates to new style ("end is behind last item")
	for (TSourceIterator it = begin(source, Standard()); it < end(source, Standard()); ++it)
	{
		for (unsigned int i = 0; i < dim; ++i)
		{
			_setLeftPosition(*it, i, leftPosition(*it, i) - 1);
		}
	}
	for (TDestIterator it = begin(dest, Standard()) + 1; it < end(dest, Standard()); ++it)
	{
		for (unsigned int i = 0; i < dim; ++i)
		{
			_setLeftPosition(*it, i, leftPosition(*it, i) - 1);
		}
	}

	//adjust right positions of the end fragment
	TDestIterator it2 = end(dest, Standard());
	--it2;
	for (unsigned int i = 0; i < dim; ++i)
	{
		_setRightPosition(*it2, i, rightPosition(*it2, i) - 1);
	}

	//return score
	return ret_value;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////





//implementation for Default
template< typename TSource, typename TDest, typename TScoring>
inline typename Value<TScoring>::Type
globalChaining(TSource & source, 
	  TDest & dest, 
	  TScoring const & scoring,
	  Default)
{
	return globalChaining(source, dest, scoring, GenericChaining()); //default is GenericChaining
}
template< typename TSource, typename TDest, typename TValue>
inline TValue
globalChaining(TSource & source, 
	  TDest & dest, 
	  Score<TValue, Zero> const & scoring,
	  Default)
{
	return globalChaining(source, dest, scoring, RangetreeChaining());
}
template< typename TSource, typename TDest, typename TValue>
inline TValue
globalChaining(TSource & source, 
	  TDest & dest, 
	  Score<TValue, Manhattan> const & scoring,
	  Default)
{
	return globalChaining(source, dest, scoring, RangetreeChaining());
}
template< typename TSource, typename TDest, typename TValue>
inline TValue
globalChaining(TSource & source, 
	  TDest & dest, 
	  Score<TValue, ChainSoP> const & scoring,
	  Default)
{
	return globalChaining(source, dest, scoring, RangetreeChaining());
}


//chain(3) => chain(4)
template< typename TSource, typename TDest, typename TScoring>
inline typename Value<TScoring>::Type
globalChaining(TSource & source, 
	  TDest & dest, 
	  TScoring const & scoring)
{
	return globalChaining(source, dest, scoring, Default());
}


///////////////////////////////////////////////////////////////////////////////////////////////////////

}

#endif	// SEQAN_HEADER_CHAIN_BASE_H
