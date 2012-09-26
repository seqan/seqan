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

#ifndef SEQAN_HEADER_META_FRAGMENT
#define SEQAN_HEADER_META_FRAGMENT

namespace seqan{

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//			Class MetaFragment_
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*DISABLED
.Class.MetaFragment_:
..summary:Basic data which associates fragments with prededing fragments and stores chain score informations
..cat:Chaining
..signature:MetaFragment_< TFragType >
..param.TFragType:Type of the fragment
..include:seqan/chaining.h
*/

		// get/set the weight of the related fragment
	template< typename TFragType > inline
	typename Weight< TFragType >::Type
	weight( MetaFragment_< TFragType > & me )
	{
		return weight( *me._frag );
	}

	template< typename TFragType, typename TWeight> inline
	void
	setWeight( MetaFragment_< TFragType > & me,
				TWeight weight )
	{
		setWeight( *me._frag, weight );
	}

		// get/set the score of the chain
	template< typename TFragType > inline
	typename Weight< TFragType >::Type 
	score( MetaFragment_< TFragType > & me )
	{
		return me._score;
	}

	template< typename TFragType, typename TWeight > inline
	void
	setScore( MetaFragment_< TFragType > & me,
						TWeight score )
	{
		me._score = score;
	}

		// get/set the priority
	template< typename TFragType > inline
	typename Weight< TFragType >::Type 
	priority( MetaFragment_< TFragType > & me )
	{
		return me._priority;
	}

	template< typename TFragType, typename TWeight > inline
	void
	setPriority( MetaFragment_< TFragType > & me,
					TWeight prio )
	{
		me._priority = prio;
	}

		// get the associated fragment
	template< typename TFragType > inline
	TFragType & 
	_getFrag( MetaFragment_< TFragType > & me )
	{
		return *me._frag;
	}

	template< typename TFragType > inline
	TFragType & 
	_getFrag( const MetaFragment_< TFragType > & me )
	{
		return *me._frag;
	}

		// get preceding fragment
	template< typename TFragType > inline
	MetaFragment_< TFragType > & 
	_getPred( MetaFragment_< TFragType > & me )
	{
		return *me._pred;
	}

	template< typename TFragType > inline 
	MetaFragment_< TFragType > & 
	_getPred( const MetaFragment_< TFragType > & me )
	{
		return *me._pred;
	}

		// set preceding fragment
	template< typename TFragType > inline
	void
	_setPred( MetaFragment_< TFragType > & me, 
				MetaFragment_< TFragType > & pred )
	{
		me._pred = &pred;
	}

	template< typename TFragType > inline
	void
	_setPred( const MetaFragment_< TFragType > & me, 
				MetaFragment_< TFragType > & pred )
	{
		me._pred = &pred;
	}

	template< typename TFragType > inline
	void 
	dump( MetaFragment_< TFragType > & me )
	{
		if( me._frag )
			dump( *me._frag );
		std::cout << me._priority << " " << me._score << std::endl;
	}


	template< typename TFragType >
	struct MetaFragment_
	{
		TFragType * _frag;
			// preceding element in a chain
		typename Weight< TFragType >::Type _priority;
		typename Weight< TFragType >::Type _score;
		MetaFragment_< TFragType > * _pred;

		MetaFragment_()		
			: _frag( NULL )
			, _priority( minValue< typename Weight< TFragType >::Type >() )
			, _score( minValue< typename Weight< TFragType >::Type >() )
			, _pred( NULL )
		{}

		MetaFragment_( TFragType & frag )
			: _frag( &frag )
			, _priority( minValue< typename Weight< TFragType >::Type >() )
			, _score( minValue< typename Weight< TFragType >::Type >() )
			, _pred( NULL )
		{}

		MetaFragment_( const MetaFragment_ & old )
			: _frag( old._frag)
			, _priority( old._priority )
			, _score( old._score )
			, _pred( old._pred )
		{}

		MetaFragment_ & operator=( const MetaFragment_ & old )
		{
			if ( this == &old ) 
				return *this;
			_frag = old._frag;
			_priority = old._priority;
			_score = old._score;
			_pred = old._pred;
			return *this;
		}

	};


}

#endif // SEQAN_HEADER_META_FRAGMENT
