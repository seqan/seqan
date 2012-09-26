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

#ifndef SEQAN_H_CHAIN_WRAPPER_POINT
#define SEQAN_H_CHAIN_WRAPPER_POINT

namespace seqan{

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//			Class WrapperPoint_
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
.Class.WrapperPoint_:
..summary:Data structure to represent a begin or end point of a fragment
..cat:Chaining
..signature:WrapperPoint_<TBorder>
..param.TBorder:Type of the class that represents the limits (multidimensional point) of an fragment.
..include:seqan/chaining.h
*/

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//			general functions
//
////////////////////////////////////////////////////////////////////////////////////////////////////////

			// get the Key of the bording points of the related fragment depends on _end 
	template< typename TFragType > inline 
	typename Key< TFragType >::Type
	key( WrapperPoint_< TFragType > & me )
	{
		return me._key;
	}

	template< typename TFragType > inline 
	typename Key< TFragType >::Type
	key( const WrapperPoint_< TFragType > & me )
	{
		return me._key;
	}


		// get the related fragment
	
	template< typename TFragType > inline
	TFragType & 
	_getFrag( WrapperPoint_< TFragType > & me )
	{
		SEQAN_ASSERT( me._frag != NULL );
		SEQAN_ASSERT( me._meta != NULL );
		SEQAN_ASSERT( &_getFrag( *me._meta ) == me._frag );
		return *me._frag;
	}

	template< typename TFragType > inline
	TFragType & 
	_getFrag( const WrapperPoint_< TFragType > & me )
	{
		SEQAN_ASSERT( me._frag != NULL );
		SEQAN_ASSERT( me._meta != NULL );
		SEQAN_ASSERT( _getFrag( *me._meta ) == me._frag );
		return *me._frag;
	}

	template< typename TFragType > inline
	MetaFragment_< TFragType > & 
	_meta( WrapperPoint_< TFragType > & me )
	{
		SEQAN_ASSERT( me._frag != NULL );
		SEQAN_ASSERT( me._meta != NULL );
		return *me._meta;
	}

	template< typename TFragType > inline
	void
	_setMeta( WrapperPoint_< TFragType > & me,
				MetaFragment_< TFragType > & meta )
	{
		me._meta = &meta;
		me._frag = _getFrag( meta );
	}

	
	template< typename TFragType > inline
	bool 
	_isEnd( WrapperPoint_< TFragType > & me )
	{
//IOREV _notio_
		return me._end;
	}

	template< typename TFragType > inline 
	bool 
	_isEnd( const WrapperPoint_< TFragType > & me )
	{
//IOREV _notio_
		return me._end;
	}
	
	template< typename TFragType > inline 
	bool
	_isBegin( WrapperPoint_< TFragType > & me )
	{
//IOREV _notio_
		return !me._end;
	}

	template< typename TFragType > inline 
	bool
	_isBegin( const WrapperPoint_< TFragType > & me )
	{
//IOREV _notio_
		return !me._end;
	}

	template< typename TFragType > inline 
	typename Size< TFragType >::Type
	dimension( const WrapperPoint_< TFragType > & me )
	{
		return dimension( *me._frag );
	}

	
	template< typename TFragType >
	struct WrapperPoint_
	{
			// related triple, which stores the preceding frgament, chain value...
		TFragType * _frag;
			// the related key
		typename Key< TFragType >::Type _key;
			// meta information struct for that point
		MetaFragment_< TFragType > * _meta;
			// the point is either the end ( _end == true ) or the beginning of a fragment
		bool _end;
			
		
	#ifdef SEQAN_CHAIN_DEBUG_
		friend inline
		void 
		dump( WrapperPoint_< TFragType > & me )
		{
			std::cout << "[ ";
			typename Size< TFragType >::Type dim = 0;
			if( me._end )
				std::cout << rightPosition( *me._frag, dim );
			else
				std::cout << leftPosition( *me._frag, dim );
			++dim;
			while( dim != dimension( *me._frag ) )
			{
				if( me._end )
					std::cout << " , " << rightPosition( *me._frag, dim );
				else
					std::cout << " , " << leftPosition( *me._frag, dim );
				++dim;
			}
			std::cout << " ]" << std::endl;
		}

	#endif // SEQAN_CHAIN_DEBUG_

		WrapperPoint_( )
			: _frag( NULL )
			, _key( 0 )
			, _meta( NULL )
			, _end( false )
		{
		}


		WrapperPoint_( const WrapperPoint_ & old )
			: _frag( old._frag )
			, _key( old._key )
			, _meta( old._meta )
			, _end( old._end )
		{
		}


		WrapperPoint_( MetaFragment_< TFragType > & meta, 
						bool end )
			: _frag( &_getFrag( meta ) )
			, _key( end ? rightPosition( _getFrag( meta ), dimension( _getFrag( meta ) ) - 1 ) : leftPosition( _getFrag( meta ), dimension( _getFrag( meta ) ) - 1 ) )
			, _meta( & meta )
			, _end( end )
		{}

		template< typename TKey >
		WrapperPoint_( MetaFragment_< TFragType > & meta,
						TKey key,
						bool end )
			: _frag( &_getFrag( meta ) )
			, _key( key )
			, _meta( & meta )
			, _end( end )
		{}


		WrapperPoint_ & operator=( const WrapperPoint_ & old )
		{
			if ( this == &old ) 
				return *this;
			_key = old._key;
			_end = old._end;
			_frag = old._frag;
			_meta = old._meta;
			return *this;
		}

	};

}
#endif // SEQAN_H_...
