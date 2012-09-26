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

#ifndef SEQAN_HEADER_CHAINPOINT_H
#define SEQAN_HEADER_CHAINPOINT_H

namespace seqan
{


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//			Class ChainPoint_
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*DISABLED
.Class.ChainPoint_:
..summary:Basic data structure to represent the end-point of a fragment for use in a RMT
..cat:Chaining
..signature:ChainPoint_< TFragType, TSpec >
..param.TFragType:Type of the class which represents the limits (multidimensional point) of an fragment.
..param.TSpec:Spec of the ChainPoint_.
..include:seqan/chaining.h
*/
	template< typename TFragType, typename TSpec >
	struct Spec< ChainPoint_< TFragType, TSpec > >
	{
		typedef TSpec Type;
	};

	template< typename TFragType, typename TSpec, typename TSize > inline
	typename Key< TFragType >::Type
	key( ChainPoint_< TFragType, TSpec > & me,
			TSize dim )
	{
        SEQAN_ASSERT_MSG(me._coords != NULL, "point not initialized");
        SEQAN_ASSERT_LT_MSG((int)dim, (int)me._dim, "dim corrupted");
		return me._coords[ dim ];
	}

	template< typename TFragType, typename TSpec, typename TSize > inline
	typename Key< TFragType >::Type
	key( const ChainPoint_< TFragType, TSpec > & me,
			TSize dim )
	{
        SEQAN_ASSERT_NEQ_MSG(me._coords, NULL, "point not initialized");
        SEQAN_ASSERT_LT((int)dim, (int)me._dim);
		return me._coords[ dim ];
	}

	template< typename TFragType, unsigned int ISize, typename TSize > inline
	typename Key< TFragType >::Type
	key( ChainPoint_< TFragType, Array< ISize> > & me,
			TSize dim )
	{
	    SEQAN_ASSERT_LT((int)dim, (int)ISize);
	    return me._coords[ dim ];
    }

	template< typename TFragType, unsigned ISize, typename TSize > inline
	typename Key< TFragType >::Type
	key( ChainPoint_< TFragType, Array< ISize> > const & me,
			TSize dim )
	{
	    SEQAN_ASSERT_LT((int)dim, (int)ISize);
	    return me._coords[ dim ];
    }

	template< typename TFragType, typename TSpec, typename TSize, typename TKey > inline
	void
	setKey( ChainPoint_< TFragType, TSpec > & me,
			TSize dim,
			TKey val )
	{
        //SEQAN_ASSERT_LEQ(dim, me._dim);
		me._coords[ dim ] = val;
	}

	template< typename TFragType, typename TSpec, typename TSize, typename TKey > inline 
	void 
	setKey( const ChainPoint_< TFragType, TSpec > & me,
				TSize dim,
				TKey val )
	{
        //SEQAN_ASSERT_LEQ(dim, me._dim);
		me._coords[ dim ] = val;
	}

	template< typename TFragType, typename TSpec, typename TSize > inline 
	void 
	_incKey(  ChainPoint_< TFragType, TSpec > & me,
				TSize dim )
	{
        //SEQAN_ASSERT_LEQ(dim, me._dim);
		++me._coords[ dim ];
	}

		// get the related fragment
	template< typename TFragType, typename TSpec > inline
	TFragType & 
	_getFrag( ChainPoint_< TFragType, TSpec > & me )
	{
        SEQAN_ASSERT_NEQ(me._meta, NULL);
		return _getFrag( *me._meta );
	}

	template< typename TFragType, typename TSpec > inline
	TFragType & 
	_getFrag( const ChainPoint_< TFragType, TSpec > & me )
	{
        SEQAN_ASSERT_NEQ(me._meta, NULL);
		return _getFrag( *me._meta );
	}

		// get the related metainformation struct
	template< typename TFragType, typename TSpec > inline
	MetaFragment_< TFragType > & 
	_meta( ChainPoint_< TFragType, TSpec > & me )
	{
		return *me._meta;
	}

	template< typename TFragType, typename TSpec > inline
	MetaFragment_< TFragType > & 
	_meta( const ChainPoint_< TFragType, TSpec > & me )
	{
		return *me._meta;
	}

		// set the related metainformation struct
	template< typename TFragType, typename TSpec > inline
	void 
	_setMeta( ChainPoint_< TFragType, TSpec > & me,
				MetaFragment_< TFragType > & meta )
	{
		me._meta = meta;
	}

	template< typename TFragType, typename TSpec > inline
	void 
	_setMeta( const ChainPoint_< TFragType, TSpec > & me,
				MetaFragment_< TFragType > & meta )
	{
		me._meta = meta;
	}

	template< typename TFragType, typename TSpec > inline
	typename Size< TFragType >::Type
	dimension( ChainPoint_< TFragType, TSpec > & me )
	{
		return me._dim;
	}

	template< typename TFragType, unsigned int ISize > inline
	typename Size< TFragType >::Type
	dimension( ChainPoint_< TFragType, Array< ISize > > &)
	{
		return ISize;
	}

	template< typename TFragType, typename TSpec > inline
	typename Weight< TFragType >::Type
	priority( ChainPoint_< TFragType, TSpec > & me )
	{
		return me._prio;
	}

	template< typename TFragType, typename TSpec, typename TPrio > inline
	void
	setPriority( ChainPoint_< TFragType, TSpec > & me,
					TPrio prio )
	{
		me._prio = prio;
	}


	template< typename TFragType, typename TSpec >
	struct ChainPoint_
	{
		
		MetaFragment_< TFragType > * _meta;
		typename Key< TFragType >::Type * _coords;
		typename Size< TFragType >::Type _dim;
		typename Weight< TFragType >::Type _prio;

	public:

			// standard constructor for use in skiplist
		ChainPoint_(  )
		: _meta( NULL )
		, _coords( NULL )
		, _dim( 1 )
		, _prio( minValue< typename Weight< TFragType >::Type >() )
		{
			allocate( *this, _coords, _dim );
		}

		ChainPoint_( typename Size< TFragType >::Type dim )
		: _meta( NULL )
		, _coords( NULL )
		, _dim( dim )
		, _prio( minValue< typename Weight< TFragType >::Type >() )
		{
			allocate( *this, _coords, _dim );
		}

		ChainPoint_( MetaFragment_< TFragType > & meta, 
						bool begin = false )
		: _meta( &meta )
		, _coords( NULL )
		, _prio( minValue< typename Weight< TFragType >::Type >() )
		{
			TFragType * frag = &_getFrag( meta );
			_dim = dimension( *frag );
			allocate( *this, _coords, _dim );
			if( begin ){
				for( typename Size< TFragType >::Type i = 0; i < _dim; ++i )
				{
					_coords[ i ] = leftPosition( *frag, i );
				}
			}
			else{
				for( typename Size< TFragType >::Type i = 0; i < _dim; ++i )
				{
					_coords[ i ] = rightPosition( *frag, i );
				}
			}
		}

		template< typename TSize >
		ChainPoint_( MetaFragment_< TFragType > & meta,
						TSize dim, 
						bool begin = false )
		: _meta( &meta )
		, _coords( NULL )
		, _dim( dim )
		, _prio( minValue< typename Weight< TFragType >::Type >() )
		{
			TFragType * frag = &_getFrag( meta );
			allocate( *this, _coords, _dim );
			if( begin ){
				for( typename Size< TFragType >::Type i = 0; i < _dim; ++i )
				{
					_coords[ i ] = leftPosition( *frag, i );
				}
			}
			else{
				for( typename Size< TFragType >::Type i = 0; i < _dim; ++i )
				{
					_coords[ i ] = rightPosition( *frag, i );
				}
			}
		}

		ChainPoint_( typename Key< TFragType >::Type * coords,
						typename Size< TFragType >::Type dim,
						MetaFragment_< TFragType > * meta, 
						bool /*begin = false*/ )
		: _meta( meta )
		, _coords( NULL )
		, _dim( dim )
		, _prio( minValue< typename Weight< TFragType >::Type >() )
		{
			allocate( *this, _coords, _dim );
			for( typename Size< TFragType >::Type i = 0; i < _dim; ++i )
			{
				_coords[ i ] = coords[ i ];
			}
		}

		ChainPoint_( typename Key< TFragType >::Type * coords,
						typename Size< TFragType >::Type dim )
		: _meta( NULL )
		, _coords( NULL )
		, _dim( dim )
		, _prio( minValue< typename Weight< TFragType >::Type >() )
		{
			allocate( *this, _coords, _dim );
			for( typename Size< TFragType >::Type i = 0; i < _dim; ++i )
			{
				_coords[ i ] = coords[ i ];
			}
		}


		~ChainPoint_()
		{
			deallocate( *this, _coords, _dim );
			_coords = NULL;
			_meta = NULL;
			_dim = 0;
		}

		ChainPoint_( const ChainPoint_ & old )
		{
			_meta = old._meta;
			_dim = old._dim;
			_prio = old._prio;
			allocate( *this, _coords, _dim );
			for( typename Size< TFragType >::Type i = 0; i < _dim; ++i )
			{
				_coords[ i ] = old._coords[ i ];
			}
		}


		ChainPoint_ & operator=( const ChainPoint_ & old )
		{
			if ( this == &old ) 
				return *this;
			_meta = old._meta;
			if( _coords )
				deallocate( *this, _coords, _dim );
			_dim = old._dim;
			_prio = old._prio;
			allocate( *this, _coords, _dim );
			for( typename Size< TFragType >::Type i = 0; i < _dim; ++i )
			{
				_coords[ i ] = old._coords[ i ];
			}
			return *this;
		}

		friend inline
		void
		dump( ChainPoint_ & me )
		{
			std::cout << "[ ";
			typename Size< Seed< TFragType > >::Type dim = 0;
			std::cout << key( me, dim );
			++dim;
			while( dim != me._dim )
			{
				std::cout << " , " << key( me, dim );
				++dim;
			}
			std::cout << " ] "<< me._prio << std::endl;
		}


	};


	template< typename TFragType, unsigned int ISize>
	struct ChainPoint_< TFragType, Array< ISize > >
	{
		
		typename Key< TFragType >::Type _coords[ISize];
		MetaFragment_< TFragType > * _meta;
		typename Weight< TFragType >::Type _prio;

	public:

			// standard constructor for use in skiplist
		ChainPoint_(  )
		: _meta( NULL )
		, _prio( minValue< typename Weight< TFragType >::Type >() )
		{
		}

		ChainPoint_( typename Size< TFragType >::Type )
		: _meta( NULL )
		, _prio( minValue< typename Weight< TFragType >::Type >() )
		{
		}

		ChainPoint_( MetaFragment_< TFragType > & meta, 
						bool begin = false )
		: _meta( &meta )
		, _prio( minValue< typename Weight< TFragType >::Type >() )
		{
			TFragType * frag = &_getFrag( meta );
			if( begin ){
				for( typename Size< TFragType >::Type i = 0; i < ISize; ++i )
				{
					_coords[ i ] = leftPosition( *frag, i );
				}
			}
			else{
				for( typename Size< TFragType >::Type i = 0; i < ISize; ++i )
				{
					_coords[ i ] = rightPosition( *frag, i );
				}
			}
		}

		template< typename TSize >
		ChainPoint_( MetaFragment_< TFragType > & meta,
						TSize, 
						bool begin = false )
		: _meta( &meta )
		, _prio( minValue< typename Weight< TFragType >::Type >() )
		{
			TFragType * frag = &_getFrag( meta );
			if( begin ){
				for( typename Size< TFragType >::Type i = 0; i < ISize; ++i )
				{
					_coords[ i ] = leftPosition( *frag, i );
				}
			}
			else{
				for( typename Size< TFragType >::Type i = 0; i < ISize; ++i )
				{
					_coords[ i ] = rightPosition( *frag, i );
				}
			}
		}

		ChainPoint_( typename Key< TFragType >::Type * coords,
						typename Size< TFragType >::Type /*dim*/,
						MetaFragment_< TFragType > * meta, 
						bool /*begin = false*/ )
		: _meta( meta )
		, _prio( minValue< typename Weight< TFragType >::Type >() )
		{
			for( typename Size< TFragType >::Type i = 0; i < ISize; ++i )
			{
				_coords[ i ] = coords[ i ];
			}
		}

		ChainPoint_( typename Key< TFragType >::Type * coords,
						typename Size< TFragType >::Type /*dim*/ )
		: _meta( NULL )
		, _prio( minValue< typename Weight< TFragType >::Type >() )
		{
			for( typename Size< TFragType >::Type i = 0; i < ISize; ++i )
			{
				_coords[ i ] = coords[ i ];
			}
		}


		~ChainPoint_()
		{
			_meta = NULL;
		}

		ChainPoint_( const ChainPoint_ & old )
		{
			_meta = old._meta;
			_prio = old._prio;
			for( typename Size< TFragType >::Type i = 0; i < ISize; ++i )
			{
				_coords[ i ] = old._coords[ i ];
			}
		}


		ChainPoint_ & operator=( const ChainPoint_ & old )
		{
			if ( this == &old ) 
				return *this;
			_meta = old._meta;
			_prio = old._prio;
			for( typename Size< TFragType >::Type i = 0; i < ISize; ++i )
			{
				_coords[ i ] = old._coords[ i ];
			}
			return *this;
		}

		friend inline
		void
		dump( ChainPoint_ & me )
		{
			std::cout << "[ ";
			typename Size< Seed< TFragType > >::Type dim = 0;
			std::cout << key( me, dim );
			++dim;
			while( dim != ISize )
			{
				std::cout << " , " << key( me, dim );
				++dim;
			}
			std::cout << " ] "<< me._prio << std::endl;
		}


	};


}	// namespace

#endif // SEQAN_HEADER_CHAINPOINT_H
