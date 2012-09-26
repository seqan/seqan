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


/*	Hendrik Woehrle
*
*	Deferred Skip List Datastructure
*
*	ClassPool -
*
*	Memory pool allocator for multiple objects at once
*
*/

#ifndef SEQAN_HEADER_SKIP_POOL_ALLOC_H
#define SEQAN_HEADER_SKIP_POOL_ALLOC_H

namespace seqan
{

struct Unlimited
{};

struct Limited
{};

/*DISABLED
.Spec.Class Pool Allocator:
..cat:Allocators
..general:Class.Allocator
..summary:Allocator that pools memory blocks for a specific class.
..signature:Allocator< ClassPool<Class, Type, ParentAllocator> >
..param.Class:The class.
..param.Type:A specialization. The Class Pool Allocator can manage either slices of a fixed size or of 
..param.ParentAllocator:An allocator that is by the pool allocator used to allocate memory.
...default:@Spec.Simple Allocator@
...note:The multi pool allocator only supports @Function.clear@ if this function is also implemented for $ParentAllocator$.
..remarks:A pool allocator allocates several memory blocks at once. 
Freed blocks are not immediately deallocated but recycled in subsequential allocations.
This way, the number of calls to the heap allocator is reduced and that speeds up memory management. The Class Pool Allocator
pools only memory blocks of a fixed size, namely $sizeof(Class)$. The Class Pool Allocator only pools the memory, constructor and destructor have to be called manually.
..include:seqan/chaining.h
*/

template< typename TClass, typename TType, typename TParentAllocator = SimpleAllocator >
struct ClassPool;

template< typename TClass, typename TSpec, typename TParentAlloc >
void clear( Allocator< ClassPool< TClass, TSpec, TParentAlloc > > & me );

template< typename TClass, typename TSpec, typename TParentAlloc, typename TSize >
void allocate( Allocator< ClassPool< TClass, TSpec, TParentAlloc > > & me, 
				TClass *& dest, 
				TSize number );

template< typename TClass, typename TSpec, typename TParentAlloc, typename TSize >
void deallocate( Allocator< ClassPool< TClass, TSpec, TParentAlloc > > & me, 
				 TClass * location );

template< typename TClass, typename TSpec, typename TParentAlloc >
TParentAlloc & 
parentAlloc( Allocator< ClassPool< TClass, TSpec, TParentAlloc > > & me );

template< typename TClass, typename TSpec, typename TParentAlloc >
void
setParentAlloc( Allocator< ClassPool< TClass, TSpec, TParentAlloc > > & me );


template< typename TClass, typename TSpec, typename TParentAlloc > inline
TClass * 
_getNextBlock( Allocator< ClassPool< TClass, TSpec, TParentAlloc > > & /*me*/,
				TClass & block )
{
	return _getNext( block );
}

template< typename TClass, typename TSpec, typename TParentAlloc > inline
void 
_setNextBlock( Allocator< ClassPool< TClass, TSpec, TParentAlloc > > & /*me*/,
				TClass & dest,
				TClass * block )
{
	_setNext( dest, block );
}


//********************************** Partially specialised Template Class of the memory allocator for SkipBaseElement
template< typename TClass, typename TParentAllocator  >
struct Allocator< ClassPool< TClass, Unlimited, TParentAllocator > >
{
		// Last free block in queue, which can be used again
	TClass * _freeBlock;
		// Pointer to end of "used" memory Block
	TClass * _end;
		// Pointer to single end of reserved memory
	TClass * _terminal_end;
		// size of blocks
	typename Size< TClass >::Type _blockSize;
		// parent Allocator
	Holder<TParentAllocator> data_parent_allocator;
		

	Allocator(Allocator const &):
		_freeBlock(0),
		_end(0),
		_terminal_end(0),
		_blockSize(0)
	{
	}

	inline Allocator &
	operator = (Allocator const &)
	{
		clear(*this);
		return *this;
	}


	Allocator( typename Size< TClass >::Type numElements = 1000 )
	{
		TClass * block;
		_blockSize = numElements > 20 ? numElements : 20;
		allocate( parentAllocator( *this ), block, _blockSize );
		_end = block;
		_terminal_end = _end + _blockSize;
		_freeBlock = NULL;
	}

	~Allocator(void)
	{
		clear( *this );
	}

};

	template< typename TClass, typename TParentAlloc >
	inline 
	TParentAlloc &
	parentAllocator( Allocator< ClassPool< TClass, Unlimited, TParentAlloc > > & me )
	{
		return value( me.data_parent_allocator );
	}

	template< typename TClass, typename TParentAlloc >
	inline 
	void
	setParentAllocator( Allocator< ClassPool< TClass, Unlimited, TParentAlloc > > & me,
					   TParentAlloc & alloc_)
	{
		setValue( me.data_parent_allocator, alloc_ );
	}


	template< typename TClass, typename TParentAlloc, typename TSize >
	void allocate( Allocator< ClassPool< TClass, Unlimited, TParentAlloc > > & me, 
					TClass *& dest, 
					TSize number )
	{
		SEQAN_ASSERT_LEQ_MSG(number, static_cast<TSize>(me._blockSize), "tried to allocate more elements than available in block");
		SEQAN_ASSERT_GT_MSG(number, static_cast<TSize>(0), "tried to allocate 0 elements");
		if( number == 1 )
		{
			if( me._end < me._terminal_end )
			{		// enough place available in current memory pool
				dest = me._end;
				++me._end;
				return;
			}
			else{
				allocate( parentAllocator( me ), me._end, me._blockSize );
				me._terminal_end = me._end + me._blockSize;
				dest = me._end;
				++me._end;
				return;
			}
		}
		else
		{
			if( me._end + number < me._terminal_end )
			{		// enough place available in current memory pool
				dest = me._end;
				me._end += number;
				return;
			}
			else{		
					// not enough memory in current pool avaiable: allocate new, adjust pointers
				if( me._end < me._terminal_end ){
					_setNext( * me._end, me._freeBlock );
					++me._end;
					while( me._end < me._terminal_end){
						_setNext( *me._end, me._end );
						++me._end;
					}
					me._freeBlock = me._end-1;
				}
				allocate( parentAllocator( me ), me._end, me._blockSize );
				me._terminal_end = me._end + me._blockSize;
				dest = me._end;
				me._end+=number;
				return;
			}
		}	
	}

	
	template< typename TClass, typename TParentAlloc, typename TSize > inline
	void deallocate( Allocator< ClassPool< TClass, Unlimited, TParentAlloc > > & me, 
					 TClass * location,
					 TSize /*count*/ )
	{
		SEQAN_ASSERT_MSG(location != NULL, "Tried to free NULL-pointer");
		_setNext( *location, me._freeBlock );
		me._freeBlock = location;
	}


	template< typename TClass, typename TParentAlloc > inline
	void
	clear( Allocator< ClassPool< TClass, Unlimited, TParentAlloc > > & me )
	{
	SEQAN_CHECKPOINT

		me._freeBlock = NULL;
		me._end = NULL;
		me._terminal_end = NULL;
		me._blockSize = 0;

		clear( parentAllocator( me ) );
	}


//***************************** Partially specialised Template Class of the memory allocator for SkipElement
template< typename TClass, typename TParentAllocator  >
struct Allocator< ClassPool< TClass, Limited, TParentAllocator > >
{
		// Map for freed Blocks, which can be used again
	TClass * _freeBlocks[ sizeof( typename Size< TClass >::Type ) * 8];
		// Pointer to end of "used" memory Block
	TClass * _end;
		// Pointer to single end of reserved memory
	TClass * _terminal_end;
		// size of Elements
	typename Size< TClass >::Type _blockSize;
		// parent Allocator
	Holder<TParentAllocator> data_parent_allocator;


public:

	Allocator( typename Size< TClass >::Type numElements = 1000 )
	{
		_blockSize = numElements > 32 ? numElements : 32;
		allocate( parentAllocator( *this ), _end, _blockSize );
		_terminal_end = _end + _blockSize;
		for( typename Size< TClass >::Type i = 0; i < sizeof( typename Size< TClass >::Type ) * 8; i++ )
			_freeBlocks[i] = NULL;
	}

	~Allocator(void)
	{
		clear( *this );
	}

};

	template< typename TClass, typename TParentAlloc >
	inline 
	TParentAlloc &
	parentAllocator( Allocator< ClassPool< TClass, Limited, TParentAlloc > > & me )
	{
		return value(me.data_parent_allocator);
	}

	template< typename TClass, typename TParentAlloc >
	inline 
	void
	setParentAllocator( Allocator< ClassPool< TClass, Limited, TParentAlloc > > & me,
					   TParentAlloc & alloc_)
	{
		setValue( me.data_parent_allocator, alloc_ );
	}

	template< typename TClass, typename TParentAlloc, typename TSize >
	void allocate( Allocator< ClassPool< TClass, Limited, TParentAlloc > > & me, 
					TClass *& dest, 
					TSize number )
	{
		SEQAN_ASSERT_LEQ_MSG(number, static_cast<TSize>(me._blockSize), "tried to allocate more elements than available in block");
		SEQAN_ASSERT_NEQ(number, static_cast<TSize>(0));
			// recycle old memory block
		if( me._freeBlocks[ number-1 ] != NULL )
		{
			dest = me._freeBlocks[number-1];
			me._freeBlocks[number-1] = _getNext( *dest );
			return;
		}
			// use in memory from current block
		else if( me._end + number < me._terminal_end )
		{
			dest = me._end;
			me._end += number;
			return;
		}
		else {	// allocate new memory block
			TClass * block;
			allocate( parentAllocator( me ), block, me._blockSize );
			typename Size< TClass >::Type rest_space = me._terminal_end - me._end;
			if( rest_space != 0 )
			{
				_setNextBlock( me, *me._end, me._freeBlocks[ rest_space - 1 ] );
				me._freeBlocks[ rest_space - 1 ] = me._end;
			}
			me._end = block;
			me._terminal_end = me._end + me._blockSize;
			dest = me._end;
			me._end += number;
		}
	}

	template< typename TClass, typename TParentAlloc, typename TSize >
	void 
	deallocate( Allocator< ClassPool< TClass, Limited, TParentAlloc > > & me,
				TClass * location,
				TSize count )
	{
		SEQAN_ASSERT_NEQ(count, static_cast<TSize>(0));
		_setNextBlock( me, *location, me._freeBlocks[ count - 1 ] );
		me._freeBlocks[ count - 1 ] = location;
	}

	template< typename TClass, typename TParentAlloc >
	void
	clear( Allocator< ClassPool< TClass, Limited, TParentAlloc > > & me )
	{
	SEQAN_CHECKPOINT

		for( size_t i = 0; i < sizeof( typename Size< TClass >::Type ) * 8; ++i )
			me._freeBlocks[i] = NULL;
		me._end = NULL;
		me._terminal_end = NULL;
		me._blockSize = 0;

		clear( parentAllocator( me ) );
	}

} // namespace seqan

#endif // SKIPPOOLALLOC_H
