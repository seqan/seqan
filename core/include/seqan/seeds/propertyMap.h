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

#ifndef SEQAN_HEADER_PropertyMap_H
#define SEQAN_HEADER_PropertyMap_H


//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN
{

/**
.Class.PropertyMap:
..cat:Seed Handling
..summary:Class used to save additional data together with a MemoryManager.
..signature:PropertyMap<TValue, Block<size> >
..param.TValue:The value type, that is the type of the items/characters stored in the manager.
...remarks:Use @Metafunction.Value@ to get the value type for a given class.
..param.size:A positive integer that specifies the number of values in each allocated block.
...remarks: Size has to be a power of 2, e.g., 1024.
..include:seqan/seeds.h
*/
template<typename TValue, typename TSPec> 
class PropertyMap;

template<typename TValue, unsigned int SPACE>
class PropertyMap<TValue, Block<SPACE> > 
{
	typedef String<TValue, Array<SPACE> >				TBlock;
	typedef TBlock*										PBlock;
	typedef Allocator< SinglePool<sizeof(TBlock)> >	TAllocator;

public:

	typedef typename Iterator<TBlock, Standard>::Type	TBlockIter;
	typedef String<PBlock>								TBlockTable;

	TBlockTable		blocks;
	TBlockIter		blockFirst, blockLast;	// current block boundaries
	TBlockIter		lastValue;				// pointer to top value
	TAllocator		alloc;
	unsigned int exponent;
	//____________________________________________________________________________
	      
	public:
	PropertyMap():
		blockFirst(TBlockIter()),
		blockLast(TBlockIter()),
		lastValue(TBlockIter()) 
	{
		SEQAN_CHECKPOINT
		unsigned int x = SPACE;
		exponent = 0;
		while (x != 1){
			x >>=1;
			++exponent;
		}
	}

	~PropertyMap() 
	{
		clear(*this);
	}

	//____________________________________________________________________________



	public:
		template<typename TPos>
		inline typename Reference<PropertyMap>::Type 
			operator[] (TPos pos) 
		{
		SEQAN_CHECKPOINT	
			return value(*this, pos);
		}

		template<typename TPos>
		inline typename Reference<PropertyMap const>::Type 
			operator[] (TPos pos) const 
		{
		SEQAN_CHECKPOINT
			return value(*this, pos);
		}
	};


	template<typename TValue, unsigned int SPACE>
	struct DefaultOverflowImplicit< PropertyMap<TValue, Block<SPACE> > >
	{
		typedef Generous Type;
	};


///.Metafunction.Value.param.T.type:Class.PropertyMap
///.Metafunction.Value.class:Class.PropertyMap

template <typename TValue, unsigned int SPACE>
struct Value<PropertyMap<TValue,Block<SPACE> > >
{
	typedef TValue Type;
};


template<typename TValue, unsigned int SPACE, typename TSource>
inline void 
assign(
	PropertyMap<TValue, Block<SPACE> >& target, 
	TSource const& source) 
{
	SEQAN_CHECKPOINT
	clear(target);
	raiseMemory(target,length(source));
	for (int i = 0; i < length(source);++i){
		target[i] = source[i];
	}

	target.lastValue = &target[length(source)-1];

	target.pEmptySpace = 0;
	const void* pTmpSource = source.pEmptySpace;
	void* pTmpTarget = &target.pEmptySpace;
	while (pTmpSource != 0){
		typename Size<String<TValue, Block<SPACE> > >::Type pos = _getPosition(source,pTmpSource);
		*((void**)pTmpTarget) = &target[pos];
		pTmpSource = (TValue *)source[pos];
		pTmpTarget = &target[pos];
	}
}


template<typename TValue, unsigned int SPACE, typename TPos>
inline typename Reference<String<TValue, Block<SPACE> > >::Type 
value(
	PropertyMap<TValue, Block<SPACE> >& stack, 
	TPos const pos) 
{
	SEQAN_CHECKPOINT
	return value(*(stack.blocks[pos>>stack.exponent]), pos % SPACE);
}


template<typename TValue, unsigned int SPACE, typename TPos>
inline typename Reference<String<TValue, Block<SPACE> > >::Type 
value(
	PropertyMap<TValue, Block<SPACE> > const& stack, 
	TPos const pos) 
{
	SEQAN_CHECKPOINT
	return value(*(stack.blocks[pos>>stack.exponent]), pos % SPACE);
}


template<typename TValue, unsigned int SPACE>
inline void 
clear(PropertyMap<TValue, Block<SPACE> >& me)
{
SEQAN_CHECKPOINT
	typedef String<TValue, Block<SPACE>	>			TBlockString;
	typedef typename TBlockString::TBlockTable		TBlockTable;
	typedef typename Iterator<TBlockTable, Standard>::Type	TIter;
	
	TIter it = begin(me.blocks), itEnd = end(me.blocks);
	while (it != itEnd) {
		deallocate(me.alloc, *it, 1);
		++it;
	}
	clear(me.blocks);
	me.lastValue = me.blockLast = typename TBlockString::TBlockIter();
}


template<typename TValue, unsigned int SPACE>
inline typename Size<String<TValue, Block<SPACE> > >::Type
length(PropertyMap<TValue, Block<SPACE> > const & me) 
{
	SEQAN_CHECKPOINT
	if (length(me.blocks))
		return (length(me.blocks) - 1) * SPACE + (me.lastValue - me.blockFirst) + 1;
	else
		return 0;
}


template<typename TValue, unsigned int SPACE>
inline typename Size<String<TValue, Block<SPACE> > >::Type
capacity(PropertyMap<TValue, Block<SPACE> > const & me) 
{
SEQAN_CHECKPOINT
	if (length(me.blocks))
		return length(me.blocks) * SPACE;
	else
		return 0;
}

template<typename TValue, unsigned int SPACE> 
void
raiseMemory(PropertyMap<TValue,Block<SPACE> > &manager)
{
	SEQAN_CHECKPOINT
	typename Size< String<TValue, Block<SPACE> > >::Type last = length(manager.blocks);
	resize(manager.blocks, last + 1, typename DefaultOverflowImplicit<String<TValue, Block<SPACE> > >::Type());
	allocate(manager.alloc, manager.blocks[last], 1);
	manager.lastValue = manager.blockFirst = begin(*manager.blocks[last]);
	manager.blockLast = (manager.blockFirst + (SPACE - 1));
}


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
