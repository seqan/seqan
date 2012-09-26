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

#ifndef SEQAN_HEADER_SEEDSET_ITERATOR_H
#define SEQAN_HEADER_SEEDSET_ITERATOR_H

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
// Iter
//////////////////////////////////////////////////////////////////////////////



template <typename TSeedSet>
class Iter<TSeedSet, SeedIterator>
{
public:
	typedef typename Size<TSeedSet>::Type TSize;
	typename std::set<TSize>::iterator data_ptr;
    typedef TSize TValue;
	TSeedSet *set;

	Iter()
	{
		this->set = 0;
	}

	Iter(TSeedSet &set)
	{
		this->set = &set;
	}

	Iter(TSeedSet const &set)
	{
		this->set = &set;
	}

	Iter(TSeedSet &set, Iter & other_):
		data_ptr(other_.data_ptr)
	{
		this->set = &set;
	}
		
	Iter(TSeedSet &set, typename std::set<TSize>::iterator other_)//:data_ptr(other_)
	{
		data_ptr = other_;
		this->set = &set;
	}
        
        
        
	Iter(TSeedSet &/*set*/, TSize &other_data_ptr):
		data_ptr(other_data_ptr)
	{
	}
	template <typename TSeedSet2>
	Iter(TSeedSet &set, Iter<TSeedSet2, SeedIterator> & other_):
		data_ptr(other_.data_ptr)
	{
		this->set = &set;
	}
	~Iter()
	{
	}

	Iter const &
	operator = (Iter const & other_)
	{
		this->data_ptr = other_.data_ptr;
                this->set = other_.set;
		return *this;
	}
	
        /*
        Iter const &
	operator = (TValue * other_data_ptr)
	{
		data_ptr = other_data_ptr;
		return *this;
	}*/
        
	template <typename TContainer2>
	Iter const &
	operator = (Iter<TContainer2, SeedIterator> const & other_)
	{
		this->data_ptr = other_.data_ptr;
                this->set = other_.set;
		return *this;
	}

	operator TValue * ()
	{
		return set->manager[*data_ptr];
	}

};


template <typename TSeedSet>
class Iter<TSeedSet const, SeedIterator>
{
public:
	typedef typename Size<TSeedSet>::Type TSize;
	typename std::set<TSize>::const_iterator data_ptr;
	TSeedSet const *set;

	Iter()
	{
		this->set = 0;
	}

	Iter(TSeedSet const &set)
	{
		this->set = &set;
	}

	Iter(TSeedSet &set, Iter const & other_):
		data_ptr(other_.data_ptr)
	{
		this->set = &set;
	}


	Iter(TSeedSet const &set, typename std::set<TSize>::const_iterator const & other_):
		data_ptr(other_)
	{
		//data_ptr = other_;
		this->set = &set;
	}

	Iter(TSeedSet const  &/*set*/, TSize * other_data_ptr):
		data_ptr(other_data_ptr)
	{
	}
	template <typename TSeedSet2>
	Iter(TSeedSet const &set, Iter<TSeedSet2 const, SeedIterator> const & other_):
		data_ptr(other_.data_ptr)
	{
		this->set = &set;
	}
		
		
	~Iter()
	{
	}

	Iter const &
	operator = (Iter const & other_)
	{
		this->data_ptr = other_.data_ptr;
		this->set = other_.set;
		return *this;
	}
	Iter const &
	operator = (TSize * other_data_ptr)
	{
		data_ptr = other_data_ptr;
		return *this;
	}
	template <typename TSeedSet2>
	Iter const &
	operator = (Iter<TSeedSet2 const, SeedIterator> const & other_)
	{
		this->data_ptr = other_.data_ptr;
		this->set = other_.set;
		return *this;
	}

};


template <typename TSeedSet>
struct Value<Iter<TSeedSet const, SeedIterator> >
{
	typedef typename Value<TSeedSet>::Type Type;
};



template<typename TSeedSet>
inline Iter<TSeedSet, SeedIterator >&
operator ++(Iter<TSeedSet, SeedIterator > &me)
{
SEQAN_CHECKPOINT
	++me.data_ptr;
	return me;
}

template<typename TSeedSet>
inline Iter<TSeedSet, SeedIterator >&
goNext(Iter<TSeedSet, SeedIterator > &it)
{
SEQAN_CHECKPOINT
	++(it.data_ptr);
	return it;
}

template<typename TSeedSet>
inline Iter<TSeedSet, SeedIterator >&
operator --(Iter<TSeedSet, SeedIterator > &it)
{
SEQAN_CHECKPOINT
	--(it.data_ptr);
	return it;
}

template<typename TSeedSet>
inline typename Reference<Iter<TSeedSet, SeedIterator> >::Type
value(Iter<TSeedSet, SeedIterator > &it)
{
SEQAN_CHECKPOINT
	return it.set->manager[*it.data_ptr];
}

template <typename TSeedSet>
inline typename Reference<Iter<TSeedSet, SeedIterator> const>::Type 
 value(Iter<TSeedSet, SeedIterator> const &me)
{
SEQAN_CHECKPOINT
	return me.set->manager[*me.data_ptr];
}



template <typename TSeedSet>
inline typename Reference<Iter<TSeedSet, SeedIterator> const>::Type 
operator * (Iter<TSeedSet, SeedIterator> & me)
{
SEQAN_CHECKPOINT
	return me.set->manager[*me.data_ptr];
}

template <typename TSeedSet>
inline typename Reference<Iter<TSeedSet const, SeedIterator> const>::Type 
operator * (Iter<TSeedSet const, SeedIterator> & me)
{
SEQAN_CHECKPOINT
	return me.set->manager[*me.data_ptr];
}

template<typename TSeedSet>
inline bool
operator !=(Iter<TSeedSet, SeedIterator > it1, Iter<TSeedSet, SeedIterator > it2)
{
	return it1.data_ptr != it2.data_ptr;
}

/**
.Function.seedScore:
..summary:Returns the score of a seed.
..cat:Seed Handling
..signature:seedScore(it);
..param.it: The seedSet iterator.
..return:Score of the seed.
..include:seqan/seeds.h
*/
template<typename TSeedSet>
inline typename ScoreType<typename ScoringScheme<TSeedSet>::Type>::Type &
seedScore(Iter<TSeedSet, SeedIterator > it)
{
SEQAN_CHECKPOINT
    return it.set->scoreMap[*it.data_ptr];
}



template<typename TSeedSet>
inline typename ScoreType<typename ScoringScheme<TSeedSet>::Type>::Type const&
seedScore(Iter<TSeedSet const, SeedIterator > &it)
{
SEQAN_CHECKPOINT
	return it.set->scoreMap[*it.data_ptr];
}

template<typename TSeedSet, typename TSize>
void
setScore(Iter<TSeedSet, SeedIterator > &it, TSize score)
{
SEQAN_CHECKPOINT
	it.set->scoreMap[*it.data_ptr] = score;
}
//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
