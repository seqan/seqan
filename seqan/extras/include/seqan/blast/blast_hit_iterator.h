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

#ifndef SEQAN_HEADER_BLAST_HIT_ITERATOR_H
#define SEQAN_HEADER_BLAST_HIT_ITERATOR_H


namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Blast Hit Iterator
//////////////////////////////////////////////////////////////////////////////




/**
.Spec.HitIterator:
..cat:Blast
..summary:Hit iterator for @Class.BlastReport@.
..signature:Iterator<TBlastReport, HitIterator>
..param.TBlastReport:A Blast report.
...type:Class.BlastReport
..general:Class.Iter
..include:seqan/blast.h
*/
template<typename TBlastReport>
class Iter<TBlastReport, SimpleBlastIterator<HitIterator> > 
{
public:

	TBlastReport * data_host;
	unsigned int data_pos;

	Iter()	
	{
	}
	
	Iter(TBlastReport & blast) : 
		data_host(&blast), 
		data_pos(0) 
	{
	SEQAN_CHECKPOINT
	}

	Iter(Iter const& it) : 
		data_host(it.data_host), 
		data_pos(it.data_pos) 
	{
	SEQAN_CHECKPOINT
	}

	~Iter() 
	{
	SEQAN_CHECKPOINT
	}

	Iter const&	operator = (Iter const & other) 
	{
	SEQAN_CHECKPOINT
		if (this == &other) return *this;
		data_host = other.data_host;
		data_pos = other.data_pos;
		return *this;
	}
//____________________________________________________________________________
};


//////////////////////////////////////////////////////////////////////////////
// SimpleBlastIterator<TSpec> - Metafunctions
//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Host.param.T.type:Spec.HitIterator
///.Metafunction.Host.param.T.type:Spec.HspIterator
template<typename TBlastObject, typename TIteratorSpec>
struct Host<Iter<TBlastObject, SimpleBlastIterator<TIteratorSpec> > >
{	
	typedef TBlastObject Type;
};



//////////////////////////////////////////////////////////////////////////////

template<typename TBlastObject, typename TIteratorSpec>
struct Reference<Iter<TBlastObject, SimpleBlastIterator<TIteratorSpec> > >
{
	typedef typename Value<Iter<TBlastObject, SimpleBlastIterator<TIteratorSpec> > >::Type& Type;
};

template<typename TBlastObject, typename TIteratorSpec>
struct Reference<Iter<TBlastObject const, SimpleBlastIterator<TIteratorSpec> > >
{
	typedef typename Value<Iter<TBlastObject const, SimpleBlastIterator<TIteratorSpec> > >::Type& Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TBlastObject, typename TIteratorSpec>
struct GetValue<Iter<TBlastObject, SimpleBlastIterator<TIteratorSpec> > >
{
	typedef typename Value<Iter<TBlastObject, SimpleBlastIterator<TIteratorSpec> > >::Type Type;
};

template<typename TBlastObject, typename TIteratorSpec>
struct GetValue<Iter<TBlastObject const, SimpleBlastIterator<TIteratorSpec> > >
{
	typedef typename Value<Iter<TBlastObject const, SimpleBlastIterator<TIteratorSpec> > >::Type Type;
};



//////////////////////////////////////////////////////////////////////////////
// SimpleBlastIterator<HitIterator> - Metafunctions
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Iterator.param.T.type:Class.BlastReport

template<typename TBlastHsp, typename TInfoSpec>
struct Iterator<BlastReport<TBlastHsp,StoreReport<TInfoSpec> >, HitIterator>
{	
	typedef Iter<BlastReport<TBlastHsp,StoreReport<TInfoSpec> >, SimpleBlastIterator<HitIterator> > Type;
};

template<typename TBlastHsp, typename TInfoSpec>
struct Iterator<BlastReport<TBlastHsp,StoreReport<TInfoSpec> > const, HitIterator>
{	
	typedef Iter<BlastReport<TBlastHsp,StoreReport<TInfoSpec> > const, SimpleBlastIterator<HitIterator> > Type;
};




//////////////////////////////////////////////////////////////////////////////

template<typename TBlastHsp, typename TSpec>
struct Value<Iter<BlastReport<TBlastHsp,TSpec>, SimpleBlastIterator<HitIterator> > >
{
	typedef typename Hit<BlastReport<TBlastHsp,TSpec> >::Type Type;
};

template<typename TBlastHsp, typename TSpec>
struct Value<Iter<BlastReport<TBlastHsp,TSpec> const, SimpleBlastIterator<HitIterator> > >
{
	typedef typename Hit<BlastReport<TBlastHsp,TSpec> const>::Type Type;
};


//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
// Blast SimpleBlastIterator<TSpec> - FUNCTIONS
// TSpecs: HspIterator and HitIterator
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TBlastObject, typename TIteratorSpec>
inline typename Reference<Iter<TBlastObject, SimpleBlastIterator<TIteratorSpec> > >::Type
operator * (Iter<TBlastObject, SimpleBlastIterator<TIteratorSpec> >& it)
{
	SEQAN_CHECKPOINT
	return value(it);
}

//////////////////////////////////////////////////////////////////////////////


/**
.Function.atBegin:
..cat:Blast
..param.iterator:
...type:Spec.HitIterator
...type:Spec.HspIterator
..include:seqan/blast.h
*/
template<typename TBlastObject, typename TIteratorSpec>
inline bool
atBegin(Iter<TBlastObject, SimpleBlastIterator<TIteratorSpec> >& it)
{
SEQAN_CHECKPOINT
	return (it.data_pos == 0);	
}

//////////////////////////////////////////////////////////////////////////////



/**
.Function.goBegin:
..cat:Blast
..param.iterator:
...type:Spec.HitIterator
...type:Spec.HspIterator
..include:seqan/blast.h
*/
template<typename TBlastObject, typename TIteratorSpec>
inline void
goBegin(Iter<TBlastObject, SimpleBlastIterator<TIteratorSpec> >& it)
{
	SEQAN_CHECKPOINT
	it.data_pos = 0;
}



//////////////////////////////////////////////////////////////////////////////

/**
.Function.goNext:
..cat:Blast
..param.iterator
...type:Spec.HitIterator
...type:Spec.HspIterator
..include:seqan/blast.h
*/
template<typename TBlastObject, typename TIteratorSpec>
inline void
goNext(Iter<TBlastObject, SimpleBlastIterator<TIteratorSpec> >& it)
{
	SEQAN_CHECKPOINT
	if (!atEnd(it)) ++it.data_pos;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TBlastObject, typename TIteratorSpec>
inline Iter<TBlastObject, SimpleBlastIterator<TIteratorSpec> >&
operator ++(Iter<TBlastObject, SimpleBlastIterator<TIteratorSpec> >& it)
{
SEQAN_CHECKPOINT
	goNext(it);
	return it;
}



//////////////////////////////////////////////////////////////////////////////

/**
.Function.goPrevious:
..cat:Blast
..param.iterator:
...type:Spec.HitIterator
...type:Spec.HspIterator
..include:seqan/blast.h
*/

template<typename TBlastObject, typename TIteratorSpec>
inline void
goPrevious(Iter<TBlastObject, SimpleBlastIterator<TIteratorSpec> >& it)
{
	SEQAN_CHECKPOINT
	if (!atBegin(it)) --it.data_pos;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TBlastObject, typename TIteratorSpec>
inline Iter<TBlastObject, SimpleBlastIterator<TIteratorSpec> >&
operator --(Iter<TBlastObject, SimpleBlastIterator<TIteratorSpec> >& it)
{
SEQAN_CHECKPOINT
	goPrevious(it);
	return it;
}

//////////////////////////////////////////////////////////////////////////////

//template<typename TBlastObject, typename TIteratorSpec>
//inline Iter<TBlastObject, SimpleBlastIterator<TIteratorSpec> >
//operator --(Iter<TBlastObject, SimpleBlastIterator<TIteratorSpec> >& it, int)
//{
//	Iter<TBlastObject, SimpleBlastIterator<TIteratorSpec> > ret = it;
//	goPrevious(it);
//	return ret;
//}

//////////////////////////////////////////////////////////////////////////////

template<typename TBlastObject, typename TIteratorSpec>
inline bool
operator ==(Iter<TBlastObject, SimpleBlastIterator<TIteratorSpec> >& it1,
			Iter<TBlastObject, SimpleBlastIterator<TIteratorSpec> >& it2)
{
SEQAN_CHECKPOINT
	return (it1.data_pos==it2.data_pos);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TBlastObject, typename TIteratorSpec>
inline bool
operator !=(Iter<TBlastObject, SimpleBlastIterator<TIteratorSpec> >& it1,
			Iter<TBlastObject, SimpleBlastIterator<TIteratorSpec> >& it2)
{
SEQAN_CHECKPOINT
	return (it1.data_pos!=it2.data_pos);
}

//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
// Blast SimpleBlastIterator<HitIterator> - FUNCTIONS
//////////////////////////////////////////////////////////////////////////////

/**
.Function.getValue:
..cat:Blast
..param.object:
...type:Spec.HitIterator
...type:Spec.HspIterator
..include:seqan/blast.h
*/
template<typename TBlastReport>
inline typename GetValue<Iter<TBlastReport, SimpleBlastIterator<HitIterator> > >::Type
getValue(Iter<TBlastReport, SimpleBlastIterator<HitIterator> >& it)
{
SEQAN_CHECKPOINT
	return it.data_host->hits[it.data_pos];
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.value:
..cat:Blast
..param.object:
...type:Spec.HitIterator
...type:Spec.HspIterator
..include:seqan/blast.h
*/
template<typename TBlastReport>
inline typename Reference<Iter<TBlastReport, SimpleBlastIterator<HitIterator> > >::Type 
value(Iter<TBlastReport, SimpleBlastIterator<HitIterator> >& it)
{
SEQAN_CHECKPOINT
	return it.data_host->hits[it.data_pos];
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.hostReport:
..cat:Blast
..summary:The BlastReport this iterator is working on.
..signature:hostReport(it)
..param.it:An iterator.
...type:Spec.HitIterator
..returns:A pointer to the host BlastReport.
..include:seqan/blast.h
*/
template<typename TBlastReport>
inline typename Host<Iter<TBlastReport, SimpleBlastIterator<HitIterator> > >::Type/* const*/ &
hostReport(Iter<TBlastReport, SimpleBlastIterator<HitIterator> > & it)
{
	SEQAN_CHECKPOINT
	return *it.data_host;
} 

//////////////////////////////////////////////////////////////////////////////


/**
.Function.atEnd:
..cat:Blast
..param.iterator:
...type:Spec.HitIterator
...type:Spec.HspIterator
..include:seqan/blast.h
*/
template<typename TBlastReport>
inline bool
atEnd(Iter<TBlastReport, SimpleBlastIterator<HitIterator> >& it)
{
SEQAN_CHECKPOINT
	return (it.data_pos == length(it.data_host->hits));	
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.goEnd:
..cat:Blast
..param.iterator:
...type:Spec.HitIterator
...type:Spec.HspIterator
..include:seqan/blast.h
*/
template<typename TBlastReport>
inline void
goEnd(Iter<TBlastReport, SimpleBlastIterator<HitIterator> >& it)
{
	SEQAN_CHECKPOINT
	it.data_pos = length(it.data_host->hits);
}


//////////////////////////////////////////////////////////////////////////////




}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
