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

#ifndef SEQAN_HEADER_BLAST_STREAM_HIT_ITERATOR_H
#define SEQAN_HEADER_BLAST_STREAM_HIT_ITERATOR_H


namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Blast Stream Hit Iterator
//////////////////////////////////////////////////////////////////////////////




////**
//.Spec.HitIterator:
//..cat:Blast
//..summary:Hit iterator for @Class.BlastReport@.
//..signature:Iterator<TBlastReport, HitIterator>
//..param.TBlastReport:A Blast report container.
//...type:Class.BlastReport
//..general:Class.Iter
//*/
template<typename TBlastHsp, typename TFile>
class Iter<BlastReport<TBlastHsp, StreamReport<TFile> >, StreamBlastIterator<HitIterator> > 
{
public:
	typedef BlastReport<TBlastHsp, StreamReport<TFile> > TBlastReport_;
	typedef typename Hit<TBlastReport_>::Type TBlastHit_;
	typedef typename Position<TFile>::Type TPosition_;

	TBlastHit_ data_hit;
	TBlastReport_* data_host;
	TPosition_ data_pos, data_next_pos;
	bool data_at_end;

	Iter()	
	{
		data_at_end = false;
	}
	
	Iter(TBlastReport_ & blast)  
	{
	SEQAN_CHECKPOINT
		data_host = &blast; 
		data_pos = blast.first_hit_pos;
		data_next_pos = data_pos;
		if(blast.hits_found)
			data_at_end = false;
		else
			data_at_end = true;
		data_hit.data_host = &blast;

	}

	Iter(Iter const& other) : 
		data_host(other.data_host), 
		data_pos(other.data_pos), 
		data_next_pos(other.data_next_pos), 
		data_at_end(other.data_at_end),
		data_hit(other.data_hit)
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
		data_next_pos = other.data_next_pos;
		data_at_end = other.data_at_end;
		data_hit = other.data_hit;
		return *this;
	}
//____________________________________________________________________________
};

//////////////////////////////////////////////////////////////////////////////
// Blast StreamHitIterator - Metafunctions
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Host.param.T.type:Class.BlastReport
template<typename TBlastReport>
struct Host<Iter<TBlastReport, StreamBlastIterator<HitIterator> > >
{	
	typedef TBlastReport Type;
};




template<typename TBlastHsp, typename TFile>
struct Iterator<BlastReport<TBlastHsp,StreamReport<TFile> >, HitIterator>
{	
	typedef Iter<BlastReport<TBlastHsp,StreamReport<TFile> >, StreamBlastIterator<HitIterator> > Type;
};

template<typename TBlastHsp, typename TFile>
struct Iterator<BlastReport<TBlastHsp,StreamReport<TFile> > const, HitIterator>
{	
	typedef Iter<BlastReport<TBlastHsp,StreamReport<TFile> > const, StreamBlastIterator<HitIterator> > Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TBlastHsp, typename TSpec>
struct Value<Iter<BlastReport<TBlastHsp,TSpec>, StreamBlastIterator<HitIterator> > >
{
	typedef typename Hit<BlastReport<TBlastHsp,TSpec> >::Type Type;
};

template<typename TBlastHsp, typename TSpec>
struct Value<Iter<BlastReport<TBlastHsp,TSpec> const, StreamBlastIterator<HitIterator> > >
{
	typedef typename Hit<BlastReport<TBlastHsp,TSpec> const>::Type Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TBlastReport>
struct Reference<Iter<TBlastReport, StreamBlastIterator<HitIterator> > >
{
	typedef typename Value<Iter<TBlastReport, StreamBlastIterator<HitIterator> > >::Type& Type;
};

template<typename TBlastReport>
struct Reference<Iter<TBlastReport const, StreamBlastIterator<HitIterator> > >
{
	typedef typename Value<Iter<TBlastReport const, StreamBlastIterator<HitIterator> > >::Type& Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TBlastReport>
struct GetValue<Iter<TBlastReport, StreamBlastIterator<HitIterator> > >
{
	typedef typename Value<Iter<TBlastReport, StreamBlastIterator<HitIterator> > >::Type Type;
};

template<typename TBlastReport>
struct GetValue<Iter<TBlastReport const, StreamBlastIterator<HitIterator> > >
{
	typedef typename Value<Iter<TBlastReport const, StreamBlastIterator<HitIterator> > >::Type Type;
};

//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
// Blast StreamBlastIterator<HitIterator> - FUNCTIONS
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////


/**
.Function.atBegin:
..cat:Blast
..signature:atBegin(file,it)
..param.file:A file stream.
..param.it:An iterator.
...type:Spec.HitIterator
...type:Spec.HspIterator
..include:seqan/blast.h
*/
template<typename TBlastReport, typename TFile>
inline bool
atBegin(TFile &,
		Iter<TBlastReport, StreamBlastIterator<HitIterator> >& it)
{
SEQAN_CHECKPOINT
	return (it.data_pos == it.data_host->first_hit_pos);	
}

//////////////////////////////////////////////////////////////////////////////



/**
.Function.goBegin:
..cat:Blast
..signature:goBegin(file,it)
..param.file:A file stream.
..param.it:An iterator.
...type:Spec.HitIterator
...type:Spec.HspIterator
..include:seqan/blast.h
*/
template<typename TBlastReport, typename TFile>
inline void
goBegin(TFile &,
		Iter<TBlastReport, StreamBlastIterator<HitIterator> >& it)
{
	SEQAN_CHECKPOINT
	it.data_pos = it.data_host->first_hit_pos;
	it.data_next_pos = it.data_pos;
	if(it.data_host->hits_found)
		it.data_at_end = false;
	else
		it.data_at_end = true;
}



//////////////////////////////////////////////////////////////////////////////

/**
.Function.goNext:
..cat:Blast
..signature:goNext(file,it)
..param.file:A file stream.
..param.it:An iterator.
...type:Spec.HitIterator
...type:Spec.HspIterator
..include:seqan/blast.h
*/
template<typename TBlastReport, typename TFile>
inline void
goNext(TFile & file, 
	   Iter<TBlastReport, StreamBlastIterator<HitIterator> >& it)
{
	SEQAN_CHECKPOINT
	if (!atEnd(file,it)) 
	{
		if(it.data_pos == it.data_next_pos)
			_getNextHitFilePos(file,it);
		if(it.data_pos == it.data_next_pos)
			it.data_at_end = true;
		else
            it.data_pos = it.data_next_pos;
	}
}

///////////////////////////////////////////////////////////////////////

//template<typename TBlastReport>
//inline void
//goPrevious(Iter<TBlastReport, StreamBlastIterator<HitIterator> >& it)
//{
//	if (!atBegin(it)) --it.data_pos;
//}
//
////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
//
//template<typename TBlastReport>
//inline Iter<TBlastReport, StreamBlastIterator<HitIterator> >
//operator --(Iter<TBlastReport, StreamBlastIterator<HitIterator> >& it, int)
//{
//	Iter<TBlastReport, StreamBlastIterator<HitIterator> > ret = it;
//	goPrevious(it);
//	return ret;
//}

//////////////////////////////////////////////////////////////////////////////

template<typename TBlastReport>
inline bool
operator ==(Iter<TBlastReport, StreamBlastIterator<HitIterator> >& it1,
			Iter<TBlastReport, StreamBlastIterator<HitIterator> >& it2)
{
SEQAN_CHECKPOINT
	return (it1.data_pos==it2.data_pos && it1.data_host==it2.data_host);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TBlastReport>
inline bool
operator !=(Iter<TBlastReport, StreamBlastIterator<HitIterator> >& it1,
			Iter<TBlastReport, StreamBlastIterator<HitIterator> >& it2)
{
SEQAN_CHECKPOINT
	return (it1.data_pos!=it2.data_pos || it1.data_host!=it2.data_host);
}

//////////////////////////////////////////////////////////////////////////////


/**
.Function.getValue:
..cat:Blast
..signature:getValue(file,it)
..param.file:A file stream.
..param.it:An iterator.
...type:Spec.HitIterator
...type:Spec.HspIterator
..include:seqan/blast.h
*/
template<typename TBlastReport, typename TFile>
inline typename GetValue<Iter<TBlastReport, StreamBlastIterator<HitIterator> > >::Type
getValue(TFile & file,
		 Iter<TBlastReport, StreamBlastIterator<HitIterator> >& it)
{
SEQAN_CHECKPOINT
	if(it.data_pos != (it.data_hit).begin_pos)
	{
		_streamSeekG(file,it.data_pos);
		it.data_host->act_c = '>';
		_parseBlastHit(file,it.data_host->act_c,it.data_hit);
		//if(_parseBlastHit(file,it.data_host->act_c,it.data_hit) == it.data_pos)
		//	it.data_at_end = true;
	}
	return it.data_hit;
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.value:
..cat:Blast
..signature:value(file,it)
..param.file:A file stream.
..param.it:An iterator.
...type:Spec.HitIterator
...type:Spec.HspIterator
..include:seqan/blast.h
*/
template<typename TBlastReport, typename TFile>
inline typename Reference<Iter<TBlastReport, StreamBlastIterator<HitIterator> > >::Type
value(TFile & file,
	  Iter<TBlastReport, StreamBlastIterator<HitIterator> >& it)
{
SEQAN_CHECKPOINT
	if(it.data_pos != (it.data_hit).begin_pos)
		it.data_hit = getValue(file,it);
	return it.data_hit;
}

//////////////////////////////////////////////////////////////////////////////


template<typename TBlastReport>
inline typename Host<Iter<TBlastReport, StreamBlastIterator<HitIterator> > >::Type &
hostReport(Iter<TBlastReport, StreamBlastIterator<HitIterator> >& it)
{
	SEQAN_CHECKPOINT
	return *it.data_host;
} 

//////////////////////////////////////////////////////////////////////////////


/**
.Function.atEnd:
..cat:Blast
..signature:atEnd(file,it)
..param.file:A file stream.
..param.it:An iterator.
...type:Spec.HitIterator
...type:Spec.HspIterator
..include:seqan/blast.h
*/

template<typename TBlastReport, typename TFile>
inline bool
atEnd(TFile &,
	  Iter<TBlastReport, StreamBlastIterator<HitIterator> >& it)
{
SEQAN_CHECKPOINT
//	return (it.data_last_pos != it.data_pos);	
	return it.data_at_end;	
}

//////////////////////////////////////////////////////////////////////////////

//template<typename TBlastReport>
//inline void
//goEnd(Iter<TBlastReport, StreamBlastIterator<HitIterator> >& it)
//{
//	it.data_pos = doof;
//}


//////////////////////////////////////////////////////////////////////////////

template<typename TBlastHsp, typename TFile>
inline void
_getNextHitFilePos(TFile & file,
				  Iter<BlastReport<TBlastHsp, StreamReport<TFile> >, StreamBlastIterator<HitIterator> >& it)
{
//IOREV maybe think about extra blast file format in file/
	typedef typename Position<TFile>::Type TPosition;

	_streamSeekG(file,it.data_pos);
	char c = '>';
	it.data_host->act_c = c;

	_parseSkipWhitespace(file,c);
	_parseSkipLine(file,c);

//	String<char> delim = ">";
	if(_parseUntilBeginLine(file,c,'>'))
	{
		TPosition event_pos = _streamTellG(file);
		if(c=='>' && ( !it.data_host->next_report || (it.data_host->next_report && event_pos < it.data_host->next_report_pos)))
			it.data_next_pos = event_pos;
		else
	        _streamSeekG(file,it.data_host->next_report_pos);

	}//end hit
	else
        _streamSeekG(file,it.data_pos);

}



}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
