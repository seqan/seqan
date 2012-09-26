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

#ifndef SEQAN_HEADER_BLAST_STREAM_HSP_ITERATOR_H
#define SEQAN_HEADER_BLAST_STREAM_HSP_ITERATOR_H


namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Blast Hsp Iterator
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////



////**
//.Spec.HspIterator:
//..cat:Blast
//..summary:Hsp iterator for @Class.BlastHit@.
//..signature:Iterator<TBlastHit, HspIterator>
//..param.TBlastHit:A Blast hit.
//...type:Class.BlastHit
//..general:Class.Iter
//*/
template<typename TBlastHsp, typename TFile>
class Iter<BlastHit<TBlastHsp,StreamReport<TFile> >, StreamBlastIterator<HspIterator> > 
{
public:

	typedef BlastHit<TBlastHsp,StreamReport<TFile> > TBlastHit_;
	typedef typename Position<TFile>::Type TPosition_;

	TBlastHsp data_hsp;
	TBlastHit_* data_host;
	TPosition_ data_pos, data_next_pos, data_hsp_begin_pos;
	bool data_at_end;


	Iter()	
	{
	data_at_end = false;
	}
	
	Iter(TBlastHit_ & blast) 
	{
	SEQAN_CHECKPOINT
		data_host = &blast; 
		data_pos = blast.first_hsp_pos;
		data_next_pos = data_pos;
		data_hsp_begin_pos = (TPosition_) 0;
		data_at_end = false;
	}

	Iter(Iter const& other): 
		data_host(other.data_host), 
		data_pos(other.data_pos), 
		data_next_pos(other.data_next_pos), 
		data_at_end(other.data_at_end),
		data_hsp(other.data_hsp),
		data_hsp_begin_pos(other.data_hsp_begin_pos) 
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
		data_at_end = other.data_at_end;
		data_hsp = other.data_hsp;
		data_next_pos = other.data_next_pos;
		return *this;
	}
//____________________________________________________________________________
};

//////////////////////////////////////////////////////////////////////////////
// Blast StreamHspIterator - Metafunctions
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////


template<typename TBlastHsp, typename TFile>
struct Iterator<BlastHit<TBlastHsp,StreamReport<TFile> >, HspIterator>
{	
	typedef Iter<BlastHit<TBlastHsp,StreamReport<TFile> >, StreamBlastIterator<HspIterator> > Type;
};

template<typename TBlastHsp, typename TFile>
struct Iterator<BlastHit<TBlastHsp,StreamReport<TFile> > const, HspIterator>
{	
	typedef Iter<BlastHit<TBlastHsp,StreamReport<TFile> >, StreamBlastIterator<HspIterator> > Type;
};


//////////////////////////////////////////////////////////////////////////////

template<typename TBlastHsp, typename TFile>
struct Iterator<BlastReport<TBlastHsp,StreamReport<TFile> >, HspIterator>
{	
	typedef Iter<typename Hit<BlastReport<TBlastHsp,StreamReport<TFile> > >::Type, StreamBlastIterator<HspIterator> > Type;
};

template<typename TBlastHsp, typename TFile>
struct Iterator<BlastReport<TBlastHsp,StreamReport<TFile> > const, HspIterator>
{	
	typedef Iter<typename Hit<BlastReport<TBlastHsp,StreamReport<TFile> > >::Type, StreamBlastIterator<HspIterator> > Type;
};


//////////////////////////////////////////////////////////////////////////////
template<typename TBlastHsp, typename TFile>
struct Value<Iter<BlastHit<TBlastHsp,StreamReport<TFile> >, StreamBlastIterator<HspIterator> > >
{
	typedef TBlastHsp Type;
};

template<typename TBlastHsp, typename TFile>
struct Value<Iter<BlastHit<TBlastHsp,StreamReport<TFile> > const, StreamBlastIterator<HspIterator> > >
{
	typedef TBlastHsp Type;
};









///.Metafunction.Host.param.T.type:Class.BlastHit
template<typename TBlastHit>
struct Host<Iter<TBlastHit, StreamBlastIterator<HspIterator> > >
{	
	typedef TBlastHit Type;
};


//////////////////////////////////////////////////////////////////////////////

template<typename TBlastHit>
struct Reference<Iter<TBlastHit, StreamBlastIterator<HspIterator> > >
{
	typedef typename Value<Iter<TBlastHit, StreamBlastIterator<HspIterator> > >::Type& Type;
};

template<typename TBlastHit>
struct Reference<Iter<TBlastHit const, StreamBlastIterator<HspIterator> > >
{
	typedef typename Value<Iter<TBlastHit, StreamBlastIterator<HspIterator> > >::Type& Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TBlastHit>
struct GetValue<Iter<TBlastHit, StreamBlastIterator<HspIterator> > >
{
	typedef typename Value<Iter<TBlastHit, StreamBlastIterator<HspIterator> > >::Type Type;
};

template<typename TBlastHit>
struct GetValue<Iter<TBlastHit const, StreamBlastIterator<HspIterator> > >
{
	typedef typename Value<Iter<TBlastHit const, StreamBlastIterator<HspIterator> > >::Type Type;
};

//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
// Blast StreamBlastIterator<HspIterator> - FUNCTIONS
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

template<typename TBlastHit, typename TFile>
inline bool
atBegin(TFile &,
		Iter<TBlastHit, StreamBlastIterator<HspIterator> >& it)
{
SEQAN_CHECKPOINT
	return (it.data_pos == it.data_host->first_hsp_pos);	
}

//////////////////////////////////////////////////////////////////////////////



template<typename TBlastHit, typename TFile>
inline void
goBegin(TFile &,
		Iter<TBlastHit, StreamBlastIterator<HspIterator> >& it)
{
	SEQAN_CHECKPOINT
	it.data_pos = it.data_host->first_hsp_pos;
	it.data_next_pos = it.data_pos;
	it.data_at_end = false;
}



//////////////////////////////////////////////////////////////////////////////

template<typename TBlastHit, typename TFile>
inline void
goNext(TFile & file,
	   Iter<TBlastHit, StreamBlastIterator<HspIterator> >& it)
{
	SEQAN_CHECKPOINT
	if (!atEnd(file,it)) 
	{
		if(it.data_pos == it.data_next_pos)
			_getNextHspFilePos(file,it);
		if(it.data_pos == it.data_next_pos)
			it.data_at_end = true;
		else
            it.data_pos = it.data_next_pos;
	}
}

//////////////////////////////////////////////////////////////////////////////


//template<typename TBlastHit>
//inline void
//goPrevious(Iter<TBlastHit, StreamBlastIterator<HspIterator> >& it)
//{
//	if (!atBegin(it)) --it.data_pos;
//}
//

//////////////////////////////////////////////////////////////////////////////

template<typename TBlastHit>
inline bool
operator ==(Iter<TBlastHit, StreamBlastIterator<HspIterator> >& it1,
			Iter<TBlastHit, StreamBlastIterator<HspIterator> >& it2)
{
SEQAN_CHECKPOINT
	return (it1.data_pos==it2.data_pos && it1.data_host==it2.data_host);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TBlastHit>
inline bool
operator !=(Iter<TBlastHit, StreamBlastIterator<HspIterator> >& it1,
			Iter<TBlastHit, StreamBlastIterator<HspIterator> >& it2)
{
SEQAN_CHECKPOINT
	return (it1.data_pos!=it2.data_pos || it1.data_host!=it2.data_host);
}

//////////////////////////////////////////////////////////////////////////////


template<typename TBlastHit, typename TFile>
inline typename GetValue<Iter<TBlastHit, StreamBlastIterator<HspIterator> > >::Type
getValue(TFile & file,
		 Iter<TBlastHit, StreamBlastIterator<HspIterator> >& it)
{
SEQAN_CHECKPOINT
	if(it.data_pos != it.data_hsp_begin_pos)
	{
		_streamSeekG(file,it.data_pos);
		(it.data_host->data_host)->act_c = ' ';
		it.data_hsp_begin_pos = it.data_pos;
		typename Position<TFile>::Type pot_next_pos = _parseBlastHsp(file,(it.data_host->data_host)->act_c,it.data_hsp);
		if(pot_next_pos > it.data_pos)
			it.data_next_pos = pot_next_pos;
		//if(_parseBlastHit(it.data_host->strm,it.data_host->act_c,it.data_hsp) == it.data_pos)
		//	it.data_at_end = true;
	}
	return it.data_hsp;
}

//////////////////////////////////////////////////////////////////////////////


template<typename TBlastHit, typename TFile>
inline typename Reference<Iter<TBlastHit, StreamBlastIterator<HspIterator> > >::Type
value(TFile & file,
	  Iter<TBlastHit, StreamBlastIterator<HspIterator> >& it)
{
SEQAN_CHECKPOINT
	if(it.data_pos != it.data_hsp_begin_pos)
		it.data_hsp = getValue(file,it);
	return it.data_hsp;
}


//////////////////////////////////////////////////////////////////////////////
//
//
//template<typename TBlastHit>
//inline typename Host< typename Iterator< typename Host< Iter<TBlastHit, StreamBlastIterator<HspIterator> >::Type>::Type, StreamBlastIterator<HitIterator> >::Type >::Type const&
//hostReport(Iter<TBlastHit, StreamBlastIterator<HspIterator> >& it)
//{
//	return *(it.data_host->data_host);
//} 

//////////////////////////////////////////////////////////////////////////////



/**
.Function.hostHit:
..cat:Blast
..summary:The BlastHit this iterator is working on.
..signature:hostHit(it)
..param.it:An iterator.
...type:Spec.HspIterator
..returns:A pointer to the host BlastHit.
..include:seqan/blast.h
*/
template<typename TBlastHit>
inline typename Host<Iter<TBlastHit, StreamBlastIterator<HspIterator> > >::Type const&
hostHit(Iter<TBlastHit, StreamBlastIterator<HspIterator> >& it)
{
	SEQAN_CHECKPOINT
	return *it.data_host;
} 

//////////////////////////////////////////////////////////////////////////////


template<typename TBlastHit, typename TFile>
inline bool
atEnd(TFile &,
	  Iter<TBlastHit, StreamBlastIterator<HspIterator> >& it)
{
SEQAN_CHECKPOINT
//	return (it.data_last_pos != it.data_pos);	
	return it.data_at_end;	
}

//////////////////////////////////////////////////////////////////////////////

//template<typename TBlastHit>
//inline void
//goEnd(Iter<TBlastHit, StreamBlastIterator<HspIterator> >& it)
//{
//	it.data_pos = doof;
//}


//////////////////////////////////////////////////////////////////////////////

template<typename TBlastHsp, typename TFile>
inline void
_getNextHspFilePos(TFile & file,
				  Iter<BlastHit<TBlastHsp, StreamReport<TFile> >, StreamBlastIterator<HspIterator> >& it)
{
//IOREV
	typedef typename Position<TFile>::Type TPosition;

	_streamSeekG(file,it.data_pos);
	char c = 'e';
	(it.data_host->data_host)->act_c = c;

	_parseSkipWhitespace(file,c);
	_parseSkipLine(file,c);

	TPosition next_event_pos;
	bool last_hit = true;
	String<char> delim = ">";
	if(_parseUntilBeginLine(file,c,'>'))
	{
		last_hit = false;
		next_event_pos = _streamTellG(file);
		if((it.data_host->data_host)->next_report && (next_event_pos > (it.data_host->data_host)->next_report_pos))
			next_event_pos = (it.data_host->data_host)->next_report_pos;
	}
	_streamSeekG(file,it.data_pos);
	c = 'e';

	String<char> search = "Score";
	if(_parseUntilBeginLine(file,c,search,5))
	{
		if(!last_hit && (_streamTellG(file) > next_event_pos))
	        _streamSeekG(file,it.data_pos);
		else
            it.data_next_pos = _streamTellG(file);
	}//end hsp
	else
        _streamSeekG(file,it.data_pos);

}



}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
