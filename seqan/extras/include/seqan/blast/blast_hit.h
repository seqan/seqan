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

#ifndef SEQAN_HEADER_BLAST_HIT_H
#define SEQAN_HEADER_BLAST_HIT_H


namespace SEQAN_NAMESPACE_MAIN
{

template<typename TBlastHsp, typename TStoreSpec>
class BlastHit;


/**
.Class.BlastHit:
..cat:Blast
..summary:Object for storing Blast hits. 
..signature:BlastHit<TBlastHsp, TSpec>  
..param.TBlastHsp:The type of HSPs that are stored.
..param.TSpec:The specializing type.
...type:Spec.StreamReport
...type:Spec.StoreReport
..remarks:Use Metafunction.Hit to get the BlastHit type used in a BlastReport object.
..include:seqan/blast.h
*/
template<typename TBlastHsp, typename TSpec>
class BlastHit<TBlastHsp, StoreReport<TSpec> > 
{
	public:
		String<char> name;
		unsigned int length; //length of whole sequence 

		String<TBlastHsp> hsps;
		
		BlastHit()
		{
		SEQAN_CHECKPOINT
		}

		//BlastHit(String<char> name, unsigned int len)
		//{
		//	name = name;
		//	length = length;
		//	clear(hsps);
		//}

		BlastHit(BlastHit const& other)
		{
		SEQAN_CHECKPOINT
			assign(hsps,other.hsps);
			name = other.name;
			length = other.length;
		}

		BlastHit & operator = (BlastHit const & other)
		{
		SEQAN_CHECKPOINT
			assign(hsps,other.hsps);
			name = other.name;
			length = other.length;
			return *this;
		}

		~BlastHit()
		{
		}

};



template<typename TBlastHsp, typename TSpec>
inline void
clear(BlastHit<TBlastHsp, StoreReport<TSpec> >& blastHit)
{
SEQAN_CHECKPOINT
	
	for(unsigned int i = 0; i < length(blastHit.hsps); ++i)
		clear(blastHit.hsps[i]);
	resize(blastHit.hsps,0);
	resize(blastHit.name,0);
}


template<typename TBlastHsp, typename TStoreSpec>
inline String<char> &
name(BlastHit<TBlastHsp, TStoreSpec>& blastHit)
{
SEQAN_CHECKPOINT
	return blastHit.name;
}

template<typename TBlastHsp, typename TStoreSpec>
inline String<char> 
getName(BlastHit<TBlastHsp, TStoreSpec>& blastHit)
{
SEQAN_CHECKPOINT
	return blastHit.name;
}

template<typename TBlastHsp, typename TStoreSpec>
inline unsigned int &
length(BlastHit<TBlastHsp, TStoreSpec>& blastHit)
{
SEQAN_CHECKPOINT
	return blastHit.length;
}

template<typename TBlastHsp, typename TStoreSpec>
inline unsigned int 
getLength(BlastHit<TBlastHsp, TStoreSpec>& blastHit)
{
SEQAN_CHECKPOINT
	return blastHit.length;
}

// for StoreReport only
template<typename TBlastHsp, typename TSpec>
inline unsigned int 
numHsps(BlastHit<TBlastHsp, StoreReport<TSpec> >& blastHit)
{
SEQAN_CHECKPOINT
	return length(blastHit.hsps);
}		


/////////////////////////////////////////////////////////////////////

//parse BlastHit
template<typename TFile, typename TChar, typename TBlastHit>
inline typename Position<TFile>::Type
_parseBlastHit(TFile & file,
			TChar & c, 
			TBlastHit & hit)
{
//IOREV _nodoc_ _hasCRef_
	typedef typename Position<TFile>::Type TPosition;
	typedef typename Hsp<TBlastHit>::Type TBlastHsp;

	String<char> pword;
	int pint;
	TPosition start_pos,act_pos;
	act_pos = _streamTellG(file);

	if(_parseUntilBeginLine(file,c,'>'))
	{
		start_pos = _streamTellG(file);
		c = _streamGet(file);
		pword = _parseReadWord(file, c);
		while (!_streamEOF(file) && c != '\n' && c != '\r')
			pword += _parseReadWord(file, c);
		if(pword[length(pword)-1] == ' ')
			resize(pword,length(pword)-1);
		hit.name = pword;
		_parseSkipWhitespace(file,c);
		String<char> search = "Length";
		if(_parseUntilBeginLine(file,c,search,6))
		{
			_parseSkipWhitespace(file,c);
			if(c == '=')
				c = _streamGet(file);
			_parseSkipWhitespace(file,c);
			pint = _parseReadNumber(file, c);
			hit.length = pint;
		}
//		TPosition temp = _streamTellG(file);
		//foreach Hsp
		//if(_parseUntilBeginLine(file,c,'S') && _parseReadWord(file,c)=="Score")
		search = "Score";
		if(_parseUntilBeginLine(file,c,search,5))
		{
			//c = _streamGet(file);
			bool in_hit = true;
			TPosition act_hsp_pos,next_hsp_pos;
			while(in_hit){
				act_hsp_pos = _streamTellG(file);
				TBlastHsp hsp;
				next_hsp_pos = _parseBlastHsp(file,c,hsp);
				//resize(hit.hsps, length(hit.hsps)+1);
				//hit.hsps[length(hit.hsps)-1] = hsp;
				appendValue(hit.hsps,hsp);
				//append(hit.hsps,hsp);
				if(next_hsp_pos == act_hsp_pos)
					in_hit = false;
				if(next_hsp_pos == (TPosition)0)
					return (TPosition) 0;
			}
			_streamSeekG(file,next_hsp_pos);
			c = _streamGet(file);
			if(_parseUntilBeginLine(file,c,'>'))
				return _streamTellG(file);
		}

		
	}//end hit
	_streamSeekG(file,act_pos);
	return act_pos;
}






////////////////////////// MetaFunctions ////////////////////////////////

template<typename TBlastHsp, typename TStoreSpec>
struct Value<BlastHit<TBlastHsp, TStoreSpec> > 
{
	typedef BlastHit<TBlastHsp,TStoreSpec> Type;
};

template<typename TBlastHsp, typename TStoreSpec>
struct Value<BlastHit<TBlastHsp, TStoreSpec> const> 
{
	typedef BlastHit<TBlastHsp,TStoreSpec> Type;
};

template<typename TBlastHsp, typename TStoreSpec>
struct Hit<BlastHit<TBlastHsp, TStoreSpec> > 
{
	typedef BlastHit<TBlastHsp,TStoreSpec> Type;
};

template<typename TBlastHsp, typename TStoreSpec>
struct Hit<BlastHit<TBlastHsp, TStoreSpec> const> 
{
	typedef BlastHit<TBlastHsp,TStoreSpec> Type;
};



template<typename TBlastHsp, typename TStoreSpec>
struct Hsp<BlastHit<TBlastHsp, TStoreSpec> > 
{
	typedef TBlastHsp Type;
};

template<typename TBlastHsp, typename TStoreSpec>
struct Hsp<BlastHit<TBlastHsp, TStoreSpec> const> 
{
	typedef TBlastHsp Type;
};


/////////////komisch///////////

	template<typename TBlastHsp, typename TStoreSpec>
	struct Hsp<String<BlastHit<TBlastHsp, TStoreSpec> > > 
	{
		typedef TBlastHsp Type;
	};

	template<typename TBlastHsp, typename TStoreSpec>
	struct Hsp<String<BlastHit<TBlastHsp, TStoreSpec> const> > 
	{
		typedef TBlastHsp Type;
	};

	template<typename TBlastSpec, typename TStoreSpec>
	struct Hsp<String<BlastHsp<TBlastSpec, TStoreSpec> > > 
	{
		typedef BlastHsp<TBlastSpec, TStoreSpec> Type;
	};

	template<typename TBlastSpec, typename TStoreSpec>
	struct Hsp<String<BlastHsp<TBlastSpec, TStoreSpec> const> > 
	{
		typedef BlastHsp<TBlastSpec, TStoreSpec> Type;
	};




}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
