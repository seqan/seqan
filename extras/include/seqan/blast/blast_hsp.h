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

#ifndef SEQAN_HEADER_BLAST_HSP_H
#define SEQAN_HEADER_BLAST_HSP_H


namespace SEQAN_NAMESPACE_MAIN
{







template<typename TBlastSpec = NucleotideBlast<>, typename TInfoSpec = BasicInfo>
class BlastHsp;




/////////////////////////////////////////////////////////////////////////////
//////////////////////// FullInfo Nucleotide Spec //////////////////////////////
/////////////////////////////////////////////////////////////////////////////

//Info types
/**
.Spec.FullInfo:
..cat:Blast
..general:Class.BlastHsp
..summary:Stores all pieces of information delivered with an alignment in a Blast report.
..signature:FullInfo
..include:blast.h
*/

/**
.Spec.BasicInfo:
..cat:Blast
..general:Class.BlastHsp
..summary:Stores only the basic pieces of information delivered with an alignment in a Blast report.
..signature:BasicInfo
..include:blast.h
*/



/**
.Class.BlastHsp:
..cat:Blast
..summary:Object for storing Blast HSPs. 
..signature:BlastHsp<TBlastSpec, TInfoSpec>  
..param.TBlastSpec:The type of Blast report to be parsed.
...type:Spec.BlastN
...type:Spec.BlastP
...default:BlastN
..param.TInfoSpec:The specializing type determining the amount of information to be stored.
...type:Spec.BasicInfo
...type:Spec.FullInfo
...default:BasicInfo
...remarks:BasicInfo stores begin and end positions on query and database sequence, as well as the alignment. FullInfo stores additional information such as score, e-value...
..include:blast.h
*/
template<typename TSpec>
class BlastHsp<NucleotideBlast<TSpec>, FullInfo> 
{

	public:
		float score;
		float bits;
		double expect;
		unsigned int identity;
		unsigned int gaps; 
		unsigned int abs_gaps; 
		bool query_strand;  //false=>minus, true=>plus
		bool db_strand;
		unsigned int query_begin;
		unsigned int db_begin;
		unsigned int query_end;	//included in the alignment
		unsigned int db_end;	//included in the alignment

		String<char> query_string;	
		String<char> db_string;
 

		BlastHsp()
		{
		SEQAN_CHECKPOINT
			//defaults for those values that are not necessarily present in each HSP in the Blast Report
			score = 0.0;
			bits = 0.0;
			expect = 0.0;
			identity = 0;
			gaps = 0; 
			abs_gaps = 0; 
			query_strand = true;  //false=>minus, true=>plus
			db_strand = true;
		}

		BlastHsp(BlastHsp const& other)
		{
		SEQAN_CHECKPOINT
			score = other.score;
			bits = other.bits;
			expect = other.expect;
			identity = other.identity;
			gaps = other.gaps;
			abs_gaps = other.abs_gaps;
			query_strand = other.query_strand;
			db_strand = other.db_strand;
			query_begin = other.query_begin;
			db_begin = other.db_begin;
			query_end = other.query_end;
			db_end = other.db_end;	
			query_string = other.query_string;	
			db_string = other.db_string;

		}

		BlastHsp & operator = (BlastHsp const & other)
		{
		SEQAN_CHECKPOINT
			score = other.score;
			bits = other.bits;
			expect = other.expect;
			identity = other.identity;
			gaps = other.gaps;
			abs_gaps = other.abs_gaps;
			query_strand = other.query_strand;
			db_strand = other.db_strand;
			query_begin = other.query_begin;
			db_begin = other.db_begin;
			query_end = other.query_end;
			db_end = other.db_end;	
			query_string = other.query_string;	
			db_string = other.db_string;

			return *this;
		}


		~BlastHsp()
		{
		SEQAN_CHECKPOINT
		}
};


//parse BlastHsp
template<typename TFile, typename TChar, typename TSpec>
inline typename Position<TFile>::Type
_parseBlastHsp(TFile & file,
			TChar & c, 
			BlastHsp<NucleotideBlast<TSpec>,FullInfo > & hsp)
{
//IOREV _hasCRef_
SEQAN_CHECKPOINT
	typedef typename Position<TFile>::Type TPosition;
	
	clear(hsp);
	String<char> pword;
	int pint;
	float pfloat;
	double pdouble;
	TPosition start_pos,act_pos;
	act_pos = _streamTellG(file);

    String<char> query_string;
	String<char> db_string;

	_parseSkipWhitespace(file,c);
	c = _streamGet(file);
	_parseSkipWhitespace(file,c);
	pfloat = (float)_parseReadEValue(file, c);
	hsp.bits = pfloat;
	if(_parseLineUntil(file,c,'('))
	{
		c = _streamGet(file);
		pfloat = _parseReadFloat(file, c);
		hsp.score = pfloat;
	}

	_parseSkipWhitespace(file,c);
	String<char> query = "Expect";
	if(_parseLineUntil(file,c,query,6))
	{
		c = _streamGet(file);
		_parseSkipWhitespace(file,c);
		if(c == '=')
			c = _streamGet(file);
		_parseSkipWhitespace(file,c);
		pdouble = _parseReadEValue(file, c);
		hsp.expect = pdouble;
	}

	query = "Identities";
	if(_parseUntilBeginLine(file,c,query,10,3))
	{
		_parseLineUntil(file,c,'(');
		c = _streamGet(file);
		pint = _parseReadNumber(file, c);
		hsp.identity = pint;
	}

	query = "Gaps";
	if(_parseLineUntil(file,c,query,4))
	{
		_parseLineUntil(file,c,'=');
		c = _streamGet(file);
		pint = _parseReadNumber(file, c);
		hsp.abs_gaps = pint;
		_parseLineUntil(file,c,'(');
		c = _streamGet(file);
		pint = _parseReadNumber(file, c);
		hsp.gaps = pint;
	}
	//else
	//{
	//	hsp.abs_gaps = 0;
	//	hsp.gaps = 0;
	//}

	query = "Strand";
	if(_parseUntilBeginLine(file,c,query,6,3))
	{
		c = _streamGet(file);
		_parseSkipWhitespace(file,c);
		if(c == '=')
			c = _streamGet(file);
		_parseSkipWhitespace(file,c);
		pword = _parseReadWord(file, c);
		if(pword == "Plus")
			hsp.query_strand = true;
		else
			hsp.query_strand = false;
		c = _streamGet(file);
		_parseSkipWhitespace(file,c);
		if(c == '/')
			c = _streamGet(file);
		_parseSkipWhitespace(file,c);
		pword = _parseReadWord(file, c);
		if(pword == "Plus")
			hsp.db_strand = true;
		else
			hsp.db_strand = false;
	}

	bool first = true;
	query = "Query";
	if(_parseUntilBeginLine(file,c,query,5))
	{
		bool in_hsp = true;
		while(in_hsp)
		{
			int end_query,end_db = -1;
			_parseSkipWhitespace(file,c);
			if(c == ':')
				c = _streamGet(file);
			_parseSkipWhitespace(file,c);
			pint = _parseReadNumber(file,c);
			_parseSkipWhitespace(file,c);
			//TPosition act_posQ = _streamTellG(file);
			if(first)
				hsp.query_begin = pint;
			query_string += _parseReadAlignmentString(file,c);
			_parseSkipWhitespace(file,c);
			end_query = _parseReadNumber(file,c);
			query = "Sbjct";
			_parseUntilBeginLine(file,c,query,5);
			_parseSkipWhitespace(file,c);
			if(c == ':')
				c = _streamGet(file);
			_parseSkipWhitespace(file,c);
			pint = _parseReadNumber(file,c);
			_parseSkipWhitespace(file,c);
			//TPosition act_posS = _streamTellG(file);
			if(first)
				hsp.db_begin = pint;
			db_string += _parseReadAlignmentString(file,c);
			_parseSkipWhitespace(file,c);
			end_db = _parseReadNumber(file,c);
			first = false;
			String<TChar> delim = "QS>R";
			if(_parseUntilBeginLineOneOf(file,c,delim,4))
			{
				// TPosition pos = _streamTellG(file);
				pword = _parseReadWord(file,c);
				
				//still in the same HSP
				if(pword == "Query")
					in_hsp = true;
				else
				{
					hsp.query_end = end_query;
					hsp.db_end = end_db;
					hsp.query_string = query_string;
					hsp.db_string = db_string;
					//next HSP
					if(pword == "Score")
						return _streamTellG(file);
					if(pword[0] == '>')
						return act_pos;
					if(pword[0] == 'R')
						return (TPosition) 0;
					in_hsp = false;
				}
			}
			else
				return act_pos;
		}
	}

	return act_pos;

}


template<typename TSpec>
inline void
clear(BlastHsp<NucleotideBlast<TSpec>, FullInfo>& blastHsp)
{
SEQAN_CHECKPOINT

	blastHsp.identity = 0;
	blastHsp.score = 0;
	blastHsp.bits = 0.0;
	blastHsp.expect = 0.0;
	blastHsp.gaps = 0; 
	blastHsp.abs_gaps = 0; 
	blastHsp.query_begin = 0;  
	blastHsp.db_end = 0;
	blastHsp.query_end = 0;  
	blastHsp.db_begin = 0;
	blastHsp.query_strand = true;  
	blastHsp.db_strand = true;
	resize(blastHsp.query_string,0);
	resize(blastHsp.db_string,0);
}

//////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////
//////////////////////// FullInfo Protein Spec //////////////////////////////
/////////////////////////////////////////////////////////////////////////////

//Protein BlastHsp
template<typename TSpec>
class BlastHsp<ProteinBlast<TSpec>, FullInfo> 
{

	public:
		float score;
		float bits;
		double expect;
		unsigned int identity;
		unsigned int positives; 
		unsigned int gaps; 
		unsigned int abs_gaps; 
		bool query_strand;  //false=>minus, true=>plus
		bool db_strand;
		unsigned int query_begin;
		unsigned int db_begin;
		unsigned int query_end;	//included in the alignment
		unsigned int db_end;	//included in the alignment
		int query_frame;
		int db_frame;
		String<char> query_string;	
		String<char> db_string;
	

		BlastHsp()
		{
		SEQAN_CHECKPOINT
			score = 0.0;
			bits = 0.0 ;
			expect = 0.0;
			identity = 0;
			positives = 0; 
			gaps = 0; 
			abs_gaps = 0;
			query_strand = true;  //false=>minus, true=>plus
			db_strand = true;
			query_frame = 1;
			db_frame = 1;
		}

		BlastHsp(BlastHsp const& other)
		{
		SEQAN_CHECKPOINT
			score = other.score;
			bits = other.bits;
			expect = other.expect;
			identity = other.identity;
			gaps = other.gaps;
			abs_gaps = other.abs_gaps;
			query_strand = other.query_strand;
			db_strand = other.db_strand;
			query_begin = other.query_begin;
			db_begin = other.db_begin;
			query_end = other.query_end;
			db_end = other.db_end;	
			query_string = other.query_string;	
			db_string = other.db_string;
			positives = other.positives; 
			query_frame = other.query_frame;
			db_frame = other.db_frame;
		}

		BlastHsp & operator = (BlastHsp const & other)
		{
		SEQAN_CHECKPOINT
			score = other.score;
			bits = other.bits;
			expect = other.expect;
			identity = other.identity;
			gaps = other.gaps;
			abs_gaps = other.abs_gaps;
			query_strand = other.query_strand;
			db_strand = other.db_strand;
			query_begin = other.query_begin;
			db_begin = other.db_begin;
			query_end = other.query_end;
			db_end = other.db_end;	
			query_string = other.query_string;	
			db_string = other.db_string;
			positives = other.positives; 
			query_frame = other.query_frame;
			db_frame = other.db_frame;
			return *this;
		}

		~BlastHsp()
		{
		SEQAN_CHECKPOINT
		}
};



//parse BlastHsp
template<typename TFile, typename TChar, typename TSpec>
inline typename Position<TFile>::Type
_parseBlastHsp(TFile & file,
			TChar & c, 
			BlastHsp<ProteinBlast<TSpec>,FullInfo> & hsp)
{
//IOREV _nodoc_ _hasCRef_
SEQAN_CHECKPOINT
	typedef typename Position<TFile>::Type TPosition;
	
	clear(hsp);
	String<char> pword;
	int pint;
	float pfloat;
	double pdouble;
	TPosition start_pos,act_pos;
	act_pos = _streamTellG(file);

	//String<char> query = "Score";
	//if(_parseUntilBeginLine(file,c,query,6))
	//{
	String<char> query_string;
	String<char> db_string;

	_parseSkipWhitespace(file,c);
	if(c == '=')
		c = _streamGet(file);
	_parseSkipWhitespace(file,c);
	pfloat = (float)_parseReadEValue(file, c);
	hsp.bits = pfloat;
	if(_parseLineUntil(file,c,'('))
	{
		c = _streamGet(file);
		pfloat = _parseReadFloat(file, c);
		hsp.score = pfloat;
	}

	_parseSkipWhitespace(file,c);
	String<char> query = "Expect";
	if(_parseLineUntil(file,c,query,6))
	{
		c = _streamGet(file);
		_parseSkipWhitespace(file,c);
		if(c == '=')
			c = _streamGet(file);
		_parseSkipWhitespace(file,c);
		pdouble = _parseReadEValue(file, c);
		hsp.expect = pdouble;
	}

	query = "Identities";
	if(_parseUntilBeginLine(file,c,query,10,3))
	{
		_parseLineUntil(file,c,'(');
		c = _streamGet(file);
		pint = _parseReadNumber(file, c);
		hsp.identity = pint;
	}

	query = "Positives";
	if(_parseLineUntil(file,c,query,9))
	{
		_parseUntil(file,c,'(');
		c = _streamGet(file);
		pint = _parseReadNumber(file, c);
		hsp.positives = pint;
	}

	query = "Gaps";
	if(_parseLineUntil(file,c,query,4))
	{
		_parseLineUntil(file,c,'=');
		c = _streamGet(file);
		pint = _parseReadNumber(file, c);
		hsp.abs_gaps = pint;
		_parseLineUntil(file,c,'(');
		c = _streamGet(file);
		pint = _parseReadNumber(file, c);
		hsp.gaps = pint;
	}

	query = "Strand";
	if(_parseUntilBeginLine(file,c,query,6,3))
	{
		c = _streamGet(file);
		_parseSkipWhitespace(file,c);
		if(c == '=')
			c = _streamGet(file);
		_parseSkipWhitespace(file,c);
		pword = _parseReadWord(file, c);
		if(pword == "Plus")
			hsp.query_strand = true;
		else
			hsp.query_strand = false;
		c = _streamGet(file);
		_parseSkipWhitespace(file,c);
		if(c == '/')
			c = _streamGet(file);
		_parseSkipWhitespace(file,c);
		pword = _parseReadWord(file, c);
		if(pword == "Plus")
			hsp.db_strand = true;
		else
			hsp.db_strand = false;
	}

	bool first = true;
	query = "Query";
	if(_parseUntilBeginLine(file,c,query,5))
	{
		bool in_hsp = true;
		while(in_hsp)
		{
			int end_query,end_db = -1;
			_parseSkipWhitespace(file,c);
			if(c == ':')
				c = _streamGet(file);
			_parseSkipWhitespace(file,c);
			pint = _parseReadNumber(file,c);
			_parseSkipWhitespace(file,c);
			//TPosition act_posQ = _streamTellG(file);
			if(first)
				hsp.query_begin = pint;
			query_string += _parseReadAlignmentString(file,c);
			_parseSkipWhitespace(file,c);
			end_query = _parseReadNumber(file,c);
			query = "Sbjct";
			_parseUntilBeginLine(file,c,query,5);
			_parseSkipWhitespace(file,c);
			if(c == ':')
				c = _streamGet(file);
			_parseSkipWhitespace(file,c);
			pint = _parseReadNumber(file,c);
			_parseSkipWhitespace(file,c);
			//TPosition act_posS = _streamTellG(file);
			if(first)
				hsp.db_begin = pint;
			db_string += _parseReadAlignmentString(file,c);
			_parseSkipWhitespace(file,c);
			end_db = _parseReadNumber(file,c);
			first = false;
			String<TChar> delim = "QS>R";
			if(_parseUntilBeginLineOneOf(file,c,delim,4))
			{
				// TPosition pos = _streamTellG(file);
				pword = _parseReadWord(file,c);
			//	_streamSeekG(file,pos);
				
				//still in the same HSP
				if(pword == "Query")
					in_hsp = true;
				else
				{
					hsp.query_end = end_query;
					hsp.db_end = end_db;
					hsp.query_string = query_string;
					hsp.db_string = db_string;
					//next HSP
					if(pword == "Score")
						return _streamTellG(file);
					if(pword[0] == '>')
						return act_pos;
					if(pword[0] == 'R')
						return (TPosition) 0;
					in_hsp = false;
				}
			}
			else
				return act_pos;
		}
	}
	return act_pos;

}


template<typename TSpec>
inline void
clear(BlastHsp<ProteinBlast<TSpec>, FullInfo>& blastHsp)
{
SEQAN_CHECKPOINT
	blastHsp.score = 0;
	blastHsp.bits = 0.0;
	blastHsp.expect = 0.0;
	blastHsp.abs_gaps = 0; 
	blastHsp.query_begin = 0;  
	blastHsp.db_end = 0;
	blastHsp.query_end = 0;  
	blastHsp.db_begin = 0;
	blastHsp.positives = 0; 
	blastHsp.query_frame = 1;
	blastHsp.db_frame = 1;
	blastHsp.identity = 0;
	blastHsp.gaps = 0; 
	blastHsp.query_strand = true;  
	blastHsp.db_strand = true;
	resize(blastHsp.query_string,0);
	resize(blastHsp.db_string,0);
}






/////////////////////////////////////////////////////////////////////////////
/////////////////////////// BasicInfo Spec //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////



template<typename TBlastSpec>
class BlastHsp<TBlastSpec, BasicInfo> 
{

	public:
		double expect;
		unsigned int query_begin;
		unsigned int db_begin;
		unsigned int query_end;	//included in the alignment
		unsigned int db_end;	//included in the alignment

		String<char> query_string;	
		String<char> db_string;


		BlastHsp()
		{
		SEQAN_CHECKPOINT
			expect = 0;	
		}

		BlastHsp(BlastHsp const& other)
		{
		SEQAN_CHECKPOINT
			expect = other.expect;
			query_begin = other.query_begin;
			db_begin = other.db_begin;
			query_end = other.query_end;
			db_end = other.db_end;	
			query_string = other.query_string;	
			db_string = other.db_string;
		}

		BlastHsp & operator = (BlastHsp const & other)
		{
		SEQAN_CHECKPOINT
			expect = other.expect;
			query_begin = other.query_begin;
			db_begin = other.db_begin;
			query_end = other.query_end;
			db_end = other.db_end;	
			query_string = other.query_string;	
			db_string = other.db_string;
			return *this;
		}
 
		~BlastHsp()
		{
		SEQAN_CHECKPOINT
		}
};


template<typename TBlastSpec>
inline void
clear(BlastHsp<TBlastSpec, BasicInfo>& blastHsp)
{
SEQAN_CHECKPOINT

	blastHsp.expect = 0.0;
	blastHsp.query_begin = 0;  
	blastHsp.db_end = 0;
	blastHsp.query_end = 0;  
	blastHsp.db_begin = 0;
	resize(blastHsp.query_string,0);
	resize(blastHsp.db_string,0);
}



//parse BlastHsp
template<typename TFile, typename TChar, typename TBlastSpec, typename TInfoSpec>
inline typename Position<TFile>::Type
_parseBlastHsp(TFile & file,
			TChar & c, 
			BlastHsp<TBlastSpec,TInfoSpec> & hsp)
{
//IOREV _nodoc_ _hasCRef_
SEQAN_CHECKPOINT
	typedef typename Position<TFile>::Type TPosition;
	
	clear(hsp);
	
	String<char> pword;
	int pint;
	double pdouble;
	TPosition start_pos,act_pos;
	act_pos = _streamTellG(file);

	String<char> query_string;
	String<char> db_string;

	_parseSkipWhitespace(file,c);
	String<char> query = "Expect";
	if(_parseLineUntil(file,c,query,6))
	{
		c = _streamGet(file);
		_parseSkipWhitespace(file,c);
		c = _streamGet(file); // = 
		_parseSkipWhitespace(file,c);
		pdouble = _parseReadEValue(file, c);
		hsp.expect = pdouble;
	}

	bool first = true;
	query = "Query";
	if(_parseUntilBeginLine(file,c,query,5))
	{
		bool in_hsp = true;
		while(in_hsp)
		{
			int end_query,end_db = -1;
			_parseSkipWhitespace(file,c);
			c = _streamGet(file);
			_parseSkipWhitespace(file,c);
			pint = _parseReadNumber(file,c);
			_parseSkipWhitespace(file,c);
			if(first)
				hsp.query_begin = pint;
			query_string += _parseReadAlignmentString(file,c);
			_parseSkipWhitespace(file,c);
			end_query = _parseReadNumber(file,c);
			query = "Sbjct";
			_parseUntilBeginLine(file,c,query,5);
			_parseSkipWhitespace(file,c);
			c = _streamGet(file);
			_parseSkipWhitespace(file,c);
			pint = _parseReadNumber(file,c);
			_parseSkipWhitespace(file,c);
			if(first)
				hsp.db_begin = pint;
			db_string += _parseReadAlignmentString(file,c);
			_parseSkipWhitespace(file,c);
			end_db = _parseReadNumber(file,c);
			first = false;
			String<TChar> delim = "QS>R";
			if(_parseUntilBeginLineOneOf(file,c,delim,4))
			{
				// TPosition pos = _streamTellG(file);
				pword = _parseReadWord(file,c);
				
				//still in the same HSP
				if(pword == "Query")
					in_hsp = true;
				else
				{
					hsp.query_end = end_query;
					hsp.db_end = end_db;
					hsp.query_string = query_string;
					hsp.db_string = db_string;
					//next HSP
					if(pword == "Score")
						return _streamTellG(file);
					if(pword[0] == '>')
						return act_pos;
					if(pword[0] == 'R')
						return (TPosition) 0;
					in_hsp = false;
				}
			}
			else
				return act_pos;
		}
	}
	return act_pos;

}




/////////////////////////////////////////////////////////////
/////////////////////// get functions ///////////////////////
/////////////////////////////////////////////////////////////

struct TagUnknownSource_;
typedef Tag<TagUnknownSource_> const UnknownSource;

struct TagKnownSource_;
typedef Tag<TagKnownSource_> const KnownSource;

/////////////////////////////////////////////////////////////
// get Alignment for Align<TSource,TSpec>
template<typename TBlastHsp, typename TSpec, typename TSource>
inline unsigned int
getAlignment(TBlastHsp & hsp,
			 Align<TSource,TSpec> & ali, 
			 UnknownSource)
{
SEQAN_CHECKPOINT

	typedef Align<TSource, TSpec> TAlign;
	typedef typename Row<TAlign>::Type TRow;

	typedef Pair<unsigned int, unsigned int> TPosLen;

	//for storing node info
	String<TPosLen> seq0info;
	String<TPosLen> seq1info;

	TSource str0,str1;

	typedef typename Iterator<String<char>, Standard >::Type TStringIterator;
	TStringIterator it_0 = begin(queryAlignmentString(hsp),Standard());
	TStringIterator it_0_end = end(queryAlignmentString(hsp),Standard());
	TStringIterator it_1 = begin(databaseAlignmentString(hsp),Standard());
	TStringIterator it_1_end = end(databaseAlignmentString(hsp),Standard());

	unsigned int act0_pos = 0;
	unsigned int act1_pos = 0;
	
	while ((it_0 != it_0_end) && (it_1 != it_1_end))
	{
		if(*it_0 == '-')
		{
			unsigned int gap_len = 0;
			while(it_0 != it_0_end && *it_0 == '-')
			{
				append(str1,getValue(it_1));
				++it_0;
				++it_1;
				++act1_pos;
				++gap_len;
			}
			append(seq0info,TPosLen(act0_pos,gap_len));
			act0_pos += gap_len;
			continue;
		}
		if(*it_1 == '-')
		{
			unsigned int gap_len = 0;
			while(it_1 != it_1_end && *it_1 == '-')
			{
				append(str0,getValue(it_0));
				++it_1;
				++it_0;
				++act0_pos;
				++gap_len;
			}
			append(seq1info,TPosLen(act1_pos,gap_len));
			act1_pos += gap_len;
			continue;
		}
		append(str0,getValue(it_0));
		append(str1,getValue(it_1));
		++it_0;
		++it_1;
		++act0_pos;
		++act1_pos;

	}
	
	resize(rows(ali),2);
	setSource(row(ali,0),str0);
	setSource(row(ali,1),str1);

	typedef typename Iterator<String<TPosLen>, Standard >::Type TNodeIterator;
	TNodeIterator nodes0it = begin(seq0info,Standard());
	TNodeIterator nodes0end = end(seq0info,Standard());

	while(nodes0it != nodes0end)
	{
		insertGaps(row(ali,0),nodes0it->i1,nodes0it->i2);
		++nodes0it;
	}


	typedef typename Iterator<String<TPosLen>, Standard >::Type TNodeIterator;
	TNodeIterator nodes1it = begin(seq1info,Standard());
	TNodeIterator nodes1end = end(seq1info,Standard());

	while(nodes1it != nodes1end)
	{
		insertGaps(row(ali,1),nodes1it->i1,nodes1it->i2);
		++nodes1it;
	}
	detach(ali);

	return length(queryAlignmentString(hsp));

}

/**
.Function.BlastHsp#getAlignment
..cat:Blast
..summary:Turns a HSP from a Blast search into an Alignment object.
..signature:getAlignment(hsp,alignment)
..param.hsp:A Blast HSP.
...type:Class.BlastHsp
..param.alignment:An Alignment object to be filled.
...type:Class.Align
...type:Spec.Alignment Graph
..include:seqan/blast.h
*/
template<typename TBlastHsp, typename TSpec, typename TSource>
inline unsigned int
getAlignment(TBlastHsp & hsp,
			 Align<TSource,TSpec> & ali)
{
SEQAN_CHECKPOINT

	if(length(rows(ali))>1)
        if(!empty(source(row(ali,0))) && !!empty(source(row(ali,1))))
			return getAlignment(hsp,ali,KnownSource());
	
	return getAlignment(hsp,ali,UnknownSource());

}



/////////////////////////////////////////////////////////////
// get Alignment for Align<TSource,TSpec>
template<typename TBlastHsp, typename TSpec, typename TSource>
inline unsigned int
getAlignment(TBlastHsp & hsp,
			 Align<TSource,TSpec> & ali, 
			 KnownSource)
{
SEQAN_CHECKPOINT

	typedef Align<TSource, TSpec> TAlign;
	typedef typename Row<TAlign>::Type TRow;

	typedef typename Iterator<TRow, Standard>::Type TAliIterator;
	TAliIterator ali_0 = begin(row(ali, 0), Standard());
	TAliIterator ali_1 = begin(row(ali, 1), Standard());

	typedef typename Iterator<String<char>, Standard >::Type TStringIterator;
	TStringIterator it_0 = begin(queryAlignmentString(hsp),Standard());
	TStringIterator it_0_end = end(queryAlignmentString(hsp),Standard());

	TStringIterator it_1 = begin(databaseAlignmentString(hsp),Standard());
	TStringIterator it_1_end = end(databaseAlignmentString(hsp),Standard());

	TSource str0,str1;

	while ((it_0 != it_0_end) && (it_1 != it_1_end))
	{
		if (*it_0 == '-')
		{
			insertGap(ali_0);
		}
		else
		{
			append(str0,getValue(it_0));
		}
		if (*it_1 == '-')
		{
			insertGap(ali_1);
		}
		else
		{
			append(str1,getValue(it_1));
		}
		
		++it_0;
		++it_1;
		++ali_0;
		++ali_1;
	}

	setClippedBeginPosition(row(ali,0),queryBegin(hsp)-1);
	setClippedBeginPosition(row(ali,1),databaseBegin(hsp)-1);
	setClippedEndPosition(row(ali,0),queryEnd(hsp));
	setClippedEndPosition(row(ali,1),databaseEnd(hsp));

	return length(databaseAlignmentString(hsp));

}





/**
.Function.BlastHsp#getAlignment
..cat:Blast
..signature:getAlignment(hsp,alignment,id0,id1)
..param.id0:The Id of the query sequence in the StringSet of the Alignment Graph.
..param.id1:The Id of the hit sequence in the StringSet of the Alignment Graph.
..include:seqan/blast.h
*/
/////////////////////////////////////////////////////////////
// get Alignment for Graph<TAlign>
template<typename TBlastHsp, typename TAlign, typename TId>
inline unsigned int
getAlignment(TBlastHsp & hsp,
			 Graph<TAlign> & ali,
			 TId id0, //query ID 
			 TId id1) //hit ID
{
SEQAN_CHECKPOINT

	typedef Graph<TAlign> TAliGraph;
	typedef typename Iterator<String<char>, Standard >::Type TStringIterator;
	typedef typename VertexDescriptor<TAliGraph>::Type TVertexDescriptor;

	TStringIterator it_0 = begin(queryAlignmentString(hsp),Standard());
	TStringIterator it_0_end = end(queryAlignmentString(hsp),Standard());

	TStringIterator it_1 = begin(databaseAlignmentString(hsp),Standard());
	TStringIterator it_1_end = end(databaseAlignmentString(hsp),Standard());
	unsigned int act0_pos = getQueryBegin(hsp)-1;
	unsigned int act1_pos = getDatabaseBegin(hsp)-1;

	unsigned int act_len = 0;
	while ((it_0 != it_0_end) && (it_1 != it_1_end))
	{
		if(*it_0 == '-')
		{
			if(act_len > 0)
			{
				TVertexDescriptor vd0 = addVertex(ali,id0,act0_pos,act_len);
				TVertexDescriptor vd1 = addVertex(ali,id1,act1_pos,act_len);
				addEdge(ali,vd0,vd1);
				act0_pos += act_len;
				act1_pos += act_len;
				act_len = 0;
			}
			unsigned int gap_len = 0;
			while(it_0 != it_0_end && *it_0 == '-')
			{
				++it_0;
				++it_1;
				++gap_len;
			}
			addVertex(ali,id1,act1_pos,gap_len);
			act1_pos += gap_len;
			continue;
		}
		if(*it_1 == '-')
		{
			if(act_len > 0)
			{
				TVertexDescriptor vd0 = addVertex(ali,id0,act0_pos,act_len);
				TVertexDescriptor vd1 = addVertex(ali,id1,act1_pos,act_len);
				addEdge(ali,vd0,vd1);
				act0_pos += act_len;
				act1_pos += act_len;
				act_len = 0;
			}
			unsigned int gap_len = 0;
			while(it_1 != it_1_end && *it_1 == '-')
			{
				++it_1;
				++it_0;
				++gap_len;
			}
			addVertex(ali,id0,act0_pos,gap_len);
			act0_pos += gap_len;
			continue;
		}
		++it_0;
		++it_1;
		++act_len;

	}
	SEQAN_ASSERT(act0_pos+act_len == getQueryEnd(hsp));
	SEQAN_ASSERT(act1_pos+act_len == getDatabaseEnd(hsp));
	if(act_len>0)
	{
		TVertexDescriptor vd0 = addVertex(ali,id0,act0_pos,act_len);
		TVertexDescriptor vd1 = addVertex(ali,id1,act1_pos,act_len);
		addEdge(ali,vd0,vd1);
	}

	return length(queryAlignmentString(hsp));

}



/////////////////////////////////////////////////////////////
// get Alignment for Graph<TAlign> (stringSet and sequence ids not given)
template<typename TBlastHsp, typename TStringSet, typename TCargo, typename TSpec>
inline unsigned int //returns the length of the alignment
getAlignment(TBlastHsp & hsp,
			 Graph<Alignment<TStringSet,TCargo,TSpec> > & ali) //hit ID
{
SEQAN_CHECKPOINT

	typedef Graph<Alignment<TStringSet,TCargo,TSpec> > TAliGraph;
	typedef typename Iterator<String<char>, Standard >::Type TStringIterator;
	typedef typename VertexDescriptor<TAliGraph>::Type TVertexDescriptor;
	typedef Pair<unsigned int, unsigned int> TPosLen;
	typedef Triple<unsigned int, unsigned int, int> TPosLenLink;
	typedef typename Value<TStringSet>::Type TString;

	//for storing node info
	String<TPosLen> seq0nodes;
	//for storing node + edge info
	String<TPosLenLink> seq1nodes;

	//for storing the sequences
	TString str0,str1;

	TStringIterator it_0 = begin(queryAlignmentString(hsp),Standard());
	TStringIterator it_0_end = end(queryAlignmentString(hsp),Standard());

	TStringIterator it_1 = begin(databaseAlignmentString(hsp),Standard());
	TStringIterator it_1_end = end(databaseAlignmentString(hsp),Standard());
	unsigned int act0_pos = 0;
	unsigned int act1_pos = 0;
	unsigned int act_len = 0;

	//walk over (gapped) sequences, retrieve ungapped sequences and store node+edge info for alignment graph
	while ((it_0 != it_0_end) && (it_1 != it_1_end))
	{
		if(*it_0 == '-')
		{
			if(act_len > 0)
			{
				append(seq0nodes,TPosLen(act0_pos,act_len));
				append(seq1nodes,TPosLenLink(act1_pos,act_len,act0_pos));
				act0_pos += act_len;
				act1_pos += act_len;
				act_len = 0;
			}
			unsigned int gap_len = 0;
			while(it_0 != it_0_end && *it_0 == '-')
			{
				append(str1,getValue(it_1));
				++it_0;
				++it_1;
				++gap_len;
			}
			append(seq1nodes,TPosLenLink(act1_pos,gap_len,-1));
			act1_pos += gap_len;
			continue;
		}
		if(*it_1 == '-')
		{
			if(act_len > 0)
			{
				append(seq0nodes,TPosLen(act0_pos,act_len));
				append(seq1nodes,TPosLenLink(act1_pos,act_len,act0_pos));
				act0_pos += act_len;
				act1_pos += act_len;
				act_len = 0;
			}
			unsigned int gap_len = 0;
			while(it_1 != it_1_end && *it_1 == '-')
			{
				append(str0,getValue(it_0));
				++it_1;
				++it_0;
				++gap_len;
			}
			append(seq0nodes,TPosLen(act0_pos,gap_len));
			act0_pos += gap_len;
			continue;
		}
		append(str0,getValue(it_0));
		append(str1,getValue(it_1));
		++it_0;
		++it_1;
		++act_len;

	}

	if(act_len>0)
	{
		append(seq0nodes,TPosLen(act0_pos,act_len));
		append(seq1nodes,TPosLenLink(act1_pos,act_len,act0_pos));
	}

	TStringSet str;
	unsigned int id0 = assignValueById(str,str0);
	unsigned int id1 = assignValueById(str,str1);

	assignStringSet(ali,str);

	typedef typename Iterator<String<TPosLen>, Standard >::Type TNodeIterator;
	TNodeIterator nodes0it = begin(seq0nodes,Standard());
	TNodeIterator nodes0end = end(seq0nodes,Standard());

	while(nodes0it != nodes0end)
	{
		addVertex(ali,id0,nodes0it->i1,nodes0it->i2);
		++nodes0it;
	}

	typedef typename Iterator<String<TPosLenLink>, Standard >::Type TNodeEdgeIterator;
	TNodeEdgeIterator nodes1it = begin(seq1nodes,Standard());
	TNodeEdgeIterator nodes1end = end(seq1nodes,Standard());
	
	while(nodes1it != nodes1end)
	{
		TVertexDescriptor vd1 = addVertex(ali,id1,nodes1it->i1,nodes1it->i2);
		if(nodes1it->i3 >= 0)
		{
			TVertexDescriptor vd2 = findVertex(ali,id0,nodes1it->i3);
			addEdge(ali,vd1,vd2);
		}
		++nodes1it;
	}

	return length(databaseAlignmentString(hsp));

}



// general (for all specs)

template<typename TBlastSpec, typename TInfoSpec>
inline unsigned int &
queryBegin(BlastHsp<TBlastSpec, TInfoSpec>& blastHsp)
{
SEQAN_CHECKPOINT
	return blastHsp.query_begin;
}




/**
.Function.getQueryBegin:
..cat:Blast
..summary:The begin position of the HSP on the query sequence.
..signature:getQueryBegin(object);
..param.object:A Blast HSP object.
...type:Class.BlastHsp
..returns:The begin position.
...type:nolink:unsigned
..include:seqan/blast.h
*/
template<typename TBlastSpec, typename TInfoSpec>
inline unsigned int
getQueryBegin(BlastHsp<TBlastSpec, TInfoSpec>& blastHsp)
{
SEQAN_CHECKPOINT
	return blastHsp.query_begin;
}
	


template<typename TBlastSpec, typename TInfoSpec>
inline unsigned int &
databaseBegin(BlastHsp<TBlastSpec, TInfoSpec>& blastHsp)
{
SEQAN_CHECKPOINT
	return blastHsp.db_begin;
}

/**
.Function.getDatabaseBegin:
..cat:Blast
..summary:The begin position of the HSP on the database sequence.
..signature:getDatabaseBegin(object);
..param.object:A Blast HSP object.
...type:Class.BlastHsp
..returns:The begin position.
...type:nolink:unsigned
..include:seqan/blast.h
*/
template<typename TBlastSpec, typename TInfoSpec>
inline unsigned int
getDatabaseBegin(BlastHsp<TBlastSpec, TInfoSpec>& blastHsp)
{
SEQAN_CHECKPOINT
	return blastHsp.db_begin;
}
	


template<typename TBlastSpec, typename TInfoSpec>
inline unsigned int &
queryEnd(BlastHsp<TBlastSpec, TInfoSpec>& blastHsp)
{
SEQAN_CHECKPOINT
	return blastHsp.query_end;
}

/**
.Function.getQueryEnd:
..cat:Blast
..summary:The end position of the HSP on the query sequence.
..signature:getQueryEnd(object);
..param.object:A Blast HSP object.
...type:Class.BlastHsp
..returns:The end position.
...type:nolink:unsigned
..include:seqan/blast.h
*/
template<typename TBlastSpec, typename TInfoSpec>
inline unsigned int
getQueryEnd(BlastHsp<TBlastSpec, TInfoSpec>& blastHsp)
{
SEQAN_CHECKPOINT
	return blastHsp.query_end;
}
	


/**
.Function.getDatabaseEnd:
..cat:Blast
..summary:The end position of the HSP on the database sequence.
..signature:getDatabaseEnd(object);
..param.object:A Blast HSP object.
...type:Class.BlastHsp
..returns:The end position.
...type:nolink:unsigned
..include:seqan/blast.h
*/
template<typename TBlastSpec, typename TInfoSpec>
inline unsigned int &
databaseEnd(BlastHsp<TBlastSpec, TInfoSpec>& blastHsp)
{
SEQAN_CHECKPOINT
	return blastHsp.db_end;
}

template<typename TBlastSpec, typename TInfoSpec>
inline unsigned int
getDatabaseEnd(BlastHsp<TBlastSpec, TInfoSpec>& blastHsp)
{
SEQAN_CHECKPOINT
	return blastHsp.db_end;
}
	

template<typename TBlastSpec, typename TInfoSpec>
inline String<char> &
queryAlignmentString(BlastHsp<TBlastSpec, TInfoSpec>& blastHsp)
{
SEQAN_CHECKPOINT
	return blastHsp.query_string;
}


template<typename TBlastSpec, typename TInfoSpec>
inline String<char>
getQueryAlignmentString(BlastHsp<TBlastSpec, TInfoSpec>& blastHsp)
{
SEQAN_CHECKPOINT
	return blastHsp.query_string;
}

	

template<typename TBlastSpec, typename TInfoSpec>
inline String<char> &
databaseAlignmentString(BlastHsp<TBlastSpec, TInfoSpec>& blastHsp)
{
SEQAN_CHECKPOINT
	return blastHsp.db_string;
}

template<typename TBlastSpec, typename TInfoSpec>
inline String<char>
getDatabaseAlignmentString(BlastHsp<TBlastSpec, TInfoSpec>& blastHsp)
{
SEQAN_CHECKPOINT
	return blastHsp.db_string;
}


template<typename TBlastSpec, typename TInfoSpec>
inline double &
eValue(BlastHsp<TBlastSpec, TInfoSpec>& blastHsp)
{
SEQAN_CHECKPOINT
	return blastHsp.expect;
}

/**
.Function.getEValue:
..cat:Blast
..summary:The e-value associated with a Blast HSP.
..signature:getEValue(object);
..param.object:A Blast HSP object.
...type:Class.BlastHsp
..returns:The e-value.
...type:nolink:double
..include:seqan/blast.h
*/
template<typename TBlastSpec, typename TInfoSpec>
inline double
getEValue(BlastHsp<TBlastSpec, TInfoSpec>& blastHsp)
{
SEQAN_CHECKPOINT
	return blastHsp.expect;
}


// for FullInfo specs

template<typename TBlastSpec>
inline float &
score(BlastHsp<TBlastSpec, FullInfo>& blastHsp)
{
SEQAN_CHECKPOINT
	return blastHsp.score;
}

/**
.Function.BlastHsp#getBlastMatchScore:
..cat:Blast
..summary:The Smith-Waterman score associated with a Blast HSP.
..signature:getBlastMatchScore(object);
..param.object:A Blast HSP object.
...type:Spec.FullInfo
..returns:The score.
...type:nolink:float
..include:seqan/blast.h
*/
template<typename TBlastSpec>
inline float 
getBlastMatchScore(BlastHsp<TBlastSpec, FullInfo>& blastHsp)
{
SEQAN_CHECKPOINT
	return blastHsp.score;
}

template<typename TBlastSpec>
inline float &
bitScore(BlastHsp<TBlastSpec, FullInfo>& blastHsp)
{
SEQAN_CHECKPOINT
	return blastHsp.bits;
}


/**
.Function.getBitScore:
..cat:Blast
..summary:The bit score associated with a Blast HSP.
..signature:getBitScore(object);
..param.object:A Blast HSP object.
...type:Spec.FullInfo
..returns:The bit score.
...type:nolink:float
..include:seqan/blast.h
*/
template<typename TBlastSpec>
inline float 
getBitScore(BlastHsp<TBlastSpec, FullInfo>& blastHsp)
{
SEQAN_CHECKPOINT
	return blastHsp.bits;
}


template<typename TBlastSpec>
inline unsigned int &
percentIdentity(BlastHsp<TBlastSpec, FullInfo>& blastHsp)
{
SEQAN_CHECKPOINT
	return blastHsp.identity;
}



template<typename TBlastSpec>
inline unsigned int
getPercentIdentity(BlastHsp<TBlastSpec, FullInfo>& blastHsp)
{
SEQAN_CHECKPOINT
	return blastHsp.identity;
}

template<typename TBlastSpec>
inline unsigned int &
percentGaps(BlastHsp<TBlastSpec, FullInfo>& blastHsp)
{
SEQAN_CHECKPOINT
	return blastHsp.gaps;
}

template<typename TBlastSpec>
inline unsigned int
getPercentGaps(BlastHsp<TBlastSpec, FullInfo>& blastHsp)
{
SEQAN_CHECKPOINT
	return blastHsp.gaps;
}
	
template<typename TBlastSpec>
inline unsigned int &
numGaps(BlastHsp<TBlastSpec, FullInfo>& blastHsp)
{
SEQAN_CHECKPOINT
	return blastHsp.abs_gaps;
}



/**
.Function.getNumGaps:
..cat:Blast
..summary:The number of gaps within a Blast HSP alignment.
..signature:getNumGaps(object);
..param.object:A Blast HSP object.
...type:Spec.FullInfo
..returns:The number of gaps.
...type:nolink:unsigned
..include:seqan/blast.h
*/
template<typename TBlastSpec>
inline unsigned int
getNumGaps(BlastHsp<TBlastSpec, FullInfo>& blastHsp)
{
SEQAN_CHECKPOINT
	return blastHsp.abs_gaps;
}



/**
.Function.queryOrientationPlus:
..cat:Blast
..summary:Orientation of the query sequence within a Blast HSP alignment.
..signature:queryOrientationPlus(object);
..param.object:A Blast HSP object.
...type:Spec.FullInfo
..returns:True if the query is in forward orientation.
...type:nolink:bool
..include:seqan/blast.h
*/
template<typename TBlastSpec>
inline bool
queryOrientationPlus(BlastHsp<TBlastSpec, FullInfo>& blastHsp)
{
SEQAN_CHECKPOINT
	return blastHsp.query_strand;
}

/**
.Function.databaseOrientationPlus:
..cat:Blast
..summary:Orientation of the database sequence within a Blast HSP alignment.
..signature:databaseOrientationPlus(object);
..param.object:A Blast HSP object.
...type:Spec.FullInfo
..returns:True if the database sequence is in forward orientation.
...type:nolink:bool
..include:seqan/blast.h
*/
template<typename TBlastSpec>
inline bool
databaseOrientationPlus(BlastHsp<TBlastSpec, FullInfo>& blastHsp)
{
SEQAN_CHECKPOINT
	return blastHsp.db_strand;
}



// only for fullinfo protein spec
template<typename TSpec>
inline int &
queryFrame(BlastHsp<ProteinBlast<TSpec>, FullInfo>& blastHsp)
{
SEQAN_CHECKPOINT
	return blastHsp.query_frame;
}

template<typename TSpec>
inline int
getQueryFrame(BlastHsp<ProteinBlast<TSpec>, FullInfo>& blastHsp)
{
SEQAN_CHECKPOINT
	return blastHsp.query_frame;
}


template<typename TSpec>
inline int &
databaseFrame(BlastHsp<ProteinBlast<TSpec>, FullInfo>& blastHsp)
{
SEQAN_CHECKPOINT
	return blastHsp.db_frame;
}

template<typename TSpec>
inline int
getDatabaseFrame(BlastHsp<ProteinBlast<TSpec>, FullInfo>& blastHsp)
{
SEQAN_CHECKPOINT
	return blastHsp.db_frame;
}

template<typename TSpec>
inline unsigned int &
percentPositives(BlastHsp<ProteinBlast<TSpec>, FullInfo>& blastHsp)
{
SEQAN_CHECKPOINT
	return blastHsp.positives;
}

template<typename TSpec>
inline unsigned int
getPercentPositives(BlastHsp<ProteinBlast<TSpec>, FullInfo>& blastHsp)
{
SEQAN_CHECKPOINT
	return blastHsp.positives;
}

/**
.Function.length:
..cat:Blast
..param.object:
...type:Class.BlastHsp
..include:seqan/blast.h
*/
template<typename TBlast, typename TSpec>
inline unsigned int
length(BlastHsp<TBlast,TSpec >& blastHsp)
{
	return length(databaseAlignmentString(blastHsp));
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
