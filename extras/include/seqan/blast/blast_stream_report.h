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

#ifndef SEQAN_HEADER_BLAST_STREAM_REPORT_H
#define SEQAN_HEADER_BLAST_STREAM_REPORT_H


namespace SEQAN_NAMESPACE_MAIN
{


//TODO macht noch nicht so richtig sinn mit dem stringSet 
//+vielleicht w�r ne map von fasta id -> stringset id gut
template<typename TBlastHsp, typename TFile>
class BlastReport<TBlastHsp, StreamReport<TFile> > 
{
//IOREV maybe think about making this a proper file format specialization
	public:
		typedef BlastHit<TBlastHsp,StreamReport<TFile> > TBlastHit_;
		typedef typename Position<TFile>::Type TPosition_;
	
		String<char> query_name;
		String<char> db_name;
	
	
		TPosition_ first_hit_pos;
		char act_c; 
		bool hits_found;
		bool next_report;
		TPosition_ next_report_pos;


		BlastReport()
		{
		SEQAN_CHECKPOINT
			next_report = true;
			next_report_pos = 0;
		}

		BlastReport(BlastReport const& other)
		{
		SEQAN_CHECKPOINT

			query_name = other.query_name;
			db_name = other.db_name;
			act_c = other.act_c;
			hits_found = other.hits_found;
			first_hit_pos = other.first_hit_pos;
			next_report = other.next_report;
			next_report_pos = other.next_report_pos;
		}

		//BlastReport(TFile file)
		//{
		//	read(file,*this, Blast());
		//}

		~BlastReport()
		{
		SEQAN_CHECKPOINT
		}


};





//read a  blast report and set the filestream to the first hit position (first '>')
template<typename TBlastHsp, typename TFile>
void 
read(TFile & file,
	 BlastReport<TBlastHsp, StreamReport<TFile> >& blastObj,	 
	 Tag<TagBlast_>) 
{
//IOREV _nodoc_ specialization not documented
SEQAN_CHECKPOINT
 


	typedef typename Position<TFile>::Type TPosition;
	typedef typename Value<TFile>::Type TValue;

	TValue c = blastObj.act_c;
	if(blastObj.next_report)
		_streamSeekG(file,blastObj.next_report_pos);

	blastObj.next_report = false;

	if (_streamEOF(file)) return;
	else c = _streamGet(file);

	String<char> query_name, db_name;
	
	//get query and database names
	_parseReadQueryAndDBName(file,c,query_name,db_name);
	blastObj.query_name = query_name;
	blastObj.db_name = db_name;

	TPosition after_dbquery_pos = _streamTellG(file);
	TValue c_before = c;

	blastObj.hits_found = false;
	TPosition next_event_pos = after_dbquery_pos;

	String<char> delim = "Reference";
	if(_parseUntilBeginLine(file,c,delim,9))
	{
		blastObj.next_report_pos = _streamTellG(file);
		blastObj.next_report = true;
	}

	_streamSeekG(file,after_dbquery_pos);
	if(_parseUntilBeginLine(file,c,'>'))
	{
		next_event_pos = _streamTellG(file);
		if(!blastObj.next_report || next_event_pos < blastObj.next_report_pos)
		{
			blastObj.hits_found = true;
			blastObj.first_hit_pos = next_event_pos;
		}
		else
			next_event_pos = after_dbquery_pos;
	}
	//get some more values aber erst sp�ter
	//_readParameters(file,c,blastObj) ;

	if(blastObj.hits_found)
	{
		_streamSeekG(file,next_event_pos);
		c = '>';
	}
	else
	{
		if(blastObj.next_report)
		{
			c = ':';
			_streamSeekG(file,blastObj.next_report_pos);
		}
		else
		{
			_streamSeekG(file,next_event_pos);
			c = c_before;
		}
	}

	blastObj.act_c = c;

}



//////////////////// Metafunctions /////////////////////////////

template<typename TBlastHsp, typename TFile>
struct Value<BlastReport<TBlastHsp, StreamReport<TFile> > > 
{
	typedef BlastHit<TBlastHsp,StreamReport<TFile> > Type;
};

template<typename TBlastHsp, typename TFile>
struct Value<BlastReport<TBlastHsp, StreamReport<TFile> > const> 
{
	typedef BlastHit<TBlastHsp,StreamReport<TFile> > Type;
};


template<typename TBlastHsp, typename TFile>
struct Hit<BlastReport<TBlastHsp, StreamReport<TFile> > > 
{
	typedef BlastHit<TBlastHsp,StreamReport<TFile> > Type;
};

template<typename TBlastHsp, typename TFile>
struct Hit<BlastReport<TBlastHsp, StreamReport<TFile> > const> 
{
	typedef BlastHit<TBlastHsp,StreamReport<TFile> > Type;
};


/////////////////////////////////////////////
//todo

	template<typename TBlastHsp, typename TFile>
	struct Hit<String<BlastReport<TBlastHsp, StreamReport<TFile> > > > 
	{
		typedef BlastHit<TBlastHsp,StreamReport<TFile> > Type;
	};

	template<typename TBlastHsp, typename TFile>
	struct Hit<String<BlastReport<TBlastHsp, StreamReport<TFile> > const> > 
	{
		typedef BlastHit<TBlastHsp,StreamReport<TFile> > Type;
	};

	template<typename TBlastHsp, typename TFile>
	struct Hit<String<BlastHit<TBlastHsp, StreamReport<TFile> > > > 
	{
		typedef BlastHit<TBlastHsp,StreamReport<TFile> > Type;
	};

	template<typename TBlastHsp, typename TFile>
	struct Hit<String<BlastHit<TBlastHsp, StreamReport<TFile> > const> > 
	{
		typedef BlastHit<TBlastHsp,StreamReport<TFile> > Type;
	};

///////////////////////////////////////////////////////////////////


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
