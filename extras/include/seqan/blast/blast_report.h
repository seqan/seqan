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

#ifndef SEQAN_HEADER_BLAST_REPORT_H
#define SEQAN_HEADER_BLAST_REPORT_H


namespace SEQAN_NAMESPACE_MAIN
{

// 

//template<typename TBlastHsp = BlastHsp<>, typename TStoreSpec = StoreReport<BasicInfo> >
template<typename TBlastHsp = BlastHsp<>, typename TStoreSpec = StreamReport<FullInfo> >
class BlastReport;


/**
.Class.BlastReport:
..cat:Blast
..summary:Object for storing Blast report information. 
..signature:BlastReport<TBlastHsp, TSpec>  
..param.TBlastHsp:The type of HSPs that are to be stored.
...remarks:see @Class.BlastHsp@.
...default:$BlastHsp<BlastN,BasicInfo>$
..param.TSpec:The specializing type.
...type:Spec.StreamReport
...type:Spec.StoreReport
...default:StreamReport 
...remarks:StoreReport requires more space while StreamReport will be more time consuming if hits are iterated over more than once.
..metafunction:Metafunction.Hsp
..metafunction:Metafunction.Hit
...remarks:Use @Metafunction.Hit@ to get the type of BlastHit.
..remarks:Hits and HSPs can be iterated over by using @Spec.HitIterator@ and @Spec.HspIterator@.
..include:blast.h
*/

template<typename TBlastHsp>
class BlastReport<TBlastHsp, StoreReport<FullInfo> > 
{
	public:

		String<BlastHit<TBlastHsp, StoreReport<FullInfo> > > hits;
		String<char> query_name;
		String<char> db_name;
		
		bool next_report;

//Search Parameters

		double min_expect;        //todo
		bool allow_gaps;          
		String<char> matrix;
		float gap_open;
		float gap_extension;



//Search Statistics

		//S' (bits) =  [lambda * S (raw)  -  ln K] / ln 2
		float lambda;
		float k; //kappa
		float h; //entropy

		float gapped_lambda;
		float gapped_k;
		float gapped_h;



//die hier auch noch?		
		//Hits_to_DB	
		//S1	
		//S1_bits	
		//S2
		//S2_bits
		//X1
		//X1_bits
		//X2
		//X2_bits
		//X3
		//X3_bits
		//dbentries	
		//dbletters	
		//num_extensions
		//num_successful_extensions
		//posted_date	
		//seqs_better_than_cutoff



		BlastReport()
		{
		SEQAN_CHECKPOINT
			allow_gaps = false;
			next_report = true;
		}

		BlastReport(BlastReport const& other)
		{
		SEQAN_CHECKPOINT
			assign(hits,other.hits);
			query_name = other.query_name;
			db_name = other.db_name;
			next_report = other.next_report;

			allow_gaps = other.allow_gaps;
			min_expect = other.min_expect;
			matrix = other.matrix;
			gap_open = other.gap_open;
			gap_extension = other.gap_extension;

			lambda = other.lambda;
			k = other.k;			
			h = other.h;
			gapped_lambda = other.gapped_lambda;
			gapped_k = other.gapped_k;
			gapped_h = other.gapped_h;
		}

		BlastReport & operator = (BlastReport const& other)
		{
		SEQAN_CHECKPOINT
			assign(hits,other.hits);
			query_name = other.query_name;
			db_name = other.db_name;

			allow_gaps = other.allow_gaps;
			min_expect = other.min_expect;
			matrix = other.matrix;
			gap_open = other.gap_open;
			gap_extension = other.gap_extension;

			lambda = other.lambda;
			k = other.k;			
			h = other.h;
			gapped_lambda = other.gapped_lambda;
			gapped_k = other.gapped_k;
			gapped_h = other.gapped_h;
			return *this;
		}

		//BlastReport(const char * file_name,
		//			double mineval = -1.0) // -1 => alle
		//{
		//	min_expect = mineval;
		//	std::fstream strm; 
		//	strm.open(file_name,::std::ios_base::in | ::std::ios_base::binary);
		//	read(strm, *this, Blast());
		//	strm.close();

		//}

//		template<typename TFile>
//		BlastReport(TFile & file,
//					double mineval = -1.0) // -1 => alle
//		{
//			min_expect = mineval;
////			strm.open(file_name.c_str(),ios_base::in | ios_base::binary);
//			read(file, *this, Blast());
////			strm.close();
//
//		}

		~BlastReport()
		{
		SEQAN_CHECKPOINT
		}


};


template<typename TBlastHsp>
inline void
clear(BlastReport<TBlastHsp, StoreReport<FullInfo> >& blastObj)
{
SEQAN_CHECKPOINT
	
	for(unsigned int i = 0; i < length(blastObj.hits); ++i)
		clear(blastObj.hits[i]);
	resize(blastObj.hits,0);
	resize(blastObj.query_name,0);
	resize(blastObj.db_name,0);
	blastObj.allow_gaps = false;
	blastObj.next_report = true;
	resize(blastObj.matrix,0);
}





//////////////////// Metafunctions /////////////////////////////

template<typename TBlastHsp, typename TSpec>
struct Value<BlastReport<TBlastHsp, TSpec> > 
{
	typedef BlastHit<TBlastHsp,TSpec> Type;
};

template<typename TBlastHsp, typename TSpec>
struct Value<BlastReport<TBlastHsp, TSpec> const> 
{
	typedef BlastHit<TBlastHsp,TSpec> const Type;
};


template<typename TBlastHsp, typename TSpec>
struct Hit<BlastReport<TBlastHsp, TSpec> > 
{
	typedef BlastHit<TBlastHsp,TSpec> Type;
};

template<typename TBlastHsp, typename TSpec>
struct Hit<BlastReport<TBlastHsp, TSpec> const> 
{
	typedef BlastHit<TBlastHsp,TSpec> const Type;
};


///////////////////////////////////////////////////
//todo
	template<typename TBlastHsp, typename TSpec>
	struct Hit<String<BlastReport<TBlastHsp, TSpec> > > 
	{
		typedef BlastHit<TBlastHsp,TSpec> Type;
	};

	template<typename TBlastHsp, typename TSpec>
	struct Hit<String<BlastReport<TBlastHsp, TSpec> const> > 
	{
		typedef BlastHit<TBlastHsp,TSpec> Type;
	};

	template<typename TBlastHsp, typename TSpec>
	struct Hit<String<BlastHit<TBlastHsp, TSpec> > > 
	{
		typedef BlastHit<TBlastHsp,TSpec> Type;
	};

	template<typename TBlastHsp, typename TSpec>
	struct Hit<String<BlastHit<TBlastHsp, TSpec> const> > 
	{
		typedef BlastHit<TBlastHsp,TSpec> Type;
	};

	template<typename TBlastSpec, typename TInfoSpec>
	struct Hit<String<BlastHsp<TBlastSpec, TInfoSpec> > > 
	{
		typedef BlastHit<BlastHsp<TBlastSpec,TInfoSpec>,StoreReport<TInfoSpec> > Type;
	};

	template<typename TBlastSpec, typename TInfoSpec>
	struct Hit<String<BlastHsp<TBlastSpec, TInfoSpec> const> > 
	{
		typedef BlastHit<BlastHsp<TBlastSpec,TInfoSpec>,StoreReport<TInfoSpec> > Type;
	};
///////////////////////////////////////////////////////////////////


template<typename TBlastHsp, typename TSpec>
struct Hsp<BlastReport<TBlastHsp, TSpec> > 
{
	typedef TBlastHsp Type;
};

template<typename TBlastHsp, typename TSpec>
struct Hsp<BlastReport<TBlastHsp, TSpec> const> 
{
	typedef TBlastHsp Type;
};



///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////



template<typename TBlastHsp>
class BlastReport<TBlastHsp, StoreReport<BasicInfo> > 
{
	public:
		typedef BlastHit<TBlastHsp, StoreReport<BasicInfo> > TBlastHit;
	
		String<TBlastHit> hits;
		String<char> query_name;
		String<char> db_name;
			
		bool next_report;

		BlastReport()
		{
			next_report = true;
		}

		BlastReport(BlastReport const& other)
		{
		SEQAN_CHECKPOINT
			assign(hits,other.hits);
			query_name = other.query_name;
			db_name = other.db_name;
			next_report = other.next_report;
		}

		BlastReport & operator = (BlastReport const& other)
		{
		SEQAN_CHECKPOINT
			assign(hits,other.hits);
			query_name = other.query_name;
			db_name = other.db_name;
			next_report = other.next_report;
			return *this;
		}

//		BlastReport(const char * file_name)
//		{
//			std::fstream strm; 
//		strm.open(file_name,::std::ios_base::in |::std::ios_base::binary);
//			read(strm, *this, Blast());
//			strm.close();
//
//		}
//
//		template<typename TFile>
//		BlastReport(TFile & file)
//		{
////			std::fstream strm; 
//			//strm.open(file_name,ios_base::in |ios_base::binary);
//			read(file, *this, Blast());
////			strm.close();
//
//		}

		~BlastReport()
		{
		SEQAN_CHECKPOINT
		}


};



template<typename TBlastHsp>
inline void
clear(BlastReport<TBlastHsp, StoreReport<BasicInfo> >& blastObj)
{
SEQAN_CHECKPOINT
	
	for(unsigned int i = 0; i < length(blastObj.hits); ++i)
		clear(blastObj.hits[i]);
	resize(blastObj.hits,0);
	resize(blastObj.query_name,0);
	resize(blastObj.db_name,0);
	blastObj.next_report = true;
}


/////////////////////////// reading //////////////////////////


template<typename TFile, typename TChar, typename TBlastHsp>
void
_readParameters(TFile & file,
				TChar & c,
                BlastReport<TBlastHsp, StoreReport<FullInfo> >& blastObj) 
{
//IOREV _nodoc_ _hasCRef_
SEQAN_CHECKPOINT
	float pfloat;

	String<char> search = "Lambda";
	if(_parseUntilBeginLine(file,c,search,6))
	{
		_parseSkipLine(file, c);
		_parseSkipWhitespace(file,c);
		pfloat = _parseReadFloat(file,c);
		blastObj.lambda = pfloat;
		_parseSkipWhitespace(file,c);
		pfloat = _parseReadFloat(file,c);
		blastObj.k = pfloat;
		_parseSkipWhitespace(file,c);
		pfloat = _parseReadFloat(file,c);
		blastObj.h = pfloat;
		if(_parseUntilBeginLine(file,c,search,6,4))
		{
			_parseSkipLine(file, c);
			_parseSkipWhitespace(file,c);
			pfloat = _parseReadFloat(file,c);
			blastObj.gapped_lambda = pfloat;
			_parseSkipWhitespace(file,c);
			pfloat = _parseReadFloat(file,c);
			blastObj.gapped_k = pfloat;
			_parseSkipWhitespace(file,c);
			pfloat = _parseReadFloat(file,c);
			blastObj.gapped_h = pfloat;
		
		}
		
	}

//	Matrix: blastn matrix:1 -3
	search = "Matrix";
	if(_parseUntilBeginLine(file,c,search,6))
	{
		c = _streamGet(file);
		_parseSkipWhitespace(file,c);
		String<char> matrix_string = _parseReadWord(file, c);
		while (!_streamEOF(file) && c != '\n' && c != '\r')
			matrix_string += _parseReadWord(file, c);
		blastObj.matrix = matrix_string;
	}


	search = "Gap";
	if(_parseUntilBeginLine(file,c,search,3))
	{
		_parseSkipWhitespace(file,c);
		String<char> pword = _parseReadWord(file,c);
		if(pword == "Penalties"){
			blastObj.allow_gaps = true;
			c = _streamGet(file);
			_parseSkipWhitespace(file,c);
			pword = _parseReadWord(file,c);
			if(pword == "Existence"){
				c = _streamGet(file);
				_parseSkipWhitespace(file,c);
				pfloat = _parseReadFloat(file,c);
				blastObj.gap_open = pfloat;
			}
			c = _streamGet(file);
			_parseSkipWhitespace(file,c);
			pword = _parseReadWord(file,c);
			if(pword == "Extension"){
				c = _streamGet(file);
				_parseSkipWhitespace(file,c);
				pfloat = _parseReadFloat(file,c);
				blastObj.gap_extension = pfloat;
			}
		}
	}

	search = "than";
	if(_parseUntil(file,c,search,4))
	{
		_parseSkipWhitespace(file,c);
		blastObj.min_expect = _parseReadEValue(file,c);
	}


}



template<typename TFile, typename TChar, typename TBlastHsp>
void
_readParameters(TFile & ,
				TChar & ,
                BlastReport<TBlastHsp, StoreReport<BasicInfo> >& ) 
{
//IOREV _nodoc_ _hasCRef_
SEQAN_CHECKPOINT
	return;
}



/**
.Function.read:
..cat:Blast
..include:seqan/blast.h
*/
template<typename TFile, typename TBlastHsp, typename TInfoSpec>
void 
read(TFile & file,
	 BlastReport<TBlastHsp, StoreReport<TInfoSpec> >& blastObj,
	 Tag<TagBlast_>) 
{
//IOREV _nodoc_ 
SEQAN_CHECKPOINT
	typedef BlastReport<TBlastHsp, StoreReport<TInfoSpec> > TBlastReport;
//	typedef BlastReport<TBlastHsp, StoreReport<FullInfo> > TBlastReport;
//	typedef TBlastReport::TBlastHit TBlastHit;
//	typedef BlastHit<TBlastHsp,StoreReport> TBlastHit;
	typedef typename Hit<TBlastReport>::Type TBlastHit;
	typedef typename Position<TFile>::Type TPosition;
	typedef typename Value<TFile>::Type TValue;

	clear(blastObj);
	blastObj.next_report = false;
	TValue c;

	if (_streamEOF(file)) return;
	else c = _streamGet(file);

	String<char> query_name, db_name;

	//get query and database names
	_parseReadQueryAndDBName(file,c,query_name,db_name);
	blastObj.query_name = query_name;
	blastObj.db_name = db_name;

	TPosition after_dbquery_pos = _streamTellG(file);
	TPosition next_report_pos = 0;

	String<char> delim = "Reference";
	if(_parseUntilBeginLine(file,c,delim,9))
	{
		next_report_pos = _streamTellG(file);
		blastObj.next_report = true;
	}

	_streamSeekG(file,after_dbquery_pos);

	//get all hits
	if(_parseUntilBeginLine(file,c,'>'))
	{
		if(!blastObj.next_report || _streamTellG(file) < next_report_pos)
		{
			bool hits_left = true;
			TPosition act_hit_pos;
			TPosition next_hit_pos;
			while(hits_left){
				act_hit_pos = _streamTellG(file);
				TBlastHit hit;
				next_hit_pos = _parseBlastHit(file,c,hit);
				//resize(blastObj.hits,length(blastObj.hits)+1);
				//blastObj.hits[length(blastObj.hits)-1] = hit;
				appendValue(blastObj.hits,hit);
				if(next_hit_pos == act_hit_pos)
					hits_left = false;
				if(next_hit_pos == (TPosition)0)
				{
					hits_left = false;
					act_hit_pos = next_report_pos;
				}
			}
			_streamSeekG(file,act_hit_pos);
		}
		else
			_streamSeekG(file,next_report_pos);
	}
	
	//get some more values
	if(!blastObj.next_report)
		_readParameters(file,c,blastObj) ;


}


/////////////////////////// writing //////////////////////////
//


//template <typename TStream>
//inline void
//_streamPutFloatBlast(TStream & target,
//			float number, 
//			char const * format_string)
//{
//	char str[32];
//	sprintf(str, format_string, number);
//	_streamWrite(target, str);
//}
//
//template <typename TStream>
//inline void
//_streamPutFloatBlast(TStream & target,
//			float number)
//{
//	_streamPutFloatBlast(target, number, "%f");
//}
//
//
//template <typename TFile, typename TBlastSpec>
//inline void
//_streamPutHspInfo(TFile & target,
//				  BlastHsp<TBlastSpec,BasicInfo> const & hsp)
//{
//			_streamWrite(target,"\n\te-val: ");
//			//todo: e value ausgabe weil doof
//			_streamPutFloatBlast(target,hsp.expect);
//			_streamPut(target,'\n');
//}
//
//
//
//
//template <typename TFile, typename TBlastSpec>
//inline void
//_streamPutHspInfo(TFile & target,
//				  BlastHsp<TBlastSpec,FullInfo> const & hsp)
//{
//			_streamWrite(target,"\n\tScore: ");
//			_streamPutFloatBlast(target,(float)hsp.bits, "%.2f");
//			_streamWrite(target," bits (");
//			_streamPutFloatBlast(target,hsp.score, "%.1f");
//			_streamPut(target,')');
//			_streamWrite(target,"\n\te-val: ");
//			//todo: e value ausgabe weil doof
//			_streamPutFloatBlast(target,hsp.expect);
//			_streamPut(target,'\n');
//}
//
//
//template <typename TFile, typename TInfoSpec, typename TBlastHsp>
//inline void
//write(TFile & target,
//	  BlastReport<TBlastHsp,StoreReport<TInfoSpec> > & blastObj, 
//	  Blast)
//{
//
//	typedef BlastReport<TBlastHsp,StoreReport<TInfoSpec> > TBlastReport;
//	typedef typename Hit<TBlastReport>::Type TBlastHit;
//	//typedef BlastHit<TBlastHsp,StoreReport> TBlastHit;
//	typedef typename Iterator<TBlastReport,HitIterator>::Type THitIterator;
//	typedef typename Iterator<TBlastHit,HspIterator>::Type THspIterator;
//
//	_streamWrite(target,"Query:    ");
//	_streamWrite(target,blastObj.query_name);
//	_streamPut(target, '\n');
//
//	_streamWrite(target,"Database: ");
//	_streamWrite(target,blastObj.db_name);
//	_streamWrite(target, "\n\n");
//
//	int hit_count = 1;
//	THitIterator hit_it(blastObj); 
//	for(goBegin(hit_it); !atEnd(hit_it); goNext(hit_it)) 
//	{
//		TBlastHit hit = getValue(hit_it);
//		int hsp_count = 1;
//		_streamPutInt(target, hit_count);
//		_streamWrite(target,". Hit: ");
//		_streamWrite(target,hit.name);
//		//_streamWrite(target,(*hit_it).name);
//		_streamPut(target,'\n');
//		THspIterator hsp_it(hit);
//		for(goBegin(hsp_it); !atEnd(hsp_it); goNext(hsp_it)) 
//		{
//			TBlastHsp hsp = getValue(hsp_it);
//			_streamPut(target,'\t');
//			_streamPut(target,'\n');
//			_streamWrite(target,"    HSP ");
//			_streamPutInt(target, hsp_count);
//			_streamWrite(target,":\n\tQuery: ");
//			_streamPutInt(target,hsp.query_begin);
//			_streamWrite(target,"...");
//			_streamPutInt(target,hsp.query_end);
//			_streamWrite(target,"\n\tSbjct: ");
//			_streamPutInt(target,hsp.db_begin);
//			_streamWrite(target,"...");
//			_streamPutInt(target,hsp.db_end);
//			_streamPut(target,'\n');
//
//			_streamPutHspInfo(target,hsp);
//
//			_streamPut(target,'\n');
//			_streamWrite(target,"\tAlignment:\n\t");
//			_streamWrite(target,hsp.query_string);
//			_streamWrite(target,"\n\t");
//			_streamWrite(target,hsp.db_string);
//			_streamPut(target,'\n');
//			_streamPut(target,'\n');
//			++hsp_count;
//		}
//		_streamPut(target,'\n');
//		++hit_count;
//	}
//	
//}
//
//
//template <typename TStream, typename TBlastHsp, typename TInfoSpec>
//inline TStream &
//operator << (TStream & target, 
//			 BlastReport<TBlastHsp,StoreReport<TInfoSpec> >& source)
//{
//	write(target, source, Blast());
//	return target;
//}


///////////////////////// loads of get functions ////////////////////

// for all BlastReports

/**
.Function.queryName:
..cat:Blast
..summary:Reference to the name (identifier) of the query in a Blast report.
..signature:queryName(blastReport);
..param.blastReport:A Blast Report.
...type:Class.BlastReport
..returns:The name (identifier) of the query.
...type:Shortcut.CharString
..see:Function.getQueryName
..see:Function.databaseName
..include:seqan/blast.h
*/
template<typename TBlastHsp, typename TSpec>
inline String<char> & 
queryName(BlastReport<TBlastHsp,TSpec> & blastObj)
{
SEQAN_CHECKPOINT
	return blastObj.query_name;
}

/**
.Function.getQueryName:
..cat:Blast
..summary:The name (identifier) of the query in a Blast report.
..signature:getQueryName(blastReport);
..param.blastReport:A Blast Report.
...type:Class.BlastReport
..returns:The name (identifier) of the query.
...type:Shortcut.CharString
..see:Function.queryName
..see:Function.getDatabaseName
..include:seqan/blast.h
*/
template<typename TBlastHsp, typename TSpec>
inline String<char> 
getQueryName(BlastReport<TBlastHsp,TSpec> & blastObj)
{
SEQAN_CHECKPOINT
	return blastObj.query_name;
}		

/**
.Function.databaseName:
..cat:Blast
..summary:Reference to the name (identifier) of the database in a Blast report.
..signature:databaseName(blastReport);
..param.blastReport:A Blast Report.
...type:Class.BlastReport
..returns:The name of the database.
...type:Shortcut.CharString
..see:Function.queryName
..see:Function.getDatabaseName
..include:seqan/blast.h
*/
template<typename TBlastHsp, typename TSpec>
inline String<char> & 
databaseName(BlastReport<TBlastHsp,TSpec> & blastObj)
{
SEQAN_CHECKPOINT
	return blastObj.db_name;
}

/**
.Function.getDatabaseName:
..cat:Blast
..summary:The name (identifier) of the database in a Blast report.
..signature:getDatabaseName(blastReport);
..param.blastReport:A Blast Report.
...type:Class.BlastReport
..returns:The name of the database.
...type:Shortcut.CharString
..see:Function.queryName
..see:Function.getDatabaseName
..include:seqan/blast.h
*/
template<typename TBlastHsp, typename TSpec>
inline String<char> 
getDatabaseName(BlastReport<TBlastHsp,TSpec> & blastObj)
{
SEQAN_CHECKPOINT
	return blastObj.db_name;
}		


//for StoreReport only
/**
.Function.numHits:
..cat:Blast
..summary:The total number of hits in a Blast report.
..signature:numHits(blastReport);
..param.blastReport:A Blast Report.
...type:Spec.StoreReport
..returns:The number of hits (unsigned int).
..see:Function.numHsps
..include:seqan/blast.h
*/
template<typename TBlastHsp, typename TInfoSpec>
inline unsigned int
numHits(BlastReport<TBlastHsp,StoreReport<TInfoSpec> > & blastObj)
{
SEQAN_CHECKPOINT
	return length(blastObj.hits);
}		


/**
.Function.numHsps:
..cat:Blast
..summary:The number of HSPs for an entire Blast report or for one Blast hit.
..signature:numHits(blastObj);
..param.blastObj:A Blast report or a Blast hit.
...type:Spec.StoreReport
..returns:The number of hsps (unsigned int)
..see:Function.numHits
..include:seqan/blast.h
*/
template<typename TBlastHsp, typename TInfoSpec>
inline unsigned int
numHsps(BlastReport<TBlastHsp,StoreReport<TInfoSpec> > & blastObj)
{
SEQAN_CHECKPOINT
	unsigned int count = 0;
	for(unsigned int i = 0; i < length(blastObj.hits); ++i)
		count += numHsps(blastObj.hits[i]);
	return count;
}		
		



// /**
//.Function.eValueCutoff:
//..cat:Blast
//..summary:The Expect-Value cutoff of a Blast report.
//..signature:eValueCutoff(blastReport);
//..param.blastReport:A Blast Report.
//...type:Spec.StoreReport
//..returns:The E-value cutoff (read from the Blast report file).
//...type:double
//*/
//for StoreReport FullInfo only
template<typename TBlastHsp>
inline double & 
eValueCutoff(BlastReport<TBlastHsp,StoreReport<FullInfo> > & blastObj)
{
SEQAN_CHECKPOINT
	return blastObj.min_expect;
}


/**
.Function.getEValueCutoff:
..cat:Blast
..summary:The Expect-Value cutoff parsed from a Blast report.
..signature:getEValueCutoff(blastReport);
..param.blastReport:A Blast Report.
...type:Spec.StoreReport
..returns:The E-value cutoff (read from the Blast report file).
...type:nolink:double
..include:seqan/blast.h
*/
template<typename TBlastHsp>
inline double
getEValueCutoff(BlastReport<TBlastHsp,StoreReport<FullInfo> > & blastObj)
{
SEQAN_CHECKPOINT
	return blastObj.min_expect;
}		

template<typename TBlastHsp>
inline String<char> & 
matrixName(BlastReport<TBlastHsp,StoreReport<FullInfo> > & blastObj)
{
SEQAN_CHECKPOINT
	return blastObj.matrix;
}

template<typename TBlastHsp>
inline String<char> 
getMatrixName(BlastReport<TBlastHsp,StoreReport<FullInfo> > & blastObj)
{
SEQAN_CHECKPOINT
	return blastObj.matrix;
}		

template<typename TBlastHsp>
inline float &
gapOpen(BlastReport<TBlastHsp,StoreReport<FullInfo> > & blastObj)
{
SEQAN_CHECKPOINT
    return blastObj.gap_open;
}		

/**
.Function.getGetOpen:
..cat:Blast
..summary:The gap open penalty parsed from a Blast report.
..signature:getGapOpen(blastReport);
..param.blastReport:A Blast Report.
...type:Spec.StoreReport
..returns:The gap open penalty.
...type:nolink:float
..include:seqan/blast.h
*/
template<typename TBlastHsp>
inline float 
getGapOpen(BlastReport<TBlastHsp,StoreReport<FullInfo> > & blastObj)
{
SEQAN_CHECKPOINT
    return blastObj.gap_open;
}		

template<typename TBlastHsp>
inline float &
gapExtension(BlastReport<TBlastHsp,StoreReport<FullInfo> > & blastObj)
{
SEQAN_CHECKPOINT
	return blastObj.gap_extension;
}		


/**
.Function.getGetExtension:
..cat:Blast
..summary:The gap extension penalty parsed from a Blast report.
..signature:getGapExtension(blastReport);
..param.blastReport:A Blast Report.
...type:Spec.StoreReport
..returns:The gap open penalty.
...type:nolink:float
..include:seqan/blast.h
*/
template<typename TBlastHsp>
inline float 
getGapExtension(BlastReport<TBlastHsp,StoreReport<FullInfo> > & blastObj)
{
SEQAN_CHECKPOINT
    return blastObj.gap_extension;
}		


template<typename TBlastHsp>
inline bool
gapsAllowed(BlastReport<TBlastHsp,StoreReport<FullInfo> > & blastObj)
{
SEQAN_CHECKPOINT
	return blastObj.allow_gaps;
}		


template<typename TBlastHsp>
inline float &
lambda(BlastReport<TBlastHsp,StoreReport<FullInfo> > & blastObj)
{
SEQAN_CHECKPOINT
	return blastObj.lambda;
}		

/**
.Function.getLambda:
..cat:Blast
..summary:The lambda value parsed from a Blast report.
..signature:getLambda(blastReport);
..param.blastReport:A Blast Report.
...type:Spec.StoreReport
..returns:The lambda value.
...type:nolink:float
..include:seqan/blast.h
*/
template<typename TBlastHsp>
inline float
getLambda(BlastReport<TBlastHsp,StoreReport<FullInfo> > & blastObj)
{
SEQAN_CHECKPOINT
	return blastObj.lambda;
}		


template<typename TBlastHsp>
inline float &
gappedLambda(BlastReport<TBlastHsp,StoreReport<FullInfo> > & blastObj)
{
SEQAN_CHECKPOINT
	return blastObj.gapped_lambda;
}		

/**
.Function.getGappedLambda:
..cat:Blast
..summary:The lambda value parsed from a Blast report.
..signature:getGappedLambda(blastReport);
..param.blastReport:A Blast Report.
...type:Spec.StoreReport
..returns:The gapped lambda value.
...type:nolink:float
..include:seqan/blast.h
*/
template<typename TBlastHsp>
inline float
getGappedLambda(BlastReport<TBlastHsp,StoreReport<FullInfo> > & blastObj)
{
SEQAN_CHECKPOINT
	return blastObj.gapped_lambda;
}		



template<typename TBlastHsp>
inline float &
kappa(BlastReport<TBlastHsp,StoreReport<FullInfo> > & blastObj)
{
SEQAN_CHECKPOINT
	return blastObj.k;
}		


/**
.Function.getKappa:
..cat:Blast
..summary:The kappa value parsed from a Blast report.
..signature:getKappa(blastReport);
..param.blastReport:A Blast Report.
...type:Spec.StoreReport
..returns:The kappa value.
...type:nolink:float
..include:seqan/blast.h
*/
template<typename TBlastHsp>
inline float
getKappa(BlastReport<TBlastHsp,StoreReport<FullInfo> > & blastObj)
{
SEQAN_CHECKPOINT
	return blastObj.k;
}		


template<typename TBlastHsp>
inline float &
gappedKappa(BlastReport<TBlastHsp,StoreReport<FullInfo> > & blastObj)
{
SEQAN_CHECKPOINT
	return blastObj.gapped_k;
}		

/**
.Function.getGappedKappa:
..cat:Blast
..summary:The gapped kappa value parsed from a Blast report.
..signature:getGappedKappa(blastReport);
..param.blastReport:A Blast Report.
...type:Spec.StoreReport
..returns:The gapped kappa value.
...type:nolink:float
..include:seqan/blast.h
*/
template<typename TBlastHsp>
inline float
getGappedKappa(BlastReport<TBlastHsp,StoreReport<FullInfo> > & blastObj)
{
SEQAN_CHECKPOINT
	return blastObj.gapped_k;
}		


template<typename TBlastHsp>
inline float &
entropy(BlastReport<TBlastHsp,StoreReport<FullInfo> > & blastObj)
{
SEQAN_CHECKPOINT
	return blastObj.h;
}		

template<typename TBlastHsp>
inline float
getEntropy(BlastReport<TBlastHsp,StoreReport<FullInfo> > & blastObj)
{
SEQAN_CHECKPOINT
	return blastObj.h;
}		


template<typename TBlastHsp>
inline float &
gappedEntropy(BlastReport<TBlastHsp,StoreReport<FullInfo> > & blastObj)
{
SEQAN_CHECKPOINT
	return blastObj.gapped_h;
}		

template<typename TBlastHsp>
inline float
getGappedEntropy(BlastReport<TBlastHsp,StoreReport<FullInfo> > & blastObj)
{
SEQAN_CHECKPOINT
	return blastObj.gapped_h;
}		





/**
.Function.atEnd:
..cat:Blast
..signature:atEnd(file,object);
..param.file:A file stream.
..param.it:A Blast report.
...type:Class.BlastReport
..see:Function.getNext
..include:seqan/blast.h
*/
template<typename TBlastHsp, typename TSpec, typename TFile>
inline bool
atEnd(TFile & file,
	  BlastReport<TBlastHsp,TSpec> & blast)
{
SEQAN_CHECKPOINT

	return !blast.next_report || _streamEOF(file);
}


/**
.Function.getNext:
..cat:Blast
..summary:Get next Blast report from a file containing multiple Blast reports.
..signature:getNext(file,blastReport);
..param.file:A file stream.
..param.blastReport:A Blast report.
...type:Class.BlastReport
..returns:The result is stored in parameter blastReport.
..see:Function.atEnd
..include:seqan/blast.h
*/
template<typename TBlastHsp, typename TSpec, typename TFile>
inline void
getNext(TFile & file,
        BlastReport<TBlastHsp,TSpec> & blast)
{
SEQAN_CHECKPOINT
	read(file,blast,Blast());
}







}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
