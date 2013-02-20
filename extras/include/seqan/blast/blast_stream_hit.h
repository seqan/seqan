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

#ifndef SEQAN_HEADER_BLAST_STREAM_HIT_H
#define SEQAN_HEADER_BLAST_STREAM_HIT_H


namespace SEQAN_NAMESPACE_MAIN
{


////////////////////////////////////////////////////////////////////////////////////////////
//  Blast Hit storing only one hsp at a time
////////////////////////////////////////////////////////////////////////////////////////////

template<typename TBlastHsp, typename TFile>
class BlastHit<TBlastHsp, StreamReport<TFile> > 
{
	public:
		typedef typename Position<TFile>::Type TPosition_;

		String<char> name;
		unsigned int length; //length of whole sequence  
		TBlastHsp act_hsp;
		TPosition_ begin_pos, first_hsp_pos;
		
		BlastReport<TBlastHsp,StreamReport<TFile> >* data_host;

	
		BlastHit()
		{
		}

		~BlastHit()
		{
		}

};




//parse BlastHit
template<typename TFile, typename TChar, typename TBlastSpec>
inline typename Position<TFile>::Type
_parseBlastHit(TFile & file,
			TChar & c, 
			BlastHit<TBlastSpec,StreamReport<TFile> > & hit)
{
//IOREV
	typedef typename Position<TFile>::Type TPosition;
	//typedef BlastHit<TBlastSpec,StreamReport<TFile> > TBlastHit;
	//typedef typename Hsp<TBlastHit>::Type TBlastHsp;

	String<char> pword;
	int pint;
	TPosition start_pos,act_pos;
	act_pos = _streamTellG(file);

	if(_parseUntilBeginLine(file,c,'>'))
	{
		hit.begin_pos = _streamTellG(file);
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
			hit.first_hsp_pos = _streamTellG(file);
			//if(_parseUntilBeginLine(file,c,'>'))
			//	return _streamTellG(file);
		}

		
	}//end hit
	_streamSeekG(file,act_pos);
	c = '>';
	return act_pos;
}







}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
