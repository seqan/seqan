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

#ifndef SEQAN_HEADER_BLAST_PARSING_H
#define SEQAN_HEADER_BLAST_PARSING_H


namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Internal Blast Parsing Functions
// remark: uses a lot of functions from graph_utility_parsing.h
//////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////////////
// read alignment string (letters and gaps)
// steht am ende dahinter 
template<typename TFile, typename TChar>
inline String<char>
_parseReadAlignmentString(TFile & file, TChar& c)
{
//IOREV _nodoc_ _hasCRef_ 
SEQAN_CHECKPOINT
	// Read word
	String<char> str(c);
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if (!_parseIsLetter(c) && !(c=='-')) break;
		append(str, c);
	}
	return str;
}





////////////////////////////////////////////////////////////////////////////
// parses query and database name (file should be pointing to the beginning of the file)
template<typename TFile, typename TChar>
typename Position<TFile>::Type
_parseReadQueryAndDBName(TFile & file,
						  TChar & c,
						  String<char> & query_name,
						  String<char> & db_name)
{
//IOREV
SEQAN_CHECKPOINT
	typedef typename Position<TFile>::Type TPosition;
	//typedef typename Value<TFile>::Type TValue;

	TChar c_before = c;
	TPosition act_pos = _streamTellG(file);
	TPosition query_pos,db_pos;

	//String<char> delim = "DQ";
	//_parseUntilBeginLineOneOf(file,c,delim,2);


	String<char> query = "Query";
	if(_parseUntilBeginLine(file,c,query,6))
		query_pos = _streamTellG(file);
	else
		return act_pos;
	_streamSeekG(file,act_pos);
	c = c_before;
	String<char> database = "Database";
	if(_parseUntilBeginLine(file,c,database,8))
		db_pos = _streamTellG(file);
	else
		return act_pos;
	
	
	//getQueryName
	_streamSeekG(file,query_pos);
	_parseSkipWhitespace(file,c);
	c = _streamGet(file);
	_parseSkipWhitespace(file,c);
	query_name = _parseReadWord(file, c);
	while (!_streamEOF(file) && c != '\n' && c != '\r')
		query_name += _parseReadWord(file, c);
	
	//getDBName
	_streamSeekG(file,db_pos);
	c = _streamGet(file);
	_parseSkipWhitespace(file,c);
	db_name = _parseReadWord(file, c);
	while (!_streamEOF(file) && c != '\n' && c != '\r')
		db_name += _parseReadWord(file, c);
	_parseSkipWhitespace(file,c);
		
	c = _streamGet(file);

	return _streamTellG(file); 
	
}








}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
