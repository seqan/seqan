 /*==========================================================================
                SeqAn - The Library for Sequence Analysis
                          http://www.seqan.de 
 ============================================================================
  Copyright (C) 2007

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

 ============================================================================
  $Id$
 ==========================================================================*/

#ifndef SEQAN_HEADER_READ_GTF_H
#define SEQAN_HEADER_READ_GTF_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Read Gtf
//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TChar>
inline bool
_parseSkipEntryUntilWhitespace(TFile& file, TChar& c)
{
    //IOREV _duplicate_ _hasCRef_ there is similar functions; return value unclear
    if (c== ' ' || c== '\t' || c == '\n' || (c == '\r' && _streamPeek(file) != '\n')) return false;
		
    while (!_streamEOF(file)) {
        c = _streamGet(file);
        if (c== ' ' || c== '\t' || c == '\n' || (c == '\r' && _streamPeek(file) != '\n')) break;
    }
    return true; 
}

template<typename TSpec, typename TConfig>
inline void
readAnnotationsFromGTF(FragmentStore<TSpec, TConfig> & me,
		       char const * fileName)
{
//IOREV _nodoc_ unclear how this relates to the other GTF implementation in store_io
	typedef typename FragmentStore<TSpec, TConfig>::TContigPos 		TContigPos;	
	typedef typename FragmentStore<TSpec, TConfig>::TAnnotationStore 	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type 				TAnnotationStoreElement;
	typedef typename TAnnotationStoreElement::TId 				TId;
	typedef typename Iterator<StringSet<CharString> >::Type 		TContigNameIter;
	typedef  	 StringSet<CharString>					TAnnotationNameStore;
	
	typedef 	 Pair<CharString, TId> 					TPair;
	typedef 	 Map<TPair > 						TMap;			
	
	static const TId  INVALID_ID  = TAnnotationStoreElement::INVALID_ID; 
	static const TContigPos  INVALID_POS = TAnnotationStoreElement::INVALID_POS;
	
	std::ifstream file;
	file.open(fileName, std::ios_base::in | std::ios_base::binary);
	if(!file.is_open())
	{
		std::cerr << "ERROR: Gff-File could not be opened!"; 
	}
	
	char c;
	c = _streamGet(file);
	unsigned count = 0;
	while (!_streamEOF(file))
	{
		c = _streamGet(file);
		_parseSkipLine(file, c);
		++count;
	}
	file.close();
	
	resize(me.annotationNameStore, count);
	
	
	//////////////////////////////////////////////////////////////////////////////
	// build annotationNameStore:
	//////////////////////////////////////////////////////////////////////////////
	
	TMap map;
	clear(map);
	
	file.open(fileName, std::ios_base::in | std::ios_base::binary);
	if(!file.is_open())
	{ 
		std::cerr << "ERROR: Gtf-File could not be opened!";
	}
	c = _streamGet(file);
	
	CharString name;
	CharString parentName;
	TPair pair;
	
	unsigned currentLine = 0;
	while (!_streamEOF(file))
	{
		// skip first entries:
		_parseSkipWhitespace(file, c);
		_parseSkipEntryUntilWhitespace(file, c);
		_parseSkipWhitespace(file, c);
		_parseSkipEntryUntilWhitespace(file, c);
		_parseSkipWhitespace(file, c);
		_parseSkipEntryUntilWhitespace(file, c);
		_parseSkipWhitespace(file, c);
		_parseSkipEntryUntilWhitespace(file, c);
		_parseSkipWhitespace(file, c);
		_parseSkipEntryUntilWhitespace(file, c);
		_parseSkipWhitespace(file, c);
		_parseSkipEntryUntilWhitespace(file, c);
		_parseSkipWhitespace(file, c);
		_parseSkipEntryUntilWhitespace(file, c);
		_parseSkipWhitespace(file, c); 
		_parseSkipEntryUntilWhitespace(file, c);
		
		// read column 9: name & parentName
		_parseSkipWhitespace(file, c);
		if (_parseReadIdentifier(file, c) != "gene_id") 
		{
			std::cout << "first feature field should be 'gene_id'"<< std::endl; 
			_parseSkipLine(file, c);
			resize(me.annotationNameStore, length(me.annotationNameStore) - 1, Generous());
			continue;
		}
		_parseSkipWhitespace(file, c); 
		name = _parseReadIdentifier(file, c);
		
		if (c == ';') c = _streamGet(file);
		_parseSkipWhitespace(file, c); 
		
		clear(parentName);       
		if (!_streamEOF(file) && !(c == '\n' || (c == '\r' && _streamPeek(file) != '\n')))
		{
			if(_parseReadIdentifier(file,c) != "transcript_id") 
			{
				std::cout << "second feature field should be 'transcript_id'"<< std::endl;
				_parseSkipLine(file, c);
				resize(me.annotationNameStore, length(me.annotationNameStore) - 1, Generous());	
				continue;
			}
			_parseSkipWhitespace(file, c); 
			parentName = _parseReadIdentifier(file, c);
		}
	
			
		// insert in AnnotationNameStore and map:
		
		if (!empty(parentName) )
		{
			if (!hasKey(map, parentName))
			{
				resize(me.annotationNameStore, length(me.annotationNameStore) + 1, Generous());
				assignValue(me.annotationNameStore, currentLine, parentName);
				
				pair.i1 = parentName;
				pair.i2 = currentLine;
				insert(map, pair);
				++currentLine;
			}
			assignValue(me.annotationNameStore, currentLine, name);
			++currentLine;
		}
		else
		{
			if (!hasKey(map, name)) 
			{
				pair.i1 = name;
				pair.i2 = currentLine;
				insert(map, pair);
				assignValue(me.annotationNameStore, currentLine, name);
				++currentLine;
			}
			else  		// for the case that parent entry in gff is behind children entries -> no insert in annoStore, because is already there
			{
				resize(me.annotationNameStore, length(me.annotationNameStore) - 1, Generous());
			}
		}
		
		_parseSkipLine(file, c);	
	}
	file.close();

	
	//////////////////////////////////////////////////////////////////////////////
	// build annotationStore:
	//////////////////////////////////////////////////////////////////////////////
	typedef typename Iterator<TAnnotationNameStore>::Type 	TNameStoreIter;
	typedef typename Iterator<TAnnotationStore>::Type	TStoreIter;
	
	resize(me.annotationStore, length(me.annotationNameStore), Generous());
	
	TNameStoreIter itName = begin(me.annotationNameStore);
	// TNameStoreIter itNameEnd = end(me.annotationNameStore);
	
	TStoreIter itStore = begin(me.annotationStore);
	
	file.open(fileName, std::ios_base::in | std::ios_base::binary);
	if(!file.is_open())
	{ 
		std::cerr << "ERROR: Gtf-File could not be opened!";
	}
	c = _streamGet(file);

	CharString contigName;
	TId contigId;
	TContigNameIter itContigName;
	TContigNameIter itContigNameEnd;
	TContigPos beginPos;
	TContigPos endPos;
	bool orientation;

	while (!_streamEOF(file))
	{
		_parseSkipWhitespace(file, c);
		
		// read 1. column: contig name
		contigName = _parseReadWordUntilWhitespace(file, c);
		
		itContigName = begin(me.contigNameStore);
		itContigNameEnd = end(me.contigNameStore);
		while ( (itContigName != itContigNameEnd) && (getValue(itContigName) !=contigName) )
		{
			goNext(itContigName);
		}
		if (itContigName != itContigNameEnd)
		{
			contigId = position(itContigName, me.contigNameStore);
		}
		else
		{
			contigId = INVALID_ID;
		}
		
		// skip column 2 and 3
		_parseSkipWhitespace(file, c);
		_parseSkipEntryUntilWhitespace(file, c);
		_parseSkipWhitespace(file, c);
		_parseSkipEntryUntilWhitespace(file, c);
		
		// read column 4: begin position
		_parseSkipWhitespace(file, c);
		if (_parseIsDigit(c)) beginPos = _parseReadNumber(file, c);
		else 
		{
			_parseSkipEntryUntilWhitespace(file, c);
			beginPos = INVALID_POS;
		}
	
		// read column 5: end position
		_parseSkipWhitespace(file, c);
		if (_parseIsDigit(c)) endPos = _parseReadNumber(file, c);
		else 
		{
			_parseSkipEntryUntilWhitespace(file, c);
			endPos = INVALID_POS;
		}	
		
	
		// skip column 6
		_parseSkipWhitespace(file, c);
		_parseSkipEntryUntilWhitespace(file, c);
	
		// read column 7: orientation
		_parseSkipWhitespace(file, c);
		if (c == '+')
		{	
			orientation = 0;	
		}
		else
		{
			orientation = 1;
		}
		c = _streamGet(file);
	
		// skip column 8
		_parseSkipWhitespace(file, c);
		_parseSkipEntryUntilWhitespace(file, c);
		
		// read column 9: name
		_parseSkipWhitespace(file, c);
		if (_parseReadIdentifier(file, c) != "gene_id") 
		{
			_parseSkipLine(file, c);
			continue;
		}
		_parseSkipWhitespace(file, c); 
		name = _parseReadIdentifier(file, c);
		
		if (c == ';') c = _streamGet(file);
		_parseSkipWhitespace(file, c); 
		
		clear(parentName);       
		if (!_streamEOF(file) && !(c == '\n' || (c == '\r' && _streamPeek(file) != '\n')))
		{
			if(_parseReadIdentifier(file,c) != "transcript_id") 
			{
				_parseSkipLine(file, c);	
				continue;
			}
			_parseSkipWhitespace(file, c); 
			parentName = _parseReadIdentifier(file, c);
		}
		
		// insert in AnnotationStore:
		
		if (parentName == getValue(itName) )
		{
			value(itStore).contigId = contigId;
			value(itStore).parentId = INVALID_ID;
			value(itStore).beginPos = INVALID_POS;
			value(itStore).endPos = INVALID_POS;
			goNext(itName);
			goNext(itStore);
			
			value(itStore).contigId = contigId;
			value(itStore).parentId = position(itStore, me.annotationStore) - 1;
	
			if (orientation == 0)
			{
				value(itStore).beginPos = beginPos;
				value(itStore).endPos =	endPos;
			}
			else
			{
				value(itStore).beginPos = endPos;
				value(itStore).endPos =	beginPos;
			}
	
			goNext(itName);
			goNext(itStore);
		}
		
		else if (name == getValue(itName) )
		{
			value(itStore).contigId = contigId;
			
			if (!empty(parentName) )
			{
				value(itStore).parentId = cargo(map, parentName);
				
				if (getValue(me.annotationStore, value(itStore).parentId).beginPos != INVALID_POS )
				{
					value(me.annotationStore, value(itStore).parentId).beginPos = INVALID_POS;
					value(me.annotationStore, value(itStore).parentId).endPos = INVALID_POS;
				}
			}
			else value(itStore).parentId = INVALID_ID;
		
			if (orientation == 0)
			{
				value(itStore).beginPos = beginPos;
				value(itStore).endPos =	endPos;
			}
			else
			{
				value(itStore).beginPos = endPos;
				value(itStore).endPos =	beginPos;
			}
		
			goNext(itName);
			goNext(itStore);
		}
		_parseSkipLine(file, c);		
	}
	file.close();
}







}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
