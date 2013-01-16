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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_HEADER_STORE_IO_H
#define SEQAN_HEADER_STORE_IO_H

#include <seqan/misc/misc_parsing.h>

/* IOREV
 *
 * _doc_
 *
 *
 * if this file is about the amos file format why isn't it named accordingly?
 *
 * altogether it is unclear why sequence io is in file/ but store io is in
 * store/
 *
 */



namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// File tags
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.File Format.tag.Amos message file:
	Amos message file.
..include:seqan/store.h
*/
struct TagAmos_;
typedef Tag<TagAmos_> const Amos;


//////////////////////////////////////////////////////////////////////////////
// Auxillary functions
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////


/**
.Function.getClrRange
..cat:Fragment Store
..summary:Get the "clear" range of a read alignment.
..remarks:The clear range of a read alignment is the range of the part of the alignment that is not clipped.
..signature:getClrRange(fragStore, alignEl, begClr, endClr)
..param.fragStore:Fragment Store to work on.
...type:Class.FragmentStore
..param.alignEl:Read alignment element.
...type:Class.AlignedReadStoreElement
..param.begClr:Start of the clear range.
..param.endClr:End of the clear range.
..include:seqan/store.h
*/

template <typename TSpec, typename TConfig, typename TPos, typename TGapAnchor, typename TSpecAlign, typename TBeginClr, typename TEndClr>
inline void
getClrRange(FragmentStore<TSpec, TConfig> const& fragStore,
			AlignedReadStoreElement<TPos, TGapAnchor, TSpecAlign> const& alignEl,
			TBeginClr& begClr,		// Out-parameter: left / begin position of the clear range
			TEndClr& endClr)		// Out-parameter: right / end position of the clear range
{
	typedef FragmentStore<TSpec, TConfig> TFragmentStore;
	typedef typename Size<TFragmentStore>::Type TSize;
	typedef typename Iterator<String<TGapAnchor>, Standard>::Type TGapIter;
	
	TSize lenRead = length(fragStore.readSeqStore[alignEl.readId]);
	TGapIter itGap = begin(alignEl.gaps, Standard() );
	TGapIter itGapEnd = end(alignEl.gaps, Standard() );
	
	// Any gaps or clipped characters?
	if (itGap == itGapEnd) {
		begClr = 0;
		endClr = lenRead;
	} else {
		// Begin clear range
		begClr = (itGap->gapPos == 0) ? itGap->seqPos : 0;
		// End clear range
		--itGapEnd;
		if (static_cast<TSize>(itGapEnd->seqPos) != lenRead) endClr = lenRead;
		else {
			int diff = (itGap != itGapEnd) ? (*(itGapEnd - 1)).gapPos - (*(itGapEnd-1)).seqPos : 0;
			int newDiff = itGapEnd->gapPos - itGapEnd->seqPos;
			endClr = (newDiff < diff) ? lenRead - (diff - newDiff) : lenRead;	
		}
	}

	// For reverse reads adapt clear ranges
	if (alignEl.beginPos > alignEl.endPos) {
		TBeginClr tmp = begClr;
		begClr = lenRead - endClr;
		endClr = lenRead - tmp;
	}
}



//////////////////////////////////////////////////////////////////////////////
// Read / Write of AMOS message files (*.afg)
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

/**
.Function.read
..cat:Fragment Store
..signature:read(file, fragStore, tag)
..param.fragStore:A fragment store. Currently @Tag.File Format.tag.Amos message file@ and @Tag.File Format.tag.Sam@ formats are supported.
...type:Class.FragmentStore
..param.tag.type:Tag.File Format.tag.Amos message file
..include:seqan/store.h
*/

template<typename TFile, typename TSpec, typename TConfig>
inline void
read(TFile & file,
	 FragmentStore<TSpec, TConfig>& fragStore,
	 Amos)
{
//IOREV
	SEQAN_CHECKPOINT
	// Basic types
	typedef FragmentStore<TSpec, TConfig> TFragmentStore;
	typedef typename Id<TFragmentStore>::Type TId;
	typedef typename Size<TFragmentStore>::Type TSize;
	typedef typename Value<TFile>::Type TValue;
	typedef typename TFragmentStore::TReadSeq TReadSeq;

	// All fragment store element types
	typedef typename Value<typename TFragmentStore::TContigStore>::Type TContigElement;
	typedef typename Value<typename TFragmentStore::TLibraryStore>::Type TLibraryStoreElement;
	typedef typename Value<typename TFragmentStore::TMatePairStore>::Type TMatePairElement;
	typedef typename Value<typename TFragmentStore::TReadStore>::Type TReadStoreElement;
	typedef typename Value<typename TFragmentStore::TAlignedReadStore>::Type TAlignedElement;

	// All maps to mirror file ids to our ids
	typedef std::map<TId, TSize> TIdMap;
	// The following maps the library/fragment/read id from AMOS into the fragment store's ids.
	TIdMap libIdMap;
	TIdMap frgIdMap;
	TIdMap readIdMap;
	// For all paired reads (inferred from FRG), a mapping from the AMOS read id to the paired match id of the one
	// alignment read from the AMOS file.  Note that this has the assumption that there is only one alignment per read
	// in the AMOS file.
	TIdMap readToPairMatchId;
	// The id of the next pair match.
	unsigned nextPairMatchId = 0;

	// Parse the file and convert the internal ids
	TValue c;
	if (_streamEOF(file)) return;
	else c = _streamGet(file);
	while (!_streamEOF(file)) {
		if (_streamEOF(file)) break;

		// New block?
		if (c == '{') {
			c = _streamGet(file);
			String<char> blockIdentifier;
			_parseReadAlnum(file, blockIdentifier, c);
			_parseSkipLine(file, c);

			// Library block
			if (blockIdentifier == "LIB") {
				TLibraryStoreElement libEl;
				TId _id = 0;
				String<char> fieldIdentifier;
				String<char> eid;
				while (c != '}') {
					clear(fieldIdentifier);
					_parseReadAlnum(file, fieldIdentifier, c);
					if (fieldIdentifier == "iid") {
						c = _streamGet(file);
						_id = _parseReadNumber(file, c);
						_parseSkipLine(file, c);
					} else if (fieldIdentifier == "eid") {
						c = _streamGet(file);
						while ((c != '\n') && (c != '\r')) {
							appendValue(eid, c, Generous());
							c = _streamGet(file);
						}
						_parseSkipLine(file, c);
					} else if (fieldIdentifier == "mea") {
						c = _streamGet(file);
						libEl.mean = _parseReadDouble(file, c);
						_parseSkipLine(file, c);
					} else if (fieldIdentifier == "std") {
						c = _streamGet(file);
						libEl.std = _parseReadDouble(file, c);
						_parseSkipLine(file, c);
					} else {
						_parseSkipLine(file, c);
					}
				}
				libIdMap.insert(std::make_pair(_id, length(fragStore.libraryStore)));
				appendValue(fragStore.libraryStore, libEl, Generous() );
				appendValue(fragStore.libraryNameStore, eid, Generous() );
			} else if (blockIdentifier == "FRG") {  // Fragment block
				TMatePairElement matePairEl;
				TId _id = 0;
				String<char> fieldIdentifier;
				String<char> eid;
				bool foundRds = false;
				while (c != '}') {
					clear(fieldIdentifier);
					_parseReadAlnum(file, fieldIdentifier, c);
					if (fieldIdentifier == "iid") {
						c = _streamGet(file);
						_id = _parseReadNumber(file, c);
						_parseSkipLine(file, c);
					} else if (fieldIdentifier == "eid") {
						c = _streamGet(file);
						while ((c != '\n') && (c != '\r')) {
							appendValue(eid, c, Generous() );
							c = _streamGet(file);
						}
						_parseSkipLine(file, c);
					} else if (fieldIdentifier == "lib") {
						c = _streamGet(file);
						matePairEl.libId = _parseReadNumber(file, c);
						_parseSkipLine(file, c);
					} else if (fieldIdentifier == "rds") {
						foundRds = true;
						c = _streamGet(file);
						matePairEl.readId[0] = _parseReadNumber(file, c);
						c = _streamGet(file);  // Skip comma.
						matePairEl.readId[1] = _parseReadNumber(file, c);
						_parseSkipLine(file, c);
						// Store mapping to pair match id.
						readToPairMatchId[matePairEl.readId[0]] = nextPairMatchId;
						readToPairMatchId[matePairEl.readId[1]] = nextPairMatchId;
						nextPairMatchId += 1;
					} else {
						_parseSkipLine(file, c);
					}
				}
				// Only insert valid mate pairs
				if (foundRds) {
					frgIdMap.insert(std::make_pair(_id, length(fragStore.matePairStore)));
					appendValue(fragStore.matePairStore, matePairEl, Generous() );
					appendValue(fragStore.matePairNameStore, eid, Generous() );
				}
			} else if (blockIdentifier == "RED") {   // Read block
				TId _id = 0;
				String<char> fieldIdentifier;
				String<char> eid;
				String<char> qual;
				// If matePairId is not updated, this yields to a singleton read below.
				TId matePairId = TReadStoreElement::INVALID_ID;
				TReadSeq seq;
				while (c != '}') {
					clear(fieldIdentifier);
					_parseReadAlnum(file, fieldIdentifier, c);
					if (fieldIdentifier == "iid") {
						c = _streamGet(file);
						_id = _parseReadNumber(file, c);
						_parseSkipLine(file, c);
					} else if (fieldIdentifier == "eid") {
						c = _streamGet(file);
						while ((c != '\n') && (c != '\r')) {
							appendValue(eid, c, Generous() );
							c = _streamGet(file);
						}
						_parseSkipLine(file, c);
					} else if (fieldIdentifier == "frg") {
						c = _streamGet(file);
						matePairId = _parseReadNumber(file, c);
						_parseSkipLine(file, c);
					} else if (fieldIdentifier == "seq") {
						c = _streamGet(file);
						_parseSkipWhitespace(file, c);
						while (c != '.') {
							_parseReadSequenceData(file,c,seq);
							_parseSkipWhitespace(file, c);
						}
					} else if (fieldIdentifier == "qlt") {
						clear(qual);
						c = _streamGet(file);
						_parseSkipWhitespace(file, c);
						while (c != '.') {
							if ((c!=' ') && (c != '\t') && (c != '\n') && (c != '\r')) appendValue(qual, c, Generous() );
							c = _streamGet(file);
						}
					} else {
						_parseSkipLine(file, c);
					}
				}
				// Set quality
				assignQualities(seq, qual);

				// Insert the read
				readIdMap.insert(std::make_pair(_id, length(fragStore.readStore)));
				appendRead(fragStore, seq, matePairId);
				appendValue(fragStore.readNameStore, eid, Generous() );
			} else if (blockIdentifier == "CTG") {   // Contig block
				TContigElement contigEl;
				TSize fromAligned = length(fragStore.alignedReadStore);
				// TId _id = 0;
				String<char> fieldIdentifier;
				String<char> eid;
				String<char> contigSeq;
				String<char> contigQual;
				while (c != '}') {
					// Are we entering a TLE block
					if (c == '{') {
						TAlignedElement alignEl;
						String<char> fdIdentifier;
						typedef typename TFragmentStore::TContigPos TContigPos;
						TContigPos offsetPos = 0;
						TContigPos clr1 = 0;
						TContigPos clr2 = 0;
						String<TContigPos> gaps;
						while (c != '}') {
							clear(fdIdentifier);
							_parseReadAlnum(file, fdIdentifier, c);
							if (fdIdentifier == "src") {
								c = _streamGet(file);
								alignEl.readId = _parseReadNumber(file, c);
								_parseSkipLine(file, c);
							} else if (fdIdentifier == "off") {
								c = _streamGet(file);
								if (c != '-') offsetPos = _parseReadNumber(file, c);
								else offsetPos = 0;
								_parseSkipLine(file, c);
							} else if (fdIdentifier == "clr") {
								c = _streamGet(file);
								clr1 = _parseReadNumber(file, c);
								c = _streamGet(file);
								clr2 = _parseReadNumber(file, c);
								_parseSkipLine(file, c);
							} else if (fdIdentifier == "gap") {
								c = _streamGet(file);
								_parseSkipWhitespace(file, c);
								while (c != '.') {
									if ((c!=' ') && (c != '\t') && (c != '\n') && (c != '\r')) {
										TSize nextGap = _parseReadNumber(file, c);
										appendValue(gaps, nextGap, Generous() );
									}
									c = _streamGet(file);
								}
							} else {
								_parseSkipLine(file, c);
							}
						}
						_parseSkipLine(file, c);

						// Get the length of the read
						TId readId = (readIdMap.find(alignEl.readId))->second;
						TSize lenRead = length(value(fragStore.readSeqStore, readId));

						// Create the gap anchors
						typedef typename TFragmentStore::TContigGapAnchor TContigGapAnchor;
						int offset = 0;
						if ((clr1 < clr2) && (clr1>0)) offset = clr1;
						else if ((clr1 > clr2) && (clr1 < static_cast<TContigPos>(lenRead))) offset = lenRead - clr1;
						int diff = -1 * (int) (offset);
						// Clipped begin
						if (offset != 0) appendValue(alignEl.gaps, TContigGapAnchor(offset, 0), Generous() );
						// Internal gaps
						typedef typename Iterator<String<TContigPos>, Standard>::Type TPosIter;
						TPosIter posIt = begin(gaps, Standard());
						TPosIter posItEnd = end(gaps, Standard());
						TContigPos lastGap = 0;
						TSize gapLen = 0;
						TSize totalGapLen = 0;
						for(;posIt!=posItEnd; goNext(posIt)) {
							if (gapLen == 0) {
								++gapLen; ++totalGapLen;
								++diff;
								lastGap = value(posIt);
							}
							else if (lastGap == value(posIt)) {
								++gapLen; ++totalGapLen;
								++diff;
							}
							else {
								appendValue(alignEl.gaps, TContigGapAnchor(offset + lastGap, offset + lastGap + diff), Generous() );
								gapLen = 1; ++totalGapLen;
								lastGap = value(posIt);
								++diff;
							}
						}
						if (gapLen > 0) appendValue(alignEl.gaps, TContigGapAnchor(offset + lastGap, offset + lastGap + diff), Generous() );
						// Clipped end
						if ((clr1 < clr2) && (clr2 < static_cast<TContigPos>(lenRead))) {
							diff -= (lenRead - clr2);
							appendValue(alignEl.gaps, TContigGapAnchor(lenRead, lenRead + diff), Generous() );
						} else if ((clr1 > clr2) && (clr2 > 0)) {
							diff -= clr2;
							appendValue(alignEl.gaps, TContigGapAnchor(lenRead, lenRead + diff), Generous() );
						}

						// Set begin and end position
						if (clr1 < clr2) {
							alignEl.beginPos = offsetPos;
							alignEl.endPos = offsetPos + totalGapLen + (clr2 - clr1);
						} else {
							alignEl.beginPos = offsetPos + totalGapLen + (clr1 - clr2);
							alignEl.endPos = offsetPos;
						}

						// Append new align fragment, note: contigId must still be set
						alignEl.id = length(fragStore.alignedReadStore);
						appendValue(fragStore.alignedReadStore, alignEl, Generous() );
					} else {
						clear(fieldIdentifier);
						_parseReadAlnum(file, fieldIdentifier, c);
						if (fieldIdentifier == "iid") {
							c = _streamGet(file);
							//_id = _parseReadNumber(file, c);
							_parseSkipLine(file, c);
						} else if (fieldIdentifier == "eid") {
							c = _streamGet(file);
							while ((c != '\n') && (c != '\r')) {
								appendValue(eid, c, Generous() );
								c = _streamGet(file);
							}
							_parseSkipLine(file, c);
						} else if (fieldIdentifier == "seq") {
							c = _streamGet(file);
							_parseSkipWhitespace(file, c);
							while (c != '.') {
								do {
									_parseReadSequenceData(file,c,contigSeq);
								} while (c == '-');
								_parseSkipWhitespace(file, c);
							}
						} else if (fieldIdentifier == "qlt") {
							c = _streamGet(file);
							_parseSkipWhitespace(file, c);
							while (c != '.') {
								if ((c!=' ') && (c != '\t') && (c != '\n') && (c != '\r')) {
									appendValue(contigQual, c, Generous() );
								}
								c = _streamGet(file);
							}
						} else {
							_parseSkipLine(file, c);
						}
					}
				}

				// Create the gap anchors
				char gapChar = gapValue<char>();
				typedef typename Iterator<String<char> >::Type TStringIter;
				TStringIter seqIt = begin(contigSeq);
				TStringIter seqItEnd = end(contigSeq);
				TStringIter qualIt = begin(contigQual);
				typedef typename TFragmentStore::TReadPos TPos;
				typedef typename TFragmentStore::TContigGapAnchor TContigGapAnchor;
				TPos ungappedPos = 0;
				TPos gappedPos = 0;
				bool gapOpen = false;
				for(;seqIt != seqItEnd; goNext(seqIt), goNext(qualIt), ++gappedPos) {
					if (value(seqIt) == gapChar)
                    {
					    gapOpen = true;
                    }
					else
					{
						if (gapOpen)
						{
							appendValue(contigEl.gaps, TContigGapAnchor(ungappedPos, gappedPos), Generous() );
							gapOpen = false;
						}
						Dna5Q letter = value(seqIt);
						assignQualityValue(letter, value(qualIt));
						appendValue(contigEl.seq, letter, Generous() );
						++ungappedPos;
					}
				}
				if (gapOpen)
				    appendValue(contigEl.gaps, TContigGapAnchor(ungappedPos, gappedPos), Generous() );

				// Set the contigId in all aligned reads
				TSize toAligned = length(fragStore.alignedReadStore);
				TId newContigId = length(fragStore.contigStore);
				for (; fromAligned < toAligned; ++fromAligned)
					fragStore.alignedReadStore[fromAligned].contigId = newContigId;

				// Insert the contig
				appendValue(fragStore.contigStore, contigEl, Generous() );
				appendValue(fragStore.contigNameStore, eid, Generous() );
			} else {
				_parseSkipLine(file, c);
			}
		} else {
			_parseSkipLine(file, c);
		}
	}

	// Renumber all ids
	typedef typename TIdMap::const_iterator TIdMapIter;
	typedef typename Iterator<typename TFragmentStore::TMatePairStore>::Type TMateIter;
	TMateIter mateIt = begin(fragStore.matePairStore);
	TMateIter mateItEnd = end(fragStore.matePairStore);
	for(;mateIt != mateItEnd; goNext(mateIt)) {
		if (mateIt->libId != TMatePairElement::INVALID_ID) {
			TIdMapIter libIdPos = libIdMap.find(mateIt->libId);
			if (libIdPos != libIdMap.end())
			    mateIt->libId = libIdPos->second;
			else
			    mateIt->libId = TMatePairElement::INVALID_ID;
		}
		if (mateIt->readId[0] != TMatePairElement::INVALID_ID) {
			TIdMapIter readIdPos = readIdMap.find(mateIt->readId[0]);
			if (readIdPos != readIdMap.end())
			    mateIt->readId[0] = readIdPos->second;
			else
			    mateIt->readId[0] = TMatePairElement::INVALID_ID;
		}
		if (mateIt->readId[1]!= TMatePairElement::INVALID_ID) {
			TIdMapIter readIdPos = readIdMap.find(mateIt->readId[1]);
			if (readIdPos != readIdMap.end())
			    mateIt->readId[1] = readIdPos->second;
			else
			    mateIt->readId[0] = TMatePairElement::INVALID_ID;
		}
	}

	// Copy data from frgIdMap into the matePairId members of the readStore.
	typedef typename Iterator<typename TFragmentStore::TReadStore>::Type TReadIter;
	TReadIter readIt = begin(fragStore.readStore);
	TReadIter readItEnd = end(fragStore.readStore);
	for (;readIt != readItEnd; goNext(readIt))
	{
		if (readIt->matePairId != TReadStoreElement::INVALID_ID)
		{
			TIdMapIter mateIdPos = frgIdMap.find(readIt->matePairId);
			if (mateIdPos != frgIdMap.end())
			    readIt->matePairId = mateIdPos->second;
			else
			    readIt->matePairId = TReadStoreElement::INVALID_ID;
		}
	}

	// Copy data from readIdMap into the pairMatchId entries of the alignedReadStore.
	typedef typename Iterator<typename TFragmentStore::TAlignedReadStore>::Type TAlignIter;
	TAlignIter alignIt = begin(fragStore.alignedReadStore);
	TAlignIter alignItEnd = end(fragStore.alignedReadStore);
	for (;alignIt != alignItEnd; goNext(alignIt))
	{
		if (alignIt->readId != TAlignedElement::INVALID_ID)
		{
			TIdMapIter readIdPos = readIdMap.find(alignIt->readId);
			if (readIdPos != readIdMap.end())
            {
                //SEQAN_ASSERT(readToPairMatchId.find(alignIt->readId) != readToPairMatchId.end());
                if (readToPairMatchId.find(alignIt->readId) != readToPairMatchId.end())
                    alignIt->pairMatchId = readToPairMatchId[alignIt->readId];
			    alignIt->readId = readIdPos->second;
            }
			else
            {
			    alignIt->readId = TAlignedElement::INVALID_ID;
            }
		}
	}
}


//////////////////////////////////////////////////////////////////////////////

/**
.Function.write
..cat:Fragment Store
..signature:write(file, fragStore, tag)
..param.fragStore:A fragment store.
...type:Class.FragmentStore
..param.tag.type:Tag.File Format.tag.Amos message file
..include:seqan/store.h
*/

template<typename TFile, typename TSpec, typename TConfig>
inline void 
write(TFile & target,
	  FragmentStore<TSpec, TConfig>& fragStore,
	  Amos) 
{
//IOREV
	SEQAN_CHECKPOINT
	// Basic types
	typedef FragmentStore<TSpec, TConfig> TFragmentStore;
	typedef typename Id<TFragmentStore>::Type TId;
	typedef typename Size<TFragmentStore>::Type TSize;
	typedef typename Value<TFile>::Type TValue;

	// All fragment store element types
	typedef typename Value<typename TFragmentStore::TContigStore>::Type TContigElement;
	typedef typename Value<typename TFragmentStore::TLibraryStore>::Type TLibraryStoreElement;
	typedef typename Value<typename TFragmentStore::TMatePairStore>::Type TMatePairElement;
	typedef typename Value<typename TFragmentStore::TReadStore>::Type TReadStoreElement;
	typedef typename Value<typename TFragmentStore::TAlignedReadStore>::Type TAlignedElement;

	// Write Header
	_streamWrite(target,"{UNV\niid:1\neid:seqan\ncom:\nafg file created with SeqAn\n.\n}\n");
	
	// Write Libraries
	typedef typename Iterator<typename TFragmentStore::TLibraryStore, Standard>::Type TLibIter;
	TLibIter libIt = begin(fragStore.libraryStore, Standard() );
	TLibIter libItEnd = end(fragStore.libraryStore, Standard() );
	bool noNamesPresent = (length(fragStore.libraryNameStore) == 0);
	for(TSize idCount = 0;libIt != libItEnd; goNext(libIt), ++idCount) {
		_streamWrite(target,"{LIB\n");
		_streamWrite(target,"iid:");
		_streamPutInt(target, idCount + 1);
		_streamPut(target, '\n');
		if (!noNamesPresent) {
			_streamWrite(target,"eid:");
			_streamWrite(target, value(fragStore.libraryNameStore, idCount));
			_streamPut(target, '\n');
		}
		_streamWrite(target,"{DST\n");
		_streamWrite(target,"mea:");
		_streamPutFloat(target, libIt->mean);
		_streamPut(target, '\n');
		_streamWrite(target,"std:");
		_streamPutFloat(target, libIt->std);
		_streamPut(target, '\n');
		_streamWrite(target,"}\n");	
		_streamWrite(target,"}\n");
	}

	// Write Fragments / mate pairs
	typedef typename Iterator<typename TFragmentStore::TMatePairStore, Standard>::Type TMateIter;
	TMateIter mateIt = begin(fragStore.matePairStore, Standard() );
	TMateIter mateItEnd = end(fragStore.matePairStore, Standard() );
	noNamesPresent = (length(fragStore.matePairNameStore) == 0);
	for(TSize idCount = 0;mateIt != mateItEnd; goNext(mateIt), ++idCount) {
		_streamWrite(target,"{FRG\n");
		_streamWrite(target,"iid:");
		_streamPutInt(target, idCount + 1);
		_streamPut(target, '\n');
		if (!noNamesPresent) {
			_streamWrite(target,"eid:");
			_streamWrite(target, value(fragStore.matePairNameStore, idCount));
			_streamPut(target, '\n');
		}
		_streamWrite(target,"lib:");
		_streamPutInt(target, mateIt->libId + 1);
		_streamPut(target, '\n');
		if ((mateIt->readId[0] != TMatePairElement::INVALID_ID) && (mateIt->readId[1] != TMatePairElement::INVALID_ID)) {
			_streamWrite(target,"rds:");
			_streamPutInt(target, mateIt->readId[0] + 1);
			_streamPut(target, ',');
			_streamPutInt(target, mateIt->readId[1] + 1);
			_streamPut(target, '\n');
		}
		_streamWrite(target,"}\n");
	}

	// Get clear ranges
	typedef Pair<typename TFragmentStore::TReadPos, typename TFragmentStore::TReadPos> TClrRange;
	String<TClrRange> clrRange;
	resize(clrRange, length(fragStore.readStore), TClrRange(0,0));
	typedef typename Iterator<typename TFragmentStore::TAlignedReadStore, Standard>::Type TAlignIter;
	TAlignIter alignIt = begin(fragStore.alignedReadStore, Standard() );
	TAlignIter alignItEnd = end(fragStore.alignedReadStore, Standard() );
	for(;alignIt != alignItEnd; goNext(alignIt)) {
		typename TFragmentStore::TReadPos begClr = 0;
		typename TFragmentStore::TReadPos endClr = 0;
		getClrRange(fragStore, value(alignIt), begClr, endClr);
		value(clrRange, alignIt->readId) = TClrRange(begClr, endClr);
	}

	// Write reads
	typedef typename Iterator<typename TFragmentStore::TReadStore, Standard>::Type TReadIter;
	TReadIter readIt = begin(fragStore.readStore, Standard() );
	TReadIter readItEnd = end(fragStore.readStore, Standard() );
	noNamesPresent = (length(fragStore.readNameStore) == 0);
	for(TSize idCount = 0;readIt != readItEnd; ++readIt, ++idCount) {
		_streamWrite(target,"{RED\n");
		_streamWrite(target,"iid:");
		_streamPutInt(target, idCount + 1);
		_streamPut(target, '\n');
		if (!noNamesPresent) {
			_streamWrite(target,"eid:");
			_streamWrite(target, value(fragStore.readNameStore, idCount));
			_streamPut(target, '\n');
		}
		_streamWrite(target,"seq:\n");
		typedef typename Iterator<typename TFragmentStore::TReadSeq>::Type TSeqIter;
		typedef typename Value<typename TFragmentStore::TReadSeq>::Type TAlphabet;
		TSeqIter seqIt = begin(value(fragStore.readSeqStore, idCount));
		TSeqIter seqItEnd = end(value(fragStore.readSeqStore, idCount));
		for(TSize k = 0;seqIt!=seqItEnd;goNext(seqIt), ++k) {
			if ((k % 60 == 0) && (k != 0)) _streamPut(target, '\n');
			_streamPut(target, getValue(seqIt));
		}
		_streamWrite(target, "\n.\n");
		_streamWrite(target,"qlt:\n");
		seqIt = begin(value(fragStore.readSeqStore, idCount));
		for(TSize k = 0;seqIt!=seqItEnd;goNext(seqIt), ++k) {
			if ((k % 60 == 0) && (k != 0)) _streamPut(target, '\n');
			Ascii c = ' ';
			convertQuality(c, getQualityValue(value(seqIt)));
			_streamPut(target, c);
		}
		_streamWrite(target, "\n.\n");
		if (readIt->matePairId != TReadStoreElement::INVALID_ID) {
			_streamWrite(target,"frg:");
			_streamPutInt(target, readIt->matePairId + 1);
			_streamPut(target, '\n');
		}
		if ((value(clrRange, idCount)).i1 != (value(clrRange, idCount)).i2) {
			_streamWrite(target,"clr:");
			_streamPutInt(target, (value(clrRange, idCount)).i1);
			_streamPut(target, ',');
			_streamPutInt(target, (value(clrRange, idCount)).i2);
			_streamPut(target, '\n');
		}
		_streamWrite(target,"}\n");
	}

	// Sort aligned reads according to contigId
	sortAlignedReads(fragStore.alignedReadStore, SortContigId());

	// Write Contigs
	typedef typename Iterator<typename TFragmentStore::TContigStore, Standard>::Type TContigIter;
	TContigIter contigIt = begin(fragStore.contigStore, Standard() );
	TContigIter contigItEnd = end(fragStore.contigStore, Standard() );
	alignIt = begin(fragStore.alignedReadStore);
	alignItEnd = end(fragStore.alignedReadStore);
	noNamesPresent = (length(fragStore.contigNameStore) == 0);
	for(TSize idCount = 0;contigIt != contigItEnd; goNext(contigIt), ++idCount) {
		_streamWrite(target,"{CTG\n");
		_streamWrite(target,"iid:");
		_streamPutInt(target, idCount + 1);
		_streamPut(target, '\n');
		if (!noNamesPresent) {
			_streamWrite(target,"eid:");
			_streamWrite(target, value(fragStore.contigNameStore, idCount));
			_streamPut(target, '\n');
		}
		String<char> qlt;
		_streamWrite(target,"seq:\n");
		typedef typename Iterator<typename TFragmentStore::TContigSeq>::Type TContigIter;
		TContigIter seqContigIt = begin(contigIt->seq);
		TContigIter seqContigItEnd = end(contigIt->seq);
		typedef typename Iterator<String<typename TFragmentStore::TContigGapAnchor> >::Type TGapsIter;
		TGapsIter itGaps = begin(contigIt->gaps);
		TGapsIter itGapsEnd = end(contigIt->gaps);
		int diff = 0;
		char gapChar = gapValue<char>();
		typename TFragmentStore::TContigPos mySeqPos = 0;
		TSize k = 0;
		for(;itGaps != itGapsEnd; goNext(itGaps)) {
			while (mySeqPos < itGaps->seqPos) {
				if ((k % 60 == 0) && (k != 0)) _streamPut(target, '\n');
				++k;
				_streamPut(target, value(seqContigIt));
				Ascii c = ' ';
				convertQuality(c, getQualityValue(value(seqContigIt)));
				appendValue(qlt, c, Generous() );
				goNext(seqContigIt);++mySeqPos;
			}
			for(int i = 0; i < ((int) itGaps->gapPos - (int) itGaps->seqPos) - diff; ++i) {
				if ((k % 60 == 0) && (k != 0)) _streamPut(target, '\n');
				++k;
				_streamPut(target, gapChar);
				appendValue(qlt, '0', Generous() );
			}
			diff = (itGaps->gapPos - itGaps->seqPos);
		}
		for(;seqContigIt != seqContigItEnd; goNext(seqContigIt)) {
			if ((k % 60 == 0) && (k != 0)) _streamPut(target, '\n');
			++k;
			_streamPut(target, value(seqContigIt));
			Ascii c = ' ';
			convertQuality(c, getQualityValue(value(seqContigIt)));
			appendValue(qlt, c, Generous() );
		}
		_streamWrite(target, "\n.\n");
		_streamWrite(target,"qlt:\n");
		for(TSize k = 0;k<length(qlt); k+=60) {
			TSize endK = k + 60;
			if (endK > length(qlt)) endK = length(qlt);
			_streamWrite(target, infix(qlt, k, endK));
			_streamPut(target, '\n');
		}
		_streamWrite(target, ".\n");
		
		while ((alignIt != alignItEnd) && (idCount < alignIt->contigId)) goNext(alignIt);
		for(;(alignIt != alignItEnd) && (idCount == alignIt->contigId); goNext(alignIt)) {
			_streamWrite(target,"{TLE\n");
			_streamWrite(target,"src:");
			_streamPutInt(target, alignIt->readId + 1);
			_streamPut(target, '\n');
			typedef typename Iterator<String<typename TFragmentStore::TReadGapAnchor> >::Type TReadGapsIter;
			TReadGapsIter itGaps = begin(alignIt->gaps);
			TReadGapsIter itGapsEnd = end(alignIt->gaps);

			// Create the gaps string and the clear ranges
			typename TFragmentStore::TReadPos lenRead = length(value(fragStore.readSeqStore, alignIt->readId));
			TSize clr1 = 0;
			TSize clr2 = lenRead;
			// Create first clear range
			if ((itGaps != itGapsEnd) && (itGaps->gapPos == 0)) clr1 = itGaps->seqPos;
			int diff = clr1;
			String<unsigned int> gaps;
			for(;itGaps != itGapsEnd; goNext(itGaps)) {
				for(int i = 0; i< diff - ((int) itGaps->seqPos - (int) itGaps->gapPos); ++i) {
					appendValue(gaps, itGaps->seqPos - clr1, Generous() );
				}
				// Clipped sequence
				if (diff - ((int) itGaps->seqPos - (int) itGaps->gapPos) < 0) {
					clr2 = lenRead + diff - ((int) itGaps->seqPos - (int) itGaps->gapPos);
				}
				diff = ((int) itGaps->seqPos - (int) itGaps->gapPos);
			}
			if (alignIt->beginPos > alignIt->endPos) {
				clr1 = lenRead - clr1;
				clr2 = lenRead - clr2;
			}
			_streamWrite(target,"off:");
			if (alignIt->beginPos < alignIt->endPos) _streamPutInt(target, alignIt->beginPos);
			else _streamPutInt(target, alignIt->endPos);
			_streamPut(target, '\n');
			_streamWrite(target,"clr:");
			_streamPutInt(target, clr1);
			_streamPut(target, ',');
			_streamPutInt(target, clr2);
			_streamPut(target, '\n');
			if (length(gaps)) {
				_streamWrite(target,"gap:\n");
				for(TSize z = 0;z<length(gaps); ++z) {
					_streamPutInt(target, value(gaps, z));
					_streamPut(target, '\n');
				}
				_streamWrite(target, ".\n");
			}
			_streamWrite(target,"}\n");
		}
		_streamWrite(target,"}\n");
	}
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.writeContigs
..class:Class.FragmentStore
..summary:Write contigs from fragment store into file.
..cat:Fragment Store
..signature:writeContigs(file, store, tag)
..param.store:The fragment store.
...type:Class.FragmentStore
..param.file:A file/stream.
..param.tag:Specify format to write, e.g. Fasta.
..returns:A $bool$ which is $true$ on success.
..include:seqan/store.h
*/
template <typename TStream, typename TFSSpec, typename TFSConfig, typename TFormat>
bool writeContigs(TStream & file, FragmentStore<TFSSpec, TFSConfig> & store, TFormat const &)
{
//IOREV _doc_
	for (unsigned i = 0; i < length(store.contigNameStore); ++i)
		write(file, store.contigStore[i].seq, store.contigNameStore[i], TFormat());
	return true;
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.loadContigs
..class:Class.FragmentStore
..summary:Loads contigs into fragment store.
..cat:Fragment Store
..signature:loadContigs(store, fileName[, loadSeqs])
..signature:loadContigs(store, fileNameList[, loadSeqs])
..param.store:The fragment store.
...type:Class.FragmentStore
..param.fileName:A sequence file name.
...type:Shortcut.CharString
..param.fileNameList:A @Class.StringSet@ of sequence file names.
...type:Class.StringSet
..param.loadSeqs:If $true$, sequences are loaded immediately. 
If $false$, an empty contig with a reference to the file is created. Its sequence can be loaded on-demand by @Function.lockContig@ or @Function.loadContig@.
...default:$true$
...type:nolink:bool
..returns:A $bool$ which is $true$ on success.
..include:seqan/store.h
*/

template <typename TFSSpec, typename TFSConfig>
bool loadContigs(FragmentStore<TFSSpec, TFSConfig> &store, StringSet<CharString> const &fileNameList, bool loadSeqs)
{
//IOREV _nodoc_ although there is dddoc, there is no entry in html-doc
	typedef FragmentStore<TFSSpec, TFSConfig>			TFragmentStore;
	typedef typename TFragmentStore::TContigStore		TContigStore;
	typedef typename TFragmentStore::TContigFileStore	TContigFileStore;
	typedef typename Value<TContigStore>::Type			TContig;
	typedef typename Value<TContigFileStore>::Type		TContigFile;
	
	unsigned seqOfs = length(store.contigStore);
	for (unsigned filecount = 0; filecount < length(fileNameList); ++filecount)
	{
		MultiSeqFile multiSeqFile;
		if (!open(multiSeqFile.concat, toCString(fileNameList[filecount]), OPEN_RDONLY))
			return false;

		TContigFile contigFile;
		guessFormat(multiSeqFile.concat, contigFile.format);		// guess file format
		split(multiSeqFile, contigFile.format);						// divide into single sequences

		contigFile.fileName = fileNameList[filecount];
		contigFile.firstContigId = seqOfs;
		appendValue(store.contigFileStore, contigFile, Generous());

		unsigned seqCount = length(multiSeqFile);
		resize(store.contigStore, seqOfs + seqCount, Generous());
		resize(store.contigNameStore, seqOfs + seqCount, Generous());
		for (unsigned i = 0; i < seqCount; ++i)
		{
			store.contigStore[seqOfs + i].usage = 0;
			store.contigStore[seqOfs + i].fileBeginPos = beginPosition(multiSeqFile[i]);
			store.contigStore[seqOfs + i].fileEndPos = endPosition(multiSeqFile[i]);
			store.contigStore[seqOfs + i].fileId = length(store.contigFileStore) - 1;
			if (loadSeqs)
				assignSeq(store.contigStore[seqOfs + i].seq, multiSeqFile[i], contigFile.format);	// read Genome sequence
			else
            {
                typename TContig::TContigSeq emptySeq;
                swap(store.contigStore[seqOfs + i].seq, emptySeq);
            }
			assignCroppedSeqId(store.contigNameStore[seqOfs + i], multiSeqFile[i], contigFile.format);
		}
		seqOfs += seqCount;
	}
	reserve(store.contigStore, seqOfs, Exact());
	return seqOfs > 0;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFSSpec, typename TFSConfig>
bool loadContigs(FragmentStore<TFSSpec, TFSConfig> &store, CharString const &fileName, bool loadSeqs)
{
	StringSet<CharString> fileNames;
	appendValue(fileNames, fileName);
	return loadContigs(store, fileNames, loadSeqs);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFSSpec, typename TFSConfig, typename TFileNames>
bool loadContigs(FragmentStore<TFSSpec, TFSConfig> &store, TFileNames const &fileNames)
{
	return loadContigs(store, fileNames, true);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.loadContig
..class:Class.FragmentStore
..summary:Manually loads a contig sequence.
..cat:Fragment Store
..signature:loadContig(store, contigId)
..param.store:The fragment store.
...type:Class.FragmentStore
..param.contigId:Id of the contig which was created earlier by @Function.loadContigs@.
..returns:A $bool$ which is $true$ on success.
..include:seqan/store.h
*/

template <typename TSpec, typename TConfig, typename TId>
bool loadContig(FragmentStore<TSpec, TConfig> &store, TId _id)
{
	typedef FragmentStore<TSpec, TConfig>				TFragmentStore;
	typedef typename TFragmentStore::TContigStore		TContigStore;
	typedef typename TFragmentStore::TContigFileStore	TContigFileStore;
	typedef typename Value<TContigStore>::Type			TContig;
	typedef typename Value<TContigFileStore>::Type		TContigFile;

	if ((TId)length(store.contigStore) <= _id) return false;
	TContig &contig = store.contigStore[_id];

	if (contig.fileId >= length(store.contigFileStore)) return false;
	
	TContigFile &contigFile = store.contigFileStore[contig.fileId];
	String<char, MMap<> > fileString(toCString(contigFile.fileName), OPEN_RDONLY);
	assignSeq(contig.seq, infix(fileString, contig.fileBeginPos, contig.fileEndPos), contigFile.format);			// read Read sequence

	return true;
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.lockContig
..summary:Locks a contig sequence from being removed.
..cat:Fragment Store
..signature:lockContig(store, contigId)
..param.store:The fragment store.
...type:Class.FragmentStore
..param.contigId:Id of the contig which was created earlier by @Function.loadContigs@.
..returns:A $bool$ which is $true$ on success.
..remarks:This function increases the contig usage counter by 1 and ensures that the contig sequence is loaded.
..include:seqan/store.h
*/

template <typename TSpec, typename TConfig, typename TId>
bool lockContig(FragmentStore<TSpec, TConfig> &store, TId _id)
{
	typedef FragmentStore<TSpec, TConfig>				TFragmentStore;
	typedef typename TFragmentStore::TContigStore		TContigStore;
	typedef typename TFragmentStore::TContigFileStore	TContigFileStore;
	typedef typename Value<TContigStore>::Type			TContig;
	typedef typename Value<TContigFileStore>::Type		TContigFile;
	
	if ((TId)length(store.contigStore) <= _id) return false;
	TContig &contig = store.contigStore[_id];
	
	if (contig.usage++ > 0 || !empty(contig.seq)) return true;
	return loadContig(store, _id);
}

/**
.Function.unlockContig
..summary:Removes a previous contig lock.
..cat:Fragment Store
..signature:unlockContig(store, contigId)
..param.store:The fragment store.
...type:Class.FragmentStore
..param.contigId:Id of the contig which was created earlier by @Function.loadContigs@.
..returns:A $bool$ which is $true$ on success.
..remarks:This function decreases the contig usage counter by 1.
..include:seqan/store.h
*/

template <typename TSpec, typename TConfig, typename TId>
bool unlockContig(FragmentStore<TSpec, TConfig> &store, TId _id)
{
	if ((TId)length(store.contigStore) <= _id) return false;
	--store.contigStore[_id].usage;
	return true;
}

/**
.Function.unlockAndFreeContig
..summary:Removes a previous contig lock and clears sequence no further lock exist.
..cat:Fragment Store
..signature:unlockAndFreeContig(store, contigId)
..param.store:The fragment store.
...type:Class.FragmentStore
..param.contigId:Id of the contig which was created earlier by @Function.loadContigs@.
..returns:A $bool$ which is $true$ on success.
..remarks:This function decreases contig usage counter by 1 and clears contig sequence if counter is 0.
..include:seqan/store.h
*/

template <typename TSpec, typename TConfig, typename TId>
bool unlockAndFreeContig(FragmentStore<TSpec, TConfig> &store, TId _id)
{
	typedef FragmentStore<TSpec, TConfig>				TFragmentStore;
	typedef typename TFragmentStore::TContigStore		TContigStore;
	typedef typename Value<TContigStore>::Type			TContig;

	if ((TId)length(store.contigStore) <= _id) return false;
	TContig &contig = store.contigStore[_id];

	if (--contig.usage == 0 && contig.fileId < length(store.contigFileStore))
	{
        typename TContig::TContigSeq emptySeq;
        swap(contig.seq, emptySeq);
		return true;
	}
	return false;
}

/**
.Function.lockContigs
..summary:Locks all contig sequences from being removed. 
..cat:Fragment Store
..signature:lockContigs(store)
..param.store:The fragment store.
...type:Class.FragmentStore
..returns:A $bool$ which is $true$ on success.
..remarks:Calls @Function.lockContig@ for all contigs.
..include:seqan/store.h
*/

template <typename TSpec, typename TConfig>
bool lockContigs(FragmentStore<TSpec, TConfig> &store)
{
	bool result = true;
	for (unsigned _id = 0; _id < length(store.contigStore); ++_id)
		result &= lockContig(store, _id);
	return result;
}

/**
.Function.unlockContigs
..summary:Removes a previous lock for all contigs.
..cat:Fragment Store
..signature:unlockContigs(store)
..param.store:The fragment store.
...type:Class.FragmentStore
..returns:A $bool$ which is $true$ on success.
..remarks:Calls @Function.unlockContig@ for all contigs.
..include:seqan/store.h
*/

template <typename TSpec, typename TConfig>
bool unlockContigs(FragmentStore<TSpec, TConfig> &store)
{
	bool result = true;
	for (unsigned _id = 0; _id < length(store.contigStore); ++_id)
		result &= unlockContig(store, _id);
	return result;
}

/**
.Function.unlockAndFreeContigs
..summary:Removes a previous lock for all contigs and clears sequences without lock.
..cat:Fragment Store
..signature:unlockAndFreeContigs(store)
..param.store:The fragment store.
...type:Class.FragmentStore
..returns:A $bool$ which is $true$ on success.
..remarks:Calls @Function.unlockAndFreeContigs@ for all contigs.
..include:seqan/store.h
*/

template <typename TSpec, typename TConfig>
bool unlockAndFreeContigs(FragmentStore<TSpec, TConfig> &store)
{
	bool result = true;
	for (unsigned _id = 0; _id < length(store.contigStore); ++_id)
		result &= unlockAndFreeContig(store, _id);
	return result;
}


//////////////////////////////////////////////////////////////////////////////

/**
.Function.loadReads
..class:Class.FragmentStore
..summary:Loads reads into fragment store.
..remarks:When two file names are given, the files are expected to contain the same number of reads and reads with the same index are assumed to be mate pairs.
Mate pairs are stored internally in an "interleaved" mode, i.e. a read is read from each file before reading the next one.
..cat:Fragment Store
..signature:loadReads(store, fileName)
..signature:loadReads(store, fileNameL, fileNameR)
..param.store:The fragment store.
...type:Class.FragmentStore
..param.fileName:A sequence file name.
...type:Shortcut.CharString
..returns:A $bool$ which is $true$ on success.
..include:seqan/store.h
*/

template <typename TFSSpec, typename TFSConfig, typename TFileName>
bool loadReads(FragmentStore<TFSSpec, TFSConfig> &store, TFileName &fileName)
{
	typedef FragmentStore<TFSSpec, TFSConfig>			TFragmentStore;
	typedef typename TFragmentStore::TContigFileStore	TContigFileStore;
	typedef typename Value<TContigFileStore>::Type		TContigFile;

	MultiSeqFile multiSeqFile;
	if (!open(multiSeqFile.concat, toCString(fileName), OPEN_RDONLY))
		return false;

	// guess file format and split into sequence fractions
	AutoSeqFormat format;
	guessFormat(multiSeqFile.concat, format);
	split(multiSeqFile, format);

	// reserve space in fragment store
	unsigned seqOfs = length(store.readStore);
	unsigned seqCount = length(multiSeqFile);
	reserve(store.readStore, seqOfs + seqCount);
	reserve(store.readSeqStore, seqOfs + seqCount);
	reserve(store.readNameStore, seqOfs + seqCount);

	// read sequences
	String<Dna5Q> seq;
	CharString qual;
	CharString _id;

	for (unsigned i = 0; i < seqCount; ++i)
	{
		assignSeq(seq, multiSeqFile[i], format);    // read sequence
		assignQual(qual, multiSeqFile[i], format);  // read ascii quality values
		assignSeqId(_id, multiSeqFile[i], format);  // read sequence id

		// convert ascii to values from 0..62
		// store dna and quality together in Dna5Q
		// TODO: support different ASCII represenations of quality values
		assignQualities(seq, qual);
		appendRead(store, seq, _id);
	}
    return true;
}


template <typename TFSSpec, typename TFSConfig, typename TFileName>
bool loadReads(FragmentStore<TFSSpec, TFSConfig> & store, TFileName & fileNameL, TFileName & fileNameR)
{
	typedef FragmentStore<TFSSpec, TFSConfig>			TFragmentStore;
	typedef typename TFragmentStore::TContigFileStore	TContigFileStore;
	typedef typename Value<TContigFileStore>::Type		TContigFile;

	MultiSeqFile multiSeqFileL, multiSeqFileR;
	if (!open(multiSeqFileL.concat, toCString(fileNameL), OPEN_RDONLY))
		return false;
	if (!open(multiSeqFileR.concat, toCString(fileNameR), OPEN_RDONLY))
		return false;

	// Guess file format and split into sequence fractions
	AutoSeqFormat formatL, formatR;
	guessFormat(multiSeqFileL.concat, formatL);
	split(multiSeqFileL, formatL);
	guessFormat(multiSeqFileR.concat, formatR);
	split(multiSeqFileR, formatR);

    // Check that both files have the same number of reads
	SEQAN_ASSERT_EQ(length(multiSeqFileL), length(multiSeqFileR));

	// Reserve space in fragment store
	unsigned seqOfs = length(store.readStore);
	unsigned seqCountL = length(multiSeqFileL);
	unsigned seqCountR = length(multiSeqFileR);
	reserve(store.readStore, seqOfs + seqCountL + seqCountR);
	reserve(store.readSeqStore, seqOfs + seqCountL + seqCountR);
	reserve(store.readNameStore, seqOfs + seqCountL + seqCountR);

	// Read in sequences
	String<Dna5Q> seq[2];
	CharString qual[2];
	CharString _id[2];

	for (unsigned i = 0; i < seqCountL; ++i) {
		assignSeq(seq[0], multiSeqFileL[i], formatL);    // read sequence
		assignQual(qual[0], multiSeqFileL[i], formatL);  // read ascii quality values
		assignSeqId(_id[0], multiSeqFileL[i], formatL);  // read sequence id
		assignSeq(seq[1], multiSeqFileR[i], formatR);    // read sequence
		assignQual(qual[1], multiSeqFileR[i], formatR);  // read ascii quality values
		assignSeqId(_id[1], multiSeqFileR[i], formatR);  // read sequence id

		// convert ascii to values from 0..62
		// store dna and quality together in Dna5Q
		// TODO: support different ASCII represenations of quality values
		for (int j = 0; j < 2; ++j)
			assignQualities(seq[j], qual[j]);
		
		appendMatePair(store, seq[0], seq[1], _id[0], _id[1]);
	}
	return true;
}

}  // namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
