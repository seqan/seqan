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
  fragStoreRCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

 ============================================================================
  $Id$
 ==========================================================================*/

#ifndef SEQAN_HEADER_FUSION_H
#define SEQAN_HEADER_FUSION_H
//#define DEBUG_OVERLAP_MODULE

namespace seqan
{ 

//////////////////////////////////////////////////////////////////////////////
////// buildTupleCountStore_Fusion 
//////////////////////////////////////////////////////////////////////////////
template <typename TId>
struct TupleCountStoreElement_Fusion
{
	typedef String<TId>			TTuple;				
	typedef String<TTuple > 		TTupleList;
	typedef String<unsigned>		TTupleCounts;
	typedef String<double>			TTupleNorm;

	TTupleList		readConnections;
	TTupleCounts		readConnectionCounts;
	TTupleNorm		readConnectionNorm;
	
	TTupleList		matePairConnections; 
	TTupleCounts		matePairConnectionCounts;
	TTupleNorm		matePairConnectionNorm;
};


//////////////////////////////////////////////////////////////////////////////
template<typename TTupleCountStore, typename TTupleCountStore_Fusion, typename TSpec, typename TConfig, typename TReadAnnoStore>
inline void
buildTupleCountStore_Fusion(TTupleCountStore & tupleCountStore, 
		     TTupleCountStore_Fusion & tupleCountStore_Fusion, 
		     FragmentStore<TSpec, TConfig> &  fragStore, 
		     TReadAnnoStore & readAnnoStore, 
		     unsigned n, 
		     bool exact_nTuple)
{
	typedef typename FragmentStore<TSpec, TConfig>::TAnnotationStore 	TAnnotationStore;
	typedef typename FragmentStore<TSpec, TConfig>::TContigPos		TPos;
	typedef typename FragmentStore<TSpec, TConfig>::TReadStore 		TReadStore;
	typedef typename Value<TReadStore>::Type 				TReadStoreElement;
	typedef typename TReadStoreElement::TId					TReadId;
	typedef typename Value<TAnnotationStore>::Type 				TAnnotationStoreElement;
	typedef typename TAnnotationStoreElement::TId 				TId; 
	typedef typename Iterator<TReadAnnoStore>::Type 			TReadIter;
	typedef typename Value<TReadAnnoStore>::Type				TReadAnnoStoreElement;
	typedef typename TReadAnnoStoreElement::TAnnoIds			TAnnoIds;
	typedef typename Iterator<TAnnoIds>::Type 				TAnnoIdsIter;
	typedef typename Value<TAnnoIds>::Type					TIds;
	typedef typename Iterator<TIds>::Type					TIdsIter;
	
	static const TReadId INVALID_READ_ID = TReadStoreElement::INVALID_ID;
	static const TId INVALID_ANNO_ID = TAnnotationStoreElement::INVALID_ID;
	
	resize(tupleCountStore, length(fragStore.annotationStore)); 
	resize(tupleCountStore_Fusion, length(fragStore.annotationStore)); 
	
	bool validMate;
	TReadIter itRead = begin(readAnnoStore);
	TReadIter itReadEnd = end(readAnnoStore);
	TIdsIter itP;
	TIdsIter itPEnd;
	TIdsIter itP2;
	TIdsIter itP2End;
	TAnnoIds annoIds;
	TAnnoIds tupleSet;
	TReadId readId;
	TReadId matePairId;
	TReadId secReadId;
	TAnnoIds secTupleSet;
	TAnnoIds tempSecTupleSet;
	TId firstAnnoId1;
	TId firstAnnoId2;
	TAnnoIdsIter itTuple;
	TAnnoIdsIter itTupleEnd;
	TAnnoIdsIter itSecTuple;
	TAnnoIdsIter itSecTupleEnd;
	TAnnoIdsIter itAnnoIds;
	TAnnoIdsIter itAnnoIdsEnd;
	TIds matePairTuple;
	TPos beginPos1;	
	TPos endPos1;
	TPos beginPos2;
	TPos endPos2;
	unsigned pos;

	
	for ( ; itRead != itReadEnd; goNext(itRead))
	{
		if (!empty(getValue(itRead).parentIds) )
		{ 
			itP = begin(getValue(itRead).parentIds);
			itPEnd = end(getValue(itRead).parentIds);
			for ( ; itP != itPEnd; goNext(itP) )
			{
				validMate = false;
				// create list of all possible tuples for current read:
				annoIds = getValue(itRead).annoIds;
				clear(tupleSet);
				// create all Tuple of length n:
				if (exact_nTuple && n <= length(annoIds)) create_nTuple(tupleSet, fragStore, annoIds, getValue(itP), n);
				// create all max-Tuple (whole read) for current parentId:
				else if (!exact_nTuple && n == 0) create_nTuple(tupleSet, fragStore, annoIds, getValue(itP), length(annoIds));	
				// create all tuple >= n for current parentId:
				else if (!exact_nTuple) create_Tuple(tupleSet, fragStore, annoIds, getValue(itP), n);	
				if (!empty(tupleSet))
				{
					// create if necessary list of all possible tuples for second matepair-read:
					readId = position(itRead, readAnnoStore);
					matePairId = getValue(fragStore.readStore, readId).matePairId;
					clear(secTupleSet);
					if (matePairId != INVALID_READ_ID)
					{
						if(getValue(getValue(fragStore.matePairStore, matePairId).readId, 0) == readId)
							secReadId = getValue(getValue(fragStore.matePairStore, matePairId).readId, 1);
						else
							secReadId = getValue(getValue(fragStore.matePairStore, matePairId).readId, 0);
				
						if ( secReadId != INVALID_READ_ID )	
						{
							//if (!empty(getValue(readAnnoStore, secReadId).annoIds)) 
							if ( isElement_unsorted(getValue(itP), getValue(readAnnoStore, secReadId).parentIds) )	// p in parents of matepair? -> annoIds is not empty
							{
								validMate = true;
								annoIds = getValue(readAnnoStore, secReadId).annoIds;
								firstAnnoId1 = front(front(tupleSet));	// ids necessary to check positions in aligment  
								firstAnnoId2 = front(front(annoIds));	// can't be INVALID_ID, because parents was checked
							
								// check if current read-position is smaller than the position of the second read -> tuple are ordered by position
								if ( (getValue(fragStore.annotationStore, firstAnnoId1).beginPos <= 
									getValue(fragStore.annotationStore,firstAnnoId1).endPos && 
								      getValue(fragStore.annotationStore, firstAnnoId1).beginPos < 
								      	getValue(fragStore.annotationStore, firstAnnoId2).endPos) ||
								     (getValue(fragStore.annotationStore, firstAnnoId1).beginPos > 
								     	getValue(fragStore.annotationStore, firstAnnoId1).endPos && 
								      getValue(fragStore.annotationStore, firstAnnoId1).endPos < 
								      	getValue(fragStore.annotationStore, firstAnnoId2).beginPos)  ) 
								{	
									if (exact_nTuple && n <= length(annoIds)) create_nTuple(secTupleSet, fragStore, annoIds, getValue(itP), n);
									else if (!exact_nTuple && n == 0) create_nTuple(secTupleSet, fragStore, annoIds, getValue(itP), length(annoIds));		
									else if (!exact_nTuple) create_Tuple(secTupleSet, fragStore, annoIds, getValue(itP), n);
								}
							}
							else if (!empty(getValue(readAnnoStore, secReadId).parentIds)) 	// parent of matepair different -> possible transfusion
							{
								//validMate = true;
								annoIds = getValue(readAnnoStore, secReadId).annoIds;
								firstAnnoId1 = front(front(tupleSet));	// ids necessary to check positions in aligment  
								firstAnnoId2 = front(front(annoIds));	// can't be INVALID_ID, because parents was checked
							
								// check if current read-position is smaller than the position of the second read -> tuple are ordered by position
								if ( (getValue(fragStore.annotationStore, firstAnnoId1).beginPos <= 
									getValue(fragStore.annotationStore,firstAnnoId1).endPos && 
								      getValue(fragStore.annotationStore, firstAnnoId1).beginPos < 
								      	getValue(fragStore.annotationStore, firstAnnoId2).endPos) ||
								     (getValue(fragStore.annotationStore, firstAnnoId1).beginPos > 
								     	getValue(fragStore.annotationStore, firstAnnoId1).endPos && 
								      getValue(fragStore.annotationStore, firstAnnoId1).endPos < 
								      	getValue(fragStore.annotationStore, firstAnnoId2).beginPos)  ) 
								{	
									itP2 = begin(value(readAnnoStore, secReadId).parentIds);
									itP2End = end(value(readAnnoStore, secReadId).parentIds);
									clear(tempSecTupleSet);
									for ( ; itP2 != itP2End; goNext(itP2))
									{
										if (exact_nTuple && n <= length(annoIds)) create_nTuple(tempSecTupleSet, fragStore, annoIds, getValue(itP2), n);
										else if (!exact_nTuple && n == 0) create_nTuple(tempSecTupleSet, fragStore, annoIds, getValue(itP2), length(annoIds));		
										else if (!exact_nTuple) create_Tuple(tempSecTupleSet, fragStore, annoIds, getValue(itP2), n);

										itSecTuple = begin(tempSecTupleSet);
										itSecTupleEnd = end(tempSecTupleSet);
										for ( ; itSecTuple != itSecTupleEnd; goNext(itSecTuple))
										{										
											appendValue(secTupleSet, *itSecTuple);
										}
									}
								}
							}
						}
					}
					else validMate = true;
			
					// access to tupleCountStore for all tuple of current read:
					if (validMate)
					{
						
						itTuple = begin(tupleSet);
						itTupleEnd = end(tupleSet);
						for ( ; itTuple != itTupleEnd; goNext(itTuple))	
						{
							firstAnnoId1 = front(*itTuple);
							erase(value(itTuple), 0);			// first id is not stored; is know by position in tupleCountStore

							// readConnections:
							if (!empty(*itTuple))
							{
								if (searchValue(pos, *itTuple, getValue(tupleCountStore, firstAnnoId1).readConnections))
									++value(value(tupleCountStore, firstAnnoId1).readConnectionCounts, pos);
								else 
								{
									if (pos != endPosition(getValue(tupleCountStore, firstAnnoId1).readConnections) )
									{
										resizeSpace(value(tupleCountStore, firstAnnoId1).readConnections, 1, pos, pos, Generous());
										assignValue(value(tupleCountStore, firstAnnoId1).readConnections, pos, *itTuple);
										insertValue(value(tupleCountStore, firstAnnoId1).readConnectionCounts, pos, 1, Generous());
									}
									else
									{
										appendValue(value(tupleCountStore, firstAnnoId1).readConnections, *itTuple, Generous());
										appendValue(value(tupleCountStore, firstAnnoId1).readConnectionCounts, 1, Generous());
									}
								}
							}
							// matePairConnections: 
							if (!empty(secTupleSet))
							{
								itSecTuple = begin(secTupleSet);
								itSecTupleEnd = end(secTupleSet);
								for ( ; itSecTuple != itSecTupleEnd; goNext(itSecTuple) )
								{
									matePairTuple = *itTuple;
									// INVALID_ID: sign for connection by matepair (apart from that, there are no INVALID_IDs in the list)
									appendValue(matePairTuple, INVALID_ANNO_ID, Generous());				
									if (!empty(*itTuple) && back(*itTuple) == front(*itSecTuple))		// no id 2x allowed
									{	
										if (exact_nTuple == 0 && n == 0) erase(value(itSecTuple), 0);
										else continue;							// tupel would be created double or tupel wouldn't have the length n anymore
									}
									append(matePairTuple, *itSecTuple, Generous());
					
									if (empty(*itTuple))
									{
										beginPos1 = getValue(fragStore.annotationStore, firstAnnoId1).beginPos;
										endPos1 = getValue(fragStore.annotationStore, firstAnnoId1).endPos;
									}
									else
									{
										beginPos1 = getValue(fragStore.annotationStore, back(*itTuple)).beginPos;
										endPos1 = getValue(fragStore.annotationStore, back(*itTuple)).endPos;
									}
									// begin position of first annotation in tuple of second read
									beginPos2 = getValue(fragStore.annotationStore, front(*itSecTuple)).beginPos;
									endPos2 = getValue(fragStore.annotationStore, front(*itSecTuple)).endPos;
									if ( (beginPos1 <= endPos1 && endPos1 < beginPos2) ||			// no overlapping annotations allowed
									     (endPos1 < beginPos1 && beginPos1 < endPos2) )
									{
										if (searchValue(pos, matePairTuple, getValue(tupleCountStore, firstAnnoId1).matePairConnections))
											++value(value(tupleCountStore, firstAnnoId1).matePairConnectionCounts, pos);
										else 
										{
											if (pos != endPosition(getValue(tupleCountStore, firstAnnoId1).matePairConnections) )
											{
												resizeSpace(value(tupleCountStore, firstAnnoId1).matePairConnections, 1, pos, pos, Generous());
												assignValue(value(tupleCountStore, firstAnnoId1).matePairConnections, pos, matePairTuple);
												insertValue(value(tupleCountStore, firstAnnoId1).matePairConnectionCounts, pos, 1, Generous());
											}
											else
											{
												appendValue(value(tupleCountStore, firstAnnoId1).matePairConnections, matePairTuple, Generous());
												appendValue(value(tupleCountStore, firstAnnoId1).matePairConnectionCounts, 1, Generous());
											}
										}
									}
								}
							}
						}
					}
					else if (!empty(secTupleSet))	// Transfusion: matepairs with different parents
					{
						itTuple = begin(tupleSet);
						itTupleEnd = end(tupleSet);
						for ( ; itTuple != itTupleEnd; goNext(itTuple))	
						{
							firstAnnoId1 = front(*itTuple);
							erase(value(itTuple), 0);			// first id is not stored; is know by position in tupleCountStore
			
							itSecTuple = begin(secTupleSet);
							itSecTupleEnd = end(secTupleSet);
							for ( ; itSecTuple != itSecTupleEnd; goNext(itSecTuple) )
							{
								matePairTuple = *itTuple;
								// INVALID_ID: sign for connection by matepair (apart from that, there are no INVALID_IDs in the list)
								appendValue(matePairTuple, INVALID_ANNO_ID, Generous());				
								if (!empty(*itTuple) && back(*itTuple) == front(*itSecTuple))		// no id 2x allowed
								{	
									if (exact_nTuple == 0 && n == 0) erase(value(itSecTuple), 0);
									else continue;							// tupel would be created double or tupel wouldn't have the length n anymore
								}
								append(matePairTuple, *itSecTuple, Generous());
						
								if (empty(*itTuple))
								{
									beginPos1 = getValue(fragStore.annotationStore, firstAnnoId1).beginPos;
									endPos1 = getValue(fragStore.annotationStore, firstAnnoId1).endPos;
								}
								else
								{
									beginPos1 = getValue(fragStore.annotationStore, back(*itTuple)).beginPos;
									endPos1 = getValue(fragStore.annotationStore, back(*itTuple)).endPos;
								}
	
								// begin position of first annotation in tuple of second read
								beginPos2 = getValue(fragStore.annotationStore, front(*itSecTuple)).beginPos;
								endPos2 = getValue(fragStore.annotationStore, front(*itSecTuple)).endPos;
								if ( (beginPos1 <= endPos1 && endPos1 < beginPos2) ||			// no overlapping annotations allowed
								     (endPos1 < beginPos1 && beginPos1 < endPos2) )
								{
									if (searchValue(pos, matePairTuple, getValue(tupleCountStore_Fusion, firstAnnoId1).matePairConnections))
									{
										++value(value(tupleCountStore_Fusion, firstAnnoId1).matePairConnectionCounts, pos);
									}
									else 
									{
										if (pos != endPosition(getValue(tupleCountStore_Fusion, firstAnnoId1).matePairConnections) )
										{
											resizeSpace(value(tupleCountStore_Fusion, firstAnnoId1).matePairConnections, 1, pos, pos, Generous());
											assignValue(value(tupleCountStore_Fusion, firstAnnoId1).matePairConnections, pos, matePairTuple);
											insertValue(value(tupleCountStore_Fusion, firstAnnoId1).matePairConnectionCounts, pos, 1, Generous());
										}
										else
										{
											appendValue(value(tupleCountStore_Fusion, firstAnnoId1).matePairConnections, matePairTuple, Generous());
											appendValue(value(tupleCountStore_Fusion, firstAnnoId1).matePairConnectionCounts, 1, Generous());
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}


//////////////////////////////////////////////////////////////////////////////
//create tupleCountGFF_Fusion
//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TTupleCountStore_Fusion, typename TSpec, typename TConfig>
inline void
createTupleCountGFF_Fusion(TFile & tupleOutput_Fusion, TTupleCountStore_Fusion & tupleCountStore_Fusion, FragmentStore<TSpec, TConfig> & fragStore, unsigned thresholdCount, double thresholdRPKM)
{
	typedef typename FragmentStore<TSpec, TConfig>::TAnnotationStore 	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type 				TAnnotationStoreElement;
	typedef typename TAnnotationStoreElement::TId 				TId;
	
	typedef typename Iterator<TTupleCountStore_Fusion>::Type 		TCountStoreIter;
	typedef typename Value<TTupleCountStore_Fusion>::Type			TTupleCountStoreElement_Fusion;
	typedef typename TTupleCountStoreElement_Fusion::TTupleList		TTupleList;
	typedef typename TTupleCountStoreElement_Fusion::TTupleCounts		TTupleCounts;
	typedef typename TTupleCountStoreElement_Fusion::TTupleNorm		TTupleNorm;
	typedef typename Value<TTupleList>::Type				TTupel;
	typedef typename Iterator<TTupleList>::Type				TTupleListIter;
	typedef typename Iterator<TTupleCounts>::Type				TCountIter;
	typedef typename Iterator<TTupleNorm>::Type				TNormIter;
	typedef typename Iterator<TTupel>::Type					TTupelIter;

	static const TId INVALID_ID = TAnnotationStoreElement::INVALID_ID;
	
	if (!empty(tupleCountStore_Fusion))
	{
		TCountStoreIter itCountStore = begin(tupleCountStore_Fusion);
		TCountStoreIter itCountStoreEnd = end(tupleCountStore_Fusion);
		TAnnotationStoreElement currentElement;
		TTupleListIter itT;
		TTupleListIter itTEnd;
		TCountIter itC;
		TNormIter itN;
		TTupelIter itId;
		TTupelIter itIdEnd;
		for ( ; itCountStore != itCountStoreEnd; goNext(itCountStore))
		{
			currentElement = getValue(fragStore.annotationStore, position(itCountStore, tupleCountStore_Fusion));

			/*
			itT = begin(getValue(itCountStore).readConnections);
			itTEnd = end(getValue(itCountStore).readConnections);
			itC = begin(getValue(itCountStore).readConnectionCounts);
			itN = begin(getValue(itCountStore).readConnectionNorm);
			// read connections:
			for ( ; itT != itTEnd; goNext(itT), goNext(itC), goNext(itN))
			{
				
				if (getValue(itC) >= thresholdCount && getValue(itN) >= thresholdRPKM)
				{
					// contig-name
					streamPut(tupleOutput_Fusion, getValue(fragStore.contigNameStore, currentElement.contigId));
					itId = begin(getValue(itT));
					itIdEnd = end(getValue(itT));
					for ( ; itId != itIdEnd; goNext(itId))
					{
						if (getValue(itId) == INVALID_ID)
						{
							goNext(itId);
							streamPut(tupleOutput_Fusion, "~");
							streamPut(tupleOutput_Fusion, getValue(fragStore.contigNameStore, getValue(fragStore.annotationStore, getValue(itId)).contigId));
						}
					}
					streamPut(tupleOutput_Fusion, '\t');					

					// parent-names
					streamPut(tupleOutput_Fusion, getValue(fragStore.annotationNameStore, currentElement.parentId));
					itId = begin(getValue(itT));
					for ( ; itId != itIdEnd; goNext(itId))
					{
						if (getValue(itId) == INVALID_ID)
						{
							goNext(itId);
							streamPut(tupleOutput_Fusion, "~");
							streamPut(tupleOutput_Fusion, getValue(fragStore.annotationNameStore, getValue(fragStore.annotationStore, getValue(itId)).parentId));
						}
					}
					streamPut(tupleOutput_Fusion, '\t');
					
					// orientation
					if ( currentElement.beginPos <= currentElement.endPos )
					{
						streamPut(tupleOutput_Fusion, "+");
					}
					else
					{
						streamPut(tupleOutput_Fusion, "-");
					}
					itId = begin(getValue(itT));
					for ( ; itId != itIdEnd; goNext(itId))
					{
						if (getValue(itId) == INVALID_ID)
						{
							goNext(itId);
							streamPut(tupleOutput_Fusion, "~");
							if ( getValue(fragStore.annotationStore, getValue(itId)).beginPos <= getValue(fragStore.annotationStore, getValue(itId)).endPos )
							{
								streamPut(tupleOutput_Fusion, "+");
							}
							else
							{
								streamPut(tupleOutput_Fusion, "-");
							}
						}
					}
					streamPut(tupleOutput_Fusion, "\t");

					// first annotationId of tuple (store implicit)
					streamPut(tupleOutput_Fusion, getValue(fragStore.annotationNameStore, position(itCountStore, tupleCountStore_Fusion)));
					// other annotationIds
					itId = begin(getValue(itT));
					itIdEnd = end(getValue(itT));
					for ( ; itId != itIdEnd ; goNext(itId))
					{
						if (getValue(itId) != INVALID_ID)
						{
							streamPut(tupleOutput_Fusion, ":");
							streamPut(tupleOutput_Fusion, getValue(fragStore.annotationNameStore, getValue(itId)));
						}
						else
						{
							streamPut(tupleOutput_Fusion, "~");
							goNext(itId);
							streamPut(tupleOutput_Fusion, getValue(fragStore.annotationNameStore, getValue(itId)));
						}
					}
					streamPut(tupleOutput_Fusion, '\t');
					// tuple count
					streamPut(tupleOutput_Fusion, getValue(itC));
					streamPut(tupleOutput_Fusion, '\t');
					// normalized tuple count
					_streamPutDouble(tupleOutput_Fusion, getValue(itN));
					streamPut(tupleOutput_Fusion, '\n');
				}
			}
			*/
	
			// matePairConnections:
			itT = begin(getValue(itCountStore).matePairConnections);
			itTEnd = end(getValue(itCountStore).matePairConnections);
			itC = begin(getValue(itCountStore).matePairConnectionCounts);
			itN = begin(getValue(itCountStore).matePairConnectionNorm);
			for ( ; itT != itTEnd; goNext(itT), goNext(itC), goNext(itN))
			{if (getValue(itC) >= thresholdCount && getValue(itN) >= thresholdRPKM)
				{
					// contig-name
					tupleOutput_Fusion << getValue(fragStore.contigNameStore, currentElement.contigId);
					itId = begin(getValue(itT));
					itIdEnd = end(getValue(itT));
					for ( ; itId != itIdEnd; goNext(itId))
					{
						if (getValue(itId) == INVALID_ID)
						{
							goNext(itId);
							tupleOutput_Fusion << '^'
                                               << getValue(fragStore.contigNameStore,
                                                           getValue(fragStore.annotationStore,
                                                                    getValue(itId)).contigId);
						}
					}
					tupleOutput_Fusion << '\t';
				
					// parent-name
					if (currentElement.parentId == INVALID_ID )
						tupleOutput_Fusion << "NO_PARENT";
					else
						tupleOutput_Fusion << getValue(fragStore.annotationNameStore, currentElement.parentId);
					itId = begin(getValue(itT));
					for ( ; itId != itIdEnd; goNext(itId))
					{
						if (getValue(itId) == INVALID_ID)
						{
							goNext(itId);
							if (getValue(fragStore.annotationStore, getValue(itId)).parentId == INVALID_ID)
								tupleOutput_Fusion << "^NO_PARENT";
							else
								tupleOutput_Fusion << "^"
                                                   << getValue(fragStore.annotationNameStore, getValue(fragStore.annotationStore, getValue(itId)).parentId);
						}
					}
					tupleOutput_Fusion << '\t';

					// orientation
					if ( currentElement.beginPos <= currentElement.endPos )
						tupleOutput_Fusion << "+";
					else
						tupleOutput_Fusion << "-";
					itId = begin(getValue(itT));
					for ( ; itId != itIdEnd; goNext(itId))
					{
						if (getValue(itId) == INVALID_ID)
						{
							goNext(itId);
							tupleOutput_Fusion << "^";
							if ( getValue(fragStore.annotationStore, getValue(itId)).beginPos <= getValue(fragStore.annotationStore, getValue(itId)).endPos )
								tupleOutput_Fusion << "+";
							else
								tupleOutput_Fusion << "-";
						}
					}
					tupleOutput_Fusion << "\t";
				
					// first annotationId of tuple
					tupleOutput_Fusion << getValue(fragStore.annotationNameStore, position(itCountStore, tupleCountStore_Fusion));
			
					// other annotationIds of first read
					itId = begin(getValue(itT));
					itIdEnd = end(getValue(itT));
					for ( ; itId != itIdEnd && getValue(itId) != INVALID_ID; goNext(itId))
					{
						tupleOutput_Fusion << ":" << getValue(fragStore.annotationNameStore, getValue(itId));
					}
					goNext(itId);
					tupleOutput_Fusion << "^";
		
					// annotationIds of second read
					tupleOutput_Fusion << getValue(fragStore.annotationNameStore, getValue(itId));
					goNext(itId);
					for ( ; itId != itIdEnd; goNext(itId))
						tupleOutput_Fusion << ":" << getValue(fragStore.annotationNameStore, getValue(itId));
					tupleOutput_Fusion << '\t';
					
					// tuple count
					tupleOutput_Fusion << *itC << '\t';
					// normalized tuple count
					tupleOutput_Fusion << formattedNumber("%f", *itN) << '\n';
				}
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////////
/// get normalized values for tuple of transfusion genes
//////////////////////////////////////////////////////////////////////////////
template<typename TTupleCountStore_Fusion, typename TSpec, typename TConfig>
inline void
normalizeTupleCounts_Fusion(TTupleCountStore_Fusion &tupleCountStore_Fusion, FragmentStore<TSpec, TConfig> &fragStore)
{
	typedef typename FragmentStore<TSpec, TConfig>::TAnnotationStore 	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type				TAnnotationStoreElement;
	typedef typename TAnnotationStoreElement::TId			TId;
	typedef typename Value<TTupleCountStore_Fusion>::Type			TTupleCountStoreElement_Fusion;
	typedef typename TTupleCountStoreElement_Fusion::TTupleList		TTupleList;
	typedef typename TTupleCountStoreElement_Fusion::TTupleCounts		TTupleCounts;
	typedef typename TTupleCountStoreElement_Fusion::TTupleNorm		TTupleNorm;
	typedef typename TTupleCountStoreElement_Fusion::TTuple			TTuple;
	typedef typename Iterator<TTupleCountStore_Fusion>::Type 		TStoreIter;
	typedef typename Iterator<TTupleList>::Type				TTupleListIter;
	typedef typename Iterator<TTupleCounts>::Type				TCountIter;
	typedef typename Iterator<TTupleNorm>::Type				TNormIter;
	typedef typename Iterator<TTuple>::Type					TTupleIter;
	typedef typename FragmentStore<TSpec, TConfig>::TContigPos		TPos;
	typedef typename Size<TPos>::Type					TSize;
	typedef typename FragmentStore<TSpec, TConfig>::TReadStore		TReadStore;
	typedef typename Size<TReadStore>::Type					TReadStoreSize; 
	
	static const TId INVALID_ID = TAnnotationStoreElement::INVALID_ID;
	
	TReadStoreSize readNo = length(fragStore.readStore) - length(fragStore.matePairStore);

	if (!empty(tupleCountStore_Fusion))
	{
		TStoreIter itS = begin(tupleCountStore_Fusion);
		TStoreIter itSEnd = end(tupleCountStore_Fusion);
		TTupleListIter itT;
		TTupleListIter itTEnd;
		TCountIter itC;
		TNormIter itN;
		TSize tupleLength;
		TTupleIter itId;
		TTupleIter itIdEnd;
		for ( ; itS != itSEnd; goNext(itS))
		{
			/*
			// readConnections:
			resize(value(itS).readConnectionNorm, length(getValue(itS).readConnections));
			if (!empty(getValue(itS).readConnections))
			{
				itT = begin(getValue(itS).readConnections);
				itTEnd = end(getValue(itS).readConnections);
				itC = begin(getValue(itS).readConnectionCounts);
				itN = begin(getValue(itS).readConnectionNorm);
				for ( ; itT != itTEnd; goNext(itT), goNext(itC), goNext(itN))
				{
					tupleLength = 0;
					itId = begin(getValue(itT));
					itIdEnd = end(getValue(itT));
					for ( ; itId != itIdEnd; goNext(itId))
					{	
						if (getValue(itId) != INVALID_ID)
						{
							if (getValue(fragStore.annotationStore, getValue(itId)).beginPos <= getValue(fragStore.annotationStore, getValue(itId)).endPos)
								tupleLength += getValue(fragStore.annotationStore, getValue(itId)).endPos - getValue(fragStore.annotationStore, getValue(itId)).beginPos;
							else
								tupleLength += getValue(fragStore.annotationStore, getValue(itId)).beginPos - getValue(fragStore.annotationStore, getValue(itId)).endPos;
						}
					}
					value(itN) = ((double)1000000000 * (double)getValue(itC)) / ((double)readNo * (double)tupleLength);
				}
			}
			*/
			// matePairConnections:
			resize(value(itS).matePairConnectionNorm, length(getValue(itS).matePairConnections));
			if (!empty(getValue(itS).matePairConnections))
			{
				itT = begin(getValue(itS).matePairConnections);
				itTEnd = end(getValue(itS).matePairConnections);
				itC = begin(getValue(itS).matePairConnectionCounts);
				itN = begin(getValue(itS).matePairConnectionNorm);
				for ( ; itT != itTEnd; goNext(itT), goNext(itC), goNext(itN))
				{
					tupleLength = 0;
					itId = begin(getValue(itT));
					itIdEnd = end(getValue(itT));
					for ( ; itId != itIdEnd; goNext(itId))
					{
						if (getValue(itId) != INVALID_ID)
						{
							if (getValue(fragStore.annotationStore, getValue(itId)).beginPos <= getValue(fragStore.annotationStore, getValue(itId)).endPos)
								tupleLength += getValue(fragStore.annotationStore, getValue(itId)).endPos - getValue(fragStore.annotationStore, getValue(itId)).beginPos;
							else
								tupleLength += getValue(fragStore.annotationStore, getValue(itId)).beginPos - getValue(fragStore.annotationStore, getValue(itId)).endPos;
						}
					}
					value(itN) = ((double)1000000000 * (double)getValue(itC)) / ((double)readNo * (double)tupleLength);
				}
			}
		}
	}
}


//////////////////////////////////////////////////////////////////////////////
////// Overlap Module
//////////////////////////////////////////////////////////////////////////////
template<typename TReadAnnoStore, typename TAnnoCountStore, typename TTupleCountStore, typename TTupleCountStore_Fusion, typename TSpec, typename TConfig>
inline void
getResults_Fusion(TReadAnnoStore & readAnnoStore,
	   TAnnoCountStore & annoCountStore,
	   TTupleCountStore & tupleCountStore, 
	   TTupleCountStore_Fusion & tupleCountStore_Fusion, 
	   FragmentStore<TSpec, TConfig> & fragStore, 
	   unsigned tupelSize,
	   bool exact_nTuple,
	   unsigned offsetInterval,
	   unsigned thresholdGaps,
	   bool unknownO)
{
	typedef typename FragmentStore<TSpec, TConfig>::TAnnotationStore 	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type 				TAnnotationStoreElement;
	typedef typename TAnnotationStoreElement::TId 				TId;
	typedef typename FragmentStore<TSpec, TConfig>::TIntervalTreeStore 	TIntervalTreeStore;
	typedef typename Iterator<TIntervalTreeStore>::Type			TIntervalTree;
	typedef typename Value<TReadAnnoStore>::Type 				TReadAnnoStoreElement;
	typedef typename TReadAnnoStoreElement::TAnnoIds 			TAnnoIds;
	
	typedef typename FragmentStore<TSpec, TConfig>::TAlignedReadStore	TAlignedReadStore;
	typedef typename Position<TAlignedReadStore>::Type 			TAlignPos;
	typedef 	 String<AlignIntervalsStoreElement<> > 			TAlignIntervalsStore;
	typedef typename Iterator<TAlignIntervalsStore>::Type 			TAlignIntervalsStoreIter;
	
	resize(readAnnoStore, length(fragStore.readStore));
	
	TIntervalTree intervalTree;
	
	// extract intervals from alignedReadStore and store them in AlignIntervalsStore:
	TAlignIntervalsStore alignIntervalsStore;
	buildAlignIntervalsStore(alignIntervalsStore, fragStore, thresholdGaps);

	if (!empty(alignIntervalsStore))
	{
		TAlignPos alignPos;
		TId contigId;
		TId readId;
		TAnnoIds  ids;
	
		TAlignIntervalsStoreIter it = begin(alignIntervalsStore);
		TAlignIntervalsStoreIter itEnd = end(alignIntervalsStore);
		// for each item in alignIntervalsStore:
		for ( ; it != itEnd; goNext(it))
		{
			// get ids from alignedReadStore (same position as in alignIntervalsStore):
			alignPos = position(it, alignIntervalsStore);
			contigId = getValue(fragStore.alignedReadStore, alignPos).contigId;
			readId = getValue(fragStore.alignedReadStore, alignPos).readId;
			// get respective intervalTree
			if (unknownO ||  getValue(fragStore.alignedReadStore, alignPos).beginPos <= getValue(fragStore.alignedReadStore, alignPos).endPos)
				intervalTree = begin(fragStore.intervalTreeStore_F, Standard()) + contigId; 	//getValue(fragStore.intervalTreeStore_F, contigId);
			else 
				intervalTree = begin(fragStore.intervalTreeStore_R, Standard()) + contigId;        //getValue(fragStore.intervalTreeStore_R, contigId);
			
			// get annotationStore-Ids for these intervals:
			clear(ids);
			if ((*intervalTree).interval_counter != 0)
				getIdsForRead(ids, fragStore, *intervalTree, getValue(it).intervals, offsetInterval);
			// assign Ids from mapped annotations to readAnnoStore:
			value(readAnnoStore, readId).contigId = contigId;
			assignToReadAnnoStore(readAnnoStore, fragStore, readId, ids);
		}
	}
	buildAnnoCountStore(annoCountStore, fragStore, readAnnoStore);
	buildTupleCountStore_Fusion(tupleCountStore, tupleCountStore_Fusion,  fragStore, readAnnoStore, tupelSize, exact_nTuple);
}



//////////////////////////////////////////////////////////////////////////////

}// namespace seqan

#endif //#ifndef SEQAN_HEADER_...
