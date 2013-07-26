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

#ifndef SEQAN_HEADER_CREATE_GFF_H
#define SEQAN_HEADER_CREATE_GFF_H

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
//create readCountGFF
//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TReadAnnoStore, typename TSpec, typename TConfig>
inline void
createReadCountGFF(TFile & readOutput, TReadAnnoStore & readAnnoStore, FragmentStore<TSpec, TConfig> & fragStore)
{	
	typedef typename Iterator<TReadAnnoStore>::Type 			TCountIter;
	typedef typename Value<TReadAnnoStore>::Type				TReadAnnoStoreElement;
	typedef typename TReadAnnoStoreElement::TAnnoIds			TAnnoIds;
	typedef typename Iterator<TAnnoIds>::Type				TAnnoIdsIter;
	typedef typename Value<TAnnoIds>::Type					TIds;
	typedef typename Iterator<TIds>::Type					TIdsIter;	
					
	typedef typename FragmentStore<TSpec, TConfig>::TAnnotationStore 	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type 				TAnnotationStoreElement;
	typedef typename TAnnotationStoreElement::TId 				TId;
	
	static const TId INVALID_ID = TAnnotationStoreElement::INVALID_ID;
	
	TCountIter itCountStore = begin(readAnnoStore);
	TCountIter itCountStoreEnd = end(readAnnoStore);
	TAnnoIdsIter itAnnoIds;
	TAnnoIdsIter itAnnoIdsEnd;
	TId firstId;
	TIds allParentIds;
	TIdsIter itId;
	TIdsIter itIdEnd;
	TIdsIter itP;
	TIdsIter itPEnd;
	bool help;
	bool invalid;
	for ( ; itCountStore != itCountStoreEnd; goNext(itCountStore))
	{
		// read-name:
		streamPut(readOutput, getValue(fragStore.readNameStore, position(itCountStore, readAnnoStore)) );
		streamPut(readOutput, '\t');
		
		if (empty(getValue(itCountStore).annoIds) )
		{
			streamPut(readOutput, ".\t.\t.\t.\t.");
		}
		else
		{
			// contig-name:
			streamPut(readOutput, getValue(fragStore.contigNameStore, getValue(itCountStore).contigId) );
			streamPut(readOutput, '\t');
			
			itAnnoIds = begin(getValue(itCountStore).annoIds);
			itAnnoIdsEnd = end(getValue(itCountStore).annoIds);
			while (itAnnoIds != itAnnoIdsEnd && front(getValue(itAnnoIds)) == INVALID_ID)
			{
				goNext(itAnnoIds);
			}
	
			if (itAnnoIds != itAnnoIdsEnd) // not only INVALID_IDS
			{
				firstId = front(getValue(itAnnoIds));
				
				// orientation:
				if (getValue(fragStore.annotationStore, firstId).beginPos <= getValue(fragStore.annotationStore, firstId).endPos)
				{
					streamPut(readOutput, "+\t");
				}
				else streamPut(readOutput, "-\t");
				
				// Annotation-Ids:
				
				allParentIds =  getValue(itCountStore).parentIds;
				// output for first parentId if possible; get other parentIds, in which the read doesn't map (entirely with all of his intervals)
				itAnnoIds = begin(getValue(itCountStore).annoIds);
				itAnnoIdsEnd = end(getValue(itCountStore).annoIds);
				for ( ; itAnnoIds != itAnnoIdsEnd; goNext(itAnnoIds))
				{
					help = false;
					itId = begin(*itAnnoIds);
					itIdEnd = end(*itAnnoIds);
					for ( ; itId != itIdEnd; goNext(itId))
					{
						if (!empty(getValue(itCountStore).parentIds) && getValue(itId) != INVALID_ID && 	// if current parentId == first id of parentIds (of read entry in readAnnoStore)
						    getValue(fragStore.annotationStore, getValue(itId)).parentId == front(allParentIds))
						{
							if (help) 	// not the first annotation for this read-interval -> ";" sign for overlapping annotations
							{
								streamPut(readOutput, ";");
							}
							streamPut(readOutput, getValue(fragStore.annotationNameStore, getValue(itId)) );
							help = true;
						}
						else if (getValue(itId) != INVALID_ID && getValue(fragStore.annotationStore, getValue(itId)).parentId != INVALID_ID &&	// get other parentIds 
							!isElement_unsorted(getValue(fragStore.annotationStore, getValue(itId)).parentId, allParentIds) ) //?
						{
							appendValue(allParentIds, getValue(fragStore.annotationStore, getValue(itId)).parentId, Generous() );
						}
					}
					if ( !empty(getValue(itCountStore).parentIds) && position(itAnnoIds, getValue(itCountStore).annoIds) != endPosition(getValue(itCountStore).annoIds) - 1)
					{
						streamPut(readOutput, ":");
					}
				}
				if (!empty(getValue(itCountStore).parentIds))
				{
					streamPut(readOutput, '\t');
					streamPut(readOutput, getValue(fragStore.annotationNameStore, front(allParentIds)) );
					streamPut(readOutput, "\t");
				}
				// outputs for all other parentIds
				itP = begin(allParentIds);
				itPEnd = end(allParentIds);
				if (!empty(getValue(itCountStore).parentIds)) goNext(itP);
				for ( ; itP != itPEnd; goNext(itP))
				{
					itAnnoIds = begin(getValue(itCountStore).annoIds);
					itAnnoIdsEnd = end(getValue(itCountStore).annoIds);
					for ( ; itAnnoIds != itAnnoIdsEnd; goNext(itAnnoIds))
					{
						invalid = true;				// if no annotation for the current parent in  interval -> UNKOWN_REGION
						itId = begin(*itAnnoIds);	
                        itIdEnd = end(*itAnnoIds);
						for ( ; itId != itIdEnd; goNext(itId))
						{
							if (getValue(itId) != INVALID_ID && getValue(fragStore.annotationStore, getValue(itId)).parentId == getValue(itP))
							{
								if (!invalid)	// not the first annotation for this interval -> ";" sign for overlapping annotations
								{
									streamPut(readOutput, ";");
								}
								streamPut(readOutput, getValue(fragStore.annotationNameStore, getValue(itId)) );
								invalid = false;
							}
						}
						if (invalid)
						{
							streamPut(readOutput, "UNKNOWN_REGION");
						}
						if (position(itAnnoIds, getValue(itCountStore).annoIds) != endPosition(getValue(itCountStore).annoIds) - 1)
						{
							streamPut(readOutput, ":");
						}
					}
					streamPut(readOutput, '\t');
					if (getValue(itP) != INVALID_ID) streamPut(readOutput, getValue(fragStore.annotationNameStore, getValue(itP)) );
					else streamPut(readOutput, "NO_PARENT" );
					streamPut(readOutput, "\t");
				}
			}
			else  // only INVALID_IDS
			{	
				streamPut(readOutput, ".\t");
						
				// invalid_ids for each interval
				itAnnoIds = begin(getValue(itCountStore).annoIds);
				itAnnoIdsEnd = end(getValue(itCountStore).annoIds);
				for ( ; itAnnoIds != itAnnoIdsEnd; goNext(itAnnoIds))
				{
					streamPut(readOutput, "UNKNOWN_REGION");
					if (position(itAnnoIds, getValue(itCountStore).annoIds) != endPosition(getValue(itCountStore).annoIds) - 1)
					{
						streamPut(readOutput, ":");
					}
					
				}
				streamPut(readOutput, '\t');
				streamPut(readOutput, "UNKNOWN_REGION");
			}
		}
		streamPut(readOutput, '\n');	
	}	
}



//////////////////////////////////////////////////////////////////////////////
//create AnnoCountGFF
//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TAnnoCountStore, typename TAnnoNormStore, typename TSpec, typename TConfig, typename TMap>
inline void
createAnnoCountGFF(TFile & annoOutput, TAnnoCountStore & annoCountStore, TAnnoNormStore &annoNormStore, FragmentStore<TSpec, TConfig> & fragStore, TMap &mapO)
{
	typedef typename FragmentStore<TSpec, TConfig>::TContigPos 		TContigPos;
	typedef typename FragmentStore<TSpec, TConfig>::TAnnotationStore 	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type 				TAnnotationStoreElement;
	typedef typename TAnnotationStoreElement::TId 				TId;
	typedef typename Iterator<TAnnoCountStore>::Type 			TCountIter;
	typedef typename Iterator<TAnnoNormStore>::Type				TNormIter;
	typedef typename Iterator<TAnnotationStore>::Type 			TAnnoIter;
	
	static const TId INVALID_ID = TAnnotationStoreElement::INVALID_ID;
	static const TContigPos INVALID_POS = TAnnotationStoreElement::INVALID_POS;

	TCountIter itCount = begin(annoCountStore);
	TCountIter itCountEnd = end(annoCountStore);
	TAnnoIter itAnno = begin(fragStore.annotationStore);
	TNormIter itNorm = begin(annoNormStore);
	
	for ( ; itCount != itCountEnd; goNext(itCount), goNext(itAnno), goNext(itNorm))
	{
        if (getValue(itAnno).typeId != INVALID_ID)
	        if (fragStore.annotationTypeStore[getValue(itAnno).typeId] == "<root>") continue;  
		// contig-name
		if (getValue(itAnno).contigId == INVALID_ID )
		{
			streamPut(annoOutput, "INVALID_ID");
			streamPut(annoOutput, '\t');
		}
		else
		{
			streamPut(annoOutput, getValue(fragStore.contigNameStore, getValue(itAnno).contigId));
			streamPut(annoOutput, '\t');
		}
		streamPut(annoOutput, "Annotation_Count\tregion\t");
		// startposition endposition orientation . 
		if (getValue(itAnno).beginPos == INVALID_POS)
		{
			streamPut(annoOutput, ".\t.\t");
			streamPut(annoOutput, getValue(itCount));
			if (getValue(itAnno).parentId == INVALID_ID || (fragStore.annotationStore[getValue(itAnno).parentId].typeId != INVALID_ID && fragStore.annotationTypeStore[fragStore.annotationStore[getValue(itAnno).parentId].typeId] == "<root>"))
			{
				if (mapValue(mapO, position(itAnno, fragStore.annotationStore)) == 0)
				{
					streamPut(annoOutput, "\t+\t.\t");
				}
				else
				{
					streamPut(annoOutput, "\t-\t.\t");
				}
			}
			else	streamPut(annoOutput, "\t.\t.\t");
		}
		else
		{
			if (getValue(itAnno).beginPos <= getValue(itAnno).endPos)
			{
				streamPut(annoOutput, getValue(itAnno).beginPos+1);
				streamPut(annoOutput, '\t');
				streamPut(annoOutput, getValue(itAnno).endPos);
				streamPut(annoOutput, '\t');
				streamPut(annoOutput, getValue(itCount));
				streamPut(annoOutput, "\t+\t.\t");
			}
			else
			{
				streamPut(annoOutput, getValue(itAnno).endPos+1);
				streamPut(annoOutput, '\t');
				streamPut(annoOutput, getValue(itAnno).beginPos);
				streamPut(annoOutput, '\t');
				streamPut(annoOutput, getValue(itCount));
				streamPut(annoOutput, "\t-\t.\t");
			}
		}
		// annotation-name (parent annotation-name)
		if (getValue(itAnno).parentId == INVALID_ID || (fragStore.annotationStore[getValue(itAnno).parentId].typeId != INVALID_ID && fragStore.annotationTypeStore[fragStore.annotationStore[getValue(itAnno).parentId].typeId] == "<root>"))
		{
			streamPut(annoOutput, "ID=");
			streamPut(annoOutput, getValue(fragStore.annotationNameStore, position(itAnno, fragStore.annotationStore)) );
			streamPut(annoOutput, ';');
		}
		else
		{
			streamPut(annoOutput, "ID=");
			streamPut(annoOutput, getValue(fragStore.annotationNameStore, position(itAnno, fragStore.annotationStore)));
			streamPut(annoOutput, ";ParentID=");
			streamPut(annoOutput, getValue(fragStore.annotationNameStore, getValue(itAnno).parentId));
			streamPut(annoOutput, ';');
		}
		_streamPutDouble(annoOutput, getValue(itNorm));
		streamPut(annoOutput, ";\n");
	}
}


//////////////////////////////////////////////////////////////////////////////
//create tupleCountGFF
//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TTupleCountStore, typename TSpec, typename TConfig>
inline void
createTupleCountGFF(TFile & tupleOutput, TTupleCountStore & tupleCountStore, FragmentStore<TSpec, TConfig> & fragStore, unsigned thresholdCount, double thresholdRPKM)
{
	typedef typename FragmentStore<TSpec, TConfig>::TAnnotationStore 	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type 				TAnnotationStoreElement;
	typedef typename TAnnotationStoreElement::TId 				TId;
	
	typedef typename Iterator<TTupleCountStore>::Type 			TCountStoreIter;
	typedef typename Value<TTupleCountStore>::Type				TTupleCountStoreElement;
	typedef typename TTupleCountStoreElement::TTupleList			TTupleList;
	typedef typename TTupleCountStoreElement::TTupleCounts			TTupleCounts;
	typedef typename TTupleCountStoreElement::TTupleNorm			TTupleNorm;
	typedef typename Value<TTupleList>::Type				TTupel;
	typedef typename Iterator<TTupleList>::Type				TTupleListIter;
	typedef typename Iterator<TTupleCounts>::Type				TCountIter;
	typedef typename Iterator<TTupleNorm>::Type				TNormIter;
	typedef typename Iterator<TTupel>::Type					TTupelIter;
	
	static const TId INVALID_ID = TAnnotationStoreElement::INVALID_ID;
	
	if (!empty(tupleCountStore))
	{
		TCountStoreIter itCountStore = begin(tupleCountStore);
		TCountStoreIter itCountStoreEnd = end(tupleCountStore);
		TAnnotationStoreElement currentElement;
		TTupleListIter itT;
		TTupleListIter itTEnd;
		TCountIter itC;
		TNormIter itN;
		TTupelIter itId;
		TTupelIter itIdEnd;
		for ( ; itCountStore != itCountStoreEnd; goNext(itCountStore))
		{
			currentElement = getValue(fragStore.annotationStore, position(itCountStore, tupleCountStore));
	
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
					streamPut(tupleOutput, getValue(fragStore.contigNameStore, currentElement.contigId));
					streamPut(tupleOutput, '\t');
					// parent-name
					if (currentElement.parentId == INVALID_ID )
					{
						streamPut(tupleOutput, "NO_PARENT\t");
					}
					else
					{
						streamPut(tupleOutput, getValue(fragStore.annotationNameStore, currentElement.parentId));
						streamPut(tupleOutput, '\t');
					}
					// orientation
					if ( currentElement.beginPos <= currentElement.endPos )
					{
						streamPut(tupleOutput, "+\t");
					}
					else
					{
						streamPut(tupleOutput, "-\t");
					}
					// first annotationId of tuple (store implicit)
					streamPut(tupleOutput, getValue(fragStore.annotationNameStore, position(itCountStore, tupleCountStore)));
					// other annotationIds
					itId = begin(getValue(itT));
					itIdEnd = end(getValue(itT));
					for ( ; itId != itIdEnd; goNext(itId))
					{
						streamPut(tupleOutput, ":");
						streamPut(tupleOutput, getValue(fragStore.annotationNameStore, getValue(itId)));
					}
					streamPut(tupleOutput, '\t');
					// tuple count
					streamPut(tupleOutput, getValue(itC));
					streamPut(tupleOutput, '\t');
					// normalized tuple count
					_streamPutDouble(tupleOutput, getValue(itN));
					streamPut(tupleOutput, '\n');
				}
			}
			//matepair connections:
			itT = begin(getValue(itCountStore).matePairConnections);
			itTEnd = end(getValue(itCountStore).matePairConnections);
			itC = begin(getValue(itCountStore).matePairConnectionCounts);
			itN = begin(getValue(itCountStore).matePairConnectionNorm);
			for ( ; itT != itTEnd; goNext(itT), goNext(itC), goNext(itN))
			{
				if (getValue(itC) >= thresholdCount && getValue(itN) >= thresholdRPKM)
				{
					// contig-name
					streamPut(tupleOutput, getValue(fragStore.contigNameStore, currentElement.contigId));
					streamPut(tupleOutput, '\t');
					// parent-name
					if (currentElement.parentId == INVALID_ID )
					{
						streamPut(tupleOutput, "NO_PARENT\t");
					}
					else
					{
						streamPut(tupleOutput, getValue(fragStore.annotationNameStore, currentElement.parentId));
						streamPut(tupleOutput, '\t');
					}
					// orientation
					if ( currentElement.beginPos <= currentElement.endPos )
					{
						streamPut(tupleOutput, "+\t");
					}
					else
					{
						streamPut(tupleOutput, "-\t");
					}
					// first annotationId of tuple
					streamPut(tupleOutput, getValue(fragStore.annotationNameStore, position(itCountStore, tupleCountStore)));
					// other annotationIds of first read
					itId = begin(getValue(itT));
					itIdEnd = end(getValue(itT));
					for ( ; itId != itIdEnd && getValue(itId) != INVALID_ID; goNext(itId))
					{
						streamPut(tupleOutput, ":");
						streamPut(tupleOutput, getValue(fragStore.annotationNameStore, getValue(itId)));
					}
					goNext(itId);
					streamPut(tupleOutput, "^");
					// annotationIds of second read
					streamPut(tupleOutput, getValue(fragStore.annotationNameStore, getValue(itId)));
					goNext(itId);
					for ( ; itId != itIdEnd; goNext(itId))
					{
						streamPut(tupleOutput, ":");
						streamPut(tupleOutput,getValue(fragStore.annotationNameStore, getValue(itId)) );
					}
					streamPut(tupleOutput, '\t');
					// tuple count
					streamPut(tupleOutput, getValue(itC));
					streamPut(tupleOutput, '\t');
					// normalized tuple count
					_streamPutDouble(tupleOutput, getValue(itN));
					streamPut(tupleOutput, '\n');
				}
			}
		}
	}
}


//////////////////////////////////////////////////////////////////////////////


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
