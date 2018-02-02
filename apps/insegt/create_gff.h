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

namespace seqan
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
        readOutput << getValue(fragStore.readNameStore, position(itCountStore, readAnnoStore))
                   << '\t';
		
		if (empty(getValue(itCountStore).annoIds) )
		{
            readOutput << ".\t.\t.\t.\t.";
		}
		else
		{
			// contig-name:
			readOutput << getValue(fragStore.contigNameStore, getValue(itCountStore).contigId)
                                   << '\t';
			
			itAnnoIds = begin(getValue(itCountStore).annoIds);
			itAnnoIdsEnd = end(getValue(itCountStore).annoIds);
			while (itAnnoIds != itAnnoIdsEnd && front(*itAnnoIds) == INVALID_ID)
			{
				goNext(itAnnoIds);
			}
	
			if (itAnnoIds != itAnnoIdsEnd) // not only INVALID_IDS
			{
				firstId = front(*itAnnoIds);
				
				// orientation:
				if (getValue(fragStore.annotationStore, firstId).beginPos <= getValue(fragStore.annotationStore, firstId).endPos)
				{
                    readOutput << "+\t";
				}
				else
                {
                    readOutput << "-\t";
                }
				
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
                                readOutput << ';';
                            readOutput << getValue(fragStore.annotationNameStore, getValue(itId));
							help = true;
						}
						else if (getValue(itId) != INVALID_ID && getValue(fragStore.annotationStore, getValue(itId)).parentId != INVALID_ID &&	// get other parentIds 
							!isElement_unsorted(getValue(fragStore.annotationStore, getValue(itId)).parentId, allParentIds) ) //?
						{
							appendValue(allParentIds, getValue(fragStore.annotationStore, getValue(itId)).parentId, Generous() );
						}
					}
					if ( !empty(getValue(itCountStore).parentIds) && position(itAnnoIds, getValue(itCountStore).annoIds) != endPosition(getValue(itCountStore).annoIds) - 1)
                        readOutput << ':';
				}
				if (!empty(getValue(itCountStore).parentIds))
					readOutput << '\t' << getValue(fragStore.annotationNameStore, front(allParentIds)) << '\t';
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
                                    readOutput << ';';
								readOutput << getValue(fragStore.annotationNameStore, getValue(itId));
								invalid = false;
							}
						}
						if (invalid)
                            readOutput << "UNKNOWN_REGION";
						if (position(itAnnoIds, getValue(itCountStore).annoIds) != endPosition(getValue(itCountStore).annoIds) - 1)
                            readOutput << ':';
					}
                    readOutput << '\t';
					if (getValue(itP) != INVALID_ID)
                        readOutput << getValue(fragStore.annotationNameStore, getValue(itP));
					else
                        readOutput << "NO_PARENT";
                    readOutput << '\t';
				}
			}
			else  // only INVALID_IDS
			{	
                readOutput << ".\t";
						
				// invalid_ids for each interval
				itAnnoIds = begin(getValue(itCountStore).annoIds);
				itAnnoIdsEnd = end(getValue(itCountStore).annoIds);
				for ( ; itAnnoIds != itAnnoIdsEnd; goNext(itAnnoIds))
				{
                    readOutput << "UNKNOWN_REGION";
					if (position(itAnnoIds, getValue(itCountStore).annoIds) != endPosition(getValue(itCountStore).annoIds) - 1)
                        readOutput << ':';
				}
                readOutput << "\tUNKNOWN_REGION";
			}
		}
        readOutput << "\n";
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
            annoOutput << "INVALID_ID\t";
        else
            annoOutput << getValue(fragStore.contigNameStore, getValue(itAnno).contigId) << '\t';
        annoOutput << "Annotation_Count\tregion\t";
        // startposition endposition orientation .
        if (getValue(itAnno).beginPos == INVALID_POS)
		{
			annoOutput << ".\t.\t" << getValue(itCount);
			if (getValue(itAnno).parentId == INVALID_ID || (fragStore.annotationStore[getValue(itAnno).parentId].typeId != INVALID_ID && fragStore.annotationTypeStore[fragStore.annotationStore[getValue(itAnno).parentId].typeId] == "<root>"))
			{
				if (mapValue(mapO, position(itAnno, fragStore.annotationStore)) == 0)
					annoOutput << "\t+\t.\t";
				else
					annoOutput << "\t-\t.\t";
			}
			else
            {
                annoOutput << "\t.\t.\t";
            }
		}
		else
		{
			if (getValue(itAnno).beginPos <= getValue(itAnno).endPos)
				annoOutput << getValue(itAnno).beginPos + 1
				           << '\t'
				           << getValue(itAnno).endPos
				           << '\t'
				           << getValue(itCount)
				           << "\t+\t.\t";
			else
				annoOutput << getValue(itAnno).endPos + 1
				           << '\t'
				           << getValue(itAnno).beginPos
				           << '\t'
				           << getValue(itCount)
				           << "\t-\t.\t";
		}
		// annotation-name (parent annotation-name)
		if (getValue(itAnno).parentId == INVALID_ID || (fragStore.annotationStore[getValue(itAnno).parentId].typeId != INVALID_ID && fragStore.annotationTypeStore[fragStore.annotationStore[getValue(itAnno).parentId].typeId] == "<root>"))
			annoOutput << "ID="
			           << getValue(fragStore.annotationNameStore, position(itAnno, fragStore.annotationStore))
			           << ';';
		else
			annoOutput << "ID="
			           << getValue(fragStore.annotationNameStore, position(itAnno, fragStore.annotationStore))
			           << ";ParentID="
			           << getValue(fragStore.annotationNameStore, getValue(itAnno).parentId)
			           << ';';
        annoOutput << formattedNumber("%f", *itNorm) << ";\n";
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
					tupleOutput << getValue(fragStore.contigNameStore, currentElement.contigId) << '\t';
					// parent-name
					if (currentElement.parentId == INVALID_ID )
						tupleOutput << "NO_PARENT\t";
					else
						tupleOutput << getValue(fragStore.annotationNameStore, currentElement.parentId) << '\t';
					// orientation
					if ( currentElement.beginPos <= currentElement.endPos )
						tupleOutput << "+\t";
					else
						tupleOutput << "-\t";
					// first annotationId of tuple (store implicit)
					tupleOutput << getValue(fragStore.annotationNameStore, position(itCountStore, tupleCountStore));
					// other annotationIds
					itId = begin(getValue(itT));
					itIdEnd = end(getValue(itT));
					for ( ; itId != itIdEnd; goNext(itId))
						tupleOutput << ":" << getValue(fragStore.annotationNameStore, getValue(itId));
					tupleOutput << '\t';
					// tuple count
					tupleOutput << getValue(itC) << '\t';
					// normalized tuple count
                    tupleOutput << formattedNumber("%f", *itN) << '\n';
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
					tupleOutput << getValue(fragStore.contigNameStore, currentElement.contigId) << '\t';
					// parent-name
					if (currentElement.parentId == INVALID_ID )
						tupleOutput << "NO_PARENT\t";
					else
						tupleOutput << getValue(fragStore.annotationNameStore, currentElement.parentId) << '\t';
					// orientation
					if ( currentElement.beginPos <= currentElement.endPos )
						tupleOutput << "+\t";
					else
						tupleOutput << "-\t";
					// first annotationId of tuple
					tupleOutput << getValue(fragStore.annotationNameStore, position(itCountStore, tupleCountStore));
					// other annotationIds of first read
					itId = begin(getValue(itT));
					itIdEnd = end(getValue(itT));
					for ( ; itId != itIdEnd && getValue(itId) != INVALID_ID; goNext(itId))
						tupleOutput << ":" << getValue(fragStore.annotationNameStore, getValue(itId));
					goNext(itId);
                    tupleOutput << '^';
					// annotationIds of second read
					tupleOutput << getValue(fragStore.annotationNameStore, getValue(itId));
					goNext(itId);
					for ( ; itId != itIdEnd; goNext(itId))
						tupleOutput << ':' << getValue(fragStore.annotationNameStore, getValue(itId));
                    tupleOutput << '\t';
					// tuple count
					tupleOutput << *itC << '\t';
					// normalized tuple count
                    tupleOutput << formattedNumber("%f", *itN) << '\n';
				}
			}
		}
	}
}


//////////////////////////////////////////////////////////////////////////////


}// namespace seqan

#endif //#ifndef SEQAN_HEADER_...
