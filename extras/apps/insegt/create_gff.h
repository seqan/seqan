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
createReadCountGFF(TFile & readOutput, TReadAnnoStore & readAnnoStore, FragmentStore<TSpec, TConfig> & me)
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
		_streamWrite(readOutput, getValue(me.readNameStore, position(itCountStore, readAnnoStore)) );
		_streamPut(readOutput, '\t');
		
		if (empty(getValue(itCountStore).annoIds) )
		{
			_streamWrite(readOutput, ".\t.\t.\t.\t.");
		}
		else
		{
			// contig-name:
			_streamWrite(readOutput, getValue(me.contigNameStore, getValue(itCountStore).contigId) );
			_streamPut(readOutput, '\t');
			
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
				if (getValue(me.annotationStore, firstId).beginPos <= getValue(me.annotationStore, firstId).endPos)
				{
					_streamWrite(readOutput, "+\t");
				}
				else _streamWrite(readOutput, "-\t");
				
				// Annotation-Ids:
				
				allParentIds =  getValue(itCountStore).parentIds;
				// output for first parentId if possible; get other parentIds, in which the read doesn't map (entirely with all of his intervals)
				itAnnoIds = begin(getValue(itCountStore).annoIds);
				itAnnoIdsEnd = end(getValue(itCountStore).annoIds);
				for ( ; itAnnoIds != itAnnoIdsEnd; goNext(itAnnoIds))
				{
					help = false;
					itId = begin(getValue(itAnnoIds));	
					itIdEnd = end(getValue(itAnnoIds));
					for ( ; itId != itIdEnd; goNext(itId))
					{
						if (!empty(getValue(itCountStore).parentIds) && getValue(itId) != INVALID_ID && 	// if current parentId == first id of parentIds (of read entry in readAnnoStore)
						    getValue(me.annotationStore, getValue(itId)).parentId == front(allParentIds))
						{
							if (help) 	// not the first annotation for this read-interval -> ";" sign for overlapping annotations
							{
								_streamWrite(readOutput, ";");
							}
							_streamWrite(readOutput, getValue(me.annotationNameStore, getValue(itId)) );
							help = true;
						}
						else if (getValue(itId) != INVALID_ID && getValue(me.annotationStore, getValue(itId)).parentId != INVALID_ID &&	// get other parentIds 
							!isElement_unsorted(getValue(me.annotationStore, getValue(itId)).parentId, allParentIds) ) //?
						{
							appendValue(allParentIds, getValue(me.annotationStore, getValue(itId)).parentId, Generous() );
						}
					}
					if ( !empty(getValue(itCountStore).parentIds) && position(itAnnoIds, getValue(itCountStore).annoIds) != endPosition(getValue(itCountStore).annoIds) - 1)
					{
						_streamWrite(readOutput, ":");
					}
				}
				if (!empty(getValue(itCountStore).parentIds))
				{
					_streamPut(readOutput, '\t');
					_streamWrite(readOutput, getValue(me.annotationNameStore, front(allParentIds)) );
					_streamWrite(readOutput, "\t");
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
						itId = begin(getValue(itAnnoIds));	
						itIdEnd = end(getValue(itAnnoIds));
						for ( ; itId != itIdEnd; goNext(itId))
						{
							if (getValue(itId) != INVALID_ID && getValue(me.annotationStore, getValue(itId)).parentId == getValue(itP))
							{
								if (!invalid)	// not the first annotation for this interval -> ";" sign for overlapping annotations
								{
									_streamWrite(readOutput, ";");
								}
								_streamWrite(readOutput, getValue(me.annotationNameStore, getValue(itId)) );
								invalid = false;
							}
						}
						if (invalid)
						{
							_streamWrite(readOutput, "UNKNOWN_REGION");
						}
						if (position(itAnnoIds, getValue(itCountStore).annoIds) != endPosition(getValue(itCountStore).annoIds) - 1)
						{
							_streamWrite(readOutput, ":");
						}
					}
					_streamPut(readOutput, '\t');
					if (getValue(itP) != INVALID_ID) _streamWrite(readOutput, getValue(me.annotationNameStore, getValue(itP)) );
					else _streamWrite(readOutput, "NO_PARENT" );
					_streamWrite(readOutput, "\t");
				}
			}
			else  // only INVALID_IDS
			{	
				_streamWrite(readOutput, ".\t");
						
				// invalid_ids for each interval
				itAnnoIds = begin(getValue(itCountStore).annoIds);
				itAnnoIdsEnd = end(getValue(itCountStore).annoIds);
				for ( ; itAnnoIds != itAnnoIdsEnd; goNext(itAnnoIds))
				{
					_streamWrite(readOutput, "UNKNOWN_REGION");
					if (position(itAnnoIds, getValue(itCountStore).annoIds) != endPosition(getValue(itCountStore).annoIds) - 1)
					{
						_streamWrite(readOutput, ":");
					}
					
				}
				_streamPut(readOutput, '\t');
				_streamWrite(readOutput, "UNKNOWN_REGION");
			}
		}
		_streamPut(readOutput, '\n');	
	}	
}



//////////////////////////////////////////////////////////////////////////////
//create AnnoCountGFF
//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TAnnoCountStore, typename TAnnoNormStore, typename TSpec, typename TConfig, typename TMap>
inline void
createAnnoCountGFF(TFile & annoOutput, TAnnoCountStore & annoCountStore, TAnnoNormStore &annoNormStore, FragmentStore<TSpec, TConfig> & me, TMap &mapO)
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
	TAnnoIter itAnno = begin(me.annotationStore);
	TNormIter itNorm = begin(annoNormStore);
	
	for ( ; itCount != itCountEnd; goNext(itCount), goNext(itAnno), goNext(itNorm))
	{
		// contig-name
		if ( getValue(itAnno).contigId == INVALID_ID )
		{
			_streamWrite(annoOutput, "INVALID_ID");
			_streamPut(annoOutput, '\t');
		}
		else
		{
			_streamWrite(annoOutput, getValue(me.contigNameStore, getValue(itAnno).contigId));
			_streamPut(annoOutput, '\t');
		}
		_streamWrite(annoOutput, "Annotation_Count\tregion\t");
		// startposition endposition orientation . 
		if (getValue(itAnno).beginPos == INVALID_POS)
		{
			_streamWrite(annoOutput, ".\t.\t");
			_streamPutInt(annoOutput, getValue(itCount));
			if (getValue(itAnno).parentId == INVALID_ID)
			{
				if (mapValue(mapO, position(itAnno, me.annotationStore)) == 0)
				{
					_streamWrite(annoOutput, "\t+\t.\t");
				}
				else
				{
					_streamWrite(annoOutput, "\t-\t.\t");
				}
			}
			else	_streamWrite(annoOutput, "\t.\t.\t");
		}
		else
		{
			if (getValue(itAnno).beginPos <= getValue(itAnno).endPos)
			{
				_streamPutInt(annoOutput, getValue(itAnno).beginPos +1);
				_streamPut(annoOutput, '\t');
				_streamPutInt(annoOutput, getValue(itAnno).endPos);
				_streamPut(annoOutput, '\t');
				_streamPutInt(annoOutput, getValue(itCount));
				_streamWrite(annoOutput, "\t+\t.\t");
			}
			else
			{
				_streamPutInt(annoOutput, getValue(itAnno).endPos);
				_streamPut(annoOutput, '\t');
				_streamPutInt(annoOutput, getValue(itAnno).beginPos +1);
				_streamPut(annoOutput, '\t');
				_streamPutInt(annoOutput, getValue(itCount));
				_streamWrite(annoOutput, "\t-\t.\t");
			}
		}
		// annotation-name (parent annotation-name)
		if (getValue(itAnno).parentId == INVALID_ID)
		{
			_streamWrite(annoOutput, "ID=");
			_streamWrite(annoOutput, getValue(me.annotationNameStore, position(itAnno, me.annotationStore)) );
			_streamPut(annoOutput, ';');
		}
		else
		{
			_streamWrite(annoOutput, "ID=");
			_streamWrite(annoOutput, getValue(me.annotationNameStore, position(itAnno, me.annotationStore)));
			_streamWrite(annoOutput, ";ParentID=");
			_streamWrite(annoOutput, getValue(me.annotationNameStore, getValue(itAnno).parentId));
			_streamPut(annoOutput, ';');
		}
		//_streamPutFloat(annoOutput, getValue(itNorm));
		_streamPutDouble(annoOutput, getValue(itNorm));
		_streamWrite(annoOutput, ";\n");
	}
}


//////////////////////////////////////////////////////////////////////////////
//create tupleCountGFF
//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TTupleCountStore, typename TSpec, typename TConfig>
inline void
createTupleCountGFF(TFile & tupleOutput, TTupleCountStore & tupleCountStore, FragmentStore<TSpec, TConfig> & me, unsigned thresholdCount, double thresholdRPKM)
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
			currentElement = getValue(me.annotationStore, position(itCountStore, tupleCountStore));
	
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
					_streamWrite(tupleOutput, getValue(me.contigNameStore, currentElement.contigId));
					_streamPut(tupleOutput, '\t');
					// parent-name
					if (currentElement.parentId == INVALID_ID )
					{
						_streamWrite(tupleOutput, "NO_PARENT\t");
					}
					else
					{
						_streamWrite(tupleOutput, getValue(me.annotationNameStore, currentElement.parentId));
						_streamPut(tupleOutput, '\t');
					}
					// orientation
					if ( currentElement.beginPos <= currentElement.endPos )
					{
						_streamWrite(tupleOutput, "+\t");
					}
					else
					{
						_streamWrite(tupleOutput, "-\t");
					}
					// first annotationId of tuple (store implicit)
					_streamWrite(tupleOutput, getValue(me.annotationNameStore, position(itCountStore, tupleCountStore)));
					// other annotationIds
					itId = begin(getValue(itT));
					itIdEnd = end(getValue(itT));
					for ( ; itId != itIdEnd; goNext(itId))
					{
						_streamWrite(tupleOutput, ":");
						_streamWrite(tupleOutput, getValue(me.annotationNameStore, getValue(itId)));
					}
					_streamPut(tupleOutput, '\t');
					// tuple count
					_streamPutInt(tupleOutput, getValue(itC));
					_streamPut(tupleOutput, '\t');
					// normalized tuple count
					//_streamPutFloat(tupleOutput, getValue(itN));
					_streamPutDouble(tupleOutput, getValue(itN));
					_streamPut(tupleOutput, '\n');
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
					_streamWrite(tupleOutput, getValue(me.contigNameStore, currentElement.contigId));
					_streamPut(tupleOutput, '\t');
					// parent-name
					if (currentElement.parentId == INVALID_ID )
					{
						_streamWrite(tupleOutput, "NO_PARENT\t");
					}
					else
					{
						_streamWrite(tupleOutput, getValue(me.annotationNameStore, currentElement.parentId));
						_streamPut(tupleOutput, '\t');
					}
					// orientation
					if ( currentElement.beginPos <= currentElement.endPos )
					{
						_streamWrite(tupleOutput, "+\t");
					}
					else
					{
						_streamWrite(tupleOutput, "-\t");
					}
					// first annotationId of tuple
					_streamWrite(tupleOutput, getValue(me.annotationNameStore, position(itCountStore, tupleCountStore)));
					// other annotationIds of first read
					itId = begin(getValue(itT));
					itIdEnd = end(getValue(itT));
					for ( ; itId != itIdEnd && getValue(itId) != INVALID_ID; goNext(itId))
					{
						_streamWrite(tupleOutput, ":");
						_streamWrite(tupleOutput, getValue(me.annotationNameStore, getValue(itId)));
					}
					goNext(itId);
					_streamWrite(tupleOutput, "^");
					// annotationIds of second read
					_streamWrite(tupleOutput, getValue(me.annotationNameStore, getValue(itId)));
					goNext(itId);
					for ( ; itId != itIdEnd; goNext(itId))
					{
						_streamWrite(tupleOutput, ":");
						_streamWrite(tupleOutput,getValue(me.annotationNameStore, getValue(itId)) );
					}
					_streamPut(tupleOutput, '\t');
					// tuple count
					_streamPutInt(tupleOutput, getValue(itC));
					_streamPut(tupleOutput, '\t');
					// normalized tuple count
					//_streamPutFloat(tupleOutput, getValue(itN));
					_streamPutDouble(tupleOutput, getValue(itN));
					_streamPut(tupleOutput, '\n');
				}
			}
		}
	}
}


//////////////////////////////////////////////////////////////////////////////


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
