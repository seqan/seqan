/*=========================================================================
  Copyright (C) 2009 by Stephan Aiche

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your options) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ==========================================================================
  $Id$
 ==========================================================================*/

#ifndef REPSEP_HEADER_ASSEMBLY_PARSER_H
#define REPSEP_HEADER_ASSEMBLY_PARSER_H

#include <seqan/store.h>
#include <seqan/modifier.h>

#include <map>

//#define REPSEP_DEBUG_ASSEMBLY_COL_PARSER

using namespace seqan;

template<typename TGapsAnchorArray, typename TUngappedSequence, typename TReadPos, typename TGappedSequence>
inline void 
_gap_sequence(TGapsAnchorArray const & gaps,
             TUngappedSequence const & ungapped_sequence,
             TReadPos beginClr,
             TReadPos endClr,
             TGappedSequence & seq)
{
    clear(seq);
    typedef typename Size<TGapsAnchorArray>::Type TSize;
    typedef typename Value<TGappedSequence>::Type TGappedValue;

    TSize ungapped_position = beginClr;
    TSize gapped_position = 0;

    for(TSize g = 0; g < length(gaps) ; ++g) {
        while (static_cast<TReadPos>(ungapped_position) < gaps[g].seqPos && static_cast<TReadPos>(ungapped_position) < endClr) {
            // put 
            append(seq, TGappedValue(ungapped_sequence[ungapped_position] , ungapped_position));
            ++ungapped_position;
            ++gapped_position;
        }
        while (static_cast<TReadPos>(gapped_position) < gaps[g].gapPos) {
            append(seq, TGappedValue('-', endClr + 1));        
            ++gapped_position;
        }
    }

    // fill up the rest
    while (static_cast<TReadPos>(ungapped_position) < endClr) {
        append(seq, TGappedValue(ungapped_sequence[ungapped_position] , ungapped_position));
        ++ungapped_position;
    }
}


// build alignment matrix based on Amos-Assembly -- the main worker
template<typename TSpec, typename TConfig, typename TId, typename TAnnotatedCandidateColumn, typename TScannerType>
inline void 
parseContig(FragmentStore<TSpec, TConfig> const& fragStore,
           TId const contigId,
           String<TAnnotatedCandidateColumn> & candidates,
           TScannerType const algo_spec)

{
    clear(candidates);

    typedef FragmentStore<TSpec, TConfig> TFragmentStore;
    typedef typename Size<TFragmentStore>::Type TSize;
    typedef typename TFragmentStore::TContigPos TContigPos;
    typedef typename TFragmentStore::TReadPos TReadPos;
    
    typedef typename Value<TAnnotatedCandidateColumn, 2>::Type TCandidateColumn;

    // must be triplet
    typedef typename Value<TCandidateColumn>::Type TAssignedReadChar;

	// All fragment store element types
    //typedef typename Value<typename TFragmentStore::TReadStore>::Type TReadStoreElement;
    typedef typename Value<typename TFragmentStore::TAlignedReadStore>::Type TAlignedElement;

    // 
    typedef Pair<char, TReadPos> TGappedSequenceChar;
    typedef String<TGappedSequenceChar> TGappedSequence;

    // Sort by begin position .. 
    sortAlignedReads(fragStore.alignedReadStore, SortBeginPos());
	// Sort aligned reads according to contig id
    sortAlignedReads(fragStore.alignedReadStore, SortContigId());

    TGappedSequence gapped_consensus;

    _gap_sequence(fragStore.contigStore[contigId].gaps, fragStore.contigStore[contigId].seq, (TReadPos) 0 , (TReadPos) length(fragStore.contigStore[contigId].seq), gapped_consensus);
    TSize l_gapped_consensus = length(gapped_consensus);

    typedef typename Iterator<typename TFragmentStore::TAlignedReadStore const>::Type TAlignIter;
    TAlignIter alignIt = lowerBoundAlignedReads(fragStore.alignedReadStore, contigId, SortContigId());
    TAlignIter alignItEnd = upperBoundAlignedReads(fragStore.alignedReadStore, contigId, SortContigId());

    typedef std::multiset< TAlignedElement , _LessAlignedRead< TAlignedElement, SortEndPos > > TAlignedReadSet;
    
    typedef typename TAlignedReadSet::iterator TAlignedReadSetIter;
    
    TAlignedReadSet current_read_set;
#ifdef REPSEP_DEBUG_ASSEMBLY_COL_PARSER // debug output of all reads in this contig
    TAlignIter debug_alignIt = lowerBoundAlignedReads(fragStore.alignedReadStore, contigId, SortContigId());
    TAlignIter debug_alignItEnd = upperBoundAlignedReads(fragStore.alignedReadStore, contigId, SortContigId());
    
    while(debug_alignIt != debug_alignItEnd)
    {
      cout << debug_alignIt->readId << " (" << fragStore.readNameStore[debug_alignIt->readId] << ")" << endl;
      goNext(debug_alignIt);
    }
    
#endif
    for(TSize p = 0 ; p < l_gapped_consensus ; ++p)  
    {
#ifdef REPSEP_DEBUG_ASSEMBLY_COL_PARSER
        cout << "Start inspecting column " << p << endl;
#endif    
        // check if we need to remove some of the reads
        TAlignedReadSetIter end_erase = current_read_set.begin();
        TAlignedReadSetIter set_end = current_read_set.end();
        if(end_erase != set_end && _max(end_erase->beginPos, end_erase->endPos) == static_cast<TContigPos>(p)) {
            // we need to remove some of our AlignedReads

            TAlignedReadSetIter start_erase = end_erase;

            while (end_erase != set_end && _max(end_erase->beginPos, end_erase->endPos) == static_cast<TContigPos>(p)) {
#ifdef REPSEP_DEBUG_ASSEMBLY_COL_PARSER
                cout << "-> will remove: " << end_erase->readId << " (" << fragStore.readNameStore[end_erase->readId] << ")" << endl;
#endif                
                ++end_erase;
            }
            // do we need to go to previous ???

            current_read_set.erase(start_erase, end_erase);
        }

        // add all open reads into set
        while(alignIt != alignItEnd && _min(alignIt->beginPos, alignIt->endPos) == static_cast<TContigPos>(p)) {
            // add
            current_read_set.insert(*alignIt);
#ifdef REPSEP_DEBUG_ASSEMBLY_COL_PARSER
            cout << "added" << endl;
            cout << alignIt->readId << " (" << fragStore.readNameStore[alignIt->readId] << ")" << endl;
#endif            
            goNext(alignIt);
        }
#ifdef REPSEP_DEBUG_ASSEMBLY_COL_PARSER
        cout << "after add cycle iterator points to " << endl;
        cout << alignIt->readId << " (" << fragStore.readNameStore[alignIt->readId] << ")" << endl;
#endif
        // since also a gap can be off interest we also inspect those        
        TCandidateColumn column; 
        
        // inspect consensus 
        TAlignedReadSetIter iter_cr = current_read_set.begin();
        TAlignedReadSetIter iter_cr_end = current_read_set.end();

#ifdef REPSEP_DEBUG_ASSEMBLY_COL_PARSER
        cout << "at position " << p << " following reads are active (" << current_read_set.size() << "):" << endl;
#endif             
        while(iter_cr != iter_cr_end) {
#ifdef REPSEP_DEBUG_ASSEMBLY_COL_PARSER
            cout << iter_cr->readId << " (" << fragStore.readNameStore[iter_cr->readId] << ")" << endl;
#endif          
            // get position in sequence 
            TReadPos offset = _min(iter_cr->beginPos, iter_cr->endPos);
            TGappedSequence gapped_read;        

	        TReadPos begClr = 0;
	        TReadPos endClr = 0;
	        getClrRange(fragStore, *iter_cr, begClr, endClr);        
            String<Dna5Q> ungapped_seq = fragStore.readSeqStore[iter_cr->readId];
            if(iter_cr->beginPos > iter_cr->endPos) {
                reverseComplement(ungapped_seq);
                // swap the clr's
	            TReadPos tmp = begClr;
	            begClr = length(ungapped_seq) - endClr;
	            endClr = length(ungapped_seq) - tmp;
            }
            _gap_sequence(iter_cr->gaps, ungapped_seq, begClr, endClr, gapped_read);
            
            // increase vote 
            append(column, TAssignedReadChar(_sequenceCharacter(gapped_read[p - offset]), iter_cr->readId, _positionInRead(gapped_read[p - offset])));

            ++iter_cr;
        }

        // inspect column and mark as candidate .. 
        if( isCandidate(fragStore, contigId, gapped_consensus[p], column, algo_spec ) )
        {
            appendValue(candidates, TAnnotatedCandidateColumn(p, column));
        }
    }
}




#endif
