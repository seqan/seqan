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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// I/O code for BAM files, based on samtools.  This mirrors the functions
// in store_io_sam.h.
//
// TODO(holtgrew): Only allowing to read from samfile_t should be circumvented somehow. It is possible to get a FILE* from a posix file via fdopen, maybe the same is also true for streams from the <iostream> library?
// TODO(holtgrew): Writing BAM.
// ==========================================================================

#ifndef SEQAN_STORE_IO_BAM_
#define SEQAN_STORE_IO_BAM_

#include <sstream>

#include <seqan/store.h>

namespace seqan {

// ============================================================================
// Enums, Tags, Classes
// ============================================================================

struct Bam_;
typedef Tag<Bam_> Bam;

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// Read a file in BAM format.
//
// TODO(holtgrew): How to support more file types, need to rewrite BAM reading in C++?
template <typename TSpec, typename TConfig>
inline void 
read(samfile_t * file,
     FragmentStore<TSpec, TConfig> & fragStore,
     Bam const &)
{
//IOREV contains lots of TODOs by holtgrew, state of implementation unclear, looks like batchreading, not sure though
    typedef Value<FILE>::Type TValue;
    typedef FragmentStore<TSpec, TConfig> TFragmentStore;
    typedef typename TFragmentStore::TContigPos TContigPos;
    typedef typename Id<TFragmentStore>::Type TId;
        
    // data structure to temporarily store the gaps that need to be inserted in the contig sequences
    typedef MatchMateInfo_<TContigPos, TId> TMatchMateInfo;
    typedef String<TMatchMateInfo> TMatchMateInfos;
    typedef StringSet<String<typename TFragmentStore::TContigGapAnchor> > TContigAnchorGaps;

    // data structure to temporarily store information about match mates
    TMatchMateInfos matchMateInfos;
    TContigAnchorGaps contigAnchorGaps;

    // Header is already loaded in samopen().
    
    // Read in alignments section.
    _readAlignments(file, fragStore, contigAnchorGaps, matchMateInfos, Bam());
    
    // Set the match mate IDs using the information stored in matchMateInfos.
    _generatePairMatchIds(fragStore, matchMateInfos);
    
    convertPairWiseToGlobalAlignment(fragStore, contigAnchorGaps);
}

// Read all alignments from a BAM file.
template <typename TSpec, typename TConfig, typename TContigAnchorGaps, typename TMatchMateInfos>
inline void 
_readAlignments(
        samfile_t * file,
        FragmentStore<TSpec, TConfig> & fragStore,
        TContigAnchorGaps & contigAnchorGaps,   
        TMatchMateInfos & matchMateInfos,
        Bam const &)
{
//IOREV
    // create dummy entries in Sam specific aligned read quality store and aligned read tag store
    // is needed so the ID in the aligned store can be use to access the other stores
    // even if there exists previous entries without
    typedef FragmentStore<TSpec, TConfig> TFragmentStore;
    typedef typename TFragmentStore::TAlignQualityStore TAlignQualityStore;
    typedef typename TFragmentStore::TNameStore TNameStore;
    typedef typename Value<TAlignQualityStore>::Type TAlignQuality;
        
    TAlignQuality q;
    q.score = maxValue(q.score);
    int diff = length(fragStore.alignedReadStore) - length(fragStore.alignQualityStore);
    for(int i = 0; i < diff; ++i)
        appendValue(fragStore.alignQualityStore, q, Generous());
        
    diff = length(fragStore.alignedReadStore) - length(fragStore.alignedReadTagStore);
    for(int i = 0; i < diff; ++i)
        appendValue(fragStore.alignedReadTagStore, "", Generous());
        
    // read in alignments
    Nothing contextBam;
    refresh(fragStore.contigNameStoreCache);
    refresh(fragStore.readNameStoreCache);

    bam1_t * record = bam_init1();
    while (samread(file, record) > 0) {
        // TODO(holtgrew): Allow enabling/disabling of validation.
        bam_validate1(file->header, record);

        _readOneAlignment(file, record, fragStore, contigAnchorGaps, matchMateInfos, Bam(), contextBam);
    }
    bam_destroy1(record);
}

// Read one alignment from a BAM file.
template <typename TSpec, typename TConfig, typename TContigAnchorGaps, typename TMatchMateInfos, typename TContextBAM>
inline void 
_readOneAlignment (
        samfile_t * file,
        bam1_t * record,
		FragmentStore<TSpec, TConfig> & fragStore,
		TContigAnchorGaps & contigAnchorGaps,
		TMatchMateInfos & matchMateInfos,
		Bam const &,
		TContextBAM & contextBam)
{
//IOREV
    // Basic types
    typedef FragmentStore<TSpec, TConfig>										TFragmentStore;
    typedef typename Id<TFragmentStore>::Type									TId;
    typedef typename Size<TFragmentStore>::Type									TSize;
        
    // All fragment store element types
    typedef typename Value<typename TFragmentStore::TContigStore>::Type			TContigElement;
    typedef typename Value<typename TFragmentStore::TLibraryStore>::Type		TLibraryStoreElement;
    typedef typename Value<typename TFragmentStore::TMatePairStore>::Type		TMatePairElement;
    typedef typename Value<typename TFragmentStore::TReadStore>::Type			TReadStoreElement;
    typedef typename Value<typename TFragmentStore::TAlignedReadStore>::Type	TAlignedElement;
    typedef typename Value<typename TFragmentStore::TAlignQualityStore>::Type	TAlignQualityElement;
    typedef typename TAlignedElement::TGapAnchors								TReadGapAnchors;
		
    // Type for sequence in readstore
    typedef typename TFragmentStore::TReadSeq TReadSeq2;
        
    // Type for gap anchor
    typedef typename TFragmentStore::TContigPos									TContigPos;
    typedef Gaps<TReadSeq2, AnchorGaps<TReadGapAnchors> >						TReadGaps;
    typedef Gaps<Nothing, AnchorGaps<typename Value<TContigAnchorGaps>::Type> >	TContigGapsPW;
        
    // Type to temporarily store information about match mates
    typedef typename Value<TMatchMateInfos>::Type								TMatchMateInfo;

    // Get query name.
    String<char> qname = bam1_qname(record);

    // Get flag.
    int flag = record->core.flag;
    bool reverse = (flag & (1 << 4)) == (1 << 4);

    // Get reference name.
    // TODO(holtgrew): Correct behaviour? Or should rname be empty in this case?
    String<char> rname = "*";
    if (record->core.tid >= 0)
        rname = file->header->target_name[record->core.tid];

    // Get begin position.
    TContigPos beginPos = record->core.pos;

    // Get map quality.
    TAlignQualityElement mapQ;
    mapQ.score = record->core.qual;

    // Get CIGAR.
    String<CigarElement<> > cigar;
    resize(cigar, record->core.n_cigar);
    for (unsigned i = 0; i < record->core.n_cigar; ++i) {
        cigar[i] = CigarElement<>(bam1_cigar(record)[i], FromBam());
    }
        
    // Calculate the end position.
    TContigPos endPos;
    _getClippedLength(cigar, endPos);
    endPos = beginPos + endPos;
    // If the read is on the antisense strand switch begin and end position.
    if (reverse)
    {
        TContigPos temp = beginPos;
        beginPos = endPos;
        endPos = temp;
    }

    // Generate gap anchor string for the read.
    TReadGapAnchors readGapAnchors;
        
    // Get read mate reference name.
    String<char> mrnm = "*";
    if (record->core.mtid >= 0)
        mrnm = file->header->target_name[record->core.mtid];

    // Get read mate position.
    TContigPos mPos = record->core.mpos;

    // Template length is ignored for now.

    // Get sequence and associated qualities.
    TReadSeq2 readSeq;
    static const char table[16] = {'=', 'A', 'C', 'N', 'G', 'N', 'N', 'N',
                                   'T', 'N', 'N', 'N', 'N', 'N', 'N', 'N'};
    for (int i = 0; i < record->core.l_qseq; ++i) {
        appendValue(readSeq, Dna5(table[bam1_seqi(bam1_seq(record), i)]));
        assignQualityValue(back(readSeq), bam1_qual(record)[i]);
    }
    SEQAN_ASSERT_GT(length(readSeq), 0u);
    if (reverse)
        reverseComplement(readSeq);

    // Insert alignment gaps.
    TReadGaps readGaps(readSeq, readGapAnchors);
    cigarToGapAnchorRead(cigar, readGaps);
        
    // Convert binary tags into SAM tag string.
    String<char> tags;
    {
        uint8_t * s = bam1_aux(record);
        while (s < record->data + record->data_len) {
            appendValue(tags, *s++);
            appendValue(tags, *s++);
            appendValue(tags, ':');
            uint8_t type = *s++;
            appendValue(tags, type);
            appendValue(tags, ':');

            std::stringstream ss;
            if (type == 'A') {
                append(tags, "a:");
                appendValue(tags, (char)*s);
                s += 1;
            } else if (type == 'C') {
                append(tags, "i:");
                ss << *(uint8_t*)s;
                append(tags, ss.str());
                s += 1;
            } else if (type == 'c') {
                append(tags, "i:");
                ss << *(int8_t*)s;
                append(tags, ss.str());
                s += 1;
            } else if (type == 'S') {
                append(tags, "i:");
                ss << *(uint16_t*)s;
                append(tags, ss.str());
                s += 2;
            } else if (type == 's') {
                append(tags, "i:");
                ss << *(int16_t*)s;
                append(tags, ss.str());
                s += 2;
            } else if (type == 'I') {
                append(tags, "i:");
                ss << *(uint32_t*)s;
                append(tags, ss.str());
                s += 4;
            } else if (type == 'i') {
                append(tags, "i:");
                ss << *(int32_t*)s;
                append(tags, ss.str());
                s += 4;
            } else if (type == 'f') {
                append(tags, "f:");
                ss << *(float*)s;
                append(tags, ss.str());
                s += 4;
            } else if (type == 'd') {
                append(tags, "d:");
                ss << *(float*)s;
                append(tags, ss.str());
                s += 8;
            } else if (type == 'Z' || type == 'H') {
                appendValue(tags, type);
                appendValue(tags, ':');
                while (*s)
                    appendValue(tags, *s++);
                s += 1;
            }

            if (s < record->data + record->data_len)
                appendValue(tags, ' ');
        }
    }

    if (empty(qname) || empty(rname))
        return;
        
    // check if read sequence is already in the store.
    // if so get the ID, otherwise create new entries in the
    // read, read name and mate pair store
        
    TId readId = 0;
    _storeAppendRead(fragStore, readId, qname, readSeq, flag, contextBam);
        
    // check if the contig is already in the store
    // get its ID or create a new one otherwise
    TId contigId = 0;
    _storeAppendContig(fragStore, contigId, rname);

    if (empty(cigar)) return;
		
    // create a new entry in the aligned read store
    TId pairMatchId = appendAlignment(fragStore, readId, contigId, beginPos, endPos, readGapAnchors);
    resize(contigAnchorGaps, length(fragStore.alignedReadStore), Generous());
    TContigGapsPW contigGaps(back(contigAnchorGaps));
    cigarToGapAnchorContig(cigar, contigGaps);
		
    // create entries in Sam specific stores
    appendValue(fragStore.alignQualityStore, mapQ, Generous());
    appendValue(fragStore.alignedReadTagStore, tags, Generous());
        
    // store additional data about match mate temporarily
    // used in the end of the read function to generate match mate IDs
    TId mcontigId = contigId;
    if (mrnm != "*")
    {
        if (mrnm != "=")
            _storeAppendContig(fragStore, mcontigId, mrnm);

        if (flag & 0x40)	// store mate info only for the first read in the pair
        {
            TMatchMateInfo matchMateInfo = {readId, mcontigId, pairMatchId, mPos};
            appendValue(matchMateInfos, matchMateInfo);
            back(fragStore.alignedReadStore).pairMatchId = pairMatchId;
        }
    }
}

// Write a file in BAM format.
//
// TODO(holtgrew): Accept file, not file pointer.
// TODO(holtgrew): How to support more file types, need to rewrite BAM reading in C++?
template <typename TSpec, typename TConfig>
inline void 
write(char const * fileName,
     FragmentStore<TSpec, TConfig> & fragStore,
     Bam const &)
{
//IOREV see above, BAM support apperently needs to be rewritten
    // -----------------------------------------------------------------------
    // Initialize Header.
    // -----------------------------------------------------------------------
    bam_header_t bamHeader;
    bamHeader.n_targets = length(fragStore.contigStore);
    bamHeader.target_name = new char*[bamHeader.n_targets];
    bamHeader.target_len = new uint32_t[bamHeader.n_targets];
    for (int i = 0; i < bamHeader.n_targets; ++i) {
        bamHeader.target_name[i] = new char[length(fragStore.contigNameStore[i]) + 1];
        strncpy(bamHeader.target_name[i], toCString(fragStore.contigNameStore[i]), length(fragStore.contigNameStore[i]) + 1);
        bamHeader.target_len[i] = length(fragStore.contigStore[i].seq);
    }
    bamHeader.hash = 0;
    bamHeader.dict = 0;
    bamHeader.rg2lib = 0;
    bamHeader.l_text = 0;
    bamHeader.n_text = 0;
    bamHeader.text = 0;

    // -----------------------------------------------------------------------
    // Perform I/O.
    // -----------------------------------------------------------------------
    std::cerr << "fileName == " << fileName << std::endl;
    samfile_t * bamFp = samopen(fileName, "wbz", &bamHeader);
    _writeAlignments(bamFp, fragStore, Bam());
    samclose(bamFp);

    // -----------------------------------------------------------------------
    // Cleanup Header.
    // -----------------------------------------------------------------------
    for (int i = 0; i < bamHeader.n_targets; ++i)
        delete [] bamHeader.target_name[i];
    delete [] bamHeader.target_name;
    delete [] bamHeader.target_len;
}

template <typename TSpec, typename TConfig>
inline void _writeAlignments(samfile_t * /*samfile*/,
                             FragmentStore<TSpec, TConfig> & store,
                             Bam const &)
{
//IOREV _stub_ doesn't actually write anything
    typedef FragmentStore<TSpec, TConfig>							TFragmentStore;

    typedef typename TFragmentStore::TReadStore						TReadStore;
    typedef typename TFragmentStore::TReadSeqStore					TReadSeqStore;
    typedef typename TFragmentStore::TAlignedReadStore				TAlignedReadStore;
    typedef typename TFragmentStore::TContigStore					TContigStore;
    typedef typename TFragmentStore::TReadSeq						TReadSeq;

    typedef typename Value<TReadStore>::Type						TRead;
    typedef typename Value<TReadSeqStore>::Type						TReadSeqStored;
    typedef typename Value<TContigStore>::Type						TContig;
    typedef typename Value<TAlignedReadStore>::Type					TAlignedRead;

    typedef typename TContig::TContigSeq							TContigSeq;
    typedef typename Iterator<TAlignedReadStore, Standard>::Type	TAlignIter;
    typedef typename Iterator<TReadSeqStored, Standard>::Type		TReadSeqIter;
    typedef typename Id<TAlignedRead>::Type							TId;

    typedef Gaps<TReadSeq, AnchorGaps<typename TAlignedRead::TGapAnchors> >	TReadGaps;
    typedef Gaps<Nothing, AnchorGaps<typename TContig::TGapAnchors> >	TContigGaps;

    String<int> mateIndex;	// store outer library size for each pair match (indexed by pairMatchId)
    calculateMateIndices(mateIndex, store);
		
    TAlignIter it = begin(store.alignedReadStore, Standard());
    TAlignIter itEnd = end(store.alignedReadStore, Standard());
    TAlignIter mit = it;
    String<CigarElement<> > cigar;
    TReadSeq readSeq;

    // The plan is to build BAM records and write them out to the file using libbam.
    bam1_t bamRecord;
    bamRecord.l_aux = 0;  // No auxiliary data / tags for now.
    
    for(; it != itEnd; ++it)
    {
        TId alignedId = (*it).id;
        TId readId = (*it).readId;
        TId mateIdx = TRead::INVALID_ID;

        if ((*it).pairMatchId != TRead::INVALID_ID)
            mateIdx = mateIndex[2*(*it).pairMatchId + getMateNo(store, (*it).readId)];

        TContigGaps	contigGaps(/*store.contigStore[(*it).contigId].seq, */store.contigStore[(*it).contigId].gaps);
        __int64 pos = positionGapToSeq(contigGaps, _min((*it).beginPos, (*it).endPos)) + 1;
        __int64 mpos = 0;
        int isize = 0;
        unsigned short flag = 0;

        if ((*it).beginPos > (*it).endPos)
            flag |= 0x0010;			

        // calculate flags, mpos, isize
        if (mateIdx < length(store.alignedReadStore))
        {
            mit = begin(store.alignedReadStore, Standard()) + mateIdx;
            if ((*it).contigId == (*mit).contigId)
            {
                mpos = positionGapToSeq(contigGaps, _min((*mit).beginPos, (*mit).endPos)) + 1;
                if ((*it).beginPos < (*mit).beginPos)
                    isize = positionGapToSeq(contigGaps, _max((*mit).beginPos, (*mit).endPos) - 1) + 2 - pos;
                else
                    isize = mpos - positionGapToSeq(contigGaps, _max((*it).beginPos, (*it).endPos) - 1) - 2;
            }
            flag |= 0x0002;
            if ((*mit).beginPos > (*mit).endPos)
                flag |= 0x0020;				
        }
        else
            flag |= 0x0008;					// mate is unmapped (actually we should check if the mate has no match at all)
			
        signed char mateNo = getMateNo(store, readId);
        if (mateNo == 0) flag |= 0x0040;	// this read is the first in the pair
        if (mateNo == 1) flag |= 0x0080;	// this read is the second in the pair

        if (readId < length(store.readStore))
        {
            TRead &read = store.readStore[readId];
            if (read.matePairId != TRead::INVALID_ID)
                flag |= 0x0001;
        }

        // Compute length of data field of the BAM record.  The contents will
        // be filled below.
        bamRecord.core.n_cigar = length(cigar);
        bamRecord.core.l_qname = length(store.readNameStore[readId]) + 1;
        bamRecord.core.l_qseq = length(store.readSeqStore[readId]);
        bamRecord.data_len = bamRecord.core.l_qname + bamRecord.core.n_cigar + bamRecord.core.l_qseq + (bamRecord.core.l_qseq + 1) / 2 + bamRecord.l_aux;
        bamRecord.m_data = bamRecord.data_len;
        bamRecord.data = new uint8_t[bamRecord.data_len];

        // <qname>
        strncpy(bam1_qname(&bamRecord), toCString(store.readNameStore[readId]), bamRecord.core.l_qname);
            
        // <flag>
        bamRecord.core.flag = flag;

        // <rname>
        bamRecord.core.tid = (*it).contigId < length(store.contigNameStore) ? (*it).contigId : -1;
            
        // <pos>
        bamRecord.core.pos = pos - 1;
            
        // <mapq>
        if (alignedId < length(store.alignQualityStore))
            bamRecord.core.qual = store.alignQualityStore[alignedId].score;
        else
            bamRecord.core.qual = 255;
            
        // get read sequence
        if (readId < length(store.readSeqStore)) {
            readSeq = store.readSeqStore[readId];
            if ((*it).beginPos <= (*it).endPos) {
                setBeginPosition(contigGaps, (*it).beginPos);
                setEndPosition(contigGaps, (*it).endPos);
            } else {
                setBeginPosition(contigGaps, (*it).endPos);
                setEndPosition(contigGaps, (*it).beginPos);
                reverseComplement(readSeq);
            }
        } else {
            clear(readSeq);
        }
			
        // <cigar>
        TReadGaps readGaps(readSeq, (*it).gaps);
        getCigarString(cigar, contigGaps, readGaps);
        for (unsigned i = 0; i < length(cigar); ++i)
            bam1_cigar(&bamRecord)[i] = toBamCigarElement(cigar[i]);
            
        // <mrnm>
        if ((mateIdx < length(store.alignedReadStore)))
            bamRecord.core.mtid = ((*mit).contigId < length(store.contigNameStore)) ? (*mit).contigId : -1;
        else
            bamRecord.core.mtid = -1;
            
        // <mpos>
        bamRecord.core.mpos = mpos;
            
        // <isize>
        bamRecord.core.isize = isize;

        // <seq>
        // TODO(holtgrew): Compressed sequence.
            
        // <qual>
        TReadSeqIter readIt = begin(store.readSeqStore[readId], Standard());
        TReadSeqIter readItEnd = end(store.readSeqStore[readId], Standard());
        for (int i = 0; readIt != readItEnd; ++readIt, ++i)
            bam1_qual(&bamRecord)[i] = (char)(getQualityValue(*readIt) + 33);
            
        // <tags>
        // TODO(holtgrew): No tags yet.

        delete [] bamRecord.data;  // TODO(holtgrew): Do not reallocate in each iteration.
    }   
}

}  // namespace seqan

#endif  // SEQAN_STORE_IO_BAM_
