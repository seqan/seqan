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

#ifndef SEQAN_HEADER_STORE_IO_GFF_H
#define SEQAN_HEADER_STORE_IO_GFF_H

namespace SEQAN_NAMESPACE_MAIN
{

/**
.Tag.File Format.tag.Gff:
	Gff annotation file.
..include:seqan/store.h
*/
struct TagGff_;
typedef Tag<TagGff_> const Gff;

/**
.Tag.File Format.tag.Gtf:
	Gtf annotation file.
..include:seqan/store.h
*/
struct TagGtf_;
typedef Tag<TagGtf_> const Gtf;

//////////////////////////////////////////////////////////////////////////////
// _parseReadGffIdentifier
    
    template<typename TFile, typename TString, typename TChar>
    inline void
    _parseReadGffIdentifier(TFile & file, TString & str, TChar& c)
    {
//IOREV _nodoc_ _hasCRef_
        if (c == ' ' || c == '\t' || c == '\n') return;
        appendValue(str, c);
        while (!_streamEOF(file)) 
		{
            c = _streamGet(file);
            if (c == ' ' || c == '\t' || c == '\n') return;
            appendValue(str, c);
        }
    }
	
//////////////////////////////////////////////////////////////////////////////
// skip entry until whitespace
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

    template<typename TFile, typename TKeyString, typename TValueString, typename TChar>
    inline bool
    _parseReadGffKeyValue(TFile & file, TKeyString & key, TValueString & value, TChar& c)
    {
//IOREV _nodoc_ _hasCRef_
		if (c == ' ' || c == '\t' || c == '\n' || c == '=') return false;
        appendValue(key, c);
        while (!_streamEOF(file)) 
		{
            c = _streamGet(file);
            if (c == ' ' || c == '\t' || c == '\n' || c == '=') break;
            appendValue(key, c);
        }
		_parseSkipSpace(file, c);
		if (c == '=')
		{
			c = _streamGet(file);
			_parseSkipSpace(file, c);
		}
		
		if (c == '"')
		{
			c = _streamGet(file);
			appendValue(value, c);
			while (!_streamEOF(file)) 
			{
				c = _streamGet(file);
				if (c == '\n') return true;
				if (c == '"')
				{
					if (!_streamEOF(file)) c = _streamGet(file);
					break;
				}
				appendValue(value, c);
			}
			if (c == ';')
			{
				if (!_streamEOF(file)) c = _streamGet(file);
			}
		}
		else
		{
			do {
				if (c == ' ' || c == '\t' || c == '\n') return true;
				if (c == ';')
				{
					if (!_streamEOF(file)) c = _streamGet(file);
					return true;
				}
				appendValue(value, c);
				if (_streamEOF(file)) break;
				c = _streamGet(file);
			} while (true);
		}
		return true;
	}


//////////////////////////////////////////////////////////////////////////////
// Read Gff
//////////////////////////////////////////////////////////////////////////////

template <typename TFragmentStore, typename TSpec = void>
struct IOContextGff_
{
	typedef typename TFragmentStore::TAnnotationStore   TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type      TAnnotation;
	typedef typename TAnnotation::TId                   TId;

	CharString contigName;
	CharString typeName;
	CharString annotationName;
	CharString parentKey;
	CharString parentName;
	
	CharString _key;
	CharString _value;
	StringSet<CharString> keys;
	StringSet<CharString> values;
	
	CharString gtfGeneId;
	CharString gtfGeneName;
	CharString gtfTranscriptName;		// transcipt_id is stored in parentName

	TId annotationId;
	TAnnotation annotation;
};

template <typename TFragmentStore, typename TSpec>
inline void clear(IOContextGff_<TFragmentStore, TSpec> &ctx)
{
	typedef typename TFragmentStore::TAnnotationStore   TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type      TAnnotation;

	clear(ctx.contigName);
	clear(ctx.typeName);
	clear(ctx.annotationName);
	clear(ctx.parentKey);
	clear(ctx.parentName);
	clear(ctx._key);
	clear(ctx._value);
	clear(ctx.gtfGeneId);
	clear(ctx.gtfGeneName);
	clear(ctx.gtfTranscriptName);
	clear(ctx.keys);
	clear(ctx.values);
	ctx.annotationId = TAnnotation::INVALID_ID;
	clear(ctx.annotation.values);
}

//////////////////////////////////////////////////////////////////////////////
// _readOneAnnotation
//
// reads in one annotation line from a Gff file

template <typename TFile, typename TChar, typename TFragmentStore, typename TSpec>
inline bool 
_readOneAnnotation (
	TFile & file,
	TChar & c,
	IOContextGff_<TFragmentStore, TSpec> & ctx)
{
//IOREV _nodoc_ _hasCRef_
	typedef typename TFragmentStore::TContigPos         TContigPos;	
	typedef typename TFragmentStore::TAnnotationStore   TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type      TAnnotation;
	typedef typename TAnnotation::TId                   TId;
		
	clear(ctx);

	// read fields of annotation line        
	_parseSkipWhitespace(file, c);
	
	// read column 1: contig name
	// The letters until the first whitespace will be read.
	// Then, we skip until we hit the first tab character.
	_parseReadGffIdentifier(file, ctx.contigName, c);
	if (!empty(ctx.contigName) && ctx.contigName[0] == '#')
	{
		_parseSkipLine(file, c);
		return false;
	}
	_parseSkipUntilChar(file, '\t', c);
	c = _streamGet(file);
	
	// skip column 2
	_parseSkipEntryUntilWhitespace(file, c);
	_parseSkipWhitespace(file, c);
	
	// read column 3: type
	_parseReadGffIdentifier(file, ctx.typeName, c);
	_parseSkipWhitespace(file, c);
	
	// read column 4: begin position
	if (_parseIsDigit(c))
		ctx.annotation.beginPos = _parseReadNumber(file, c) - 1;
	else
	{
		ctx.annotation.beginPos = TAnnotation::INVALID_POS;
		_parseSkipEntryUntilWhitespace(file, c);
	}
	_parseSkipWhitespace(file, c);

	// read column 5: end position
	if (_parseIsDigit(c))
		ctx.annotation.endPos = _parseReadNumber(file, c);
	else 
	{
		ctx.annotation.endPos = TAnnotation::INVALID_POS;
		_parseSkipEntryUntilWhitespace(file, c);
	}
	_parseSkipWhitespace(file, c);	

	// skip column 6
	_parseSkipEntryUntilWhitespace(file, c);
	_parseSkipWhitespace(file, c);

	// read column 7: orientation
	if (c == '-')
	{
		TContigPos tmp = ctx.annotation.beginPos;
		ctx.annotation.beginPos = ctx.annotation.endPos;
		ctx.annotation.endPos = tmp;
	}
	c = _streamGet(file);
	_parseSkipWhitespace(file, c);

	// skip column 8
	_parseSkipEntryUntilWhitespace(file, c);
	_parseSkipSpace(file, c);
	
	// read column 9: name
	while (!_streamEOF(file) &&	_parseReadGffKeyValue(file, ctx._key, ctx._value, c))
	{
		if (ctx._key == "ID") 
			ctx.annotationName = ctx._value;
		else
			if (!empty(ctx._key) && !empty(ctx._value))
			{
				appendValue(ctx.keys, ctx._key);
				appendValue(ctx.values, ctx._value);
			}
				
		if (ctx._key == "Parent" || ctx._key == "ParentID" || ctx._key == "transcript_id") 
		{
			ctx.parentKey = ctx._key;
			ctx.parentName = ctx._value;
		}
		else if (ctx._key == "transcript_name")
		{
				ctx.gtfTranscriptName = ctx._value;
		} 
		else if (ctx._key == "gene_id")
		{
			ctx.gtfGeneId = ctx._value;
		}
		else if (ctx._key == "gene_name")
		{
				ctx.gtfGeneName = ctx._value;
		} 

		clear(ctx._key);
		clear(ctx._value);
		_parseSkipSpace(file, c);
	}
	return true;
}

template <typename TAnnotation>
inline void 
_adjustParent (
	TAnnotation &parent,
	TAnnotation const &child)
{
	if (child.contigId == TAnnotation::INVALID_ID || child.beginPos == TAnnotation::INVALID_POS || child.endPos == TAnnotation::INVALID_POS)
		return;

	parent.contigId = child.contigId;	

    // Has parent an invalid begin and end position?
	if ((parent.beginPos == TAnnotation::INVALID_POS) && (parent.endPos == TAnnotation::INVALID_POS))
    {
        parent.beginPos = child.beginPos;
        parent.endPos = child.endPos;
        return;
    }

	if ((parent.beginPos == TAnnotation::INVALID_POS) || (parent.endPos == TAnnotation::INVALID_POS))
		return;
    
	typename TAnnotation::TPos childBegin, childEnd;
	if (child.beginPos < child.endPos)
	{
		childBegin = child.beginPos;
		childEnd = child.endPos;
	} else {
		childBegin = child.endPos;
		childEnd = child.beginPos;
	}

    // Keep parent's orientation and maximize begin and end using child's boundaries.
	if (parent.beginPos < parent.endPos)
	{
		if (parent.beginPos == TAnnotation::INVALID_POS || parent.beginPos > childBegin)
			parent.beginPos = childBegin;
		if (parent.endPos == TAnnotation::INVALID_POS || parent.endPos < childEnd)
			parent.endPos = childEnd;
	} else
	{
		if (parent.endPos == TAnnotation::INVALID_POS || parent.endPos > childBegin)
			parent.endPos = childBegin;
		if (parent.beginPos == TAnnotation::INVALID_POS || parent.beginPos < childEnd)
			parent.beginPos = childEnd;
	}
}

template <typename TFragmentStore, typename TSpec>
inline void 
_storeOneAnnotation (
	TFragmentStore & fragStore,
	IOContextGff_<TFragmentStore, TSpec> & ctx)
{
	typedef typename TFragmentStore::TAnnotationStore   TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type      TAnnotation;
	typedef typename TAnnotation::TId                   TId;
	
	TId maxId = 0;

	// for lines in Gtf format get/add the parent gene first
	TId geneId = TAnnotation::INVALID_ID;
	if (!empty(ctx.gtfGeneId))
	{
		_storeAppendAnnotationName(fragStore, geneId, ctx.gtfGeneId, (TId)TFragmentStore::ANNO_GENE);
		if (maxId < geneId)
			maxId = geneId;
	}	

	// if we have a parent transcript, get/add the parent transcript then
	if (!empty(ctx.parentName))
	{
// From now, we support gtf files with genes/transcripts having the same name.
//
//        // if gene and transcript names are equal (like in some strange gtf files)
//        // try to make the transcript name unique
//        if (ctx.gtfGeneId == ctx.parentName)
//            append(ctx.parentName, "_1");

        if (ctx.parentKey == "transcript_id")
        {
            // type is implicitly given (mRNA)
            _storeAppendAnnotationName(fragStore, ctx.annotation.parentId, ctx.parentName, (TId)TFragmentStore::ANNO_MRNA);
        }
        else
        {
            // type is unknown
            _storeAppendAnnotationName(fragStore, ctx.annotation.parentId, ctx.parentName);
        }
		if (maxId < ctx.annotation.parentId)
			maxId = ctx.annotation.parentId;
	}
	else
		ctx.annotation.parentId = 0;	// if we have no parent, we are a child of the root

	// add contig and type name
	_storeAppendContig(fragStore, ctx.annotation.contigId, ctx.contigName);	
	_storeAppendType(fragStore, ctx.annotation.typeId, ctx.typeName);

	// add annotation name of the current line
	_storeAppendAnnotationName(fragStore, ctx.annotationId, ctx.annotationName, ctx.annotation.typeId);
	if (maxId < ctx.annotationId)
		maxId = ctx.annotationId;
	
	for (unsigned i = 0; i < length(ctx.keys); ++i)
    {
        // don't store gene_name as key/value pair unless it is a gene
        if (ctx.keys[i] == "gene_name" && ctx.annotation.typeId != TFragmentStore::ANNO_GENE)
            continue;

        // don't store transcript_name as key/value pair unless it is a transcript
        if (ctx.keys[i] == "transcript_name" && ctx.annotation.typeId != TFragmentStore::ANNO_MRNA)
            continue;

        // don't store Parent, transcript_id or gene_id as key/value pair (the are used to link annotations)
        if (ctx.keys[i] != ctx.parentKey && ctx.keys[i] != "gene_id")
            annotationAssignValueByKey(fragStore, ctx.annotation, ctx.keys[i], ctx.values[i]);
    }

	if (length(fragStore.annotationStore) <= maxId)
		resize(fragStore.annotationStore, maxId + 1, Generous());
	fragStore.annotationStore[ctx.annotationId] = ctx.annotation;
	
	if (geneId != TAnnotation::INVALID_ID)
	{
		// link and adjust our gtf ancestors
		TAnnotation &gene = fragStore.annotationStore[geneId];
		TAnnotation &transcript = fragStore.annotationStore[ctx.annotation.parentId];

		gene.parentId = 0;
		gene.typeId = TFragmentStore::ANNO_GENE;
		_adjustParent(gene, ctx.annotation);
//		std::cout<<"gene_name "<<ctx.gtfGeneName<<"  transcript_name  " << ctx.gtfTranscriptName<<std::endl;

		if (!empty(ctx.gtfGeneName))
			annotationAssignValueByKey(fragStore, gene, "gene_name", ctx.gtfGeneName);

		transcript.parentId = geneId;
		transcript.typeId = TFragmentStore::ANNO_MRNA;
		_adjustParent(transcript, ctx.annotation);
		if (!empty(ctx.gtfTranscriptName))
			annotationAssignValueByKey(fragStore, transcript, "transcript_name", ctx.gtfTranscriptName);
	}
}

template<typename TFile, typename TSpec, typename TConfig>
inline void 
read (
	TFile & file,
	FragmentStore<TSpec, TConfig> & fragStore,
	Gff)
{
//IOREV _nodoc_
	typedef FragmentStore<TSpec, TConfig> TFragmentStore;
	
	if (_streamEOF(file)) return;

	// get first character from the stream
	char c = _streamGet(file);
	IOContextGff_<TFragmentStore> ctx;
	
	refresh(fragStore.contigNameStoreCache);
	refresh(fragStore.annotationNameStoreCache);
	refresh(fragStore.annotationTypeStoreCache);
	
	while (!_streamEOF(file))
	{
		if (_readOneAnnotation(file, c, ctx))
			_storeOneAnnotation(fragStore, ctx);
	}
	_storeClearAnnoBackLinks(fragStore.annotationStore);
	_storeCreateAnnoBackLinks(fragStore.annotationStore);
	_storeRemoveTempAnnoNames(fragStore);
}

template<typename TFile, typename TSpec, typename TConfig>
inline void 
read (
	TFile & file,
	FragmentStore<TSpec, TConfig> & fragStore,
	Gtf)
{
//IOREV _nodoc_ how do Gtf and Gff compare? nodoc for gtf
	read (file, fragStore, Gff());
}

//////////////////////////////////////////////////////////////////////////////
// Write Gff
//////////////////////////////////////////////////////////////////////////////

template<typename TTargetStream, typename TSpec, typename TConfig, typename TAnnotation, typename TId>
inline void 
_writeOneAnnotation (
	TTargetStream & target,
	FragmentStore<TSpec, TConfig> & store,
	TAnnotation &annotation,
	TId id,
	Gff)
{
//IOREV _nodoc_
	typedef FragmentStore<TSpec, TConfig>       TFragmentStore;
	typedef typename TFragmentStore::TContigPos TContigPos;
	
	if (id == 0) return;
	
	// write column 1: contig name
	if (annotation.contigId < length(store.contigNameStore))
		_streamWrite(target, store.contigNameStore[annotation.contigId]);
	_streamPut(target, '\t');
	
	// skip column 2: source
	_streamWrite(target, ".\t");
	
	// write column 3: type
	if (annotation.typeId < length(store.annotationTypeStore))
		_streamWrite(target, store.annotationTypeStore[annotation.typeId]);
	_streamPut(target, '\t');
	
	TContigPos beginPos = annotation.beginPos;
	TContigPos endPos = annotation.endPos;
	char orienation = '+';
	if (endPos < beginPos)
	{
		TContigPos tmp = beginPos;
		beginPos = endPos;
		endPos = tmp;
		orienation = '-';
	}
	
	// write column 4: begin position
	if (beginPos != TAnnotation::INVALID_POS)
		_streamPutInt(target, beginPos + 1);
	else
		_streamPut(target, '.');
	_streamPut(target, '\t');

	// write column 5: end position
	if (endPos != TAnnotation::INVALID_POS)
		_streamPutInt(target, endPos);
	else
		_streamPut(target, '.');
	_streamPut(target, '\t');

	// skip column 6: score
	_streamWrite(target, "0\t");

	// write column 7: orientation
	_streamPut(target, orienation);
	_streamPut(target, '\t');

	// skip column 8: frame
	_streamWrite(target, ".\t");
	
	// write column 9: group
	// write column 9.1: annotation id
	bool semicolon = false;
	if (id < length(store.annotationNameStore) && !empty(getAnnoName(store, id)))
	{
		_streamWrite(target, "ID=");
		_streamWrite(target, getAnnoName(store, id));
		semicolon = true;
	} 
	else if (annotation.lastChildId != TAnnotation::INVALID_ID)
	{
		_streamWrite(target, "ID=");
		_streamWrite(target, getAnnoUniqueName(store, id));
		semicolon = true;
	}

	// write column 9.2: parent id
	if (store.annotationStore[annotation.parentId].typeId > 1)	// ignore root/deleted nodes
	{
		if (semicolon) _streamPut(target, ';');
		_streamWrite(target, "Parent=");
		if (annotation.parentId < length(store.annotationNameStore) && !empty(getAnnoName(store, annotation.parentId)))
			_streamWrite(target, getAnnoName(store, annotation.parentId));
		else
			_streamWrite(target, getAnnoUniqueName(store, annotation.parentId));
		semicolon = true;
	}
	
	// write column 9.3-...: key, value pairs
	for (unsigned keyId = 0; keyId < length(annotation.values); ++keyId)
		if (!empty(annotation.values[keyId]))
		{
			if (semicolon) _streamPut(target, ';');
			_streamWrite(target, store.annotationKeyStore[keyId]);
			_streamPut(target, '=');
			_streamWrite(target, annotation.values[keyId]);
			semicolon = true;
		}
	
	_streamPut(target, '\n');	
}

template<typename TTargetStream, typename TSpec, typename TConfig, typename TAnnotation, typename TId>
inline void 
_writeOneAnnotation (
	TTargetStream & target,
	FragmentStore<TSpec, TConfig> & store,
	TAnnotation &annotation,
	TId id,
	Gtf)
{
//IOREV _nodoc_
	typedef FragmentStore<TSpec, TConfig>				TFragmentStore;
	typedef typename TFragmentStore::TContigPos			TContigPos;
	
	if (annotation.typeId <= TFragmentStore::ANNO_MRNA) return;
	
	// write column 1: contig name
	if (annotation.contigId < length(store.contigNameStore))
		_streamWrite(target, store.contigNameStore[annotation.contigId]);
	_streamPut(target, '\t');
	
	// skip column 2: source
	_streamWrite(target, ".\t");
	
	// write column 3: type
	if (annotation.typeId < length(store.annotationTypeStore))
		_streamWrite(target, store.annotationTypeStore[annotation.typeId]);
	_streamPut(target, '\t');
	
	TContigPos beginPos = annotation.beginPos;
	TContigPos endPos = annotation.endPos;
	char orienation = '+';
	if (endPos < beginPos)
	{
		TContigPos tmp = beginPos;
		beginPos = endPos;
		endPos = tmp;
		orienation = '-';
	}
	
	// write column 4: begin position
	if (beginPos != TAnnotation::INVALID_POS)
		_streamPutInt(target, beginPos + 1);
	else
		_streamPut(target, '.');
	_streamPut(target, '\t');

	// write column 5: end position
	if (endPos != TAnnotation::INVALID_POS)
		_streamPutInt(target, endPos);
	else
		_streamPut(target, '.');
	_streamPut(target, '\t');

	// skip column 6: score
	_streamWrite(target, "0\t");

	// write column 7: orientation
	_streamPut(target, orienation);
	_streamPut(target, '\t');

	// skip column 8: frame
	_streamWrite(target, ".\t");
	
	// write column 9: group
	bool semicolon = false;

    // step up until we reach a transcript
    TId transcriptId = annotation.parentId;
    while (transcriptId < length(store.annotationStore) && store.annotationStore[transcriptId].typeId != TFragmentStore::ANNO_MRNA)
        transcriptId = store.annotationStore[transcriptId].parentId;

    // step up until we reach a gene
    TId geneId = transcriptId;
    while (geneId < length(store.annotationStore) && store.annotationStore[geneId].typeId != TFragmentStore::ANNO_GENE)
        geneId = store.annotationStore[geneId].parentId;	
	
    CharString tmpStr;
    if (geneId < length(store.annotationStore) && annotationGetValueByKey(store, store.annotationStore[geneId], "gene_name", tmpStr))
    {
        if (semicolon) _streamWrite(target, "; ");
        _streamWrite(target, "gene_name \"");
        _streamWrite(target, tmpStr);
        _streamPut(target, '"');
        semicolon = true;
    }
    if (transcriptId < length(store.annotationStore) && annotationGetValueByKey(store, store.annotationStore[transcriptId], "transcript_name", tmpStr))
    {
        if (semicolon) _streamWrite(target, "; ");
        _streamWrite(target, "transcript_name \"");
        _streamWrite(target, tmpStr);
        _streamPut(target, '"');
        semicolon = true;
    }

	if (id < length(store.annotationNameStore) && !empty(getAnnoName(store, id)))
	{
		if (semicolon) _streamWrite(target, "; ");
		_streamWrite(target, "ID \"");
		_streamWrite(target, getAnnoName(store, id));
		_streamPut(target, '"');
		semicolon = true;
	} 
		
	// write key, value pairs
	for (unsigned keyId = 0; keyId < length(annotation.values); ++keyId)
		if (!empty(annotation.values[keyId]))
		{
			if (semicolon) _streamWrite(target, "; ");
			_streamWrite(target, store.annotationKeyStore[keyId]);
			_streamWrite(target, " \"");
			_streamWrite(target, annotation.values[keyId]);
			_streamPut(target, '"');
			semicolon = true;
		}
        
    // The GTF format version 2.2 requires the keys gene_id and transcript_id to be the last keys of line
    // read http://mblab.wustl.edu/GTF22.html and http://www.bioperl.org/wiki/GTF
    
    if (geneId < length(store.annotationStore))
    {
        if (semicolon) _streamWrite(target, "; ");
        _streamWrite(target, "gene_id \"");
        if (geneId < length(store.annotationNameStore) && !empty(getAnnoName(store, geneId)))
            _streamWrite(target, getAnnoName(store, geneId));
        else
            _streamWrite(target, getAnnoUniqueName(store, geneId));
		_streamPut(target, '"');
        semicolon = true;
    }

	if (transcriptId < length(store.annotationStore))
	{
		if (semicolon) _streamWrite(target, "; ");
		_streamWrite(target, "transcript_id \"");
		if (transcriptId < length(store.annotationNameStore) && !empty(getAnnoName(store, transcriptId)))
			_streamWrite(target, getAnnoName(store, transcriptId));
		else
			_streamWrite(target, getAnnoUniqueName(store, transcriptId));
		_streamPut(target, '"');
		semicolon = true;
	}

	if (semicolon) _streamWrite(target, ';');	
	_streamPut(target, '\n');	
}

template<typename TTargetStream, typename TSpec, typename TConfig, typename TFormat>
inline void 
_writeGffGtf (
	TTargetStream & target,
	FragmentStore<TSpec, TConfig> & store,
	TFormat format)
{
//IOREV _nodoc_
	typedef FragmentStore<TSpec, TConfig>							TFragmentStore;
	typedef typename TFragmentStore::TAnnotationStore				TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type					TAnnotation;
	typedef typename Iterator<TAnnotationStore, Standard>::Type		TAnnoIter;
	typedef typename Id<TAnnotation>::Type							TId;

	TAnnoIter it = begin(store.annotationStore, Standard());
	TAnnoIter itEnd = end(store.annotationStore, Standard());
	
	for(TId id = 0; it != itEnd; ++it, ++id)
		_writeOneAnnotation(target, store, *it, id, format);
}

template<typename TTargetStream, typename TSpec, typename TConfig>
inline void 
write (
	TTargetStream & target,
	FragmentStore<TSpec, TConfig> & store,
	Gff format)
{
//IOREV _nodoc_
	_writeGffGtf(target, store, format);
}

template<typename TTargetStream, typename TSpec, typename TConfig>
inline void 
write (
	TTargetStream & target,
	FragmentStore<TSpec, TConfig> & store,
	Gtf format)
{
//IOREV _nodoc_
	_writeGffGtf(target, store, format);
}


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
