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

#ifndef SEQAN_HEADER_STORE_IO_UCSC_H
#define SEQAN_HEADER_STORE_IO_UCSC_H

/* IOREV
 *
 * _doc_
 *
 *
 * maybe move this to file/ because its a file format
 *
 */


namespace SEQAN_NAMESPACE_MAIN
{

template <typename TSpec>
struct Ucsc_;

/**
.Tag.File Format.tag.Ucsc:
	Ucsc Genome Browser annotation file (a.k.a. knownGene format).
..include:seqan/store.h
*/

struct UcscKnownGene_;
typedef Tag<Ucsc_<UcscKnownGene_> > const Ucsc;

/**
.Tag.File Format.tag.UcscIsoforms:
	Ucsc Genome Browser isoforms file (a.k.a. knownIsoforms format).
..include:seqan/store.h
*/
struct UcscKnownIsoforms_;
typedef Tag<Ucsc_<UcscKnownIsoforms_> > const UcscIsoforms;
	
//////////////////////////////////////////////////////////////////////////////
// _parseReadUcscIdentifier
    
    template<typename TFile, typename TString, typename TChar>
    inline void
    _parseReadUcscIdentifier(TFile & file, TString & str, TChar& c)
    {
//IOREV
        if (c == ' ' || c == '\t' || c == '\n') return;
        appendValue(str, c);
        while (!_streamEOF(file)) 
		{
            c = _streamGet(file);
            if (c == ' ' || c == '\t' || c == '\n') return;
            appendValue(str, c);
        }
    }

	template<typename TFile, typename TChar>
	inline void 
	_parseSkipWhiteComma(TFile& file, TChar& c)
	{
//IOREV
		if (c != ',' && c != ' ') return;
		while (!_streamEOF(file)) {
			c = _streamGet(file);
			if (c != ',' && c != ' ') break;
		}
	}
	
//////////////////////////////////////////////////////////////////////////////
// Read Ucsc
//////////////////////////////////////////////////////////////////////////////

template <typename TFragmentStore, typename TSpec = void>
struct IOContextUcsc_
{
	typedef typename TFragmentStore::TAnnotationStore   TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type      TAnnotation;
	typedef typename TAnnotation::TId                   TId;

	CharString		transName;
	CharString		contigName;
	__int64			cdsBegin;
	__int64			cdsEnd;
	String<__int64>	exonBegin;
	String<__int64>	exonEnd;
	CharString		proteinName;
	
	enum { KNOWN_GENE, KNOWN_ISOFORMS } format;
	TAnnotation annotation;
};

template <typename TFragmentStore, typename TSpec>
inline void clear(IOContextUcsc_<TFragmentStore, TSpec> &ctx)
{
	typedef typename TFragmentStore::TAnnotationStore   TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type      TAnnotation;

	clear(ctx.transName);
	clear(ctx.contigName);
	clear(ctx.exonBegin);
	clear(ctx.exonEnd);
	clear(ctx.proteinName);
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
	IOContextUcsc_<TFragmentStore, TSpec> & ctx)
{
//IOREV
	typedef typename TFragmentStore::TContigPos         TContigPos;	
	typedef typename TFragmentStore::TAnnotationStore   TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type      TAnnotation;
	typedef typename TAnnotation::TId                   TId;
		
	clear(ctx);

	// read fields of alignments line        
	_parseSkipWhitespace(file, c);
	
	// read column 1: transcript name
	// The letters until the first whitespace will be read.
	// Then, we skip until we hit the first tab character.
	_parseReadUcscIdentifier(file, ctx.transName, c);
	if (!empty(ctx.transName) && ctx.transName[0] == '#')
	{
		_parseSkipLine(file, c);
		return false;
	}
	_parseSkipUntilChar(file, '\t', c);
	c = _streamGet(file);
	
	// read column 2: contig name
	_parseReadUcscIdentifier(file, ctx.contigName, c);
	_parseSkipSpace(file, c);
	
	// read column 3: orientation
	if (c != '+' && c != '-')
	{
		ctx.format = ctx.KNOWN_ISOFORMS;
		assign(prefix(ctx.transName, 0), "GENE");
		_parseSkipLine(file, c);
		return true;
	}
	ctx.format = ctx.KNOWN_GENE;
	char orientation = c;
	c = _streamGet(file);
	_parseSkipWhitespace(file, c);

	// read column 4: transcript begin position
	if (_parseIsDigit(c))
		ctx.annotation.beginPos = _parseReadNumber(file, c);
	else
	{
		ctx.annotation.beginPos = TAnnotation::INVALID_POS;
		_parseSkipEntryUntilWhitespace(file, c);
	}
	_parseSkipWhitespace(file, c);

	// read column 5: transcript end position
	if (_parseIsDigit(c))
		ctx.annotation.endPos = _parseReadNumber(file, c);
	else 
	{
		ctx.annotation.endPos = TAnnotation::INVALID_POS;
		_parseSkipEntryUntilWhitespace(file, c);
	}
	_parseSkipWhitespace(file, c);	

	// read column 6: CDS begin position
	if (_parseIsDigit(c))
		ctx.cdsBegin = _parseReadNumber(file, c);
	else
	{
		ctx.cdsBegin = TAnnotation::INVALID_POS;
		_parseSkipEntryUntilWhitespace(file, c);
	}
	_parseSkipWhitespace(file, c);
	
	// read column 7: CDS end position
	if (_parseIsDigit(c))
		ctx.cdsEnd = _parseReadNumber(file, c);
	else 
	{
		ctx.cdsEnd = TAnnotation::INVALID_POS;
		_parseSkipEntryUntilWhitespace(file, c);
	}
	_parseSkipWhitespace(file, c);	
	
	// read column 8: exon count
	int exons = -1;
	if (_parseIsDigit(c))
		exons = _parseReadNumber(file, c);
	_parseSkipWhitespace(file, c);

	// read column 9: exon begin positions
	for (int i = 0; i < exons; ++i)
	{
		appendValue(ctx.exonBegin, _parseReadNumber(file, c), Generous());
		_parseSkipWhiteComma(file, c);
	}
	_parseSkipUntilChar(file, '\t', c);
	c = _streamGet(file);
	
	// read column 10: exon end positions
	for (int i = 0; i < exons; ++i)
	{
		appendValue(ctx.exonEnd, _parseReadNumber(file, c));
		_parseSkipWhiteComma(file, c);
	}
	_parseSkipUntilChar(file, '\t', c);
	c = _streamGet(file);
	
	// read column 10: protein name
	_parseReadUcscIdentifier(file, ctx.proteinName, c);
	_parseSkipUntilChar(file, '\t', c);
	c = _streamGet(file);

	// skip column 11
	_parseSkipEntryUntilWhitespace(file, c);
	_parseSkipWhitespace(file, c);

	// adapt positions
	if (orientation == '-')
	{
		TContigPos tmp = ctx.annotation.beginPos;
		ctx.annotation.beginPos = ctx.annotation.endPos;
		ctx.annotation.endPos = tmp;
		tmp = ctx.cdsBegin;
		ctx.cdsBegin = ctx.cdsEnd;
		ctx.cdsEnd = tmp;
		for (int i = 0; i < exons; ++i)
		{
			tmp = ctx.exonBegin[i];
			ctx.exonBegin[i] = ctx.exonEnd[i];
			ctx.exonEnd[i] = tmp;
		}
	}

	return true;
}

template <typename TFragmentStore, typename TSpec>
inline void 
_storeOneAnnotationKnownGene (
	TFragmentStore & fragStore,
	IOContextUcsc_<TFragmentStore, TSpec> & ctx)
{
	typedef typename TFragmentStore::TAnnotationStore   TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type      TAnnotation;
	typedef typename TAnnotation::TId                   TId;
	
	SEQAN_ASSERT_EQ(length(fragStore.annotationStore), length(fragStore.annotationNameStore));

	TId annoStoreLen = length(fragStore.annotationStore);
	TId transId = TAnnotation::INVALID_ID;
	TId cdsId = TAnnotation::INVALID_ID;

	// add transcript and CDS
	_storeAppendAnnotationName(fragStore, transId, ctx.transName);
	cdsId = length(fragStore.annotationNameStore);
	appendName(fragStore.annotationNameStore, ctx.proteinName, fragStore.annotationNameStoreCache);
	
	if (annoStoreLen <= transId)
		annoStoreLen = transId + 1;
	
	if (annoStoreLen <= cdsId)
		annoStoreLen = cdsId + 1;
	
	resize(fragStore.annotationStore, annoStoreLen + length(ctx.exonBegin), Generous());
	resize(fragStore.annotationNameStore, annoStoreLen + length(ctx.exonBegin), Generous());

	// add contig name
	_storeAppendContig(fragStore, ctx.annotation.contigId, ctx.contigName);	

	TAnnotation &transcript = fragStore.annotationStore[transId];
	TId geneId = transcript.parentId;
	if (geneId == TAnnotation::INVALID_ID) geneId = 0;
	transcript = ctx.annotation;
	transcript.parentId = geneId;
	transcript.typeId = TFragmentStore::ANNO_MRNA;
	
	TAnnotation &cds = fragStore.annotationStore[cdsId];
	cds = ctx.annotation;
	cds.parentId = transId;
	cds.typeId = TFragmentStore::ANNO_CDS;
	cds.beginPos = ctx.cdsBegin;
	cds.endPos = ctx.cdsEnd;
	_adjustParent(transcript, cds);
	
	// insert exons
	ctx.annotation.parentId = transId;
	ctx.annotation.typeId = TFragmentStore::ANNO_EXON;
	for (unsigned i = 0; i < length(ctx.exonBegin); ++i)
	{
		ctx.annotation.beginPos = ctx.exonBegin[i];
		ctx.annotation.endPos = ctx.exonEnd[i];
		fragStore.annotationStore[annoStoreLen + i] = ctx.annotation;
		_adjustParent(transcript, ctx.annotation);
	}
	if (geneId != 0)
		_adjustParent(fragStore.annotationStore[geneId], transcript);
}

template <typename TFragmentStore, typename TSpec>
inline void 
_storeOneAnnotationKnownIsoforms (
	TFragmentStore & fragStore,
	IOContextUcsc_<TFragmentStore, TSpec> & ctx)
{
	typedef typename TFragmentStore::TAnnotationStore   TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type      TAnnotation;
	typedef typename TAnnotation::TId                   TId;
	
	SEQAN_ASSERT_EQ(length(fragStore.annotationStore), length(fragStore.annotationNameStore));
	
	TId annoStoreLen = length(fragStore.annotationStore);
	TId geneId = TAnnotation::INVALID_ID;
	TId transId = TAnnotation::INVALID_ID;
	
	// add transcript and CDS
	_storeAppendAnnotationName(fragStore, geneId, ctx.transName);
	_storeAppendAnnotationName(fragStore, transId, ctx.contigName);
	
	if (annoStoreLen <= geneId)
		annoStoreLen = geneId + 1;
	
	if (annoStoreLen <= transId)
		annoStoreLen = transId + 1;
	
	resize(fragStore.annotationStore, annoStoreLen, Generous());
	resize(fragStore.annotationNameStore, annoStoreLen, Generous());
	
	// set parent link locus->root
	TAnnotation &locus = fragStore.annotationStore[geneId];
	locus.parentId = 0;
	locus.typeId = TFragmentStore::ANNO_GENE;

	// set parent link transcript->locus
	TAnnotation &transcript = fragStore.annotationStore[transId];
	transcript.parentId = geneId;
	transcript.typeId = TFragmentStore::ANNO_MRNA;
	
	_adjustParent(locus, transcript);
}

template <typename TFragmentStore, typename TSpec>
inline void 
_storeOneAnnotation (
	TFragmentStore & fragStore,
	IOContextUcsc_<TFragmentStore, TSpec> & ctx)
{
	if (ctx.format == ctx.KNOWN_GENE)
		_storeOneAnnotationKnownGene(fragStore, ctx);
	else
		_storeOneAnnotationKnownIsoforms(fragStore, ctx);
}

template<typename TFile, typename TSpec, typename TConfig, typename TFormatSpec>
inline void 
read (
	TFile & file,
	FragmentStore<TSpec, TConfig> & fragStore,
	Tag<Ucsc_<TFormatSpec> > const)
{
//IOREV _nodoc_
	typedef FragmentStore<TSpec, TConfig> TFragmentStore;
	
	if (_streamEOF(file)) return;

	// get first character from the stream
	char c = _streamGet(file);
	IOContextUcsc_<TFragmentStore> ctx;
	
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

//////////////////////////////////////////////////////////////////////////////
// Write Ucsc
//////////////////////////////////////////////////////////////////////////////

template <typename TFragmentStore, typename TSpec, typename TAnnotation, typename TId>
inline bool 
_retrieveOneAnnotation (
	TFragmentStore & fragStore,
	IOContextUcsc_<TFragmentStore, TSpec> & ctx,
	TAnnotation &annotation,
	TId id,
	Ucsc)
{	
	if (annotation.typeId != TFragmentStore::ANNO_MRNA) return false;
	
	ctx.format = ctx.KNOWN_GENE;
	ctx.transName = getAnnoUniqueName(fragStore, id);
	if (annotation.contigId < length(fragStore.contigNameStore))
		ctx.contigName = fragStore.contigNameStore[annotation.contigId];
	else
		clear(ctx.contigName);
	
	ctx.annotation = annotation;
	clear(ctx.proteinName);
	clear(ctx.exonBegin);
	clear(ctx.exonEnd);
	
	TId lastChildId = annotation.lastChildId;
	TId i = lastChildId;
	do {
		i = fragStore.annotationStore[i].nextSiblingId;
		TAnnotation &anno = fragStore.annotationStore[i];
		if (anno.typeId == TFragmentStore::ANNO_CDS)
		{
			if (i < length(fragStore.annotationNameStore))
				ctx.proteinName = fragStore.annotationNameStore[i];
			ctx.cdsBegin = anno.beginPos;
			ctx.cdsEnd = anno.endPos;
		}
		if (anno.typeId == TFragmentStore::ANNO_EXON)
		{
			appendValue(ctx.exonBegin, anno.beginPos, Generous());
			appendValue(ctx.exonEnd, anno.endPos, Generous());
		}
	} while (i != lastChildId);
	return true;
}

template <typename TFragmentStore, typename TSpec, typename TAnnotation, typename TId>
inline bool 
_retrieveOneAnnotation (
	TFragmentStore & fragStore,
	IOContextUcsc_<TFragmentStore, TSpec> & ctx,
	TAnnotation &annotation,
	TId id,
	UcscIsoforms)
{	
	if (annotation.typeId != TFragmentStore::ANNO_MRNA) return false;
	if (annotation.parentId == TAnnotation::INVALID_ID || annotation.parentId == 0) return false;
	
	ctx.format = ctx.KNOWN_ISOFORMS;
	ctx.transName = getAnnoUniqueName(fragStore, annotation.parentId);
	ctx.contigName = getAnnoUniqueName(fragStore, id);
	return true;
}

template<typename TTargetStream, typename TFragmentStore, typename TSpec>
inline void 
_writeOneAnnotation (
	TTargetStream & file,
	IOContextUcsc_<TFragmentStore, TSpec> & ctx)
{
//IOREV
	typedef typename TFragmentStore::TContigPos         TContigPos;	
	
	unsigned suf = 0;
	if (ctx.format == ctx.KNOWN_ISOFORMS && length(ctx.transName) >= 4 && prefix(ctx.transName, 4) == "GENE")
		suf = 4;
	
	// read column 1: transcript name
	// The letters until the first whitespace will be read.
	// Then, we skip until we hit the first tab character.
	_streamWrite(file, suffix(ctx.transName, suf));
	_streamPut(file, '\t');
	
	// read column 2: contig name
	_streamWrite(file, ctx.contigName);
	if (ctx.format == ctx.KNOWN_ISOFORMS)
	{
		_streamPut(file, '\n');
		return;
	}
	_streamPut(file, '\t');
	
	// read column 3: orientation
	TContigPos transBeginPos, transEndPos;
	TContigPos cdsBeginPos, cdsEndPos;
	if (ctx.annotation.beginPos < ctx.annotation.endPos)
	{
		_streamPut(file, '+');
		transBeginPos = ctx.annotation.beginPos;
		transEndPos = ctx.annotation.endPos;
		cdsBeginPos = ctx.cdsBegin;
		cdsEndPos = ctx.cdsEnd;
	}
	else
	{
		_streamPut(file, '-');
		transEndPos = ctx.annotation.beginPos;
		transBeginPos = ctx.annotation.endPos;
		cdsEndPos = ctx.cdsBegin;
		cdsBeginPos = ctx.cdsEnd;
	}
	_streamPut(file, '\t');

	// read column 4: transcript begin position
	_streamPutInt(file, transBeginPos);
	_streamPut(file, '\t');

	// read column 5: transcript end position
	_streamPutInt(file, transEndPos);
	_streamPut(file, '\t');

	// read column 6: CDS begin position
	_streamPutInt(file, cdsBeginPos);
	_streamPut(file, '\t');
	
	// read column 7: CDS end position
	_streamPutInt(file, cdsEndPos);
	_streamPut(file, '\t');
	
	// read column 8: exon count
	_streamPutInt(file, length(ctx.exonBegin));
	_streamPut(file, '\t');

	// read column 9: exon begin positions
	for (unsigned i = 0; i < length(ctx.exonBegin); ++i)
	{
		_streamPutInt(file, _min(ctx.exonBegin[i], ctx.exonEnd[i]));
		_streamPut(file, ',');
	}
	_streamPut(file, '\t');
	
	// read column 10: exon end positions
	for (unsigned i = 0; i < length(ctx.exonBegin); ++i)
	{
		_streamPutInt(file, _max(ctx.exonBegin[i], ctx.exonEnd[i]));
		_streamPut(file, ',');
	}
	_streamPut(file, '\t');
	
	// read column 10: protein name
	_streamWrite(file, ctx.proteinName);
	_streamPut(file, '\t');

	// skip column 11
	_streamWrite(file, ctx.transName);
	_streamPut(file, '\n');
}

template<typename TTargetStream, typename TSpec, typename TConfig, typename TFormatSpec>
inline void 
write (
	TTargetStream & target,
	FragmentStore<TSpec, TConfig> & store,
	Tag<Ucsc_<TFormatSpec> > const format)
{
//IOREV _nodoc_
	typedef FragmentStore<TSpec, TConfig>							TFragmentStore;
	typedef typename TFragmentStore::TAnnotationStore				TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type					TAnnotation;
	typedef typename Iterator<TAnnotationStore, Standard>::Type		TAnnoIter;
	typedef typename Id<TAnnotation>::Type							TId;

	IOContextUcsc_<TFragmentStore> ctx;

	TAnnoIter it = begin(store.annotationStore, Standard());
	TAnnoIter itEnd = end(store.annotationStore, Standard());
	
	for(TId id = 0; it != itEnd; ++it, ++id)
	{
		if (_retrieveOneAnnotation(store, ctx, *it, id, format))
			_writeOneAnnotation(target, ctx);
	}
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
