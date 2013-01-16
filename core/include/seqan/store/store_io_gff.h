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

namespace SEQAN_NAMESPACE_MAIN {

/**
.Tag.File Format.tag.Gff:
    Gff annotation file.
..include:seqan/store.h
*/
struct TagGff_;
typedef Tag<TagGff_> const Gff;

/**
// .Tag.File Format.tag.Gtf:
    Gtf annotation file.
..include:seqan/store.h
*/
struct TagGtf_;
typedef Tag<TagGtf_> const Gtf;

// Returns an error code, i.e. == 0 for OK, != 0 for error.
template <typename TReader, typename TKeyString, typename TValueString>
inline int
_parseReadGffKeyValue(TReader & reader, TKeyString & key, TValueString & outValue)
{
    char c = value(reader);
    if (c == ' ' || c == '\t' || c == '\n' || c == '=')
        return 1;

    for (; !atEnd(reader); goNext(reader))
    {
        c = value(reader);
        if (c == ' ' || c == '\t' || c == '\n' || c == '=')
            break;
        appendValue(key, c);
    }
    if (skipWhitespaces(reader) != 0)
        return 1;

    if (value(reader) == '=')
    {
        goNext(reader);
        if (skipWhitespaces(reader) != 0)
            return 1;

        if (atEnd(reader))
            return 1;
    }

    // Handle the case of a string literal.
    if (value(reader) == '"')
    {
        goNext(reader);
        // Append all characters in the literal to outValue until the first line
        // break or the closing '"'.
        for (; !atEnd(reader); goNext(reader))
        {
            if (value(reader) == '\n')
                return 1;

            if (value(reader) == '"')
            {
                goNext(reader);
                break;
            }
            appendValue(outValue, value(reader));
        }
        // Go over the trailing simicolon and any trailing space.
        while (!atEnd(reader) && (value(reader) == ';' || value(reader) == ' '))
            goNext(reader);
    }
    else
    {
        // Read until the first semicolon, return at whitespace.
        for (; !atEnd(reader); goNext(reader))
        {
            if (isspace(value(reader)))
                return 0;

            if (value(reader) == ';')
                break;
            appendValue(outValue, value(reader));
        }
        // Skip semicolon and spaces if any.
        while (!atEnd(reader) && (value(reader) == ';' || value(reader) == ' '))
            goNext(reader);
    }
    return 0;
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
    CharString gtfTranscriptName;       // transcipt_id is stored in parentName

    TId annotationId;
    TAnnotation annotation;
};

template <typename TFragmentStore, typename TSpec>
inline void clear(IOContextGff_<TFragmentStore, TSpec> & ctx)
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

template <typename TRecordReader, typename TFragmentStore, typename TSpec>
inline bool
_readOneAnnotation(
    TRecordReader & reader,
    IOContextGff_<TFragmentStore, TSpec> & ctx)
{
//IOREV _nodoc_ _hasCRef_
    typedef typename TFragmentStore::TContigPos         TContigPos;
    typedef typename TFragmentStore::TAnnotationStore   TAnnotationStore;
    typedef typename Value<TAnnotationStore>::Type      TAnnotation;
    typedef typename TAnnotation::TId                   TId;

    clear(ctx);

    // read column 1: contig name
    // The letters until the first whitespace will be read.
    // Then, we skip until we hit the first tab character.
    if (readUntilWhitespace(ctx.contigName, reader))
        return false;

    if (!empty(ctx.contigName) && ctx.contigName[0] == '#')
    {
        if (skipLine(reader))
            return false;

        return false;
    }
    if (skipWhitespaces(reader))
        return false;

    // skip column 2
    if (skipUntilWhitespace(reader) || skipBlanks(reader))
        return false;

    // read column 3: type
    if (readUntilWhitespace(ctx.typeName, reader))
        return false;

    if (skipWhitespaces(reader))
        return false;

    // read column 4: begin position
    String<char> temp;
    if (readDigits(temp, reader))
        return false;

    if (length(temp) > 0u)
    {
        if (!lexicalCast2(ctx.annotation.beginPos, temp))
            return false;

        --ctx.annotation.beginPos;
    }
    else
    {
        ctx.annotation.beginPos = TAnnotation::INVALID_POS;
        if (skipUntilWhitespace(reader))
            return false;
    }
    if (skipBlanks(reader))
        return false;

    // read column 5: end position
    clear(temp);
    if (readDigits(temp, reader))
        return false;

    if (length(temp) > 0u)
    {
        if (!lexicalCast2(ctx.annotation.endPos, temp))
            return false;
    }
    else
    {
        ctx.annotation.endPos = TAnnotation::INVALID_POS;
        if (skipUntilWhitespace(reader))
            return false;
    }
    if (skipBlanks(reader))
        return false;


    // skip column 6
    if (skipUntilWhitespace(reader) || skipBlanks(reader))
        return false;

    // read column 7: orientation
    clear(temp);
    if (readUntilWhitespace(temp, reader))
        return false;

    if (temp == "-")
    {
        TContigPos tmp = ctx.annotation.beginPos;
        ctx.annotation.beginPos = ctx.annotation.endPos;
        ctx.annotation.endPos = tmp;
    }
    if (skipBlanks(reader))
        return false;

    // skip column 8
    if (skipUntilWhitespace(reader) || skipBlanks(reader))
        return false;

    // read column 9: name
    while (!atEnd(reader))
    {
        // Read next key/value pair.
        if (_parseReadGffKeyValue(reader, ctx._key, ctx._value) != 0)
            return false;

        if (ctx._key == "ID")
        {
            ctx.annotationName = ctx._value;
        }
        else if (!empty(ctx._key) && !empty(ctx._value))
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

        // At end of line:  Skip EOL and break.
        if (!atEnd(reader) && (value(reader) == '\r' || value(reader) == '\n'))
        {
            if (skipLine(reader) != 0)
                return false;

            break;
        }
    }
    return true;
}

template <typename TAnnotation>
inline void
_adjustParent(
    TAnnotation & parent,
    TAnnotation const & child)
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
    }
    else
    {
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
    }
    else
    {
        if (parent.endPos == TAnnotation::INVALID_POS || parent.endPos > childBegin)
            parent.endPos = childBegin;
        if (parent.beginPos == TAnnotation::INVALID_POS || parent.beginPos < childEnd)
            parent.beginPos = childEnd;
    }
}

template <typename TFragmentStore, typename TSpec>
inline void
_storeOneAnnotation(
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
        _storeAppendAnnotationName(fragStore, geneId, ctx.gtfGeneId, (TId) TFragmentStore::ANNO_GENE);
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
            _storeAppendAnnotationName(fragStore, ctx.annotation.parentId, ctx.parentName, (TId) TFragmentStore::ANNO_MRNA);
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
        ctx.annotation.parentId = 0;    // if we have no parent, we are a child of the root

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
        TAnnotation & gene = fragStore.annotationStore[geneId];
        TAnnotation & transcript = fragStore.annotationStore[ctx.annotation.parentId];

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

template <typename TFile, typename TSpec, typename TConfig>
inline void
read(
    TFile & file,
    FragmentStore<TSpec, TConfig> & fragStore,
    Gff)
{
//IOREV _nodoc_
    typedef FragmentStore<TSpec, TConfig> TFragmentStore;

    if (streamEof(file))
        return;

    IOContextGff_<TFragmentStore> ctx;

    refresh(fragStore.contigNameStoreCache);
    refresh(fragStore.annotationNameStoreCache);
    refresh(fragStore.annotationTypeStoreCache);

    RecordReader<TFile, SinglePass<> > reader(file);

    while (!atEnd(reader))
    {
        if (_readOneAnnotation(reader, ctx))
            _storeOneAnnotation(fragStore, ctx);
    }
    _storeClearAnnoBackLinks(fragStore.annotationStore);
    _storeCreateAnnoBackLinks(fragStore.annotationStore);
    _storeRemoveTempAnnoNames(fragStore);
}

template <typename TFile, typename TSpec, typename TConfig>
inline void
read(
    TFile & file,
    FragmentStore<TSpec, TConfig> & fragStore,
    Gtf)
{
//IOREV _nodoc_ how do Gtf and Gff compare? nodoc for gtf
    read(file, fragStore, Gff());
}

//////////////////////////////////////////////////////////////////////////////
// Write Gff
//////////////////////////////////////////////////////////////////////////////

// This function checks if the string to be written contains a semicolon. If
// this is the case parenthesis are written around the string.
// Returns false on success.
template <typename TTargetStream, typename TString>
inline bool
_writeSemicolonSensitive(TTargetStream & target, TString & temp)
{
    if (std::find(begin(temp), end(temp), ';') != end(temp))
    {
        if (streamWriteChar(target, '"') ||
            streamWriteBlock(target, &temp[0], length(temp)) < length(temp) ||
            streamWriteChar(target, '"'))
            return true;
    }
    else
    {
        if (streamWriteBlock(target, &temp[0], length(temp)) < length(temp))
            return true;
    }
    return false;
}

// This function write the information that are equal for gff and gtf files.
template <typename TTargetStream, typename TSpec, typename TConfig, typename TAnnotation, typename TId>
inline bool
_writeCommonGffGtfInfo(
    TTargetStream & target,
    FragmentStore<TSpec, TConfig> & store,
    TAnnotation & annotation,
    TId /*id*/)
{

    typedef FragmentStore<TSpec, TConfig>       TFragmentStore;
    typedef typename TFragmentStore::TContigPos TContigPos;

    // write column 1: contig name
    if (annotation.contigId < length(store.contigNameStore))
    {
        if (length(store.contigNameStore[annotation.contigId]) > 0u)
        {
            if (streamWriteBlock(target, &(store.contigNameStore[annotation.contigId])[0], length(store.contigNameStore[annotation.contigId])) < length(store.contigNameStore[annotation.contigId]))
                return false;
        }
    }
    if (streamWriteChar(target, '\t'))
        return false;

    // skip column 2: source
    if (streamWriteBlock(target, ".\t", 2) < 2u)
        return false;

    // write column 3: type
    if (annotation.typeId < length(store.annotationTypeStore))
    {
        if (length(store.annotationTypeStore[annotation.typeId]) > 0u)
        {
            if (streamWriteBlock(target, &(store.annotationTypeStore[annotation.typeId])[0], length(store.annotationTypeStore[annotation.typeId])) < length(store.annotationTypeStore[annotation.typeId]))
                return false;
        }

    }
    if (streamWriteChar(target, '\t'))
        return false;

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
    {
        if (streamPut(target, beginPos + 1))
            return false;
    }
    else
    {
        if (streamWriteChar(target, '.'))
            return false;
    }
    if (streamWriteChar(target, '\t'))
        return false;


    // write column 5: end position
    if (endPos != TAnnotation::INVALID_POS)
    {
        if (streamPut(target, endPos))
            return false;
    }
    else
    {
        if (streamWriteChar(target, '.'))
            return false;
    }
    if (streamWriteChar(target, '\t'))
        return false;

    // skip column 6: score
    if (streamWriteBlock(target, ".\t", 2) < 2u)
        return false;

    // write column 7: orientation
    if (streamWriteChar(target, orienation))
        return false;

    if (streamWriteChar(target, '\t'))
        return false;

    // skip column 8: frame
    if (streamWriteBlock(target, ".\t", 2) < 2u)
        return false;

    return true;

}

template <typename TTargetStream, typename TSpec, typename TConfig, typename TAnnotation, typename TId>
inline bool
_writeOneAnnotation(
    TTargetStream & target,
    FragmentStore<TSpec, TConfig> & store,
    TAnnotation & annotation,
    TId id,
    Gff)
{

    if (id == 0)
        return false;

    if (!_writeCommonGffGtfInfo(target, store, annotation, id))
        return false;

    // write column 9: group
    // write column 9.1: annotation id
    bool semicolon = false;
    String<char> temp;
    if (id < length(store.annotationNameStore) && !empty(getAnnoName(store, id)))
    {
        temp = getAnnoName(store, id);
    }
    else if (annotation.lastChildId != TAnnotation::INVALID_ID)
    {
        temp = getAnnoUniqueName(store, id);
    }

    if (length(temp) > 0)
    {
        if (streamWriteBlock(target, "ID=", 3u) < 3u)
            return false;

        if (_writeSemicolonSensitive(target, temp))
            return false;

        semicolon = true;
    }

    // write column 9.2: parent id
    if (store.annotationStore[annotation.parentId].typeId > 1)  // ignore root/deleted nodes
    {
        if (semicolon)
            if (streamWriteChar(target, ';'))
                return false;

        if (streamWriteBlock(target, "Parent=", 7u) < 7u)
            return false;

        String<char> temp = getAnnoUniqueName(store, annotation.parentId);
        if (_writeSemicolonSensitive(target, temp))
            return false;

        semicolon = true;
    }

    // write column 9.3-...: key, value pairs
    for (unsigned keyId = 0; keyId < length(annotation.values); ++keyId)
        if (!empty(annotation.values[keyId]))
        {
            if (semicolon)
                if (streamWriteChar(target, ';'))
                    return false;

            String<char> temp = store.annotationKeyStore[keyId];
            if (_writeSemicolonSensitive(target, temp))
                return false;

            if (streamWriteChar(target, '='))
                return false;

            temp = annotation.values[keyId];
            if (_writeSemicolonSensitive(target, temp))
                return false;

            semicolon = true;
        }

    if (streamWriteChar(target, '\n'))
        return false;

    return true;
}

template <typename TTargetStream, typename TSpec, typename TConfig, typename TAnnotation, typename TId>
inline bool
_writeOneAnnotation(
    TTargetStream & target,
    FragmentStore<TSpec, TConfig> & store,
    TAnnotation & annotation,
    TId id,
    Gtf)
{
    typedef FragmentStore<TSpec, TConfig>               TFragmentStore;

    if (annotation.typeId <= TFragmentStore::ANNO_MRNA)
        return false;

    if (!_writeCommonGffGtfInfo(target, store, annotation, id))
        return false;

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
        if (semicolon)
            if (streamWriteBlock(target, "; ", 2) < 2u)
                return false;

        if (streamWriteBlock(target, "gene_name \"", 11u) < 11u)
            return false;

        if (streamWriteBlock(target, &tmpStr[0], length(tmpStr)) < length(tmpStr))
            return false;

        if (streamWriteChar(target, '"'))
            return false;

        semicolon = true;
    }
    if (transcriptId < length(store.annotationStore) && annotationGetValueByKey(store, store.annotationStore[transcriptId], "transcript_name", tmpStr))
    {
        if (semicolon)
            if (streamWriteBlock(target, "; ", 2) < 2u)
                return false;

        if (streamWriteBlock(target, "transcript_name \"", 17u) < 11u)
            return false;

        if (streamWriteBlock(target, &tmpStr[0], length(tmpStr)) < length(tmpStr))
            return false;

        if (streamWriteChar(target, '"'))
            return false;

        semicolon = true;
    }

    if (id < length(store.annotationNameStore) && !empty(getAnnoName(store, id)))
    {
        if (semicolon)
            if (streamWriteBlock(target, "; ", 2) < 2u)
                return false;

        if (streamWriteBlock(target, "ID \"", 3u) < 11u)
            return false;

        String<char> temp = getAnnoName(store, id);
        if (streamWriteBlock(target, &temp[0], length(temp)) < length(temp))
            return false;

        if (streamWriteChar(target, '"'))
            return false;

        semicolon = true;
    }

    // write key, value pairs
    for (unsigned keyId = 0; keyId < length(annotation.values); ++keyId)
        if (!empty(annotation.values[keyId]))
        {
            if (semicolon)
                if (streamWriteBlock(target, "; ", 2) < 2u)
                    return false;

            if (streamWriteBlock(target, &(store.annotationKeyStore[keyId])[0], length(store.annotationKeyStore[keyId])) < length(store.annotationKeyStore[keyId]))
                return false;

            if (streamWriteBlock(target, " \"", 2u) < 2u)
                return false;

            if (streamWriteBlock(target, &(annotation.values[keyId])[0], length(annotation.values[keyId])) < length(annotation.values[keyId]))
                return false;

            if (streamWriteChar(target, '"'))
                return false;

            semicolon = true;
        }

    // The GTF format version 2.2 requires the keys gene_id and transcript_id to be the last keys of line
    // read http://mblab.wustl.edu/GTF22.html and http://www.bioperl.org/wiki/GTF

    if (geneId < length(store.annotationStore))
    {
        if (semicolon)
            if (streamWriteBlock(target, "; ", 2) < 2u)
                return false;

        if (streamWriteBlock(target, "gene_id \"", 9u) < 9u)
            return false;

        String<char> temp = getAnnoUniqueName(store, geneId);
        if (streamWriteBlock(target, &temp[0], length(temp)) < length(temp))
            return false;

        if (streamWriteChar(target, '"'))
            return false;

        semicolon = true;
    }

    if (transcriptId < length(store.annotationStore))
    {
        if (semicolon)
            if (streamWriteBlock(target, "; ", 2) < 2u)
                return false;

        if (streamWriteBlock(target, "transcript_id \"", 15u) < 15u)
            return false;

        String<char> temp = getAnnoUniqueName(store, transcriptId);
        if (streamWriteBlock(target, &temp[0], length(temp)) < length(temp))
            return false;

        if (streamWriteChar(target, '"'))
            return false;

        semicolon = true;
    }

    if (semicolon)
        if (streamWriteChar(target, ';'))
            return false;

    if (streamWriteChar(target, '\n'))
        return false;

    return true;
}

template <typename TTargetStream, typename TSpec, typename TConfig, typename TFormat>
inline void
_writeGffGtf(
    TTargetStream & target,
    FragmentStore<TSpec, TConfig> & store,
    TFormat format)
{
    typedef FragmentStore<TSpec, TConfig>                           TFragmentStore;
    typedef typename TFragmentStore::TAnnotationStore               TAnnotationStore;
    typedef typename Value<TAnnotationStore>::Type                  TAnnotation;
    typedef typename Iterator<TAnnotationStore, Standard>::Type     TAnnoIter;
    typedef typename Id<TAnnotation>::Type                          TId;

    TAnnoIter it = begin(store.annotationStore, Standard());
    TAnnoIter itEnd = end(store.annotationStore, Standard());

    for (TId id = 0; it != itEnd; ++it, ++id)
        _writeOneAnnotation(target, store, *it, id, format);
}

template <typename TTargetStream, typename TSpec, typename TConfig>
inline void
write(
    TTargetStream & target,
    FragmentStore<TSpec, TConfig> & store,
    Gff format)
{
    _writeGffGtf(target, store, format);
}

template <typename TTargetStream, typename TSpec, typename TConfig>
inline void
write(
    TTargetStream & target,
    FragmentStore<TSpec, TConfig> & store,
    Gtf format)
{
    _writeGffGtf(target, store, format);
}

} // namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
