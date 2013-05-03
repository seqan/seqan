// ==========================================================================
//                                   Gustaf
// ==========================================================================
// Copyright (c) 2011, Knut Reinert, FU Berlin
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
// Author: Kathrin Trappe <kathrin.trappe@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_EXTRAS_APPS_GUSTAF_MSPLAZER_H_
#define SEQAN_EXTRAS_APPS_GUSTAF_MSPLAZER_H_

#define BREAKPOINT_DEBUG

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/arg_parse.h>

using namespace seqan;

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/**
.Tag.DotDrawingMSplazer
..cat:Input/Output
..summary:Switch to trigger drawing in dot format.
..value.DotDrawingMSplazer:Graphs in dot format.
..include:seqan/msplazer.h
*/

struct DotDrawingMSplazer_;
typedef Tag<DotDrawingMSplazer_> DotDrawingMSplazer;

/**
.Tag.DotDrawingMSplazerBestChain
..cat:Input/Output
..summary:Switch to trigger drawing in dot format.
..value.DotDrawingMSplazerBestChain:Best chain graphs in dot format.
..include:seqan/msplazer.h
*/

struct DotDrawingMSplazerBestChain_;
typedef Tag<DotDrawingMSplazer_> DotDrawingMSplazerBestchain;



// ----------------------------------------------------------------------------
// Class MSplazerOptions
// ----------------------------------------------------------------------------

struct MSplazerOptions
{
    // I/O options
    CharString databaseFile;        // name of database (db) file
    CharString queryFile;           // name of query file
    // CharString queryFile2;        // name of 2nd query file (mate pairs)
    CharString emptyReadsOutFile;   // file with reads that have broken or empty chains
    CharString disabledQueriesFile; // name of result file containing disabled queries
    CharString breakpointOutFile;   // breakpoint output file name
    CharString jobName;             // job name, used for dot output
    CharString stellarInputFile;    // optional input file with stellar matches
    bool dotOut;

    // Transposition penalties and thresholds (criteria for inserting edges in read graphs)
    unsigned diffDBPen;             // Penalty for matches on different databases
    unsigned diffStrandPen;         // Pen. for matches on diff. strand of the same database
    unsigned diffOrderPen;          // Pen. for matches with diff. order wrt the order within the query
    double simThresh;               // Allowed similarity between overlapping sequences, i.e. percentage of overlap
    int gapThresh;                  // Allowed gap or distance between matches
    int initGapThresh;              // Maximal allowed start or end gap length
    unsigned breakpointPosRange;    // Allowed range of breakpoint positions
    unsigned support;
    unsigned librSize;              // Library size (mate pairs)

    MSplazerOptions() :
        databaseFile("reference.fa"),
        queryFile("reads.fa"),
        emptyReadsOutFile("emptyReads.fa"),
        disabledQueriesFile("msplazer.disabled.fasta"),
        breakpointOutFile("breakpoints.gff"),
        dotOut(false),
        diffDBPen(5),
        diffStrandPen(5),
        diffOrderPen(0),
        simThresh(0.5),
        gapThresh(10),
        initGapThresh(15),
        breakpointPosRange(6),
        support(2),
        librSize(0){}
};

// ----------------------------------------------------------------------------
// Class Breakpoint
// ----------------------------------------------------------------------------

// BreakPoint class: container for structural variant or RNA-seq breakpoint information
template <typename TSequence_, typename TId_>
struct Breakpoint
{
    typedef TSequence_                          TSequence;
    typedef TId_                                TId;
    typedef typename Position<TSequence>::Type  TPos;

    // Ids of the two sequences
    TId startSeqId;
    TId endSeqId;
    // Sequence orientation
    bool startSeqStrand;
    bool endSeqStrand;
    // Last position in start sequence and first position in end sequence
    TPos startSeqPos;
    TPos endSeqPos;
    TPos readStartPos;
    TPos readEndPos;
    // Counter of occurrences (read support)
    unsigned support;
    // Query Sequence Ids (queries/reads that support the breakpoint)
    StringSet<TId> supportIds;
    // SV type
    TId svtype;
    TSequence insertionSeq;
    bool revStrandDel;

    Breakpoint() :
        startSeqId("####"),
        endSeqId("####"),
        startSeqStrand('+'),
        endSeqStrand('-'),
        startSeqPos(0),
        endSeqPos(0),
        readStartPos(0),
        readEndPos(0),
        support(1),
        svtype("svtype"),
        insertionSeq("NNNN"),
        revStrandDel(false)
    {}

    Breakpoint(TId const & sId,
               TId const & eId,
               bool const & sStrand,
               bool const & eStrand,
               TPos const & sPos,
               TPos const & ePos,
               TPos const & rsPos,
               TPos const & rePos) :
        startSeqId(sId),
        endSeqId(eId),
        startSeqStrand(sStrand),
        endSeqStrand(eStrand),
        startSeqPos(sPos),
        endSeqPos(ePos),
        readStartPos(rsPos),
        readEndPos(rePos),
        support(1),
        svtype("svtype"),
        insertionSeq("NNNN"),
        revStrandDel(false)
    {}

    Breakpoint(TId const & sId,
               TId const & eId,
               bool const & sStrand,
               bool const & eStrand,
               TPos const & sPos,
               TPos const & ePos,
               TPos const & rsPos,
               TPos const & rePos,
               TId const & spId) :
        startSeqId(sId),
        endSeqId(eId),
        startSeqStrand(sStrand),
        endSeqStrand(eStrand),
        startSeqPos(sPos),
        endSeqPos(ePos),
        readStartPos(rsPos),
        readEndPos(rePos),
        support(1),
        svtype("svtype"),
        insertionSeq("NNNN"),
        revStrandDel(false)
    {appendValue(supportIds, spId); }

    Breakpoint(TId const & sId,
               TId const & eId,
               bool const & sStrand,
               bool const & eStrand,
               TPos const & sPos,
               TPos const & ePos,
               TPos const & rsPos,
               TPos const & rePos,
               unsigned const & s,
               TId const & spId) :
        startSeqId(sId),
        endSeqId(eId),
        startSeqStrand(sStrand),
        endSeqStrand(eStrand),
        startSeqPos(sPos),
        endSeqPos(ePos),
        readStartPos(rsPos),
        readEndPos(rePos),
        support(s),
        svtype("svtype"),
        insertionSeq("NNNN"),
        revStrandDel(false)
    {appendValue(supportIds, spId); }

    Breakpoint(TId const & sId,
               TId const & eId,
               bool const & sStrand,
               bool const & eStrand,
               TPos const & sPos,
               TPos const & ePos,
               TPos const & rsPos,
               TPos const & rePos,
               unsigned const & s,
               StringSet<TId> const & spId) :
        startSeqId(sId),
        endSeqId(eId),
        startSeqStrand(sStrand),
        endSeqStrand(eStrand),
        startSeqPos(sPos),
        endSeqPos(ePos),
        readStartPos(rsPos),
        readEndPos(rePos),
        support(s),
        supportIds(spId),
        svtype("svtype"),
        insertionSeq("NNNN"),
        revStrandDel(false)
    {}
};

// ----------------------------------------------------------------------------
// Class SparsePropertyMap
// ----------------------------------------------------------------------------

// Sparse property map class: A property map that has only a few object or where most of the object would be empty
/**
.Class.SparsePropertyMap
..summary:Stores only a partial property map, instead of one with many empty entries, and a lookup table for the
  indices to the small property map.
..cat:Allocators
..signature:SparsePropertyMap<TValue, TPos>
..param.TValue:Type of stored objects.
..param.TPos:Type to store positions.

.Memfunc.SparsePropertyMap#SparsePropertyMap
..summary:Constructor
..signature:SparsePropertyMap<TValue,TPos> ()
..signature:SparsePropertyMap<TValue,TPos>  (TValueTable vt, TSlotLookupTable slt)
..param.vt:Allocator for objects.
..param.slt:Allocator for positions.
..class:Class.SparsePropertyMap
.Memvar.SparsePropertyMap#valueTable
..summary:Allocator for objects.
..class:Class.SparsePropertyMap
.Memvar.SparsePropertyMap#slotLookupTable
..summary:Allocator for positions.
..class:Class.SparsePropertyMap
..include:msplazer.h
*/
template <typename TValue, typename TPos>
struct SparsePropertyMap
{
    typedef String<TValue>       TValueTable;
    typedef String<TPos>         TSlotLookupTable;

    TValueTable valueTable;
    TSlotLookupTable    slotLookupTable;

    SparsePropertyMap(){}
    SparsePropertyMap(TValueTable vt, TSlotLookupTable slt) :
        valueTable(vt), slotLookupTable(slt) {}
};
/*
 * Only works within seqan namespace
template <typename TValue, typename TPos>
struct Value<SparsePropertyMap<TValue, TPos> >
{
    typedef TValue Type;
};
template <typename TValue, typename TPos>
struct Value<SparsePropertyMap<TValue, TPos> const>
{
    typedef TValue Type;
};
*/

// ----------------------------------------------------------------------------
// Class MSplazerChain
// ----------------------------------------------------------------------------

// Container for storing chaining graph, matchDistanceScores, start and end vertex for one read
template <typename TGraph_, typename TVertexDescriptor_, typename TScoreAlloc_, typename TSparsePropertyMap_,
          typename TMatchAlloc_>
struct MSplazerChain
{
    typedef TGraph_                 TGraph;
    typedef TVertexDescriptor_      TVertexDescriptor;
    typedef TScoreAlloc_            TScoreAlloc;
    typedef TSparsePropertyMap_     TSparsePropertyMap;
    typedef typename Size<TGraph>::Type TGraphSize;
    typedef TMatchAlloc_        TMatchAlloc;

    TGraph graph;                       // Contains (Stellar)matches as vertices and edges between compatible matches
    TVertexDescriptor startVertex;      // Artificial start and end vertex (represents start/end of the read sequence)
    TVertexDescriptor endVertex;
    String<TVertexDescriptor> predMap;  // Predecessor and distance map for dagShortestPath
    String<TGraphSize> distMap;         // Distance map for shortest path
    TScoreAlloc matchDistanceScores;    // Distance scores of matches (edit distance of reference mapping)
    TSparsePropertyMap breakpoints;
    String<TMatchAlloc> bestChains;
    bool isEmpty;
    bool isPartial;

    MSplazerChain(TScoreAlloc & _scores) :
        matchDistanceScores(_scores), isEmpty(false), isPartial(false)
    {}
};

struct Options
{
    bool showHelp;
    bool showVersion;
    int i;
    String<CharString> texts;

    Options(bool const & h, int & it) :
        showHelp(h), i(it)
    {}

};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function assignProperty()
// ----------------------------------------------------------------------------

/**
.Function.assignProperty:
..signature:assignProperty(spm, descr)
..param.spm: SparsePropertyMap.
...type:Class.SparsePropertyMap<TObject,TPos>
..param.descr: Vertex or edge descriptor.
...type:Metafunction.VertexDescriptor
...type:Metafunction.EdgeDescriptor
..include:msplazer.h
 */

template <typename TObject, typename TPos, typename TDescriptor>
inline void assignProperty(SparsePropertyMap<TObject, TPos> & spm, TDescriptor const & descr)
{
    assignProperty(spm.slotLookupTable, descr, -1);
}

/**
.Function.assignProperty:
..signature:assignProperty(spm, descr, val)
..param.spm: SparsePropertyMap.
...type:Class.SparsePropertyMap<TObject,TPos>
..param.descr: Vertex or edge descriptor.
...type:Metafunction.VertexDescriptor
...type:Metafunction.EdgeDescriptor
..param.val:The new value.
...remarks:Type of the new value must match the value type of the properties map object table.
.include:msplazer.h
 */

template <typename TObject, typename TPos, typename TDescriptor, typename TValue>
inline void assignProperty(SparsePropertyMap<TObject, TPos> & spm,
                           TDescriptor const & descr,
                           TValue const & val)
{
    assignProperty(spm.slotLookupTable, descr, length(spm.valueTable));
    appendValue(spm.valueTable, val);
}

// ----------------------------------------------------------------------------
// Function property()
// ----------------------------------------------------------------------------

/**
.Function.property:
..summary: Returns object at the given position.
..signature:getProperty(spm, descr, obj)
..param.spm: SparsePropertyMap.
...type:Class.SparsePropertyMap<TObject,TPos>
..param.descr: Vertex or edge descriptor.
...type:Metafunction.VertexDescriptor
...type:Metafunction.EdgeDescriptor
..returns:Reference to the item in the property map.
..param.obj: Output parameter (object contained in sparse property map.
...remarks:Type must match the value type of the properties map object table. Assumes that there is an object at this
  position!
.include:msplazer.h
 */

template <typename TObject, typename TPos, typename TDescriptor>
inline typename Reference<TObject>::Type property(SparsePropertyMap<TObject, TPos> & spm,
                                                  TDescriptor const & descr)
{

    TPos index = getProperty(spm.slotLookupTable, descr);
    // obj = getValue(spm.valueTable, index);
    return value(spm.valueTable, index);
}

template <typename TObject, typename TPos, typename TDescriptor>
inline typename Reference<TObject>::Type property(SparsePropertyMap<TObject, TPos> const & spm,
                                                  TDescriptor const & descr)
{

    TPos index = getProperty(spm.slotLookupTable, descr);
    // obj = getValue(spm.valueTable, index);
    return value(spm.valueTable, index);
}

// ----------------------------------------------------------------------------
// Function getProperty()
// ----------------------------------------------------------------------------

/**
.Function.getProperty:
..summary: Returns false if there is no object at the given position. Otherwise the object is written to the output
  parameter obj.
..signature:getProperty(spm, descr, obj)
..param.spm: SparsePropertyMap.
...type:Class.SparsePropertyMap<TObject,TPos>
..param.descr: Vertex or edge descriptor.
...type:Metafunction.VertexDescriptor
...type:Metafunction.EdgeDescriptor
..param.obj: Output parameter (object contained in sparse property map.
...remarks:Type must match the value type of the properties map object table.
.include:msplazer.h
 */

template <typename TObject, typename TPos, typename TDescriptor>
inline bool getProperty(SparsePropertyMap<TObject, TPos> & spm, TDescriptor const & descr, TObject & obj)
{
    TPos index = getProperty(spm.slotLookupTable, descr);
    if (index == static_cast<TPos>(-1))
        return 0;

    // obj = getValue(spm.valueTable, index);
    obj = value(spm.valueTable, index);
    return 1;
}

template <typename TObject, typename TPos, typename TDescriptor>
inline bool getProperty(SparsePropertyMap<TObject, TPos> const & spm, TDescriptor const & descr, TObject & obj)
{
    TPos index = getProperty(spm.slotLookupTable, descr);
    if (index == static_cast<TPos>(-1))
        return 0;

    // obj = getValue(spm.valueTable, index);
    obj = value(spm.valueTable, index);
    return 1;
}

// ----------------------------------------------------------------------------
// Function appendSupportId()
// ----------------------------------------------------------------------------

/**
.Function.appendSupportId:
..signature:appendSupportId(bp, id)
..param.bp: Breakpoint.
...type:Class.Breakpoint
..param.id: Id to be appended.
...type:Class.TId
..include:msplazer.h
 */

template <typename TBreakpoint, typename TId>
inline void appendSupportId(TBreakpoint & bp, TId const & id)
{
    for (unsigned i = 0; i < length(bp.supportIds); ++i)
    {
        if (bp.supportIds[i] == id)
        {
            ++bp.support;
            return;
        }
    }
    appendValue(bp.supportIds, id);
    ++bp.support;
}

/**
.Function.appendSupportId:
..signature:appendSupportId(bp, id)
..param.bp: Breakpoint.
...type:Class.Breakpoint
..param.id: Id to be appended.
...type:TId
..include:msplazer.h
 */

template <typename TBreakpoint, typename TId>
inline void appendSupportId(TBreakpoint & bp, StringSet<TId> const & ids)
{
/*
    typedef typename Iterator<StringSet<TId> >::Type TIterator;
    TIterator it = begin(ids);
    for(;!atEnd(it);goNext(it))
        appendSupportId(bp, *it);
*/
    for (unsigned i = 0; i < length(ids); ++i)
        appendSupportId(bp, ids[i]);
}

/**
.Function.setSupport:
..signature:setSupport(bp, value)
..param.bp: Breakpoint.
...type:Class.Breakpoint
..param.value: New support value.
...type:TId
..include:msplazer.h
 */

// ----------------------------------------------------------------------------
// Function setSupport()
// ----------------------------------------------------------------------------
template <typename TBreakpoint>
inline void setSupport(TBreakpoint & bp, unsigned const & value)
{
    assignValue(bp.support, value);
}

// ----------------------------------------------------------------------------
// Function getSVType()
// ----------------------------------------------------------------------------

/**
.Function.getSVType:
..signature:getSVType(bp)
..summary:Return the breakpoints svtype.
..param.bp: Breakpoint.
...type:Class.Breakpoint
..include:msplazer.h
 */

template <typename TSequence, typename TId>
inline TId getSVType(Breakpoint<TSequence, TId> & bp)
{
    return bp.svtype;
}

// ----------------------------------------------------------------------------
// Function setSVType()
// ----------------------------------------------------------------------------

/**
.Function.setSVType:
..signature:setSVType(bp, type)
..summary:Sets the breakpoints svtype to "type".
..param.bp: Breakpoint.
...type:Class.Breakpoint
..param.type: TId.
...type:Class.Breakpoint.TId
..include:msplazer.h
 */

template <typename TBreakpoint, typename TId>
inline void setSVType(TBreakpoint & bp, TId const & type)
{
    bp.svtype = type;
}

/**
.Function.setSVType:
..signature:setSVType(bp)
..summary:Computes SV type of the beakpoint. Returns "true" for insertion and "false" otherwise.
..param.bp: Breakpoint.
...type:Class.Breakpoint
..include:msplazer.h
 */

template <typename TBreakpoint>
inline bool setSVType(TBreakpoint & bp)
{
    typedef typename TBreakpoint::TId TId;
    // if insertion return 1; else return 0;
    if (bp.startSeqId != bp.endSeqId)
    {
        bp.svtype = static_cast<TId>("translocation");
        // setSVType(bp, static_cast<TId>("translocation"));
        return false;
    }
    if (bp.startSeqStrand != bp.endSeqStrand)
    {
        bp.svtype = static_cast<TId>("inversion");
        // setSVType(bp, static_cast<TId>("inversion"));
        return false;
    }
    if (bp.startSeqPos < bp.endSeqPos)
    {
        bp.svtype = static_cast<TId>("deletion");
        // setSVType(bp, static_cast<TId>("deletion"));
        return false;
    }
    if (bp.startSeqPos > bp.endSeqPos)
    {
        std::swap(bp.startSeqPos, bp.endSeqPos);
        if (bp.startSeqStrand)
        {
            bp.svtype = static_cast<TId>("translocation");
            return false;
        }
        setSVType(bp, static_cast<TId>("deletion"));
        bp.revStrandDel = true;
        return false;
    }
    bp.svtype = static_cast<TId>("insertion");
    // setSVType(bp, static_cast<TId>("insertion"));
    return true;
}

/**
.Function.setInsertionSeq:
..signature:setInsertionSeq(bp, inSeq)
..summary:Computes SV type of the beakpoint. Returns true for insertion and false otherwise.
..param.bp: Breakpoint.
...type:Class.Breakpoint
..param.inSeq: Insertion sequence.
...type:TSequence
..include:msplazer.h
 */

template <typename TBreakpoint, typename TSequence>
inline void setInsertionSeq(TBreakpoint & bp, TSequence & inSeq)
{
    bp.insertionSeq = inSeq;
    // assignValue(bp.insertionSeq, inSeq);
}

/**
.Function.posInSameRange:
..signature:posInSameRange(pos1, pos2, range)
..param.pos1: First position to be compared.
...type:Class.Position
..param.pos2: Snd position to be compared.
...type:Class.Position
..param.range: Valid range for position difference.
...type:Class.Position
..include:msplazer.h
 */

template <typename TPos, typename TPosR>
inline bool posInSameRange(TPos const & pos1, TPos const & pos2, TPosR const & range)
{
    return abs(pos2 - pos1) < range;
}

/**
.Function.similarBreakpoints:
..summary:Tests two breakpoints for similarity, i.e. if they have the same sequence Ids and lie within a specified range.
..signature:similarBreakpoints(bp1, bp2i)
..param.bp1:First breakpoint to be compared.
...type:Class.Breakpoint
..param.bp2:Snd breakpoint to be compared.
...type:Class.Breakpoint
..include:seqan/msplazer.h
 */

template <typename TId, typename TPos>
inline bool similarBreakpoints(Breakpoint<TId, TPos> const & bp1, Breakpoint<TId, TPos> const & bp2)
{
    bool sameSeqs = (bp1.startSeqId == bp2.startSeqId) && (bp1.endSeqId == bp2.endSeqId);
    bool samePosRange = posInSameRange(bp1.startSeqPos, bp2.startSeqPos) && posInSameRange(bp1.endSeqPos,
                                                                                           bp2.endSeqPos);
    return sameSeqs && samePosRange;
}

/**
.Function.operator==:
..summary:Operator== implementation for Breakpoints. A Breakpoint has two Ids, two strands, two positions and a SV type,
and can be compared according to them (in this order of priority).
..include:seqan/msplazer.h
 */

template <typename TId, typename TPos>
inline bool operator==(Breakpoint<TId, TPos> const & bp1, Breakpoint<TId, TPos> const & bp2)
{
    if (bp1.startSeqId != bp2.startSeqId)
        return false;

    if (bp1.endSeqId != bp2.endSeqId)
        return false;

    /*
    if(bp1.startSeqStrand != bp2.startSeqStrand)
        return false;
    if(bp1.endSeqStrand != bp2.endSeqStrand)
        return false;
    */
    if (bp1.startSeqPos != bp2.startSeqPos)
        return false;

    if (bp1.endSeqPos != bp2.endSeqPos)
        return false;

    if (bp1.svtype != bp2.svtype)
        return false;

    if (bp1.svtype == "insertion" && bp2.svtype == "insertion")
        return length(bp1.insertionSeq) == length(bp2.insertionSeq);

    return true;
}

// ----------------------------------------------------------------------------
// Function operator<(Breakpoint)
// ----------------------------------------------------------------------------

/**
.Function.operator<:
..summary:Operator< implementation for Breakpoints. A Breakpoint has two Ids and two positions,
and can be sorted according to them (in this order of priority).
..include:seqan/msplazer.h
 */

template <typename TId, typename TPos>
inline bool operator<(Breakpoint<TId, TPos> const & bp1, Breakpoint<TId, TPos> const & bp2)
{
    if (bp1.startSeqId != bp2.startSeqId)
        return bp1.startSeqId < bp2.startSeqId;

    if (bp1.endSeqId != bp2.endSeqId)
        return bp1.endSeqId < bp2.endSeqId;

    if (bp1.startSeqPos != bp2.startSeqPos)
        return bp1.startSeqPos < bp2.startSeqPos;

    return bp1.endSeqPos < bp2.endSeqPos;
}

/**
 * Ostream operator << for Breakpoint class
 */
template <typename TSequence, typename TId, typename TStream>
// std::ostream & operator<<(std::ostream & out, Breakpoint<TSequence, TId> const & value)
TStream & operator<<(TStream & out, Breakpoint<TSequence, TId> const & value)
{
    out << "Breakpoint: seq1 --> seq2; posInSeq1 --> posInSeq2; readPos1 --> readPos2 :" << std::endl;
    out << value.startSeqId << " ( " << value.startSeqStrand << " ) " << " --> " << value.endSeqId << " ( " <<
    value.endSeqStrand << " ) " << std::endl;
    out << " ( " << value.startSeqPos + 1 << " ) --> ( " << value.endSeqPos + 1 << " ) " << std::endl;
    out << " ( " << value.readStartPos + 1 << " ) --> ( " << value.readEndPos + 1 << " ) " << std::endl;
    out << "SVType: " << value.svtype << " insertionSeq: " << value.insertionSeq << std::endl;
    out << "Support: " << value.support << " Ids: ";
    for (unsigned i = 0; i < length(value.supportIds); ++i)
        out << value.supportIds[i] << ", ";
    out << std::endl;
    return out;
}

/**
 * Ostream operator << for StellarMatch
 */
template <typename TSequence, typename TId, typename TStream>
TStream & operator<<(TStream & out, StellarMatch<TSequence, TId> & match)
{
    out << "DB Id: " << match.id;
    if (match.orientation)
        out << " + " << std::endl;
    else
        out << " - " << std::endl;
    out << "DB pos: " << match.begin1 << " ... " << match.end1 << std::endl;
    out << "Query pos: " << match.begin2 << " ... " << match.end2 << std::endl;

    if (!match.orientation)
        reverseComplement(infix(source(match.row1), match.begin1, match.end1));

    typedef typename StellarMatch<TSequence, TId>::TAlign TAlign;
    TAlign align;
    // assignSource makes a local copy of the source sequence (only the part of the row)
    appendValue(align.data_rows, match.row1);
    appendValue(align.data_rows, match.row2);
    out << align;
    out << std::endl;
    if (!match.orientation)
        reverseComplement(infix(source(match.row1), match.begin1, match.end1));
    return out;
}

struct BreakpointFileHeader_;
typedef Tag<BreakpointFileHeader_> BreakpointFileHeader;

template <typename TSequence, typename TId>
inline void _streamPut(::std::FILE * target, Breakpoint<TSequence, TId> const & bp, unsigned const & counter)
{
    TId sId;
    // typedef typename Breakpoint<TSequence, TId>::TPos TPos;
    typedef typename Position<TSequence>::Type  TPos;
    _getShortId(sId, bp.startSeqId);

    _streamWrite(target, sId);
    _streamPut(target, '\t');
    _streamWrite(target, "GUSTAF");
    _streamPut(target, '\t');
    _streamWrite(target, bp.svtype);
    _streamPut(target, '\t');
    _streamPutInt(target, bp.startSeqPos + 1);
    _streamPut(target, '\t');
    if (bp.svtype == "deletion" || bp.svtype == "inversion")
        _streamPutInt(target, bp.endSeqPos);
    else
        _streamPutInt(target, bp.startSeqPos + 1);
    _streamPut(target, '\t');
    _streamWrite(target, '.');
    _streamPut(target, '\t');
    if (bp.startSeqStrand)
        _streamWrite(target, '+');
    else
        _streamWrite(target, '-');
    _streamPut(target, '\t');
    _streamWrite(target, '.');
    _streamPut(target, '\t');

    _streamWrite(target, "ID=");
    _streamPutInt(target, counter);
    _streamWrite(target, ';');
    if (bp.svtype == "insertion")
    {
        _streamWrite(target, "size=-");
        _streamPutInt(target, length(bp.insertionSeq));
        _streamWrite(target, ';');
        _streamWrite(target, "seq=");
        _streamWrite(target, bp.insertionSeq);
        _streamWrite(target, ';');
    }
    else if (bp.svtype == "deletion")
    {
        _streamWrite(target, "size=");
        _streamPutInt(target, static_cast<TPos>(bp.endSeqPos - bp.startSeqPos));
        _streamWrite(target, ';');
    }
    else
    {
        _streamWrite(target, "endChr=");
        _getShortId(sId, bp.endSeqId);
        _streamWrite(target, sId);
        _streamWrite(target, ';');
        // if(bp.svtype == "translocation"){
        _streamWrite(target, "endPos=");
        _streamPutInt(target, bp.endSeqPos);
        _streamWrite(target, ';');
        // }
        _streamWrite(target, "endStrand=");
        if (bp.endSeqStrand)
            _streamWrite(target, '+');
        else
            _streamWrite(target, '-');
        _streamWrite(target, ';');
    }
    _streamWrite(target, "support=");
    _streamPutInt(target, bp.support);
    _streamWrite(target, ';');
    _streamWrite(target, "supportIds=");
    for (unsigned i = 0; i < length(bp.supportIds); ++i)
    {
        _streamWrite(target, bp.supportIds[i]);
        _streamWrite(target, ',');
    }

    _streamWrite(target, ';');
    _streamWrite(target, "breakpoint:");
    _streamPutInt(target, bp.readStartPos + 1);

    _streamPut(target, '\n');
}

template <typename TSequence, typename TId>
inline void _streamWrite(::std::FILE * target, Breakpoint<TSequence, TId> const & bp, unsigned const & counter)
{
    _streamPut(target, bp, counter);
}

inline void _streamPut(::std::FILE * target, BreakpointFileHeader const &)
{
    _streamWrite(target, "StartSeqId");
    _streamPut(target, '\t');
    _streamWrite(target, "Label");
    _streamPut(target, '\t');
    _streamWrite(target, "SV type");
    _streamPut(target, '\t');
    _streamWrite(target, "sPos");
    _streamPut(target, '\t');
    _streamWrite(target, "sPos");
    _streamPut(target, '\t');
    _streamWrite(target, ".");
    _streamPut(target, '\t');
    _streamWrite(target, "+/-");
    _streamPut(target, '\t');
    _streamWrite(target, ".");
    _streamPut(target, '\t');
    _streamWrite(target, "Tags:ID=i;size=(-)indel;EndSeqId=;EndSeqPos=;EndSeqStrand=;Support=;SupportIds=;\n");
}

inline void _streamWrite(::std::FILE * target, BreakpointFileHeader const &)
{
    _streamPut(target, BreakpointFileHeader());
}

// ----------------------------------------------------------------------------
// Function insertBestChain
// ----------------------------------------------------------------------------

/**
.Function.insertBestChain:
..signature:insertBestChain(msplChain, chain)
..param.msplChain: MSplazerChain object.
...type:Class.MSplazerChain
..param.chain: New chain.
...type:Class.MSplazerChain.TMatchAlloc
..include:seqan/msplazer.h
 */
template <typename TMSplazerChain, typename TMatchAlloc>
void insertBestChain(TMSplazerChain & mspChain, TMatchAlloc const & chain)
{
    appendValue(mspChain.bestChains, chain);
}

/**
.Function.insertBestChain:
..signature:write(msplChain, chain)
..param.msplChain: MSplazerChain object.
...type:Class.MSplazerChain
..param.chain: New chain.
...type:Class.MSplazerChain.TMatchAlloc
..include:seqan/msplazer.h
 */
template <typename TMSplazerChain, typename TMatchAlloc>
void insertBestChain(TMSplazerChain & mspChain, TMatchAlloc & chain)
{
    appendValue(mspChain.bestChains, chain);
}

template <typename T>
inline std::string toString(const T & t)
{
    std::stringstream ss;
    ss << t;
    return ss.str();
}

int mainWithOptions(Options & options)
{
    typedef Iterator<String<CharString> >::Type TIterator;
    std::cout << "Non-option Arguments:" << std::endl;
    for (TIterator it = begin(options.texts); it != end(options.texts); ++it)
    {
        std::cout << "  " << *it << std::endl;
    }

    return 0;
}

#endif  // #ifndef SANDBOX_MY_SANDBOX_APPS_MSPLAZER_MSPLAZER_H_