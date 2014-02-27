// ==========================================================================
//                                   Gustaf
// ==========================================================================
// Copyright (c) 2011-2013, Kathrin Trappe, FU Berlin
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

#include "create_stellarmatches_from_file.h"

using namespace seqan;

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================


// ----------------------------------------------------------------------------
// Class MSplazerOptions
// ----------------------------------------------------------------------------

struct MSplazerOptions
{
    // I/O options
    CharString databaseFile;        // name of database (db) file
    StringSet<CharString> queryFile;           // name of query file(s) (two in case of paired-end)
    // CharString queryFile2;        // name of 2nd query file (mate pairs)
    CharString emptyReadsOutFile;   // file with reads that have broken or empty chains
    CharString disabledQueriesFile; // name of result file containing disabled queries
    CharString vcfOutFile;   // breakpoint output file name
    CharString gffOutFile;   // breakpoint output file name
    CharString jobName;             // job name, used for dot output
    CharString stellarInputFile;    // optional input file with stellar matches
    bool dotOut;

    // Transposition penalties and thresholds (criteria for inserting edges in read graphs)
    unsigned diffDBPen;             // Penalty for matches on different databases
    unsigned diffStrandPen;         // Pen. for matches on diff. strand of the same database
    unsigned diffOrderPen;          // Pen. for matches with diff. order wrt the order within the query
    unsigned noMateMatchesPen;      // Pen. for matches with no confirming mate matches
    double simThresh;               // Allowed similarity between overlapping sequences, i.e. percentage of overlap
    int gapThresh;                  // Allowed gap or distance between matches
    int initGapThresh;              // Maximal allowed start or end gap length
    int breakendThresh;             // Maximal allowed length for a breakend
    int tandemThresh;               // Minimal length of insertion (double overlap) to be called as tandem repeat
    unsigned breakpointPosRange;    // Allowed range of breakpoint positions
    unsigned support;
    unsigned mateSupport;
    unsigned libSize;               // Library size (mate pairs)
    unsigned libError;              // Library size (mate pairs)
    bool pairedEndMode;             // Whether or not to run Gustaf in paired-end mode
    bool revCompl;                  // Whether or not to rev-compl the second input file


    MSplazerOptions() :
        databaseFile("reference.fa"),
        emptyReadsOutFile("emptyReads.fa"),
        disabledQueriesFile("msplazer.disabled.fasta"),
        vcfOutFile("breakpoints.vcf"),
        gffOutFile("breakpoints.gff"),
        dotOut(false),
        diffDBPen(5),
        diffStrandPen(5),
        diffOrderPen(0),
        noMateMatchesPen(5),
        simThresh(0.5),
        gapThresh(10),
        initGapThresh(15),
        breakendThresh(30),
        tandemThresh(50),
        breakpointPosRange(5),
        support(2),
        mateSupport(1),
        libSize(0),
        libError(0),
        pairedEndMode(false),
        revCompl(true){}
};

// ----------------------------------------------------------------------------
// Class Breakpoint
// ----------------------------------------------------------------------------

// BreakPoint class: container for structural variant or RNA-seq breakpoint information
template <typename TSequence_, typename TId_>
struct Breakpoint
{
    enum SVType
    {
        INVALID,            // 0
        INSERTION,          // 1
        DELETION,           // 2
        INVERSION,          // 3
        TANDEM,             // 4
        DISPDUPLICATION,    // 5
        INTRATRANSLOCATION, // 6
        TRANSLOCATION,      // 7
        BREAKEND            // 8
    };

    typedef TSequence_                          TSequence;
    typedef TId_                                TId;
    typedef typename Position<TSequence>::Type  TPos;

    // Ids of the two sequences
    TId startSeqId;
    TId endSeqId;
    TId midPosId;
    // Sequence orientation
    bool startSeqStrand;
    bool endSeqStrand;
    bool midPosStrand;
    // Last position in start sequence and first position in end sequence
    TPos startSeqPos;  // End position of split alignment, i.e. first position of variant
    TPos endSeqPos;    // Begin position of split alignment, i.e. position after split in reference, i.e. after the variant
    TPos dupTargetPos; // Pos. before inserted, duplicated seq (should be equal to either startSeqPos or endSeqPos)
    TPos dupMiddlePos; // TPos dupMidPos; Middle position of duplication or translocation
                       // per default position before event, has to be adjustet depending on dupTargetPos
    TPos readStartPos;
    TPos readEndPos;
    TPos cipos;       // Confidence after startSeqPos for imprecise breakpoint
    TPos ciend;       // Confidence after endSeqPos for imprecise breakpoint
    TPos cimiddle;    // Confidence after dupMiddlePos for imprecise breakpoint
    // Counter of occurrences (read support)
    unsigned support;
    // Query Sequence Ids (queries/reads that support the breakpoint)
    StringSet<TId> supportIds;
    // SV type
    SVType svtype;
    TSequence insertionSeq;
    bool revStrandDel;
    // If both of these flags are true, then we have seen two (pseudo)deletions supporting both start and end position
    // of a translocation.
    bool translSuppStartPos;
    bool translSuppEndPos;
    // bool imprecise = false;
    // Storing on which site the breakend is: 
    // 0: left breakend, i.e. sequence continues right of position
    // 1: right breakend, i.e. sequence continues left of position
    bool breakend;

    Breakpoint() :
        startSeqId("####"),
        endSeqId("####"),
        midPosId("####"),
        startSeqStrand(true),
        endSeqStrand(true),
        startSeqPos(0),
        endSeqPos(0),
        dupTargetPos(maxValue<unsigned>()),
        dupMiddlePos(maxValue<unsigned>()),
        readStartPos(0),
        readEndPos(0),
        cipos(maxValue<unsigned>()),
        ciend(maxValue<unsigned>()),
        cimiddle(maxValue<unsigned>()),
        support(1),
        svtype(INVALID),
        insertionSeq("NNNN"),
        revStrandDel(false),
        translSuppStartPos(false),
        translSuppEndPos(false),
        breakend(false)
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
        midPosId("####"),
        startSeqStrand(sStrand),
        endSeqStrand(eStrand),
        startSeqPos(sPos),
        endSeqPos(ePos),
        dupTargetPos(maxValue<unsigned>()),
        dupMiddlePos(maxValue<unsigned>()),
        readStartPos(rsPos),
        readEndPos(rePos),
        cipos(maxValue<unsigned>()),
        ciend(maxValue<unsigned>()),
        cimiddle(maxValue<unsigned>()),
        support(1),
        svtype(INVALID),
        insertionSeq("NNNN"),
        revStrandDel(false),
        translSuppStartPos(false),
        translSuppEndPos(false),
        breakend(false)
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
        midPosId("####"),
        startSeqStrand(sStrand),
        endSeqStrand(eStrand),
        startSeqPos(sPos),
        endSeqPos(ePos),
        dupTargetPos(maxValue<unsigned>()),
        dupMiddlePos(maxValue<unsigned>()),
        readStartPos(rsPos),
        readEndPos(rePos),
        cipos(maxValue<unsigned>()),
        ciend(maxValue<unsigned>()),
        cimiddle(maxValue<unsigned>()),
        support(1),
        svtype(INVALID),
        insertionSeq("NNNN"),
        revStrandDel(false),
        translSuppStartPos(false),
        translSuppEndPos(false)
    {appendValue(supportIds, spId); }
    /*
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
        svtype(INVALID),
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
        svtype(INVALID),
        insertionSeq("NNNN"),
        revStrandDel(false)
    {}
    */
};

// ----------------------------------------------------------------------------
// Class SparsePropertyMap
// ----------------------------------------------------------------------------

// Sparse property map class: A property map that has only a few objects or where most of the object would be empty
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
    // slotLookupTable: Stores at each position (descriptor value) the position of the
    // corresponding object in valueTable, or -1 if there is no object
    // for this descriptor
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
// Vertex descriptor value of each vertex in the graph corresponds to the position of the stellar match within
// container for Stellar matches (Query Matches)
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
    typedef IntervalAndCargo<unsigned, unsigned> TInterval;
    typedef IntervalTree<unsigned, unsigned> TIntervalTree;

    TGraph graph;                       // Contains (Stellar)matches as vertices and edges between compatible matches
    TVertexDescriptor startVertex;      // Artificial start and end vertex (represents start/end of the read sequence)
    TVertexDescriptor endVertex;
    String<TVertexDescriptor> predMap;  // Predecessor and distance map for dagShortestPath
    String<TGraphSize> distMap;         // Distance map for shortest path
    TScoreAlloc matchDistanceScores;    // Distance scores of matches (edit distance of reference mapping)
    TSparsePropertyMap breakpoints;
    String<TMatchAlloc> bestChains;
    TIntervalTree rightMateTree;
    TIntervalTree leftMateTree;
    unsigned mateJoinPosition;
    bool isEmpty;
    bool isPartial;
    // bool transl/dupl;

    MSplazerChain(TScoreAlloc & _scores) :
        matchDistanceScores(_scores), mateJoinPosition(0), isEmpty(false), isPartial(false)
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

    obj = value(spm.valueTable, index);
    return 1;
}

template <typename TObject, typename TPos, typename TDescriptor>
inline bool getProperty(SparsePropertyMap<TObject, TPos> const & spm, TDescriptor const & descr, TObject & obj)
{
    TPos index = getProperty(spm.slotLookupTable, descr);
    if (index == static_cast<TPos>(-1))
        return 0;

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

template <typename TSequence, typename TId, typename TSVType>
inline TSVType getSVType(Breakpoint<TSequence, TId> & bp)
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

template <typename TBreakpoint, typename TSVType>
inline void setSVType(TBreakpoint & bp, TSVType type)
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
    // if insertion return 1; else return 0;
    if (bp.startSeqId != bp.endSeqId)
    {
        bp.svtype = TBreakpoint::INTRATRANSLOCATION;
        return false;
    }
    if (bp.startSeqStrand != bp.endSeqStrand)
    {
        if (bp.startSeqPos > bp.endSeqPos)
        {
            std::swap(bp.startSeqPos, bp.endSeqPos);
            // bp.startSeqPos = bp.startSeqPos - 1;
            // bp.endSeqPos = bp.endSeqPos + 1;
        }
        setSVType(bp, TBreakpoint::INVERSION);
        return false;
    }
    if (bp.startSeqPos < bp.endSeqPos)
    {
        setSVType(bp, TBreakpoint::DELETION);
        return false;
    }
    if (bp.startSeqPos > bp.endSeqPos)
    {
        std::swap(bp.startSeqPos, bp.endSeqPos);
        // bp.startSeqPos = bp.startSeqPos - 1;
        // bp.endSeqPos = bp.endSeqPos + 1;
        if (bp.startSeqStrand)
        {
            setSVType(bp, TBreakpoint::DISPDUPLICATION);
            // setSVType(bp, TBreakpoint::TRANSLOCATION);
            return false;
        }
        setSVType(bp, TBreakpoint::DELETION);
        bp.revStrandDel = true;
        return false;
    }
    setSVType(bp, TBreakpoint::INSERTION);
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
inline bool _posInSameRange(TPos const & pos1, TPos const & pos2, TPosR const & range)
{
    return (__int64)abs(pos2 - pos1) < (__int64)(range + 1);
}

// Breakends are distinguishable by there one reference Id (startId=endId) and position
// (start and end position are the same), strand doesn't matter
template <typename TId, typename TPos, typename TPosR>
inline bool _similarBreakends(Breakpoint<TId, TPos> & be1, Breakpoint<TId, TPos> & be2, TPosR const & range)
{
    if (be1.startSeqId != be2.startSeqId)
        return false;
    return _posInSameRange(be1.startSeqPos, be2.startSeqPos, range);
}

template <typename TId, typename TPos, typename TPosR>
inline bool _breakendSupport(Breakpoint<TId, TPos> & be, Breakpoint<TId, TPos> & bp, TPosR const & range)
{
    typedef Breakpoint<TId, TPos> TBreakpoint;
    if (be.startSeqId != bp.startSeqId && be.startSeqId != bp.endSeqId)
        return false;
    // If bp is duplication or translocation, also check targetpos
    if ((bp.svtype == TBreakpoint::DISPDUPLICATION || bp.svtype == TBreakpoint::TRANSLOCATION) && bp.dupMiddlePos != maxValue<unsigned>())
        return (_posInSameRange(be.startSeqPos, bp.startSeqPos, range) ||
                _posInSameRange(be.startSeqPos, bp.endSeqPos, range)   ||
                _posInSameRange(be.startSeqPos, bp.dupMiddlePos, range) );

    return (_posInSameRange(be.startSeqPos, bp.startSeqPos, range) ||
            _posInSameRange(be.startSeqPos, bp.endSeqPos, range) );
}
/**
.Function.similarBreakpoints:
..summary:Tests two breakpoints for similarity, i.e. if they have the same sequence Ids and lie within a specified range.
..signature:similarBreakpoints(bp1, bp2)
..param.bp1:First breakpoint to be compared.
...type:Class.Breakpoint
..param.bp2:Snd breakpoint to be compared.
...type:Class.Breakpoint
..include:seqan/msplazer.h
 */

template <typename TId, typename TPos>
inline bool _similarBreakpoints(Breakpoint<TId, TPos> & bp1, Breakpoint<TId, TPos> & bp2)
{
    typedef Breakpoint<TId, TPos> TBP;
    unsigned const range = 5;
    bool sameOrient = true;
    bool sameSeqs = (bp1.startSeqId == bp2.startSeqId) && (bp1.endSeqId == bp2.endSeqId);
    // bool sameOrient = (bp1.startSeqStrand == bp2.startSeqStrand) && (bp1.endSeqStrand == bp2.endSeqStrand);
    bool samePosRange = 0;
    if (bp1.svtype == TBP::BREAKEND || bp2.svtype == TBP::BREAKEND)
        samePosRange = posInSameRange(bp1.startSeqPos, bp2.startSeqPos, range)
            || posInSameRange(bp1.endSeqPos, bp2.endSeqPos, range);
    else
        samePosRange = posInSameRange(bp1.startSeqPos, bp2.startSeqPos, range)
            && posInSameRange(bp1.endSeqPos, bp2.endSeqPos, range);
    return sameSeqs && sameOrient && samePosRange;
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

    if (bp1.svtype == 1 && bp2.svtype == 1)
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
    typedef Breakpoint<TSequence, TId> TBreakpoint;
    out << "Breakpoint: seq1 --> seq2; posInSeq1 --> posInSeq2; readPos1 --> readPos2 :" << std::endl;
    out << value.startSeqId << " ( " << value.startSeqStrand << " ) " << " --> " << value.endSeqId << " ( " <<
    value.endSeqStrand << " ) " << std::endl;
    out << " ( " << value.startSeqPos + 1 << " ) --> ( " << value.endSeqPos + 1 << " ) " << std::endl;
    if (value.dupMiddlePos != maxValue<unsigned>())
        out << "target pos " << value.dupTargetPos << std::endl;
    out << " ( " << value.readStartPos + 1 << " ) --> ( " << value.readEndPos + 1 << " ) " << std::endl;
    switch (value.svtype)
    {
        case TBreakpoint::INVALID:
        out << "SVType: invalid";
        break;
        case TBreakpoint::INSERTION:
        out << "SVType: insertion";
        break;
        case TBreakpoint::DELETION:
        out << "SVType: deletion";
        break;
        case TBreakpoint::INVERSION:
        out << "SVType: inversion";
        break;
        case TBreakpoint::TANDEM:
        out << "SVType: tandem";
        break;
        case TBreakpoint::DISPDUPLICATION:
        out << "SVType: duplication";
        break;
        case TBreakpoint::INTRATRANSLOCATION:
        out << "SVType: intra-transl";
        break;
        case TBreakpoint::TRANSLOCATION:
        out << "SVType: inter-transl";
        break;
        case TBreakpoint::BREAKEND:
        out << "SVType: breakend";
    }
    out << " insertionSeq: " << value.insertionSeq << std::endl;
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
