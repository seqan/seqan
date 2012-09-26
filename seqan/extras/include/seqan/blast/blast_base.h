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

#ifndef SEQAN_HEADER_BLAST_BASE_H
#define SEQAN_HEADER_BLAST_BASE_H


namespace SEQAN_NAMESPACE_MAIN
{


//////////////////////////////////////////////////////////////////////////////
// Blast Report types


struct FullInfo;
//...remarks:BasicInfo stores begin and end positions on query and database sequence, as well as the alignment. FullInfo stores additional information such as score, e-value...

struct BasicInfo;



/**
.Spec.StoreReport:
..cat:Blast
..general:Class.BlastReport
..summary:BlastReport specialization. Parses a Blast report and stores all hits and HSPs.
..signature:BlastReport<TBlastHsp,StoreReport<TSpec> >
..param.TBlastHsp:The type of HSPs to be stored. See @Class.BlastHsp@
...metafunction:Metafunction.Hsp
...default:BlastHsp<BlastN,BasicInfo> 
..param.TSpec:The specializing type.
...default:BasicInfo
...type:Spec.BasicInfo
...type:Spec.FullInfo 
..include:blast.h
*/
//
//...remarks:BasicInfo only stores query name, database name and a String of all hits found. FullInfo also stores the following 
//parameters: lambda, k, h, gapped_lambda, gapped_k, gapped_h, gap_open, gap_extension; String<char> matrix; double min_expect;
//
template<typename TInfoSpec = BasicInfo>
struct StoreReport;		//stores the whole report


/**
.Spec.StreamReport:
..cat:Blast
..general:Class.BlastReport
..summary:BlastReport specialization that works on a file stream (parses hits/HSPs when iterating over them).
..signature:BlastReport<TBlastHsp,StreamReport<TFile> >
..param.TBlastHsp:The type of HSPs to be stored. See @Class.BlastHsp@
...metafunction:Metafunction.Hsp
...default:BlastHsp<BlastN,BasicInfo> 
..param.TFile:The type of the stream.
...default:std::fstream
..include:blast.h
*/
template<typename TFile = std::fstream>    //works on a stream
struct StreamReport;




//////////////////////////////////////////////////////////////////////////////
//Blast Meta functions


/**
.Metafunction.Hit:
..cat:Blast
..summary:Blast Hit type of a Blast object.
..signature:Hsp<T>::Type
..param.T:A Blast report object.
...type:Class.BlastReport
..returns.param.Type:BlastHit type.
..include:seqan/blast.h
*/
template<typename T>
struct Hit;

/**
.Metafunction.Hsp:
..cat:Blast
..summary:Blast HSP type of a Blast object.
..signature:Hsp<T>::Type
..param.T:A Blast object.
...type:Class.BlastReport
...type:Class.BlastHit
..returns.param.Type:BlastHsp type.
..include:seqan/blast.h
*/
template<typename T>
struct Hsp;


//////////////////////////////////////////////////////////////////////////////
// Blast Tag

struct TagBlast_;
typedef Tag<TagBlast_> const Blast;


//////////////////////////////////////////////////////////////////////////////
// Blat Tag 
//already defined in seeds/seedHandlingTags.h (49)
//moved to basic_tag
/*
struct Blat_;
typedef Tag<Blat_> const Blat;
*/
//////////////////////////////////////////////////////////////////////////////



/**
.Spec.BlastN:
..cat:Blast
..general:Class.BlastHsp
..summary:For BlastN Blast reports.
..signature:BlastN
..include:blast.h
*/


/**
.Spec.BlastP:
..cat:Blast
..general:Class.BlastHsp
..summary:For BlastP Blast reports.
..signature:BlastP
..include:blast.h
*/



struct TagBlastN_;
struct TagMegaBlast_;
struct TagBlastP_;
struct TagBlastX_;
struct TagTBlastN_;
struct TagTBlastX_;

template <typename TSpec = TagBlastN_>
class NucleotideBlast {};

template <typename TSpec = TagBlastP_>
class ProteinBlast {};

typedef NucleotideBlast<TagBlastN_> BlastN;
typedef NucleotideBlast<TagMegaBlast_> MegaBlast;

typedef ProteinBlast<TagBlastP_> BlastP;
typedef ProteinBlast<TagBlastX_> BlastX;
typedef ProteinBlast<TagTBlastN_> TBlastN;
typedef ProteinBlast<TagTBlastX_> TBlastX;



}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
