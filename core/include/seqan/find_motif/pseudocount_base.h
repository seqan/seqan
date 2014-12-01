// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
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

#ifndef SEQAN_HEADER_PSEUDOCOUNT_BASE_H
#define SEQAN_HEADER_PSEUDOCOUNT_BASE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

/**
.Class.Pseudocount:
..summary:Holds the pseudocounts for each residue of a given sequence alphabet.
..cat:Motif Search
..signature:Pseudocount<TValue, TSpec>
..param.TValue:The type of sequence which is considered.
...metafunction:Metafunction.Value
...type:Spec.Dna
...type:Spec.AminoAcid
..param.TSpec:Specialization tag for determining the pseudocount method.
...type:Spec.CMode 
...type:Spec.PMode
..include:seqan/find_motif.h
*/

template<typename TValue, typename TSpec>
class Pseudocount;

//////////////////////////////////////////////////////////////////////////////
//Metafunctions
//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Value.param.T.type:Class.Pseudocount
///.Metafunction.Value.class:Class.Pseudocount

template<typename TValue, typename TSpec>
struct Value< Pseudocount<TValue, TSpec> >
{
	typedef TValue Type;
};
template<typename TValue, typename TSpec>
struct Value< Pseudocount<TValue, TSpec> const>
{
	typedef TValue const Type;
};

/////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
