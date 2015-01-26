// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
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

#ifndef SEQAN_HEADER_INDEX_QGRAM_H
#define SEQAN_HEADER_INDEX_QGRAM_H

namespace seqan
{

	template <typename TText, typename TShapeSpec, typename TSpec>
	inline bool indexCreate(
		Index<TText, IndexQGram<MinimizerShape<...>, TSpec> > &index,
		FibreSADir, 
		Default const) 
	{		
		resize(indexSA(index), _qgramQGramCount(index), Exact());
		resize(indexDir(index), _fullDirLength(index), Exact());
		createQGramIndex(index);
		resize(indexSA(index), back(indexDir(index)), Exact());     // shrink if some buckets were disabled
        _refineQGramIndex(index, ...);
		return true;
	}



	template < typename TText, typename TShapeSpec, typename TSpec, typename TShapeSpec2, typename TValue >
	inline typename Infix< typename Fibre< Index< TText, IndexQGram<TShapeSpec, TSpec> >, FibreSA>::Type const >::Type 
	getOccurrences(
		Index< TText, IndexQGram<TShapeSpec, TSpec> > const &index,
		Shape< TValue, TShapeSpec2 > const &shape) 
	{
		typedef typename Size<typename Fibre< Index< TText, IndexQGram<TShapeSpec, TSpec> >, FibreDir>::Type>::Type TDirSize;
		TDirSize bucket = getBucket(indexBucketMap(index), value(shape));
		TSA sa = infix(indexSA(index), indexDir(index)[bucket], indexDir(index)[bucket + 1]);

        TSAIterator saBegin = begin(sa, Standard()) + value(it).range.i1;
        TSASize saLen = isRoot(it) ? length(sa) : value(it).range.i2 - value(it).range.i1;
        TSearchTreeIterator node(saBegin, saLen);
        Pair<TSAIterator> range = _equalRangeSA(text, node, pattern, value(it).repLen);

        return range;
	}	

	template < typename TText, typename TShapeSpec, typename TSpec, typename TShapeSpec2, typename TValue >
	inline typename Infix< typename Fibre< Index< TText, IndexQGram<TShapeSpec, TSpec> >, FibreSA>::Type const >::Type 
	getOccurrences(
		Index< TText, IndexQGram<TShapeSpec, TSpec> > &index,
		Shape< TValue, TShapeSpec2 > const &shape) 
	{
		indexRequire(index, QGramSADir());
		return getOccurrences(const_cast<Index< TText, IndexQGram<TShapeSpec, TSpec> > const &>(index), shape);
	}	


    //TODO(singer): Why Size of FibreSA and not Size<Index>?
	template < typename TText, typename TShapeSpec, typename TSpec, typename TShapeSpec2, typename TValue >
	inline typename Size< typename Fibre< Index< TText, IndexQGram<TShapeSpec, TSpec> >, FibreSA>::Type const >::Type 
	countOccurrences(
		Index< TText, IndexQGram<TShapeSpec, TSpec> > const &index,
		Shape< TValue, TShapeSpec2 > const &shape) 
	{
		typedef typename Size<typename Fibre< Index< TText, IndexQGram<TShapeSpec, TSpec> >, FibreDir>::Type>::Type TDirSize;
		TDirSize bucket = getBucket(indexBucketMap(index), value(shape));
		return indexDir(index)[bucket + 1] - indexDir(index)[bucket];
	}	

	template < typename TText, typename TShapeSpec, typename TSpec, typename TShapeSpec2, typename TValue >
	inline typename Size< typename Fibre< Index< TText, IndexQGram<TShapeSpec, TSpec> >, FibreSA>::Type const >::Type 
	countOccurrences(
		Index< TText, IndexQGram<TShapeSpec, TSpec> > &index,
		Shape< TValue, TShapeSpec2 > const &shape) 
	{
		indexRequire(index, QGramDir());
		return countOccurrences(const_cast<Index< TText, IndexQGram<TShapeSpec, TSpec> > const &>(index), shape);
	}

}

#endif //#ifndef SEQAN_HEADER_...
