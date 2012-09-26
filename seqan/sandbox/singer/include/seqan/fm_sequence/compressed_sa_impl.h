// ==========================================================================
//                                  FMIndex
// ==========================================================================
// Copyright (c) 2006-2011, Knut Reinert, FU Berlin
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
// Author: Jochen Singer <jochen.singer@fu-berlin.de>
// ==========================================================================

#ifndef SANDBOX_MY_SANDBOX_APPS_FMINDEX_COMPRESSEDSAIMPL_H_
#define SANDBOX_MY_SANDBOX_APPS_FMINDEX_COMPRESSEDSAIMPL_H_

namespace seqan {

template <typename TSparseString, typename TLfTable, typename TSpec, typename TValue>
void assignBlockSize(CompressedSA<TSparseString, TLfTable, TSpec> & container, TValue value)
{
    assignBlockSize(getFibre(container, FibreSA()), value);
}

template <typename TSparseString, typename TLfTable, typename TSpec>
void assignLfTable(CompressedSA<TSparseString, TLfTable, TSpec> & container, TLfTable & lfTable)
{
    container.lfTable = &lfTable;
}

template <typename TSparseString, typename TLfTable, typename TSpec, typename TPos, typename TValue>
void assignValue(CompressedSA<TSparseString, TLfTable, TSpec> & container, TPos pos, TValue value)
{
    assignValue(container.compressedSA, pos, value);
}

template <typename TSparseString, typename TLfTable, typename TSpec>
inline typename Size<typename Fibre<CompressedSA<TSparseString, TLfTable, TSpec>, FibreSparseString>::Type>::Type
length(CompressedSA<TSparseString, TLfTable, TSpec> & compressedSA)
{
    return length(compressedSA.compressedSA);
}

template <typename TSparseString, typename TLfTable, typename TSpec>
inline void clear(CompressedSA<TSparseString, TLfTable, TSpec> & compressedSA)
{
    clear(compressedSA.compressedSA);
}

template <typename TSparseString, typename TLfTable, typename TSpec>
inline typename Size<typename Fibre<TSparseString, FibreSparseString>::Type>::Type
getBlockSize(CompressedSA<TSparseString, TLfTable, TSpec> & compressedSA)
{
    return getBlockSize(compressedSA.compressedSA);
}

template <typename TSparseString, typename TLfTable, typename TSpec, typename TSize>
inline void resize(CompressedSA<TSparseString, TLfTable, TSpec> & container, TSize size)
{
	//std::cerr << "resize: " << size << std::endl;
    resize(container.compressedSA, size);
}

template <typename TSparseString, typename TLfTable, typename TSpec, typename TSize>
inline void resize(CompressedSA<String<TSparseString>, TLfTable, TSpec> & container, TSize size)
{
	//std::cerr << "resize: " << size << std::endl;
    resize(container.compressedSA, size);
}

template <typename TSparseString, typename TLfTable, typename TSpec, typename TSize>
inline void reserve(CompressedSA<TSparseString, TLfTable, TSpec> & compressedSA,
             TSize size)
{
	//std::cerr << "reserve: " << size << std::endl;
    resize(compressedSA.compressedSA, size);
}

template <typename TSparseString, typename TLfTable, typename TSpec>
struct Reference<CompressedSA<TSparseString, TLfTable, TSpec> >
{
    typedef typename Value<CompressedSA<TSparseString, TLfTable, TSpec> >::Type Type;
};

template <typename TSparseString, typename TLfTable, typename TSpec>
struct Reference<const CompressedSA<TSparseString, TLfTable, TSpec> >
{
    typedef typename Value<CompressedSA<TSparseString, TLfTable, TSpec> >::Type const Type;
};

template <typename TSparseString, typename TLfTable, typename TSpec>
inline Iterator<CompressedSA<TSparseString, TLfTable, TSpec> >
begin(CompressedSA<TSparseString, TLfTable, TSpec> & compressedSA)
{
    return Iterator<CompressedSA<TSparseString, TLfTable, TSpec> >(compressedSA, 0);
}

template <typename TSparseString, typename TLfTable, typename TSpec>
inline Iterator<CompressedSA<TSparseString, TLfTable, TSpec> >
end(CompressedSA<TSparseString, TLfTable, TSpec> & compressedSA)
{
    return Iterator<CompressedSA<TSparseString, TLfTable, TSpec> >(compressedSA, length(compressedSA.compressedSA));
}

template <typename TSparseString, typename TLfTable, typename TSpec, typename TPos>
inline typename Value<typename Fibre<TSparseString, FibreSparseString>::Type>::Type
value(CompressedSA<TSparseString, TLfTable, TSpec> & compressedSA, TPos pos)
{
    return compressedSA[pos];
}

template <typename TSparseString, typename TLfTable, typename TSpec, typename TPos>
inline typename Value<typename Fibre<TSparseString, FibreSparseString>::Type>::Type const
value(const CompressedSA<TSparseString, TLfTable, TSpec> & compressedSA, TPos pos)
{
//		TSparseString test;
//		static_cast<Nothing>(test);
    return compressedSA[pos];
}

}

#endif // SANDBOX_MY_SANDBOX_APPS_FMINDEX_COMPRESSEDSA_H_
