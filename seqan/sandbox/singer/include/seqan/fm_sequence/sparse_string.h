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

#ifndef SANDBOX_MY_SANDBOX_APPS_FMINDEX_SPARSE_STRING_H_
#define SANDBOX_MY_SANDBOX_APPS_FMINDEX_SPARSE_STRING_H_

namespace seqan {

//forward declaration
// template <typename TSparseString, typename TSpec>
// struct SparseString;
// 
// struct FibreSparseString_;
// struct FibreIndicatorString_;
// 
// typedef Tag<FibreSparseString_>       const FibreSparseString;
// typedef Tag<FibreIndicatorString_>    const FibreIndicatorString;
// 
// template <typename TSparseString, typename TSpec>
// struct Fibre<SparseString<TSparseString, TSpec>, FibreSparseString>
// {
//     typedef TSparseString Type;
// };
// 
// template <typename TSparseString, typename TSpec>
// struct Fibre<SparseString<TSparseString, TSpec>, FibreIndicatorString>
// {
//     typedef RankSupportBitString<void> Type;
// };
// 
// template <typename TSparseString, typename TSpec>
// struct SparseString
// {
// 
//     typedef typename Fibre<SparseString, FibreIndicatorString>::Type TIndicatorString;
// 
//     TSparseString                           string;
//     TIndicatorString                        indicatorString;
//     typename Size<TSparseString>::Type    blockSize;
// 
//     SparseString() :
//     	string(),
//     	indicatorString(),
//     	blockSize(0)
//     {}
// 
//     inline SparseString & operator=(SparseString const & other)
//     {
//     	string = other.string;
//     	indicatorString = other.indicatorString;
//     	blockSize = other.blockSize;
//     	return *this;
//     }
// 
//     inline bool operator==(const SparseString & b) const
//     {
//         return string == b.string &&
//                indicatorString == b.indicatorString &&
//                //blockSize == b.blockSize);
//                1;
//     }
// 
// };
// template <typename TSparseString, typename TSpec>
// inline typename Fibre<SparseString<TSparseString, TSpec>, FibreSparseString>::Type &
// getFibre(SparseString<TSparseString, TSpec> & sparseString, FibreSparseString)
// {
//     return sparseString.string;
// }
// 
// template <typename TSparseString, typename TSpec>
// inline typename Fibre<SparseString<TSparseString, TSpec>, FibreIndicatorString>::Type &
// getFibre(SparseString<TSparseString, TSpec> & sparseString, FibreIndicatorString)
// {
//     return sparseString.indicatorString;
// }
// 
// template <typename TSparseString, typename TSpec>
// inline typename Fibre<SparseString<TSparseString, TSpec>, FibreSparseString>::Type const &
// getFibre(SparseString<TSparseString, TSpec> const & sparseString, FibreSparseString)
// {
//     return sparseString.string;
// }
// 
// template <typename TSparseString, typename TSpec>
// inline typename Fibre<SparseString<TSparseString, TSpec>, FibreIndicatorString>::Type const &
// getFibre(SparseString<TSparseString, TSpec> const & sparseString, FibreIndicatorString)
// {
//     return sparseString.indicatorString;
// }
// 
// template <typename TSparseString, typename TSpec>
// struct Iterator<SparseString<TSparseString, TSpec> const, Standard>
// {
//     typedef Iter<SparseString<TSparseString, TSpec> const, PositionIterator> Type;
// };
// 
// template <typename TSparseString, typename TSpec>
// struct Iterator<SparseString<TSparseString, TSpec>, Standard>
// {
//     typedef Iter<SparseString<TSparseString, TSpec>, PositionIterator> Type;
// };
// 
// template <typename TSparseString, typename TSpec>
// struct Iterator<SparseString<TSparseString, TSpec>, Rooted>:
//     Iterator<SparseString<TSparseString, TSpec>, Standard>{};
// 
// template <typename TSparseString, typename TSpec>
// struct Iterator<SparseString<TSparseString, TSpec> const, Rooted>:
//     Iterator<SparseString<TSparseString, TSpec> const, Standard>{};
// 
// template <typename TSparseString, typename TSpec>
// struct Value<SparseString<TSparseString, TSpec> >
// {
//     typedef typename Value<TSparseString>::Type Type;
// };
// 
// template <typename TSparseString, typename TSpec>
// struct Value<SparseString<TSparseString, TSpec> const>
// {
//     typedef typename Value<TSparseString>::Type const Type;
// };
// 
// template <typename TSparseString, typename TSpec, typename TValue>
// inline void assignBlockSize(SparseString<TSparseString, TSpec> & string, TValue value)
// {
//     //std::cerr << "sparse" << std::endl;
//     string.blockSize = value;
// }
// 
// template <typename TSparseString, typename TSpec, typename TPos, typename TValue>
// inline void assignValue(SparseString<TSparseString, TSpec> & string, TPos pos, TValue value)
// {
//     string.string[pos] = value;
// }
// 
// template <typename TSparseString, typename TSpec>
// inline typename Size<typename Fibre<SparseString<TSparseString, TSpec>, FibreSparseString>::Type>::Type
// getBlockSize(SparseString<TSparseString, TSpec> & string)
// {
//     return string.blockSize;
// }
// 
// template <typename TSparseString, typename TSpec, typename TPos>
// inline typename Value<typename Fibre<SparseString<TSparseString, TSpec>, FibreSparseString>::Type>::Type
// getValue(SparseString<TSparseString, TSpec> & string, TPos pos)
// {
//     return getValue(string.string, pos);
// }
// 
// template <typename TSparseString, typename TSpec, typename TPos>
// inline typename Value<typename Fibre<SparseString<TSparseString, TSpec> const, FibreSparseString>::Type>::Type
// getValue(SparseString<TSparseString, TSpec> const & string, TPos pos)
// {
//     return getValue(string.string, pos);
// }
// 
// template <typename TSparseString, typename TSpec>
// inline typename Size<typename Fibre<SparseString<TSparseString, TSpec>, FibreSparseString>::Type>::Type
// length(SparseString<TSparseString, TSpec> & string)
// {
//     return length(string.string);
// }
// 
// template <typename TSparseString, typename TSpec>
// inline void clear(SparseString<TSparseString, TSpec> & string)
// {
//     clear(string.string);
//     clear(string.indicatorString);
// }
// 
// template <typename TSparseString, typename TSpec, typename TSize>
// inline void resize(SparseString<TSparseString, TSpec> & string,
//             TSize size)
// {
//     resize(getFibre(string, FibreSparseString()), size / getBlockSize(string) + 1);
//     resize(getFibre(string, FibreIndicatorString()), size, 0);
// }
// 
// template <typename TSparseString, typename TSpec, typename TPos>
// inline typename Value<typename Fibre<SparseString<TSparseString, TSpec>, FibreSparseString>::Type>::Type
// value(SparseString<TSparseString, TSpec> & string, TPos pos)
// {
//     return string.string[pos];
// }
// 
// template <typename TSparseString, typename TSpec, typename TPos>
// inline typename Value<typename Fibre<SparseString<TSparseString, TSpec>, FibreSparseString>::Type>::Type const
// value(SparseString<TSparseString, TSpec> const & string, TPos pos)
// {
//     return string.string[pos];
// }
// 
// template <typename TSparseString, typename TSpec>
// inline bool open(
//     SparseString<TSparseString, TSpec> & sparseString,
//     const char * fileName,
//     int openMode)
// {
//     String<char> name;
//     name = fileName;    append(name, ".sstring");
//     if (!open(getFibre(sparseString, FibreSparseString()), toCString(name), openMode))
//     {
//         return false;
//     }
//     name = fileName;    append(name, ".istring");   open(getFibre(sparseString, FibreIndicatorString()), toCString(name), openMode);
//     return true;
// }
// 
// template <typename TSparseString, typename TSpec>
// inline bool open(
//     SparseString<TSparseString, TSpec> & sparseString,
//     const char * fileName)
// {
//     return open(sparseString, fileName, DefaultOpenMode<SparseString<TSparseString, TSpec> >::VALUE);
// }
// 
// template <typename TSparseString, typename TSpec>
// inline bool save(
//     SparseString<TSparseString, TSpec> const & sparseString,
//     const char * fileName,
//     int openMode)
// {
//     String<char> name;
//     name = fileName;    append(name, ".sstring");
//     if (!save(getFibre(sparseString, FibreSparseString()), toCString(name), openMode))
//     {
//         return false;
//     }
//     name = fileName;    append(name, ".istring");   save(getFibre(sparseString, FibreIndicatorString()), toCString(name), openMode);
//     return true;
// }
// 
// template <typename TSparseString, typename TSpec>
// inline bool save(
//     SparseString<TSparseString, TSpec> const & sparseString,
//     const char * fileName)
// {
//     return save(sparseString, fileName, DefaultOpenMode<SparseString<TSparseString, TSpec> >::VALUE);
// }

}
#endif // SANDBOX_MY_SANDBOX_APPS_FMINDEX_COMPRESSEDSA_H_
