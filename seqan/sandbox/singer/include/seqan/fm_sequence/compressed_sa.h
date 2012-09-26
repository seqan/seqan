// // ==========================================================================
// //                                  FMIndex
// // ==========================================================================
// // Copyright (c) 2006-2012, Knut Reinert, FU Berlin
// // All rights reserved.
// //
// // Redistribution and use in source and binary forms, with or without
// // modification, are permitted provided that the following conditions are met:
// //
// //     * Redistributions of source code must retain the above copyright
// //       notice, this list of conditions and the following disclaimer.
// //     * Redistributions in binary form must reproduce the above copyright
// //       notice, this list of conditions and the following disclaimer in the
// //       documentation and/or other materials provided with the distribution.
// //     * Neither the name of Knut Reinert or the FU Berlin nor the names of
// //       its contributors may be used to endorse or promote products derived
// //       from this software without specific prior written permission.
// //
// // THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// // AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// // IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// // ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// // FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// // DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// // SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// // CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// // LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// // OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// // DAMAGE.
// //
// // ==========================================================================
// // Author: Jochen Singer <jochen.singer@fu-berlin.de>
// // ==========================================================================
// 
// #ifndef SANDBOX_MY_SANDBOX_APPS_FMINDEX_COMPRESSEDSA_H_
// #define SANDBOX_MY_SANDBOX_APPS_FMINDEX_COMPRESSEDSA_H_
// 
// namespace seqan {
// 
// template <typename TOccTable, typename TPrefixSumTable>
// struct LFTable
// {
//     TOccTable occTable;
//     TPrefixSumTable prefixSumTable;
// 
//     LFTable() :
//     	occTable(),
//     	prefixSumTable()
//     {}
// 
//     LFTable(TOccTable occTable, TPrefixSumTable prefixSumTable) :
//         occTable(occTable),
//         prefixSumTable(prefixSumTable)
//     {}
// 
//     inline LFTable & operator=(LFTable const & other)
//     {
//     	occTable = other.occTable;
//     	prefixSumTable = other.prefixSumTable;
//     	return *this;
//     }
// 
//     inline bool operator==(const LFTable & b) const
//     {
//         return occTable == b.occTable &&
//                prefixSumTable == b.prefixSumTable;
//     }
// 
// };
// 
// struct FibreOccTable_;
// struct FibrePrefixSumTable_;
// 
// typedef Tag<FibreOccTable_> const FibreOccTable;
// typedef Tag<FibrePrefixSumTable_> const FibrePrefixSumTable;
// 
// template <typename TOccTable, typename TPrefixSumTable>
// struct Fibre<LFTable<TOccTable, TPrefixSumTable>, FibreOccTable>
// {
//     typedef TOccTable Type;
// };
// 
// template <typename TOccTable, typename TPrefixSumTable>
// struct Fibre<LFTable<TOccTable, TPrefixSumTable> const, FibreOccTable>
// {
//     typedef TOccTable const Type;
// };
// 
// template <typename TOccTable, typename TPrefixSumTable>
// struct Fibre<LFTable<TOccTable, TPrefixSumTable>, FibrePrefixSumTable>
// {
//     typedef TPrefixSumTable Type;
// };
// 
// template <typename TOccTable, typename TPrefixSumTable>
// inline typename Fibre<LFTable<TOccTable, TPrefixSumTable>, FibrePrefixSumTable>::Type const &
// getFibre(LFTable<TOccTable, TPrefixSumTable> const & lfTable, FibrePrefixSumTable)
// {
//     return lfTable.prefixSumTable;
// }
// 
// template <typename TOccTable, typename TPrefixSumTable>
// inline typename Fibre<LFTable<TOccTable, TPrefixSumTable>, FibrePrefixSumTable>::Type &
// getFibre(LFTable<TOccTable, TPrefixSumTable> & lfTable, FibrePrefixSumTable)
// {
//     return lfTable.prefixSumTable;
// }
// 
// template <typename TOccTable, typename TPrefixSumTable>
// inline typename Fibre<LFTable<TOccTable, TPrefixSumTable>, FibreOccTable>::Type &
// getFibre(LFTable<TOccTable, TPrefixSumTable> & lfTable, FibreOccTable)
// {
// 	return lfTable.occTable;
// }
// 
// template <typename TOccTable, typename TPrefixSumTable>
// inline typename Fibre<LFTable<TOccTable, TPrefixSumTable>, FibreOccTable>::Type const &
// getFibre(LFTable<TOccTable, TPrefixSumTable> const & lfTable, FibreOccTable)
// {
//     return lfTable.occTable;
// }
// 
// template <typename TOccTable, typename TPrefixSumTable>
// inline void clear(LFTable<TOccTable, TPrefixSumTable> & lfTable)
// {
//     clear(lfTable.occTable);
//     clear(lfTable.prefixSumTable);
// }
// 
// template <typename TOccTable, typename TPrefixSumTable>
// inline bool open(
//     LFTable<TOccTable, TPrefixSumTable> & lfTable,
//     const char * fileName,
//     int openMode)
// {
//     String<char> name;
//     name = fileName;    append(name, ".occ");
//     if (!open(getFibre(lfTable, FibreOccTable()), toCString(name), openMode))
//     {
//         return false;
//     }
//     name = fileName;    append(name, ".psum");  open(getFibre(lfTable, FibrePrefixSumTable()), toCString(name), openMode);
//     return true;
// 
// }
// 
// template <typename TOccTable, typename TPrefixSumTable>
// inline bool open(
//     LFTable<TOccTable, TPrefixSumTable> & lfTable,
//     const char * fileName)
// {
//     return open(lfTable, fileName, DefaultOpenMode<LFTable<TOccTable, TPrefixSumTable> >::VALUE);
// }
// 
// template <typename TOccTable, typename TPrefixSumTable>
// inline bool save(
//     LFTable<TOccTable, TPrefixSumTable> const & lfTable,
//     const char * fileName,
//     int openMode)
// {
//     String<char> name;
//     name = fileName;    append(name, ".occ");
//     if (!save(getFibre(lfTable, FibreOccTable()), toCString(name), openMode))
//     {
//         return false;
//     }
//     name = fileName;    append(name, ".psum");  save(getFibre(lfTable, FibrePrefixSumTable()), toCString(name), openMode);
//     return true;
// 
// }
// 
// template <typename TOccTable, typename TPrefixSumTable>
// inline bool save(
//     LFTable<TOccTable, TPrefixSumTable> const & lfTable,
//     const char * fileName)
// {
//     return save(lfTable, fileName, DefaultOpenMode<LFTable<TOccTable, TPrefixSumTable> >::VALUE);
// }
// 
// template <typename TSparseString, typename TLfTable, typename TSpec>
// struct CompressedSA
// {
//     TSparseString 	compressedSA;
//     TLfTable * 		lfTable;
// 
//     CompressedSA() :
//     	compressedSA(),
//     	lfTable()
//     {}
// 
//     inline CompressedSA & operator=(CompressedSA const & other)
//     {
//     	compressedSA = other.compressedSA;
//     	lfTable = other.lfTable;
//     	return *this;
//     }
// 
//     typedef typename Value<typename Fibre<TSparseString, FibreSparseString>::Type>::Type TCompressedSaValue;
//     typedef typename Fibre<TSparseString, FibreIndicatorString>::Type TIndicatorString;
// 
//     template <typename TPos>
//     inline TCompressedSaValue const operator[](TPos pos)
//     {
//     	TIndicatorString const & indicatorString = getFibre(compressedSA, FibreIndicatorString());
//     	TPos counter = 0;
// 
//     	while (!getBit(indicatorString, pos))
//     	{
//     		pos = lfMapping(*lfTable, pos);
//     		++counter;
//     	}
//     	return getValue(compressedSA, getRank(indicatorString, pos) - 1) + counter;
//     }
// 
//     template <typename TPos>
//     inline TCompressedSaValue operator[](TPos pos) const
//     {
//         TIndicatorString const & indicatorString = getFibre(compressedSA, FibreIndicatorString());
//         TPos counter = 0;
//         while (!getBit(indicatorString, pos))
//         {
//             pos = lfMapping(*lfTable, pos);
//             ++counter;
//         }
//         return getValue(compressedSA, getRank(indicatorString, pos) - 1) + counter;
//     }
// 
//     inline bool operator==(const CompressedSA & b) const
//     {
//         return compressedSA == b.compressedSA &&
//                *lfTable == *(b.lfTable);
//     }
// 
// };
// 
// template <typename TSpecPairI1, typename TSpecPairI2, typename TSpecPairSpec, typename TStringSpec, typename TSparseStringSpec, typename TLfTable, typename TSpec>
// struct CompressedSA<SparseString<String<Pair<TSpecPairI1, TSpecPairI2, TSpecPairSpec>, TStringSpec>, TSparseStringSpec>, TLfTable, TSpec>
// {
//     typedef SparseString<String<Pair<TSpecPairI1, TSpecPairI2, TSpecPairSpec>, TStringSpec>, TSparseStringSpec>   TSparseString;
//     typedef typename Value<typename Fibre<SparseString<TSparseString, TStringSpec>, FibreSparseString>::Type>::Type TCompressedSaValue;
//     typedef typename Fibre<SparseString<TSparseString, TStringSpec>, FibreIndicatorString>::Type TIndicatorString;
//     TSparseString                                                                                       compressedSA;
//     TLfTable * lfTable;
// 
//     CompressedSA() :
//     	compressedSA(),
//     	lfTable()
//     {}
// 
//     template <typename TPos>
//     inline TCompressedSaValue const operator[](TPos pos)
//     {
//         TIndicatorString const & indicatorString = getFibre(compressedSA, FibreIndicatorString());
//         TPos counter = 0;
//         while (!getBit(indicatorString, pos))
//         {
//             pos = lfMapping(*lfTable, pos);
//             ++counter;
//         }
//         TCompressedSaValue temp = getValue(compressedSA, getRank(indicatorString, pos) - 1);
//         temp.i2 += counter;
//         return temp;
//     }
// 
//     template <typename TPos>
//     inline TCompressedSaValue operator[](TPos pos) const
//     {
//         TIndicatorString const & indicatorString = getFibre(compressedSA, FibreIndicatorString());
//         TPos counter = 0;
//         while (!getBit(indicatorString, pos))
//         {
//             pos = lfMapping(*lfTable, pos);
//             ++counter;
//         }
//         TCompressedSaValue temp = getValue(compressedSA, getRank(indicatorString, pos) - 1);
//         temp.i2 += counter;
//         return temp;
//     }
// 
//     inline bool operator==(const CompressedSA & b) const
//     {
//         return compressedSA == b.compressedSA &&
//                *lfTable == *(b.lfTable);
//     }
// 
// };
// 
// template <typename TSparseString, typename TLfTable, typename TSpec>
// struct Fibre<CompressedSA<TSparseString, TLfTable, TSpec>, FibreSA>
// {
//     typedef TSparseString Type;
// };
// 
// 
// template <typename TSparseString, typename TLfTable, typename TSpec>
// inline typename Fibre<CompressedSA<TSparseString, TLfTable, TSpec>, FibreSA>::Type const &
// getFibre(CompressedSA<TSparseString, TLfTable, TSpec> const & compressedSA, FibreSA)
// {
//     return compressedSA.compressedSA;
// }
// 
// template <typename TSparseString, typename TLfTable, typename TSpec>
// inline typename Fibre<CompressedSA<TSparseString, TLfTable, TSpec>, FibreSA>::Type &
// getFibre(CompressedSA<TSparseString, TLfTable, TSpec> & compressedSA, FibreSA)
// {
//     return compressedSA.compressedSA;
// }
// 
// template <typename TSparseString, typename TLfTable, typename TSpec>
// inline bool open(
//     CompressedSA<TSparseString, TLfTable, TSpec> & compressedSA,
//     const char * fileName,
//     int openMode)
// {
//     String<char> name;
//     name = fileName;    append(name, ".sstring");
//     if (!open(getFibre(compressedSA, FibreSA()), toCString(name), openMode))
//     {
//         return false;
//     }
//     return true;
// 
// }
// 
// template <typename TSparseString, typename TLfTable, typename TSpec>
// inline bool open(
//     CompressedSA<TSparseString, TLfTable, TSpec> & compressedSA,
//     const char * fileName)
// {
//     return open(compressedSA, fileName, DefaultOpenMode<CompressedSA<TSparseString, TLfTable, TSpec> >::VALUE);
// }
// 
// template <typename TSparseString, typename TLfTable, typename TSpec>
// inline bool save(
//     CompressedSA<TSparseString, TLfTable, TSpec> const & compressedSA,
//     const char * fileName,
//     int openMode)
// {
//     String<char> name;
//     name = fileName;    append(name, ".sstring");
//     if (!save(getFibre(compressedSA, FibreSA()), toCString(name), openMode))
//     {
//         return false;
//     }
//     return true;
// }
// 
// template <typename TSparseString, typename TLfTable, typename TSpec, typename TPos>
// inline bool getNextPos(CompressedSA<TSparseString, TLfTable, TSpec> const & compressedSA, TPos & pos)
// {
//     typedef typename Fibre<TSparseString, FibreIndicatorString>::Type TIndicatorString;
// 	TIndicatorString const & indicatorString = compressedSA.compressedSA.indicatorString;
// 
// 	if (getBit(indicatorString, pos))
// 	{
// 		return true;
// 	}
// 	pos = lfMapping(*compressedSA.lfTable, pos);
// 	return false;
// }
// 
// template <typename TSparseString, typename TLfTable, typename TSpec>
// inline bool save(
//     CompressedSA<TSparseString, TLfTable, TSpec> const & compressedSA,
//     const char * fileName)
// {
//     return save(compressedSA, fileName, DefaultOpenMode<CompressedSA<TSparseString, TLfTable, TSpec> >::VALUE);
// }
// 
// template <typename TSparseString, typename TLfTable, typename TSpec>
// struct Iterator<CompressedSA<TSparseString, TLfTable, TSpec> const, Standard>
// {
//     typedef Iter<CompressedSA<TSparseString, TLfTable, TSpec> const, PositionIterator> Type;
// };
// 
// template <typename TSparseString, typename TLfTable, typename TSpec>
// struct Iterator<CompressedSA<TSparseString, TLfTable, TSpec>, Standard>
// {
//     typedef Iter<CompressedSA<TSparseString, TLfTable, TSpec>, PositionIterator> Type;
// };
// 
// template <typename TSparseString, typename TLfTable, typename TSpec>
// struct Iterator<CompressedSA<TSparseString, TLfTable, TSpec>, Rooted>:
//     Iterator<CompressedSA<TSparseString, TLfTable, TSpec>, Standard>{};
// 
// template <typename TSparseString, typename TLfTable, typename TSpec>
// struct Iterator<CompressedSA<TSparseString, TLfTable, TSpec> const, Rooted>:
//     Iterator<CompressedSA<TSparseString, TLfTable, TSpec> const, Standard>{};
// 
// template <typename TSparseString, typename TLfTable, typename TSpec>
// struct Value<CompressedSA<TSparseString, TLfTable, TSpec> >
// {
//     typedef typename Value<TSparseString>::Type Type;
// };
// 
// template <typename TSparseString, typename TLfTable, typename TSpec>
// struct Value<CompressedSA<TSparseString, TLfTable, TSpec> const>
// {
//     typedef typename Value<TSparseString>::Type const Type;
// };
// 
// }
// 
// 
// #endif // SANDBOX_MY_SANDBOX_APPS_FMINDEX_COMPRESSEDSA_H_
