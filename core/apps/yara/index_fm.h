// ==========================================================================
//                      Yara - Yet Another Read Aligner
// ==========================================================================
// Copyright (c) 2011-2014, Enrico Siragusa, FU Berlin
// Copyright (c) 2013 NVIDIA Corporation
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
//     * Neither the name of Enrico Siragusa or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ENRICO SIRAGUSA OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================
// This file contains Index overloads.
// ==========================================================================

#ifndef APP_YARA_INDEX_FM_H_
#define APP_YARA_INDEX_FM_H_

namespace seqan {

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function indexRequire()
// ----------------------------------------------------------------------------
// This function is overloaded to avoid compiling the code building the index.

#ifndef YARA_INDEXER
template <typename TText, typename TSSetSpec, typename TSpec, typename TConfig, typename TFibre>
inline bool indexRequire(Index<StringSet<TText, TSSetSpec>,  FMIndex<TSpec, TConfig> > & index, Tag<TFibre> const fibre)
{
    if (!indexSupplied(index, fibre))
        throw RuntimeError("The reference index file was not loaded correctly.");

    return true;
}
#endif

// ----------------------------------------------------------------------------
// Function save()
// ----------------------------------------------------------------------------
// This function is overloaded to avoid saving the text.

template <typename TText, typename TSSetSpec, typename TSpec, typename TConfig>
inline bool save(Index<StringSet<TText, TSSetSpec>, FMIndex<TSpec, TConfig> > const & index,
                 const char * fileName, int openMode)
{
    String<char> name;

    name = fileName;    append(name, ".sa");
    if (!save(getFibre(index, FibreSA()), toCString(name), openMode)) return false;

    name = fileName;    append(name, ".lf");
    if (!save(getFibre(index, FibreLF()), toCString(name), openMode)) return false;

    return true;
}

// ----------------------------------------------------------------------------
// Function open()
// ----------------------------------------------------------------------------
// This function is overloaded to avoid loading the text.

template <typename TText, typename TSSetSpec, typename TSpec, typename TConfig>
inline bool open(Index<StringSet<TText, TSSetSpec>, FMIndex<TSpec, TConfig> > & index,
                 const char * fileName, int openMode)
{
    String<char> name;

    name = fileName;    append(name, ".sa");
    if (!open(getFibre(index, FibreSA()), toCString(name), openMode)) return false;

    name = fileName;    append(name, ".lf");
    if (!open(getFibre(index, FibreLF()), toCString(name), openMode)) return false;

    setFibre(getFibre(index, FibreSA()), getFibre(index, FibreLF()), FibreLF());

    return true;
}

// ----------------------------------------------------------------------------
// Function _getNodeByChar()
// ----------------------------------------------------------------------------
// This function is overloaded for Dna5 and Dna5Q to let Ns always mismatch.

template <typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec>
SEQAN_HOST_DEVICE inline bool
_getNodeByChar(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec> >, VSTree<TopDown<TSpec> > > const & it,
               typename VertexDescriptor<Index<TText, FMIndex<TOccSpec, TIndexSpec> > >::Type const & vDesc,
               Pair<typename Size<Index<TText, FMIndex<TOccSpec, TIndexSpec> > >::Type> & _range,
               Dna5 c)
{
    return _getNodeByCharImpl(it, vDesc, _range, c);
}

template <typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec>
SEQAN_HOST_DEVICE inline bool
_getNodeByChar(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec> >, VSTree<TopDown<TSpec> > > const & it,
               typename VertexDescriptor<Index<TText, FMIndex<TOccSpec, TIndexSpec> > >::Type const & vDesc,
               Pair<typename Size<Index<TText, FMIndex<TOccSpec, TIndexSpec> > >::Type> & _range,
               Dna5Q c)
{
    return _getNodeByCharImpl(it, vDesc, _range, c);
}

template <typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec, typename TChar>
SEQAN_HOST_DEVICE inline bool
_getNodeByCharImpl(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec> >, VSTree<TopDown<TSpec> > > const & it,
               typename VertexDescriptor<Index<TText, FMIndex<TOccSpec, TIndexSpec> > >::Type const & vDesc,
               Pair<typename Size<Index<TText, FMIndex<TOccSpec, TIndexSpec> > >::Type> & _range,
               TChar c)
{
    typedef Index<TText, FMIndex<TOccSpec, TIndexSpec> >        TIndex;
    typedef typename Fibre<TIndex, FibreLF>::Type               TLF;
    typedef typename Value<TIndex>::Type                        TAlphabet;

    TIndex const & index = container(it);
    TLF const & lf = indexLF(index);

    if (ordValue(c) >= ValueSize<TAlphabet>::VALUE) return false;

    _range = range(index, vDesc);

    _range.i1 = lf(_range.i1, c);
    _range.i2 = lf(_range.i2, c);

    return _range.i1 < _range.i2;
}

}

#endif  // #ifndef APP_YARA_INDEX_FM_H_
