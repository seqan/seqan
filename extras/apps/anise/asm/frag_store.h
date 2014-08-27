// ==========================================================================
//                                   ANISE
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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

#ifndef SANDBOX_HERBARIUM_APPS_ANISE_ANISE_FRAG_STORE_H_
#define SANDBOX_HERBARIUM_APPS_ANISE_ANISE_FRAG_STORE_H_

#include <seqan/store.h>

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class MyStoreConfig
// ----------------------------------------------------------------------------

// Configuration for the FragmentStore we use for visualizing the multi-read alignments of the contigs.

struct MyStoreConfig
{
	typedef seqan::String<seqan::Dna5Q>	TReadSeq;
	typedef seqan::String<seqan::Dna5Q>	TContigSeq;
	
	typedef double			TMean;
	typedef double			TStd;
	typedef signed char		TMappingQuality;
		
	typedef void					TReadStoreElementSpec;
	typedef seqan::Owner<>                 TReadSeqStoreSpec;
	typedef void					TMatePairStoreElementSpec;
	typedef void					TLibraryStoreElementSpec;
	typedef void					TContigStoreElementSpec;
	typedef void					TContigFileSpec;
	typedef void					TAlignedReadStoreElementSpec;
	typedef seqan::Owner<seqan::ConcatDirect<> >	TAlignedReadTagStoreSpec;
	typedef void					TAnnotationStoreElementSpec;

    typedef seqan::Alloc<>					TReadNameSpec;
	typedef seqan::Owner<>	TReadNameStoreSpec;
};

// ----------------------------------------------------------------------------
// Typedef TFragmentStore and related
// ----------------------------------------------------------------------------

typedef seqan::FragmentStore<void, MyStoreConfig>                                 TFragmentStore;
typedef seqan::Value<TFragmentStore::TAlignedReadStore>::Type                     TAlignedReadStoreElement;
typedef seqan::Iterator<TFragmentStore::TAlignedReadStore, seqan::Standard>::Type TAlignedReadIter;

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

#endif  // #ifndef SANDBOX_HERBARIUM_APPS_ANISE_ANISE_FRAG_STORE_H_
