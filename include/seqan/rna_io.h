// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2024, Knut Reinert, FU Berlin
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
// Authors: Lily Shellhammer <lily.shellhammer@gmail.com>
//          Joerg Winkler <j.winkler@fu-berlin.de>
// ==========================================================================
// Facade header for module rna_io.
// ==========================================================================

#ifndef INCLUDE_SEQAN_RNA_IO_H_
#define INCLUDE_SEQAN_RNA_IO_H_

// ===========================================================================
// Prerequisites.
// ===========================================================================

#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/sequence.h>
#include <seqan/journaled_set.h>
#include <seqan/graph_types.h>
// forEach loop
#include <seqan/parallel/parallel_algorithms.h>

// ===========================================================================
// RNA structure containers and file formats
// ===========================================================================

// containers
#include <seqan/rna_io/rna_record.h>
#include <seqan/rna_io/rna_header.h>
#include <seqan/rna_io/rna_io_context.h>

// file format specific
#include <seqan/rna_io/connect_read_write.h>
#include <seqan/rna_io/dot_bracket_read_write.h>
#include <seqan/rna_io/bpseq_read_write.h>
#include <seqan/rna_io/stockholm_read_write.h>
#include <seqan/rna_io/ebpseq_read_write.h>
#include <seqan/rna_io/vienna_read_write.h>

// general file I/O
#include <seqan/rna_io/rna_struct_file.h>

#endif  // INCLUDE_SEQAN_RNA_IO_H_
