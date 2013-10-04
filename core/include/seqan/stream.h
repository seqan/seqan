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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================
// Facade header for the stream module.
// ==========================================================================

#ifndef SEQAN_STREAM_H_
#define SEQAN_STREAM_H_

/*!
 * @macro SEQAN_HAS_ZLIB
 * @headerfile <seqan/stream.h>
 * @brief Defined as 0 or 1, depending on zlib being available.
 *
 * @signature #define SEQAN_HAS_ZLIB 0  // or 1
 */

/*!
 * @macro SEQAN_HAS_BZIP2
 * @headerfile <seqan/stream.h>
 * @brief Defined as 0 or 1, depending on bzlib being available.
 *
 * @signature #define SEQAN_HAS_BZIP 0  // or 1
 */

/**
.Macro.SEQAN_HAS_ZLIB
..cat:Input/Output
..cat:From Outside
..signature:SEQAN_HAS_ZLIB
..summary:If set to 1 then zlib is available, i.e. including $<zlib.h>$ and linking against libz works.
..remarks:This flag is normally set from the outside by your build system using compiler flags.

.Macro.SEQAN_HAS_BZIP2
..cat:Input/Output
..cat:From Outside
..signature:SEQAN_HAS_BZLIB
..summary:If set to 1 then bzlib2 is available, i.e. including $<bzlib.h>$ and linking against libbzip2 works.
..remarks:This flag is normally set from the outside by your build system using compiler flags.
 */

// ===========================================================================
// Prerequisites.
// ===========================================================================

#include <iostream>
#include <fstream>
#include <sstream>

#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/sequence.h>

// ===========================================================================
// Stream Concept, Adaptions.
// ===========================================================================

#include <seqan/stream/concept_stream.h>
#include <seqan/stream/stream_base.h>

#include <seqan/stream/adapt_ios.h>

// TODO(somebody): refactor cstdio adaption.

//#include <seqan/stream/adapt_cstdio.h>
//
//#include <seqan/stream/file_stream.h>
//
//#if SEQAN_HAS_ZLIB
//// Enable Stream<GZFile> and Stream<Bgzf> if available.
//#include <seqan/stream/stream_gz_file.h>
//#include <seqan/stream/stream_bgzf.h>
//#endif  // #if SEQAN_HAS_ZLIB
//
//#if SEQAN_HAS_BZIP2  // Enable Stream<BZ2File> if available.
//#include <seqan/stream/stream_bz2_file.h>
//#endif  // #if SEQAN_HAS_BZIP2

// ===========================================================================
// Stream Iterators.
// ===========================================================================

#include <seqan/stream/iter_stream.h>

// ===========================================================================
// Tokenizing and *is
// ===========================================================================

#include <seqan/stream/is.h>
//#include <seqan/stream/tokenize.h>
#include <seqan/stream/tokenization.h>
#include <seqan/stream/lexical_cast.h>

#endif  // SEQAN_STREAM_H_
