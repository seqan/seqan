// ==========================================================================
//                                  ANISE
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

#include "file_name_tokens.h"

char const * GLOBAL_STATE_TOKEN = "global_state";
char const * GLOBAL_STATE_EXT = ".txt";

// Token and extension for site state files.
char const * SITE_STATE_TOKEN = "site_state";
char const * SITE_STATE_EXT = ".txt";

// Token and extension for orphans FASTQ file.
char const * ORPHANS_TOKEN = "orphans";
char const * ORPHANS_EXT = ".fq";

// Token and extension for orphans active map.
char const * ORPHANS_ACTIVE_TOKEN = "orphans_active";
char const * ORPHANS_ACTIVE_EXT = ".bin";

// Token and extension for reads SAM files.
char const * READS_TOKEN = "reads";
char const * READS_EXT = ".sam";

// Token and extension for orphans active map.
char const * SCAFFOLD_SEQS_TOKEN = "scaffold_seqs";
char const * SCAFFOLD_SEQS_EXT = ".fa";

// Token and extension for time log.
char const * TIME_LOG_TOKEN = "time_log";
char const * TIME_LOG_EXT = ".tsv";
