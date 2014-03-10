// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2013, Knut Reinert, FU Berlin
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
// Author: Hannes Hauswedell <hauswedell@mi.fu-berlin.de>
// ==========================================================================
// Module with reduced versions of alphabets
// ==========================================================================

#ifndef SEQAN_EXTRAS_REDUCED_AMINOACID_H_
#define SEQAN_EXTRAS_REDUCED_AMINOACID_H_

#ifdef SEQAN_CXX11_STANDARD

#include <seqan/basic.h>
#include <seqan/score.h>

#include <seqan/reduced_aminoacid/base.h>

#include <seqan/reduced_aminoacid/murphy10_base.h>
#include <seqan/reduced_aminoacid/murphy10_tables.h>

#include <seqan/reduced_aminoacid/cluster_red_base.h>
#include <seqan/reduced_aminoacid/cluster_red_tables_20_to_n_b62.h>
#include <seqan/reduced_aminoacid/cluster_red_tables_22_to_n_b62.h>
#include <seqan/reduced_aminoacid/cluster_red_tables_24_to_n_b62.h>

#include <seqan/reduced_aminoacid/base_late.h>

#else
#error Module reduced_aminoacid is only available when SEQAN_CXX11_STANDARD is \
defined and your compiler supports alias templates.
#endif

#endif // def SEQAN_EXTRAS_REDUCED_AMINOACID_H_
