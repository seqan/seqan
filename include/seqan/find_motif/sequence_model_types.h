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

#ifndef SEQAN_HEADER_SEQUENCE_MODEL_TYPES_H
#define SEQAN_HEADER_SEQUENCE_MODEL_TYPES_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Tags
//////////////////////////////////////////////////////////////////////////////


/**
.Tag.Oops:
..summary:Represents the One Occurrence Per Sequence model.
..cat:Motif Search
..remarks:The @Tag.Oops@ model, which was introduced by Lawrence and Reilly permits 
          exactly one motif occurrence in each sequence.
..include:seqan/find_motif.h
*/

struct Oops
{
	enum{VALUE=0};
};

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Omops:
..summary:Represents the One or More Occurences Per Sequence model.
..cat:Motif Search
..remarks:The @Tag.Omops@ model is comparable with the @Tag.Tcm@ model with the one difference
          that zero occurrence in a sample sequence is not permitted.
..include:seqan/find_motif.h
*/

struct Omops
{
	enum{VALUE=1};
};

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Zoops:
..summary:Represents the Zero or One Occurence Per Sequence model.
..cat:Motif Search
..remarks:The @Tag.Zoops@ model formulated by Bailey and Elkan permits at most one
          motif occurrence in each sequence.
..include:seqan/find_motif.h
*/

struct Zoops
{
	enum{VALUE=2};
	double threshold;

	Zoops():
		threshold((double)0.5)
	{
	}
	Zoops(double val):
		threshold(val)
	{
	}
	~Zoops()
	{
	}
};

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Tcm:
..summary:Represents the Two-Component-Mixture Sequence model.
..cat:Motif Search
..remarks:The @Tag.Tcm@ model formulated by Bailey and Elkan permits any number pf
          non-overlapping motif occurrences per sequence.
..include:seqan/find_motif.h
*/

struct Tcm
{
	enum{VALUE=3};
	double threshold;

	Tcm():
		threshold((double)0.5)
	{
	}
	Tcm(double val):
		threshold(val)
	{
	}
	~Tcm()
	{
	}
};

//////////////////////////////////////////////////////////////////////////////

} // namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
