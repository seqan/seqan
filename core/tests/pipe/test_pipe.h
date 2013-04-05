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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_HEADER_TEST_PIPE_H

#define SEQAN_HEADER_TEST_PIPE_H



namespace SEQAN_NAMESPACE_MAIN

{



	template < typename TBuffer >

	void permute(TBuffer buf) {

		typename Size<TBuffer>::Type i, j, s = length(buf);

//        srand( (unsigned)time( NULL ) );

		for(i = 0; i < s; i++)

			buf[i] = s-i-1;

		for(i = 0; i < s; i++) {

            if (i > 0) {

			    j = i - (rand() % i) - 1;

                assert(/*0 <= j && */j < s);

            } else

                j = 0;

			unsigned tmp = buf[i];

			buf[i] = buf[j];

			buf[j] = tmp;

		}

	}



	template < typename TBuffer >

	void randomize(TBuffer buf) {

		typename Size<TBuffer>::Type i, s = length(buf);

		for(i = 0; i < s; i++)

            buf[i] = rand() % s;

	}



	template < typename TValue >

	struct IdentityMap : public ::std::unary_function< TValue, TValue > {

		inline TValue operator() (TValue const i) { return i; }

	};



	template < typename TValue >

	struct SimpleCompare : public ::std::binary_function< TValue const, TValue const, int > {

		inline int operator() (TValue const a, TValue const b) const {

            if (a < b) return -1;

            if (a > b) return 1;

            return 0;

        }

	};

	

}



#endif



