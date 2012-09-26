// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
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

/*	
*	Random bit generator
*	
*	from Numerical Recipes in C
*	
*/


#ifndef SEQAN_HEADER_RAND_GEOM
#define SEQAN_HEADER_RAND_GEOM


#define SEQAN_RG_IB1 1
#define SEQAN_RG_IB2 2
#define SEQAN_RG_IB5 16
#define SEQAN_RG_IB18 131072
#define SEQAN_RG_MASK ( SEQAN_RG_IB1 + SEQAN_RG_IB2 + SEQAN_RG_IB5 )

namespace seqan {

	template< typename T > inline
	T
	_geomRand( )
	{
		static unsigned long seed = rand();
		T value = 0;
		while( true )
		{
			if( ( seed & SEQAN_RG_IB18 ) ){
				seed = ( ( seed ^ SEQAN_RG_MASK ) << 1 ) | SEQAN_RG_IB1;
				++value;
			}
			else {
				seed <<= 1;
				break;
			}
		}
		return value;
	}

}

#endif
