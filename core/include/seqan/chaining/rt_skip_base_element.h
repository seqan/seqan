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

//SEQAN_NO_DDDOC: do not generate documentation for this file

#ifndef RT_MAX_SKIP_BASE_ELEMENT_H_
#define RT_MAX_SKIP_BASE_ELEMENT_H_


namespace seqan
{


//___________________________ struct SkipBaseElement< TObject, TModus, RT< TSpec >, Deferred > > _______________________
// Adaption of the struct SkipBaseElement for use in a RangeTree
// Instead saving it's own theKey and value, it has a pointer to the related 
// RTEntry object, which stores the theKey and value
// To get the correct theKey, SkipBaseElement needs information about the dimension
// of the RangeTree it is part of.


	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
	void 
	dump( SkipBaseElement< TObject, TModus, RT< TSpec >, TStructuring > & me )
	{
		typename Size< TObject >::Type dim = dimension( *getObject( &me ) );
		typename Size< TObject >::Type counter = 0;
		while( counter < dim )
		{
			if( key( &getObject( &me ), counter ) == minValue< typename Key< TObject >::Type >( ) )
					std::cout << std::left << "L";
			else
				std::cout << key( getObject( &me ), counter );
			++counter;
			std::cout << " ";
		}
		std::cout << std::endl;
	}


	
} // namespace SEQAN_NAMESPACE_SKIPLIST
#endif //RT_SKIP_BASE_ELEMENT_H_
