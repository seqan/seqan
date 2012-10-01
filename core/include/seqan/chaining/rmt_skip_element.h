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

#ifndef SEQAN_RMT_SKIP_ELEMENT_H
#define SEQAN_RMT_SKIP_ELEMENT_H

namespace seqan
{

// Modifications of the struct SkipElement for use in a SkipListStatic< True, RT< MaxTag > >

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
	struct RangeMaxCargo_
	{

		SkipList< TObject, TModus, TSpec, TStructuring > * _assocStruct;
		TObject * _maxObj;

		RangeMaxCargo_()
			: _assocStruct( NULL )
			, _maxObj( NULL )
		{}

	};

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
	SkipList< TObject, TModus, TSpec, TStructuring > *
	_getAssoc( RangeMaxCargo_< TObject, TModus, TSpec, TStructuring > & me )
	{
		return me._assocStruct;
	}

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
	void
	_setAssoc( RangeMaxCargo_< TObject, TModus, TSpec, TStructuring > & me,
				SkipList< TObject, TModus, TSpec, TStructuring > * list )
	{
		me._assocStruct = list;
	}

	
		// handling for max objects

	template< typename TObject, typename TSpec, typename TStructuring > inline
	TObject *
	_getMaxObject( RangeMaxCargo_< TObject, SkipListStatic, TSpec, TStructuring > & me )
	{
		return me._maxObj;
	}

	template< typename TObject, typename TSpec, typename TStructuring > inline
	void
	_setMaxObject( RangeMaxCargo_< TObject, SkipListStatic, TSpec, TStructuring > & me,
				TObject * obj )
	{
		me._maxObj = obj;
	}

	template< typename TObject,typename TSpec, typename TStructuring > inline
	TObject *
	_getMaxObject( SkipElement< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > * me )
	{
		return me->_cargo._maxObj;
	}

	template< typename TObject, typename TSpec, typename TStructuring > inline
	void
	_setMaxObject( SkipElement< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > * me,
					TObject * maxObj )
	{
		me->_cargo._maxObj = maxObj;
	}
	
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
	struct Cargo< SkipElement< TObject, TModus, RT< MaxTree< TSpec > >, TStructuring > >
	{
		typedef RangeMaxCargo_< TObject, TModus, RT< MaxTree< TSpec > >, TStructuring > Type;
	};

		// get the score value of the related object
	template< typename TObject, typename TSpec, typename TStructuring > inline
	typename Weight< TObject >::Type
	weight( SkipElement< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > * me )
	{
		return weight( *me->_cargo._maxObj );
	}

		// get the chain score value of the related object
	template< typename TObject, typename TSpec, typename TStructuring > inline
	typename Weight< TObject >::Type
	priority( SkipElement< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > * me )
	{
		return priority( *me->_cargo._maxObj );
	}

} // namespace SEQAN_NAMESPACE_SKIPLIST
#endif // SEQAN_RMT_SKIP_ELEMENT_H
