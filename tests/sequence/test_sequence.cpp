// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
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

#define SEQAN_DEBUG


#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>

using namespace seqan;

#include "test_string.h"
#include "test_string_packed_extension.h"
#include "test_stringset.h"
#include "test_segment.h"
#include "test_sequence_std_adaptions.h"


SEQAN_BEGIN_TESTSUITE(Sequence tests)
{
    // -----------------------------------------------------------------------
    // Tests for STL adaptions.
    // -----------------------------------------------------------------------
    //
    // Test adaptions for std::string.
    SEQAN_CALL_TEST(test_sequence_adaptions_metafunctions_std_string);
    SEQAN_CALL_TEST(test_sequence_adaptions_iterators_std_string);
    SEQAN_CALL_TEST(test_sequence_adaptions_sequence_interface_std_string);
    SEQAN_CALL_TEST(test_sequence_adaptions_sequence_memory_std_string);

	// Test adaptions for std::vector.
    SEQAN_CALL_TEST(test_sequence_adaptions_metafunctions_std_vector);
    SEQAN_CALL_TEST(test_sequence_adaptions_iterators_std_vector);
    SEQAN_CALL_TEST(test_sequence_adaptions_sequence_interface_std_vector);
    SEQAN_CALL_TEST(test_sequence_adaptions_sequence_memory_std_vector);

    // Test adaptions for std::array.
    SEQAN_CALL_TEST(test_sequence_adaptions_metafunctions_std_array);
    SEQAN_CALL_TEST(test_sequence_adaptions_iterators_std_array);
    SEQAN_CALL_TEST(test_sequence_adaptions_sequence_interface_std_array);

    // Test adaptions for std::list.
    SEQAN_CALL_TEST(test_sequence_adaptions_metafunctions_std_list);
    SEQAN_CALL_TEST(test_sequence_adaptions_iterators_std_list);
    SEQAN_CALL_TEST(test_sequence_adaptions_sequence_interface_std_list);

    // Use the constant EMPTY_STRING once to get rid of linker error in Visual Studio.
    (void)seqan::String<char const, seqan::CStyle>::EMPTY_STRING;

	SEQAN_CALL_TEST(Sequence_Interface);
	SEQAN_CALL_TEST(String_Base);
	SEQAN_CALL_TEST(String_Alloc);
	SEQAN_CALL_TEST(String_Array);
	SEQAN_CALL_TEST(String_Stack);
	SEQAN_CALL_TEST(String_Pointer);
	SEQAN_CALL_TEST(String_CStyle);
	SEQAN_CALL_TEST(String_Packed);
	SEQAN_CALL_TEST(Std_String);

	SEQAN_CALL_TEST(Lexical);
	SEQAN_CALL_TEST(Combinatoric);

	SEQAN_CALL_TEST(Segment);

    SEQAN_CALL_TEST(StringSet_Owner_Default);
    SEQAN_CALL_TEST(StringSet_Concat_Owner_Default);
    SEQAN_CALL_TEST(StringSet_Concat_Owner_ConcatDirect);
    SEQAN_CALL_TEST(StringSet_Id_Dependent_Tight);
    SEQAN_CALL_TEST(StringSet_Id_Dependent_Generous);
    SEQAN_CALL_TEST(StringSetIdHolder_Char_Dependent_Tight);
    SEQAN_CALL_TEST(StringSetIdHolder_Char_Dependent_Generous);

    SEQAN_CALL_TEST(Infix);
    SEQAN_CALL_TEST(Suffix);
    SEQAN_CALL_TEST(ticket317);
    SEQAN_CALL_TEST(ticket848);
    SEQAN_CALL_TEST(test_find_motif_memory_leak_ticket_364);
    SEQAN_CALL_TEST(ticket901);
    #ifndef __OpenBSD__
    // TODO(h-2): fix this test on OpenBSD (some problem with mmap)
    SEQAN_CALL_TEST(ticket1108);
    #endif
    SEQAN_CALL_TEST(String_Packed_Extension);
}
SEQAN_END_TESTSUITE
