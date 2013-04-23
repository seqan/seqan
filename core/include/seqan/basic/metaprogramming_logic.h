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
// ==========================================================================
// Logic types and operations for metaprogramming.
//
// We can perform logical computations both type and value based.  The type
// based metafunctions are the "main" ones, the value based meta functions
// have the postfix "C" (e.g. "Or" vs "OrC"), similar to how Boost's MPL does
// this.
// ==========================================================================

#ifndef SEQAN_CORE_INCLUDE_SEQAN_BASIC_METAPROGRAMMING_LOGIC_H_
#define SEQAN_CORE_INCLUDE_SEQAN_BASIC_METAPROGRAMMING_LOGIC_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Tags True and False
// ----------------------------------------------------------------------------

/**
.Tag.Logical Values:
..cat:Metaprogramming
..summary:Tag that represents true and false.
..tag.True:The logical value "true".
..tag.False:The logical value "false".
..remarks:These tags also function as Metafunctions that return a boolean value $VALUE$ and themselves ($True$/$False$) as $Type$.
..example.text:Print the values of these tags/metafunctions.
..example.code:
std::cout << False::VALUE << std::endl;                         // 0
std::cout << True::VALUE << std::endl;                          // 1
std::cout << IsSameType<False,False::Type>::VALUE << std::endl; // 1
..include:seqan/basic.h
..see:Metafunction.Eval
..see:Metafunction.IsSameType
*/

struct True
{
    typedef True Type;
	static const bool VALUE = true;
};

struct False
{
    typedef False Type;
	static const bool VALUE = false;
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Eval
// ----------------------------------------------------------------------------

/**
.Metafunction.Eval
..cat:Metaprogramming
..summary:Convert from $bool$ values to types ($True$ and $False$).
..signature:Eval<b>::Type
..param.b:The boolean to evaluate.
...type:nolink:$bool$
..returns:Either @Tag.Logical Values.tag.True@ or @Tag.Logical Values.tag.False@, depending on $b$.
..include:seqan/basic.h
..example.text:Demonstrate the usage of $Eval$ to bridge between $bool$ values and the logical tags.
..example.code:
void printBoolType(True const &)
{
    std::cout << "true" << std::endl;
}

void printBoolType(False const &)
{
    std::cout << "false" << std::endl;
}

int main(int argc, char ** argv)
{
    using namespace seqan;

    printBoolType(Eval<true>::Type());   // => "true\n"
    printBoolType(Eval<false>::Type());  // => "false\n"
    return 0;
}
*/

template <bool B>
struct Eval : False {};

template <>
struct Eval<true> : True {};

// ----------------------------------------------------------------------------
// Metafunction Not
// ----------------------------------------------------------------------------

/**
.Metafunction.Not
..cat:Metaprogramming
..summary:Metaprogramming boolean "not" operator.
..signature:Not<B>::Type
..param.B:Argument to invert.
...type:Tag.Logical Values.tag.True
...type:Tag.Logical Values.tag.False
..returns:One of @Tag.Logical Values.tag.True@ and @Tag.Logical Values.tag.False@, the result of logical not.
The arguments $B$ can either be @Tag.Logical Values.tag.True@/@Tag.Logical Values.tag.False@
or boolean metafunctions that return @Tag.Logical Values.tag.True@/@Tag.Logical Values.tag.False@.
..example.code:
Not<False>::Type
Not<True>::Type
..include:seqan/basic.h
*/

template <typename TBool1>
struct Not : Not<typename TBool1::Type> {}; 

template <>
struct Not<False> : True {};
template <>
struct Not<True> : False {};

// ----------------------------------------------------------------------------
// Metafunction NotC
// ----------------------------------------------------------------------------

/**
.Metafunction.NotC
..cat:Metaprogramming
..summary:Metaprogramming boolean "not" operator, value variant.
..signature:NotC<B>::Type
..param.B:Argument
...type:nolink:$bool$
..returns:$!B1$.
..example.code:
NotC<false, true>::Type
NotC<false, false>::VALUE
NotC<AndC<false, true>::VALUE, false>::Type
NotC<AndC<false, true>, true>::Type
..see:Metafunction.Or
..see:Metafunction.And
..see:Metafunction.AndC
..include:seqan/basic.h
*/

template <bool B>
struct NotC;

template <>
struct NotC<false> : True {};
template <>
struct NotC<true> : False {};

// ----------------------------------------------------------------------------
// Metafunction Or
// ----------------------------------------------------------------------------

/**
.Metafunction.Or
..cat:Metaprogramming
..summary:Metaprogramming boolean "or" operator.
..signature:Or<B1, B2>::Type
..signature:Or<B1, B2>::VALUE
..param.B1:Left-hand argument.
...type:Tag.Logical Values.tag.True
...type:Tag.Logical Values.tag.False
..param.B2:Right-hand argument.
...type:Tag.Logical Values.tag.True
...type:Tag.Logical Values.tag.False
..returns:One of @Tag.Logical Values.tag.True@ and @Tag.Logical Values.tag.False@, the result of logical or.
The arguments $B1$/$B2$ can either be @Tag.Logical Values.tag.True@/@Tag.Logical Values.tag.False@
or boolean metafunctions that return @Tag.Logical Values.tag.True@/@Tag.Logical Values.tag.False@.
..example.code:
Or<False,False>::Type
Or<False,True>::Type
Or<typename And<T1,T2>::Type,T3>::Type
Or<And<T1,T2>,T3>::Type
..see:Metafunction.And
..see:Metafunction.AndC
..see:Metafunction.OrC
..include:seqan/basic.h
*/

template <typename TBool1, typename TBool2>
struct Or : Or<typename TBool1::Type, typename TBool2::Type> {}; 

template <>
struct Or<False, False> : False {};
template <>
struct Or<False, True> : True {};
template <>
struct Or<True, False> : True {};
template <>
struct Or<True, True> : True {};

// ----------------------------------------------------------------------------
// Metafunction OrC
// ----------------------------------------------------------------------------

/**
.Metafunction.OrC
..cat:Metaprogramming
..summary:Metaprogramming boolean "or" operator, value variant.
..signature:OrC<B1, B2>::Type
..signature:OrC<B1, B2>::VALUE
..param.B1:Left-hand argument.
...type:nolink:$bool$
..param.B2:Right-hand argument.
...type:nolink:$bool$
..returns:$B1 || B2$.
..example.code:
OrC<false, true>::Type
OrC<false, false>::VALUE
OrC<AndC<false, true>::VALUE, false>::Type
OrC<AndC<false, true>, true>::Type
..see:Metafunction.Or
..see:Metafunction.And
..see:Metafunction.AndC
..include:seqan/basic.h
*/

template <bool B1, bool B2>
struct OrC;

template <>
struct OrC<false, false> : False {};
template <>
struct OrC<false, true> : True {};
template <>
struct OrC<true, false> : True {};
template <>
struct OrC<true, true> : True {};

// ----------------------------------------------------------------------------
// Metafunction And
// ----------------------------------------------------------------------------

/**
.Metafunction.And
..cat:Metaprogramming
..summary:Metaprogramming boolean "And" operator.
..signature:Or<B1, B2>::Type
..param.B1:Left-hand argument.
...type:Tag.Logical Values.tag.True
...type:Tag.Logical Values.tag.False
..param.B2:Right-hand argument.
...type:Tag.Logical Values.tag.True
...type:Tag.Logical Values.tag.False
..returns:One of @Tag.Logical Values.tag.True@ and @Tag.Logical Values.tag.False@, the result of logical and.
..include:seqan/basic.h
*/

template <typename TBool1, typename TBool2>
struct And : And<typename TBool1::Type, typename TBool2::Type> {}; 

template <>
struct And<False, False> : False {};
template <>
struct And<False, True> : False {};
template <>
struct And<True, False> : False {};
template <>
struct And<True, True> : True {};

// ----------------------------------------------------------------------------
// Metafunction OrC
// ----------------------------------------------------------------------------

/**
.Metafunction.AndC
..cat:Metaprogramming
..summary:Metaprogramming boolean "and" operatand, value variant.
..signature:AndC<B1, B2>::Type
..signature:AndC<B1, B2>::VALUE
..param.B1:Left-hand argument.
...type:nolink:$bool$
..param.B2:Right-hand argument.
...type:nolink:$bool$
..returns:$B1 && B2$.
..example.code:
AndC<false, true>::Type
AndC<false, false>::VALUE
AndC<AndC<false, true>::VALUE, false>::Type
AndC<AndC<false, true>, true>::Type
..see:Metafunction.And
..see:Metafunction.And
..see:Metafunction.AndC
..include:seqan/basic.h
*/

template <bool B1, bool B2>
struct AndC;

template <>
struct AndC<false, false> : False {};
template <>
struct AndC<false, true> : False {};
template <>
struct AndC<true, false> : False {};
template <>
struct AndC<true, true> : True {};

// ----------------------------------------------------------------------------
// Metafunction If
// ----------------------------------------------------------------------------

/**
.Metafunction.If
..cat:Metaprogramming
..summary:Metaprogramming "if".
..signature:If<T, T1, T2>::Type
..param.T:The condition.
..param.T1:Result if $b$.
..param.T2:Result if not $b$.
..see:Metafunction.IfC
..returns:If $T$ is $True$ then $T1$, otherwise $T2$.
..include:seqan/basic.h
*/

// TODO(holtgrew): Inherit from T1/T2 to follow pattern for Or<>, And<>?
template <typename TCondition, typename T1, typename T2>
struct If : If<typename TCondition::Type, T1, T2>{};

template <typename T1, typename T2>
struct If<True, T1, T2>
{
    typedef T1 Type;
};

template <typename T1, typename T2>
struct If<False, T1, T2>
{
    typedef T2 Type;
};

// ----------------------------------------------------------------------------
// Metafunction IfC
// ----------------------------------------------------------------------------

/**
.Metafunction.IfC
..cat:Metaprogramming
..summary:Metaprogramming "if", value based version.
..signature:If<b, T1, T2>::Type
..param.b:The condition to evaluate.
...type:nolink:$bool$
..param.T1:Result if $b$.
..param.T2:Result if not $b$.
..returns:If $b$ is $true$ then $T1$, otherwise $T2$.
..see:Metafunction.If
..include:seqan/basic.h
*/

// TODO(holtgrew): Inherit from T1/T2 to follow pattern for Or<>, And<>?
template <bool FLAG, typename T1, typename T2>
struct IfC
{
    typedef T1 Type;
};

template <typename T1, typename T2>
struct IfC<false, T1, T2>
{
    typedef T2 Type;
};

// ============================================================================
// Functions
// ============================================================================

}  // namespace seqan

#endif  // #ifndef SEQAN_CORE_INCLUDE_SEQAN_BASIC_METAPROGRAMMING_LOGIC_H_
