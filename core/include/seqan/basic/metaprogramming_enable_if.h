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
// "enable if" functionality.
// ==========================================================================

// SEQAN_NO_GENERATED_FORWARDS

// TODO(holtgrew): Rename *2 Metafunctions to *C metafunctions.
// TODO(holtgrew): Document properly.

#ifndef SEQAN_CORE_INCLUDE_SEQAN_BASIC_METAPROGRAMMING_ENABLE_IF_H_
#define SEQAN_CORE_INCLUDE_SEQAN_BASIC_METAPROGRAMMING_ENABLE_IF_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction EnableIf
// ----------------------------------------------------------------------------


template <typename TBool, typename T = void>
struct EnableIf : EnableIf<typename TBool::Type, T> {};

template <typename T>
struct EnableIf<False, T> {};

template <typename T>
struct EnableIf<True, T>
{
    typedef T Type;
};

//template <typename T>
//struct EnableIf<False, T> {};

// ----------------------------------------------------------------------------
// Metafunction DisableIf
// ----------------------------------------------------------------------------

template <typename TBool, typename T = void>
struct DisableIf : DisableIf<typename TBool::Type, T> {};

template <typename T>
struct DisableIf<True, T> {};

template <typename T>
struct DisableIf<False, T>
{
    typedef T Type;
};

// ----------------------------------------------------------------------------
// Metafunction EnableIf2
// ----------------------------------------------------------------------------

//
// Example for enable-if function: 
//
//    template <typename TContainer>
//    typename EnableIf<
//        IsContainer<TContainer>,                  // 1st arg: enable-if condition
//        typename Size<TContainer>::Type >::Type   // 2nd arg: return type
//    length(TContainer & cont) 
//    {
//        SEQAN_CONCEPT_ASSERT((ContainerConcept<TContainer>));
//        return end(cont) - begin(cont);
//    }
    
template <bool b, typename T>
struct EnableIf2;

template <bool b, typename T = void>
struct EnableIf2
{
    typedef T Type;
};

template <typename T>
struct EnableIf2<false, T> {};

// ----------------------------------------------------------------------------
// Metafunction DisableIf2
// ----------------------------------------------------------------------------

template <bool b, typename T>
struct DisableIf2;

template <bool b, typename T = void>
struct DisableIf2
{
    typedef T Type;
};

template <typename T>
struct DisableIf2<true, T> {};

}  // namespace seqan

// ============================================================================
// Macros
// ============================================================================

/**
.Macro.SEQAN_CTOR_ENABLE_IF
..cat:Concepts
..summary:Bind the visibility of a constructor to an expression.
..signature:SEQAN_CTOR_ENABLE_IF(cond)
..param.cond:Boolean type. If @Tag.Logical Values.tag.True@ or a metafunction that returns @Tag.Logical Values.tag.True@, the following function is visible, otherwise not.
...remarks:The boolean value must be available at compile-time, e.g. $sizeof(T)>4$.
..remarks:This macro allows to bind the visibility of a constructor to a boolean expression
by using the @http://en.wikipedia.org/wiki/Substitution_failure_is_not_an_error|SFINAE@ principle for an optional argument with default value.
It can be used as the last dummy-argument of a constructor.
To avoid an unused argument warning, call $ignoreUnusedVariableWarning(dummy)$ in the constructor's body.
..example.code:
Rational(T const & n, SEQAN_CTOR_ENABLE_IF(Is<IntegerConcept<T> >)) :  // macro must be extra c'tor argument
    num(n), den(1)
{ 
    ignoreUnusedVariableWarning(dummy);		// necessary to avoid unused warning
}
..include:seqan/basic.h
 */

#define SEQAN_CTOR_ENABLE_IF(...) typename EnableIf<__VA_ARGS__>::Type * dummy = 0

/**
.Macro.SEQAN_CTOR_DISABLE_IF
..cat:Concepts
..summary:Bind the visibility of a constructor to an expression.
..signature:SEQAN_CTOR_DISABLE_IF(cond)
..param.cond:Boolean type. If @Tag.Logical Values.tag.False@ or a metafunction that returns @Tag.Logical Values.tag.False@, the following function is visible, otherwise not.
...remarks:The boolean value must be available at compile-time, e.g. $sizeof(T)>4$.
..remarks:This macro allows to bind the visibility of a constructor to a boolean expression
by using the @http://en.wikipedia.org/wiki/Substitution_failure_is_not_an_error|SFINAE@ principle for an optional argument with default value.
It can be used as the last dummy-argument of a constructor.
To avoid an unused argument warning, call $ignoreUnusedVariableWarning(dummy)$ in the constructor's body.
..example.code:
String(T const & s, SEQAN_CTOR_DISABLE_IF(IsProxy<T> >)) :  // macro must be extra c'tor argument
    str(s)
{
    ignoreUnusedVariableWarning(dummy);     // necessary to avoid unused warning
}
..include:seqan/basic.h
 */

#define SEQAN_CTOR_DISABLE_IF(...) typename DisableIf<__VA_ARGS__>::Type * dummy = 0

/**
.Macro.SEQAN_FUNC_ENABLE_IF
..cat:Concepts
..summary:Bind the visibility of a function to an expression.
..signature:SEQAN_FUNC_ENABLE_IF(cond, retType)
..param.cond:Boolean type. If @Tag.Logical Values.tag.True@ or a metafunction that returns @Tag.Logical Values.tag.True@, the following function is visible, otherwise not.
...remarks:The boolean value must be available at compile-time, e.g. $sizeof(T)>4$.
..param.retType:Function return type, e.g. $typename Size<T>::Type$ or $void$.
..remarks:This macro allows to bind the visibility of a function to a boolean expression
by using the @http://en.wikipedia.org/wiki/Substitution_failure_is_not_an_error|SFINAE@ principle.
It can be used in a function declaration/definition instead of its return type.
As constructors have no return types, they must be disabled with @Macro.SEQAN_CTOR_ENABLE_IF@.
..example.code:
template <typename TContainer>
SEQAN_FUNC_ENABLE_IF(
    IsContainer<TContainer>, 
    typename Size<TContainer>::Type)
length(TContainer & cont) 
{
    SEQAN_CONCEPT_ASSERT((ContainerConcept<TContainer>));
    return end(cont) - begin(cont);
}
..include:seqan/basic.h
 */

#define SEQAN_FUNC_ENABLE_IF(...) typename EnableIf<__VA_ARGS__>::Type

/**
.Macro.SEQAN_FUNC_DISABLE_IF
..cat:Concepts
..summary:Bind the visibility of a function to an expression.
..signature:SEQAN_FUNC_DISABLE_IF(cond, retType)
..param.cond:Boolean type. If @Tag.Logical Values.tag.False@ or a metafunction that returns @Tag.Logical Values.tag.False@, the following function is visible, otherwise not.
...remarks:The boolean value must be available at compile-time, e.g. $sizeof(T)>4$.
..param.retType:Function return type, e.g. $typename Size<T>::Type$ or $void$.
..remarks:This macro allows to bind the visibility of a function to a boolean expression
by using the @http://en.wikipedia.org/wiki/Substitution_failure_is_not_an_error|SFINAE@ principle.
It can be used in a function declaration/definition instead of its return type.
As constructors have no return types, they must be disabled with @Macro.SEQAN_CTOR_DISABLE_IF@.
..example.code:
template <typename TContainer>
SEQAN_FUNC_DISABLE_IF(
    IsIntergral<TContainer>,
    typename Size<TContainer>::Type)
length(TContainer & cont)
{
    SEQAN_CONCEPT_ASSERT((ContainerConcept<TContainer>));
    return end(cont) - begin(cont);
}
..include:seqan/basic.h
 */

#define SEQAN_FUNC_DISABLE_IF(...) typename DisableIf<__VA_ARGS__>::Type

// ============================================================================
// Functions
// ============================================================================

#endif  // #ifndef SEQAN_CORE_INCLUDE_SEQAN_BASIC_METAPROGRAMMING_ENABLE_IF_H_
