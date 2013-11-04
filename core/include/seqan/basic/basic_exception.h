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
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================
// Definition of basic exceptions.
// ==========================================================================

#ifndef SEQAN_BASIC_BASIC_EXCEPTION_H_
#define SEQAN_BASIC_BASIC_EXCEPTION_H_

// ============================================================================
// Prerequisites
// ============================================================================

#include <typeinfo>
#include <exception>
#include <stdexcept>

#ifdef PLATFORM_GCC
#include <cxxabi.h>
#endif

// ============================================================================
// Macros
// ============================================================================

/*!
 * @macro SEQAN_EXCEPTIONS
 * @headerfile <seqan/basic.h>
 * @brief Determines whether exceptions are enabled or not.
 * 
 * @signature SEQAN_EXCEPTIONS
 *
 * @see SEQAN_TRY
 * @see SEQAN_CATCH
 * @see SEQAN_THROW
 * @see Exception
 */

#define SEQAN_EXCEPTIONS    __EXCEPTIONS

/*!
 * @macro SEQAN_TRY
 * @headerfile <seqan/basic.h>
 * @brief Replaces the C++ try keyword.
 * 
 * @signature SEQAN_TRY {} SEQAN_CATCH() {}
 *
 * @section Remarks
 * 
 * When exceptions are disabled, i.e. SEQAN_EXCEPTIONS is set to false, the code inside the try block is always executed".
 * 
 * @see SEQAN_CATCH
 * @see SEQAN_THROW
 * @see Exception
 *
 * @section Examples
 *
 * @code{.cpp}
 *
 * SEQAN_TRY
 * {
 *     SEQAN_THROW(Exception)
 * }
 * SEQAN_CATCH(Exception const & e)
 * {
 *     std::cerr << e.what() << std::endl;
 * }
 *
 * @endcode
 */

/*!
 * @macro SEQAN_CATCH
 * @headerfile <seqan/basic.h>
 * @brief Replaces the C++ catch keyword.
 * 
 * @signature SEQAN_TRY {} SEQAN_CATCH() {}
 *
 * @section Remarks
 *
 * When exceptions are disabled, i.e. SEQAN_EXCEPTIONS is set to false, the code inside the catch block is never executed".
 * 
 * @see SEQAN_TRY
 * @see SEQAN_THROW
 * @see Exception
 *
 * @section Examples
 *
 * See @link SEQAN_TRY @endlink for a full example.
 *
 */

/*!
 * @macro SEQAN_THROW
 * @headerfile <seqan/basic.h>
 * @brief Replaces the C++ throw keyword.
 * 
 * @signature SEQAN_THROW(Exception);
 *
 * @section Remarks
 *
 * When exceptions are disabled, i.e. SEQAN_EXCEPTIONS is set to false, the macro turns into SEQAN_FAIL".
 * 
 * @see SEQAN_TRY
 * @see SEQAN_CATCH
 * @see SEQAN_FAIL
 * @see Exception
 *
 * @section Examples
 *
 * See @link SEQAN_TRY @endlink for a full example.
 *
 */
 
#ifdef SEQAN_EXCEPTIONS

#define SEQAN_TRY           try
#define SEQAN_CATCH(E)      catch(E)
#define SEQAN_THROW(E)      throw E
#define SEQAN_RETHROW       throw

#else

#define SEQAN_TRY           if (true)
#define SEQAN_CATCH(E)      if (false)
//#define SEQAN_CATCH(E)      for (E ; false; )
#define SEQAN_THROW(E)      SEQAN_FAIL(#E)
#define SEQAN_RETHROW

#endif // #ifdef SEQAN_EXCEPTIONS

namespace seqan {

// ============================================================================
// Exceptions
// ============================================================================

// ----------------------------------------------------------------------------
// Basic Exception
// ----------------------------------------------------------------------------

/*!
 * @class Exception
 * @headerfile <seqan/basic.h>
 * @brief Generic SeqAn exception.
 * @signature Exception();
 */

typedef std::exception          Exception;

// ----------------------------------------------------------------------------
// Exception BadAlloc
// ----------------------------------------------------------------------------

/*!
 * @class BadAlloc
 * @headerfile <seqan/basic.h>
 * @brief Bad memory allocation exception.
 * @signature BadAlloc();
 */

typedef std::bad_alloc          BadAlloc;

// ----------------------------------------------------------------------------
// Exception BadCast
// ----------------------------------------------------------------------------

/*!
 * @class BadCast
 * @headerfile <seqan/basic.h>
 * @brief Bad cast exception.
 * @signature BadCast();
 */

typedef std::bad_cast           BadCast;

// ----------------------------------------------------------------------------
// Exceptions Bad*
// ----------------------------------------------------------------------------
// NOTE(esiragusa): These exceptions can be introduced as long as we need them.

//typedef std::bad_exception      BadException;
//typedef std::bad_typeid         BadTypeId;
//typedef std::bad_function_call  BadFunctionCall;
//typedef std::bad_weak_ptr       BadWeakPtr;

// ----------------------------------------------------------------------------
// Exception RuntimeError
// ----------------------------------------------------------------------------

/*!
 * @class RuntimeError
 * @headerfile <seqan/basic.h>
 * @brief Runtime error exception.
 * @signature RuntimeError("Message");
 */

typedef std::runtime_error      RuntimeError;

// ----------------------------------------------------------------------------
// Exception LogicError
// ----------------------------------------------------------------------------
// NOTE(esiragusa): Always prefer SEQAN_ASSERT to logic error exceptions.

//typedef std::logic_error        LogicError;

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction ExceptionMessage
// ----------------------------------------------------------------------------

template <typename T, typename TSpec = void>
struct ExceptionMessage
{
    static const std::string VALUE;
};

template <typename T, typename TSpec>
const std::string ExceptionMessage<T, TSpec>::VALUE;

// ============================================================================
// Functors
// ============================================================================

// ----------------------------------------------------------------------------
// Functor Asserter
// ----------------------------------------------------------------------------

template <typename TFunctor, typename TException, typename TContext = void, bool RETURN_VALUE = false>
struct Asserter
{
    TFunctor func;

    Asserter() {}

    Asserter(TFunctor & func) :
        func(func)
    {}

    template <typename TValue>
    bool operator() (TValue const & val) const
    {
        if (SEQAN_UNLIKELY(!func(val)))
            throw TException(std::string("Value '") + val + "' produced an error. " +
                             ExceptionMessage<TFunctor, TContext>::VALUE);
        return RETURN_VALUE;
    }
};

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class Demangler
// ----------------------------------------------------------------------------
// Holds the name of a given C++ type T.
// NOTE(esiragusa): this class could become a subclass of CStyle String...

template <typename T>
struct Demangler
{
    char *data_begin;

    Demangler()
    {
        T t;
        _demangle(*this, t);
    }

    Demangler(T const & t)
    {
        _demangle(*this, t);
    }

    ~Demangler()
    {
#ifdef PLATFORM_GCC
        free(data_begin);
#endif
    }
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _demangle(Demangler)
// ----------------------------------------------------------------------------

template <typename T>
inline void _demangle(Demangler<T> & me, T const & t)
{
#ifdef PLATFORM_GCC
    int status;
    me.data_begin = abi::__cxa_demangle(typeid(t).name(), NULL, NULL, &status);
#else
    me.data_begin = typeid(t).name();
#endif
}

// ----------------------------------------------------------------------------
// Function toCString(Demangler)
// ----------------------------------------------------------------------------

template <typename T>
inline char * toCString(Demangler<T> const & me)
{
    return me.data_begin;
}

// ----------------------------------------------------------------------------
// Function toCTypeName()
// ----------------------------------------------------------------------------

template <typename T>
inline char * toCTypeName(T const & t)
{
    return toCString(Demangler<T>(t));
}

// ----------------------------------------------------------------------------
// Function globalExceptionHandler()
// ----------------------------------------------------------------------------

static void globalExceptionHandler()
{
    SEQAN_TRY
    {
        SEQAN_RETHROW;
    }
    SEQAN_CATCH(Exception & e)
    {
        SEQAN_FAIL("Uncaught exception of type %s: %s", toCTypeName(e), e.what());
    }
}

// Install global exception handler.
#ifdef SEQAN_EXCEPTIONS
static const std::terminate_handler _globalExceptionHandler = std::set_terminate(globalExceptionHandler);
#endif

}  // namespace seqan

#endif  // #ifndef SEQAN_BASIC_BASIC_EXCEPTION_H_
