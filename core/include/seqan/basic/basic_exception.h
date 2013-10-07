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
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class Exception
// ----------------------------------------------------------------------------

/*!
 * @class Exception
 * @headerfile <seqan/basic.h>
 * @brief Generic SeqAn exception.
 * @signature Exception;
 */

typedef std::exception          Exception;

// ----------------------------------------------------------------------------
// Class BadAlloc
// ----------------------------------------------------------------------------

/*!
 * @class BadAlloc
 * @headerfile <seqan/basic.h>
 * @brief Bad memory allocation exception.
 * @signature BadAlloc;
 */

typedef std::bad_alloc          BadAlloc;

// ----------------------------------------------------------------------------
// Classes Bad*
// ----------------------------------------------------------------------------
// NOTE(esiragusa): These exceptions can be introduced as long as we need them.

//typedef std::bad_exception      BadException;
//typedef std::bad_cast           BadCast;
//typedef std::bad_typeid         BadTypeId;
//typedef std::bad_function_call  BadFunctionCall;
//typedef std::bad_weak_ptr       BadWeakPtr;

// ----------------------------------------------------------------------------
// Class RuntimeError
// ----------------------------------------------------------------------------

/*!
 * @class RuntimeError
 * @headerfile <seqan/basic.h>
 * @brief Runtime error exception.
 * @signature RuntimeError("Message");
 */

typedef std::runtime_error      RuntimeError;

// ----------------------------------------------------------------------------
// Class LogicError
// ----------------------------------------------------------------------------
// NOTE(esiragusa): Always prefer SEQAN_ASSERT to logic error exceptions.

//typedef std::logic_error        LogicError;

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function globalExceptionHandler()
// ----------------------------------------------------------------------------

#ifdef SEQAN_EXCEPTIONS
static void globalExceptionHandler()
{
    SEQAN_TRY
    {
        SEQAN_RETHROW;
    }
    SEQAN_CATCH(Exception & e)
    {
        SEQAN_FAIL("Uncaught exception of type %s: %s", typeid(e).name(), e.what());
    }
}

// Install global exception handler.
static const std::terminate_handler _globalExceptionHandler = std::set_terminate(globalExceptionHandler);
#endif

}  // namespace seqan

#endif  // #ifndef SEQAN_BASIC_BASIC_EXCEPTION_H_
