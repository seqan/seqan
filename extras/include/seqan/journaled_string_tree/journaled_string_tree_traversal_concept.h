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
// Author: Rene Rahn <rene.rahn@fu-berlin.de>
// ==========================================================================
// Concept to be implemented by external algorithms that interact with the
// traversal of a Journaled String Tree.
// ==========================================================================

#ifndef EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_JOURNALED_STRING_TREE_TRAVERSAL_CONCEPT_H_
#define EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_JOURNALED_STRING_TREE_TRAVERSAL_CONCEPT_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

struct TraverseStateMaster_;
typedef Tag<TraverseStateMaster_> StateTraverseMaster;

struct TraverseStateBranch_;
typedef Tag<TraverseStateBranch_> StateTraverseBranch;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/*!
 * @concept JstTraversalConcept
 *
 * @headerfile <seqan/journaled_string_tree.h>
 *
 * @brief Concept for a type that triggers the traversal over a @link JournaledStringTree @endlink.
 *
 * @signature concept JstTraversalConcept;
 *
 * The @link JournaledStringTree @endlink offers no iterator but a global traversal function to iterate over
 * the data. To facilitate iterator-like behavior, the JstTraversalConcept offers a common interface to interrupt the
 * traversal at any position and allows communications between the caller and the callee. Using this interface algorithms
 * can generically implement the traversal over a @link JournaledStringTree @endlink.
 */

/*!
 * @mfn JstTraversalConcept#GetState
 * @brief Returns type of the state of the external caller.
 * @signature GetState<TCaller>::Type;
 * @tparam TCaller The type of the caller that triggers the traversal.
 * @return  TState The type of state objects. Default type is @link Nothing @endlink.
 */

/*!
 * @mfn JstTraversalConcept#GetJstTraverser
 * @brief Returns type of the traverser for the caller.
 * @signature GetJstTraverser<TCaller>::Type;
 * @tparam  TCaller The type of the caller that triggers the traversal.
 * @return  TTraverser The type of the traverser (@link JstTraverser @endlink).
 */

/*!
 * @fn JstTraversalConcept#deliverContext
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Interrupts the traversal and delivers the current sequence context to the caller.
 *
 * @signature TSize deliverContext(caller, delegate, callee, tag);
 *
 * @param[in,out] caller   The caller that triggered the traversal.
 * @param[in,out] delegate An additional functor that can be called.
 * @param[in]      callee   The traverser managing the traversal. Must be of type @link JstTraverser @endlink.
 * @param[in]      tag      The tag indicating the current state. One of @link JstTraversalStates @endlink.
 *
 * @return TSize The step size the window is moved by the traverser.
 */

/*!
 * @fn JstTraversalConcept#getState
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Returns the current state of the caller.
 *
 * @signature TState getState(caller);
 * @param[in] caller The caller that triggered the traversal.
 * @return TState The current state of type @link JstTraversalConcept#GetState @endlink. Returns per default
 * @link Nothing @endlink.
 */

/*!
 * @fn JstTraversalConcept#setState
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Sets the current state to the caller.
 *
 * @signature setState(caller, state);
 * @param[in,out] caller   The caller to set the state to.
 * @param[in]      state    The state to be set. Of type @link JstTraversalConcept#GetState @endlink.
 *
 * The default implementation is a <i>no-op</i> function.
 */

SEQAN_CONCEPT(JstTraversalConcept,(TCaller))
{
    typedef typename GetState<TCaller>::Type                  TState;
    typedef typename GetJstTraverser<TCaller>::Type           TTraverser;
    typedef typename Size<TTraverser>::Type                   TSize;
    typedef Nothing                                           TDelegate;

    TCaller    c;
    TState     state;
    TTraverser traverser;
    TSize      size;

    SEQAN_CONCEPT_USAGE(JstTraversalConcept)
    {
        // Concept for setState and getState.
        sameType(getState(c), state);
        state = getState(c);
        setState(c, state);

        // Concept for deliverContext.
        sameType(deliverContext(c, TDelegate(), traverser, StateTraverseMaster()), size);
        sameType(deliverContext(c, TDelegate(), traverser, StateTraverseBranch()), size);
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

}

#endif // EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_JOURNALED_STRING_TREE_TRAVERSAL_CONCEPT_H_
