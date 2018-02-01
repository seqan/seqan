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
// Author: Rene Rahn <rene.rahn@fu-berlin.de>
// ==========================================================================

#ifndef INCLUDE_SEQAN_ALIGN_ALIGN_INTERFACE_WRAPPER_H_
#define INCLUDE_SEQAN_ALIGN_ALIGN_INTERFACE_WRAPPER_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _alignWrapperSequential(); Score; StringSet vs. StringSet
// ----------------------------------------------------------------------------

template <typename TSetH,
          typename TSetV,
          typename TScoreValue, typename TScoreSpec,
          typename TAlignConfig,
          typename TGapModel,
          std::enable_if_t<And<And<Is<ContainerConcept<TSetH>>,
                                   Is<ContainerConcept<typename Value<TSetH>::Type>>>,
                               And<Is<ContainerConcept<TSetV>>,
                                   Is<ContainerConcept<typename Value<TSetV>::Type>>>
                               >::VALUE,
                           int> = 0>
inline auto
_alignWrapperSequential(TSetH const & stringsH,
                        TSetV const & stringsV,
                        Score<TScoreValue, TScoreSpec> const & scoringScheme,
                        TAlignConfig const & config,
                        TGapModel const & /*gaps*/)

{
    String<TScoreValue> results;
    resize(results, length(stringsV));

    auto zipCont = makeZipView(results, stringsH, stringsV);
    forEach(zipCont,
            [&] (auto tuple)
            {
                using namespace seqan;
                DPScoutState_<Default> dpScoutState;
                String<TraceSegment_<unsigned, unsigned> > traceSegments;  // Dummy segments.
                std::get<0>(tuple) = _setUpAndRunAlignment(traceSegments, dpScoutState,
                                                           std::get<1>(tuple), std::get<2>(tuple), scoringScheme,
                                                           config, TGapModel());
            });
    return results;
}

// ----------------------------------------------------------------------------
// Function _alignWrapperSequential(); Score; String vs. StringSet
// ----------------------------------------------------------------------------

template <typename TSeqH,
          typename TSetV,
          typename TScoreValue, typename TScoreSpec,
          typename TAlignConfig,
          typename TGapModel,
          std::enable_if_t<And<And<Is<ContainerConcept<TSeqH>>,
                                   Not<Is<ContainerConcept<typename Value<TSeqH>::Type>>>>,
                               And<Is<ContainerConcept<TSetV>>,
                                   Is<ContainerConcept<typename Value<TSetV>::Type>>>
                               >::VALUE,
                           int> = 0>
inline auto
_alignWrapperSequential(TSeqH const & stringH,
                        TSetV const & stringsV,
                        Score<TScoreValue, TScoreSpec> const & scoringScheme,
                        TAlignConfig const & config,
                        TGapModel const & /*gaps*/)

{
    String<TScoreValue> results;
    resize(results, length(stringsV));

    auto zipCont = makeZipView(results, stringsV);
    forEach(zipCont,
            [&] (auto tuple)
            {
                using namespace seqan;
                DPScoutState_<Default> dpScoutState;
                String<TraceSegment_<unsigned, unsigned> > traceSegments;  // Dummy segments.
                std::get<0>(tuple) = _setUpAndRunAlignment(traceSegments, dpScoutState, stringH, std::get<1>(tuple),
                                                           scoringScheme, config, TGapModel());
            });
    return results;
}

// ----------------------------------------------------------------------------
// Function _alignWrapperSequential(); Gaps
// ----------------------------------------------------------------------------

template <typename TSetH,
          typename TSetV,
          typename TScoreValue, typename TScoreSpec,
          typename TAlignConfig,
          typename TGapModel,
          std::enable_if_t<And<And<Is<ContainerConcept<TSetH>>,
                                   Is<AlignedSequenceConcept<typename Value<TSetH>::Type>>>,
                               And<Is<ContainerConcept<TSetV>>,
                                   Is<AlignedSequenceConcept<typename Value<TSetV>::Type>>>
                               >::VALUE,
                           int> = 0>
inline auto
_alignWrapperSequential(TSetH & gapSeqSetH,
                        TSetV & gapSeqSetV,
                        Score<TScoreValue, TScoreSpec> const & scoringScheme,
                        TAlignConfig const & config,
                        TGapModel const & /*gaps*/)

{
    typedef typename Size<TSetH>::Type TSize;
    typedef typename Position<TSetH>::Type TPosition;
    typedef TraceSegment_<TPosition, TSize> TTraceSegment;

    String<TScoreValue> results;
    resize(results, length(gapSeqSetH));
    
    auto zipCont = makeZipView(results, gapSeqSetH, gapSeqSetV);
    forEach(zipCont,
            [&] (auto tuple)
            {
                using namespace seqan;
                String<TTraceSegment> trace;
                DPScoutState_<Default> dpScoutState;
                std::get<0>(tuple) = _setUpAndRunAlignment(trace, dpScoutState, source(std::get<1>(tuple)),
                                                           source(std::get<2>(tuple)), scoringScheme, config,
                                                           TGapModel());
                _adaptTraceSegmentsTo(std::get<1>(tuple), std::get<2>(tuple), trace);
            });
    return results;
}

// ----------------------------------------------------------------------------
// Function _alignWrapper()
// ----------------------------------------------------------------------------

template <typename... TArgs>
inline auto _alignWrapper(TArgs && ...args)
{
// NOTE(marehr): ume_simd is currently not working, thus falling back to sequential case
#if defined(SEQAN_SIMD_ENABLED) && !defined(SEQAN_UMESIMD_ENABLED)
    return _alignWrapperSimd(std::forward<TArgs>(args)...);
#else
    return _alignWrapperSequential(std::forward<TArgs>(args)...);
#endif
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_ALIGN_ALIGN_INTERFACE_WRAPPER_H_
