// ==========================================================================
//                      Yara - Yet Another Read Aligner
// ==========================================================================
// Copyright (c) 2011-2018, Enrico Siragusa, FU Berlin
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
//     * Neither the name of Enrico Siragusa or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ENRICO SIRAGUSA OR THE FU BERLIN BE LIABLE
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

#ifndef APP_YARA_BITS_CONTEXT_H_
#define APP_YARA_BITS_CONTEXT_H_

using namespace seqan;

// ============================================================================
// Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class ReadsContext
// ----------------------------------------------------------------------------

template <typename TSpec = void, typename TConfig = void>
struct ReadsContext
{
    String<unsigned char>       seedErrors;
    String<unsigned char>       minErrors;
    String<bool>                mapped;
    String<bool>                paired;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void clear(ReadsContext<TSpec, TConfig> & ctx)
{
    clear(ctx.seedErrors);
    clear(ctx.minErrors);
    clear(ctx.mapped);
    clear(ctx.paired);
    shrinkToFit(ctx.seedErrors);
    shrinkToFit(ctx.minErrors);
    shrinkToFit(ctx.mapped);
    shrinkToFit(ctx.paired);
}

// ----------------------------------------------------------------------------
// Function resize()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TReadSeqs>
inline void resize(ReadsContext<TSpec, TConfig> & ctx, TReadSeqs const & readSeqs)
{
    resize(ctx.seedErrors, getReadSeqsCount(readSeqs), 0, Exact());
    resize(ctx.minErrors, getReadSeqsCount(readSeqs), std::numeric_limits<unsigned char>::max(), Exact());
    resize(ctx.mapped, getReadsCount(readSeqs), false, Exact());
    resize(ctx.paired, getReadsCount(readSeqs), false, Exact());
}

// ----------------------------------------------------------------------------
// Function getSeedErrors()
// ----------------------------------------------------------------------------

template <typename TReadsContext, typename TReadSeqId>
inline unsigned char getSeedErrors(TReadsContext const & ctx, TReadSeqId readSeqId)
{
    return ctx.seedErrors[readSeqId];
}

// ----------------------------------------------------------------------------
// Function setSeedErrors()
// ----------------------------------------------------------------------------

template <typename TReadsContext, typename TReadSeqId, typename TErrors>
inline void setSeedErrors(TReadsContext & ctx, TReadSeqId readSeqId, TErrors errors)
{
    assignValue(ctx.seedErrors, readSeqId, errors);
}

// ----------------------------------------------------------------------------
// Function getMinErrors()
// ----------------------------------------------------------------------------

template <typename TReadsContext, typename TReadId>
inline unsigned char getMinErrors(TReadsContext const & ctx, TReadId readId)
{
    return ctx.minErrors[readId];
}

// ----------------------------------------------------------------------------
// Function setMinErrors()
// ----------------------------------------------------------------------------

template <typename TReadsContext, typename TReadId, typename TErrors>
inline void setMinErrors(TReadsContext & ctx, TReadId readId, TErrors errors)
{
    if (errors < getMinErrors(ctx,readId))
        assignValue(ctx.minErrors, readId, errors);
}

// ----------------------------------------------------------------------------
// Function setMapped()
// ----------------------------------------------------------------------------

template <typename TReadsContext, typename TReadId>
inline void setMapped(TReadsContext & ctx, TReadId readId)
{
    assignValue(ctx.mapped, readId, true);
}

// ----------------------------------------------------------------------------
// Function isMapped()
// ----------------------------------------------------------------------------

template <typename TReadsContext, typename TReadId>
inline bool isMapped(TReadsContext const & ctx, TReadId readId)
{
    return ctx.mapped[readId];
}

// ----------------------------------------------------------------------------
// Function setPaired()
// ----------------------------------------------------------------------------

template <typename TReadsContext, typename TReadId>
inline void setPaired(TReadsContext & ctx, TReadId readId)
{
    assignValue(ctx.paired, readId, true);
}

// ----------------------------------------------------------------------------
// Function isPaired()
// ----------------------------------------------------------------------------

template <typename TReadsContext, typename TReadId>
inline bool isPaired(TReadsContext const & ctx, TReadId readId)
{
    return ctx.paired[readId];
}

#endif  // #ifndef APP_YARA_BITS_CONTEXT_H_
