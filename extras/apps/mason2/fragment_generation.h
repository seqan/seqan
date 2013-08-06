// ==========================================================================
//                         Mason - A Read Simulator
// ==========================================================================
// Copyright (c) 2006-2012, Knut Reinert, FU Berlin
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
// Simulation of the fragmentation step.
//
// DNA is shattered into fragments after PCR.  These fragments are then
// sequenced from one or both sides to form single/paired reads.  At the
// moment, we assume that no errors are introduced in PCR and also ignore any
// biases for PCR or fragment selection.
// ==========================================================================

#ifndef EXTRAS_APPS_MASON2_FRAGMENT_GENERATION_H_
#define EXTRAS_APPS_MASON2_FRAGMENT_GENERATION_H_

#include <sstream>
#include <fstream>
#include <vector>
#include <memory>

#include <seqan/random.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

#include "external_split_merge.h"
#include "mason_options.h"

// ============================================================================
// Forwards
// ============================================================================

typedef seqan::Rng<seqan::MersenneTwister> TRng;

inline void trimAfterSpace(seqan::CharString & s);

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class Fragment
// ----------------------------------------------------------------------------

class Fragment
{
public:
    // Index of contig to simulate on.
    int rId;
    // Begin and end position of the fragment.
    int beginPos, endPos;

    Fragment() : rId(-1), beginPos(0), endPos(0)
    {}
};

// ----------------------------------------------------------------------------
// Class FragmentSamplerImpl
// ----------------------------------------------------------------------------

// Abstract base class for the fragment generation classes.

class FragmentSamplerImpl
{
public:
    virtual void generate(Fragment & frag, int rId, unsigned contigLength) = 0;
    virtual void generateMany(std::vector<Fragment> & frags, int rId, unsigned contigLength,
                              unsigned count) = 0;
};


// ----------------------------------------------------------------------------
// Class UniformFragmentSamplerImpl
// ----------------------------------------------------------------------------

// Generation of fragments with uniform length.

class UniformFragmentSamplerImpl : public FragmentSamplerImpl
{
public:
    // Minimal and maximal fragment length.
    int minLength;
    int maxLength;

    // The random number generator to use.
    TRng & rng;
    // The probability density function for the simulation.
    seqan::Pdf<seqan::Uniform<int> > pdf;

    UniformFragmentSamplerImpl(TRng & rng, int minLength, int maxLength) :
            minLength(minLength), maxLength(maxLength), rng(rng), pdf(minLength, maxLength)
    {}

    virtual void generate(Fragment & frag, int rId, unsigned contigLength);

    virtual void generateMany(std::vector<Fragment> & frags, int rId, unsigned contigLength, unsigned count);

    void _generate(Fragment & frag, int rId, unsigned contigLength);
};

// ----------------------------------------------------------------------------
// Class NormalFragmentSamplerImpl
// ----------------------------------------------------------------------------

// Generation of fragments with normally distributed length.

class NormalFragmentSamplerImpl : public FragmentSamplerImpl
{
public:
    // Mean length and length standard deviation.
    int meanLength;
    int stdDevLength;

    // The random number generator to use.
    TRng & rng;
    // The probability density function for the simulation.
    seqan::Pdf<seqan::Normal> pdf;

    NormalFragmentSamplerImpl(TRng & rng, int meanLength, int stdDevLength) :
            meanLength(meanLength), stdDevLength(stdDevLength), rng(rng), pdf(meanLength, stdDevLength)
    {}

    virtual void generate(Fragment & frag, int rId, unsigned contigLength);

    virtual void generateMany(std::vector<Fragment> & frags, int rId, unsigned contigLength, unsigned count);

    void _generate(Fragment & frag, int rId, unsigned contigLength);
};

// ----------------------------------------------------------------------------
// Class FragmentSampler
// ----------------------------------------------------------------------------

// Generator for sequence fragments.

class FragmentSampler
{
public:
    // Configuration for the generator.
    FragmentSamplerOptions options;

    // The actual generator implementation to use.
    std::auto_ptr<FragmentSamplerImpl> impl;

    FragmentSampler(TRng & rng, FragmentSamplerOptions const & options) : options(options)
    {
        if (options.model == FragmentSamplerOptions::UNIFORM)
            impl.reset(new UniformFragmentSamplerImpl(rng, options.minFragmentSize,
                                                        options.maxFragmentSize));
        else
            impl.reset(new NormalFragmentSamplerImpl(rng, options.meanFragmentSize,
                                                       options.stdDevFragmentSize));
    }

    // Generate a fragment given a contig id and the length of the contig.
    void generate(Fragment & frag, int rId, unsigned contigLength)
    {
        impl->generate(frag, rId, contigLength);
    }

    // Generate multiple fragments.
    void generateMany(std::vector<Fragment> & frags, int rId, unsigned contigLength, unsigned count)
    {
        impl->generateMany(frags, rId, contigLength, count);
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// --------------------------------------------------------------------------
// Function trimAfterSpace()
// --------------------------------------------------------------------------

// TODO(holtgrew): Put into header or SeqAn library.

// Trim after the first whitespace.
inline
void trimAfterSpace(seqan::CharString & s)
{
    unsigned i = 0;
    for (; i < length(s); ++i)
        if (isspace(s[i]))
            break;
    resize(s, i);
}

// TODO(holtgrew): Too much redundancy here.

void UniformFragmentSamplerImpl::generate(Fragment & frag, int rId, unsigned contigLength)
{
    _generate(frag, rId, contigLength);
}

void UniformFragmentSamplerImpl::generateMany(std::vector<Fragment> & frags, int rId,
                                                unsigned contigLength, unsigned count)
{
    frags.resize(count);
    for (unsigned i = 0; i < count; ++i)
        _generate(frags[0], rId, contigLength);
}

void UniformFragmentSamplerImpl::_generate(Fragment & frag, int rId, unsigned contigLength)
{
    int fragLength = 0;
    unsigned const MAX_TRIES = 1000;
    for (unsigned tryNo = 0; tryNo < MAX_TRIES; ++tryNo)
    {
        fragLength = pickRandomNumber(rng, pdf);
        if (fragLength <= 0 || fragLength > (int)contigLength)
            continue;  // Try again

        seqan::Pdf<seqan::Uniform<int> > posPdf(0, contigLength - fragLength);
        int beginPos = pickRandomNumber(rng, posPdf);

        frag.rId = rId;
        frag.beginPos = beginPos;
        frag.endPos = beginPos + fragLength;
        SEQAN_ASSERT_LEQ(frag.endPos, (int)contigLength);
        break;
    }
}

void NormalFragmentSamplerImpl::generate(Fragment & frag, int rId, unsigned contigLength)
{
    _generate(frag, rId, contigLength);
}

void NormalFragmentSamplerImpl::generateMany(std::vector<Fragment> & frags, int rId,
                                               unsigned contigLength, unsigned count)
{
    frags.resize(count);
    for (unsigned i = 0; i < count; ++i)
        _generate(frags[i], rId, contigLength);
}

void NormalFragmentSamplerImpl::_generate(Fragment & frag, int rId, unsigned contigLength)
{
    int fragLength = 0;
    unsigned const MAX_TRIES = 1000;
    for (unsigned tryNo = 0; tryNo < MAX_TRIES; ++tryNo)
    {
        fragLength = pickRandomNumber(rng, pdf);
        if (fragLength <= 0 || fragLength > (int)contigLength)
            continue;  // Try again

        seqan::Pdf<seqan::Uniform<int> > posPdf(0, contigLength - fragLength);
        int beginPos = pickRandomNumber(rng, posPdf);

        frag.rId = rId;
        frag.beginPos = beginPos;
        frag.endPos = beginPos + fragLength;
        SEQAN_ASSERT_LEQ(frag.endPos, (int)contigLength);
        break;
    }
}

#endif  // #ifndef EXTRAS_APPS_MASON2_FRAGMENT_GENERATION_H_
