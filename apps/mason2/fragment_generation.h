// ==========================================================================
//                         Mason - A Read Simulator
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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Simulation of the fragmentation step.
//
// DNA is shattered into fragments after PCR.  These fragments are then
// sequenced from one or both sides to form single/paired reads.  At the
// moment, we assume that no errors are introduced in PCR and also ignore any
// biases for PCR or fragment selection.
// ==========================================================================

#ifndef APPS_MASON2_FRAGMENT_GENERATION_H_
#define APPS_MASON2_FRAGMENT_GENERATION_H_

#include <sstream>
#include <fstream>
#include <vector>
#include <memory>
#include <random>

#include <seqan/sequence.h>
#include <seqan/seq_io.h>

#include "external_split_merge.h"
#include "mason_options.h"

// ============================================================================
// Forwards
// ============================================================================

typedef std::mt19937 TRng;

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
    virtual void generate(Fragment & frag, int rId, unsigned contigLength,
                          std::vector<std::pair<int, int> > const & gapIntervals) = 0;
    virtual void generateMany(std::vector<Fragment> & frags, int rId, unsigned contigLength,
                              std::vector<std::pair<int, int> > const & gapIntervals,
                              unsigned count) = 0;

    virtual ~FragmentSamplerImpl() = default;
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
    // Lower bound on fragment size.
    int fragSizeLowerBound;

    // The random number generator to use.
    TRng & rng;
    // The probability density function for the simulation.
    std::uniform_int_distribution<int> dist;

    UniformFragmentSamplerImpl(TRng & rng, int minLength, int maxLength, int fragSizeLowerBound) :
            minLength(minLength), maxLength(maxLength), fragSizeLowerBound(fragSizeLowerBound),
            rng(rng), dist(minLength, maxLength)
    {}

    virtual void generate(Fragment & frag, int rId, unsigned contigLength,
                          std::vector<std::pair<int, int> > const & gapIntervals);

    virtual void generateMany(std::vector<Fragment> & frags, int rId, unsigned contigLength,
                              std::vector<std::pair<int, int> > const & gapIntervals,
                              unsigned count);

    void _generate(Fragment & frag, int rId, unsigned contigLength,
                   std::vector<std::pair<int, int> > const & gapIntervals);
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
    // Lower bound on fragment size.
    int fragSizeLowerBound;

    // The random number generator to use.
    TRng & rng;
    // The probability density function for the simulation.
    std::normal_distribution<> dist;

    NormalFragmentSamplerImpl(TRng & rng, int meanLength, int stdDevLength, int fragSizeLowerBound) :
            meanLength(meanLength), stdDevLength(stdDevLength), fragSizeLowerBound(fragSizeLowerBound),
            rng(rng), dist(meanLength, stdDevLength)
    {}

    virtual void generate(Fragment & frag, int rId, unsigned contigLength,
                          std::vector<std::pair<int, int> > const & gapIntervals);

    virtual void generateMany(std::vector<Fragment> & frags, int rId, unsigned contigLength,
                              std::vector<std::pair<int, int> > const & gapIntervals, unsigned count);

    void _generate(Fragment & frag, int rId, unsigned contigLength,
                   std::vector<std::pair<int, int> > const & gapIntervals);
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
    std::unique_ptr<FragmentSamplerImpl> impl;

    FragmentSampler(TRng & rng, FragmentSamplerOptions const & options) : options(options)
    {
        if (options.model == FragmentSamplerOptions::UNIFORM)
            impl.reset(new UniformFragmentSamplerImpl(rng, options.minFragmentSize,
                                                      options.maxFragmentSize,
                                                      options.fragSizeLowerBound));
        else
            impl.reset(new NormalFragmentSamplerImpl(rng, options.meanFragmentSize,
                                                     options.stdDevFragmentSize,
                                                     options.fragSizeLowerBound));
    }

    // Generate a fragment given a contig id and the length of the contig.
    void generate(Fragment & frag, int rId, unsigned contigLength,
                  std::vector<std::pair<int, int> > const & gapIntervals)
    {
        impl->generate(frag, rId, contigLength, gapIntervals);
    }

    // Generate multiple fragments.
    void generateMany(std::vector<Fragment> & frags, int rId, unsigned contigLength,
                  std::vector<std::pair<int, int> > const & gapIntervals, unsigned count)
    {
        impl->generateMany(frags, rId, contigLength, gapIntervals, count);
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// --------------------------------------------------------------------------
// Function overlapsWithInterval()
// --------------------------------------------------------------------------

template <class TIter, class T>
TIter _intervalLowerBound (TIter first, TIter last, const T & val)
{
    TIter it;
    typename std::iterator_traits<TIter>::difference_type count, step;
    count = std::distance(first,last);
    while (count>0)
    {
        it = first;
        step=count / 2;
        std::advance(it, step);
        if (it->first < val)
        {
            first=++it;
            count-=step+1;
        }
        else count=step;
    }
    return first;
}

template <class TIter, class T>
TIter _intervalUpperBound (TIter first, TIter last, const T & val)
{
    TIter it;
    typename std::iterator_traits<TIter>::difference_type count, step;
    count = std::distance(first,last);
    while (count > 0)
    {
        it = first;
        step = count / 2;
        std::advance(it, step);
        if (!(val < it->second))
        {
            first = ++it;
            count -= step + 1;
        }
        else count=step;
    }
    return first;
}

inline bool overlaps(int begin0, int end0, int begin1, int end1)
{ return (begin1 < end0 && begin0 < end1); }

// Returns whether the given interval overlaps with another interval.

inline bool overlapsWithInterval(std::vector<std::pair<int, int> > const & intervals,
                                 int beginPos, int endPos)
{
    typedef std::vector<std::pair<int, int> >::const_iterator TIter;
    TIter itL = _intervalLowerBound(intervals.begin(), intervals.end(), beginPos);
    if (itL != intervals.begin())
        --itL;
    TIter itU = _intervalUpperBound(intervals.begin(), intervals.end(), endPos);
    if (itU != intervals.end())
        ++itU;
    for (TIter it = itL; it != itU; ++it)
        if (overlaps(it->first, it->second, beginPos, endPos))
            return true;
    return false;
}


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

void UniformFragmentSamplerImpl::generate(Fragment & frag, int rId, unsigned contigLength,
                                          std::vector<std::pair<int, int> > const & gapIntervals)
{
    _generate(frag, rId, contigLength, gapIntervals);
}

void UniformFragmentSamplerImpl::generateMany(std::vector<Fragment> & frags, int rId, unsigned contigLength,
                                              std::vector<std::pair<int, int> > const & gapIntervals,
                                              unsigned count)
{
    frags.resize(count);
    for (unsigned i = 0; i < count; ++i)
        _generate(frags[0], rId, contigLength, gapIntervals);
}

void UniformFragmentSamplerImpl::_generate(Fragment & frag, int rId, unsigned contigLength,
                                           std::vector<std::pair<int, int> > const & gapIntervals)
{
    int fragLength = 0;
    unsigned const MAX_TRIES = 1000;
    for (unsigned tryNo = 0; tryNo < MAX_TRIES; ++tryNo)
    {
        fragLength = dist(rng);
        if (fragLength <= 0 || fragLength > (int)contigLength ||
            (fragSizeLowerBound && fragSizeLowerBound > fragLength))
            continue;  // Try again

        std::uniform_int_distribution<int> posDist(0, contigLength - fragLength);
        int beginPos = posDist(rng);

        frag.rId = rId;
        frag.beginPos = beginPos;
        frag.endPos = beginPos + fragLength;

        if (overlapsWithInterval(gapIntervals, frag.beginPos, frag.endPos))
            continue;  // overlaps with a gap

        SEQAN_ASSERT_LEQ(frag.endPos, (int)contigLength);
        break;
    }
}

void NormalFragmentSamplerImpl::generate(Fragment & frag, int rId, unsigned contigLength,
                                         std::vector<std::pair<int, int> > const & gapIntervals)
{
    _generate(frag, rId, contigLength, gapIntervals);
}

void NormalFragmentSamplerImpl::generateMany(std::vector<Fragment> & frags, int rId, unsigned contigLength,
                                             std::vector<std::pair<int, int> > const & gapIntervals, unsigned count)
{
    frags.resize(count);
    for (unsigned i = 0; i < count; ++i)
        _generate(frags[i], rId, contigLength, gapIntervals);
}

void NormalFragmentSamplerImpl::_generate(Fragment & frag, int rId, unsigned contigLength,
                                          std::vector<std::pair<int, int> > const & gapIntervals)
{
    int fragLength = 0;
    unsigned const MAX_TRIES = 1000;
    unsigned tryNo = 0;
    for (; tryNo < MAX_TRIES; ++tryNo)
    {
        fragLength = static_cast<int>(dist(rng));
        if (fragLength <= 0 || fragLength > (int)contigLength ||
            (fragSizeLowerBound && fragSizeLowerBound > fragLength))
            continue;  // Try again

        std::uniform_int_distribution<int> posDist(0, contigLength - fragLength);
        int beginPos = posDist(rng);

        frag.rId = rId;
        frag.beginPos = beginPos;
        frag.endPos = beginPos + fragLength;

        if (overlapsWithInterval(gapIntervals, frag.beginPos, frag.endPos))
            continue;  // overlaps with a gap

        SEQAN_ASSERT_LEQ(frag.endPos, (int)contigLength);
        break;
    }

    if (tryNo == MAX_TRIES)
    {
        std::cerr << "WARNING: Tried " << MAX_TRIES << " times to sample fragment from contig of size "
                  << contigLength << ".  Giving up.\n";
    }
}

#endif  // #ifndef APPS_MASON2_FRAGMENT_GENERATION_H_
