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
// Management of fragment/read-to-content distribution.
//
// Simulating contig-wise saves memory to a comfortable amount and is easy.
// We first simulate for the fragments 0..(n-1) from which contig they come
// from write their ids to one file per contig.  In a second step, we read
// contig-wise through the id files and simulate the fragments/reads.  In a
// final step, we merge the fragments/reads by their id and write them to
// the final file.
//
// This header provides the data structures and routines to manage this.
// ==========================================================================

#ifndef APPS_MASON2_EXTERNAL_SPLIT_MERGE_H_
#define APPS_MASON2_EXTERNAL_SPLIT_MERGE_H_

#include <vector>
#include <iostream>

#include <seqan/bam_io.h>
#include <seqan/seq_io.h>

#include "mason_types.h"

// ============================================================================
// Forwards
// ============================================================================

inline bool ltBamAlignmentRecord(seqan::BamAlignmentRecord const & lhs,
                                 seqan::BamAlignmentRecord const & rhs);
inline int strnum_cmp(const char *a, const char *b);

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// --------------------------------------------------------------------------
// Class ContigPicker
// --------------------------------------------------------------------------

// Distribute to contig and haplotypes.
//
// Contigs are picked with a probability proportional to their length and haplotypes are picked uniformly at random.

class ContigPicker
{
public:
    // The random number generator to use.
    TRng & rng;

    // The length of the contigs.
    std::vector<int64_t> lengthSums;
    // The number of haplotypes.
    int numHaplotypes;

    ContigPicker(TRng & rng) : rng(rng)
    {}

    // Return a position (contig, haplotype) to distribute the read to.
    std::pair<int, int> pick();

    // Convert a position (contig, haplotype) to an integer.
    int toId(std::pair<int, int> pos) const
    {
        return pos.first * numHaplotypes + pos.second;
    }
};

// ----------------------------------------------------------------------------
// Class IdSplitter
// ----------------------------------------------------------------------------

// Allows distributing ids from/to files.
//
// General protocol:
//
// * construct
// * open()
// * write ids, splitting
// * reset()
// * read ids, contig-wise
// * close()

// TODO(holtgrew): Name bogus, FileBundle would be better.

class IdSplitter
{
public:
    // The number of contigs to split to.
    unsigned numContigs;

    // The file pointers for each contig.
    std::vector<std::fstream *> files;
    // The names of the temporary files (required on Windows).
    std::vector<std::string> fileNames;

    IdSplitter() : numContigs(0)
    {}

    IdSplitter(unsigned numContigs) : numContigs(numContigs)
    {}

    ~IdSplitter()
    {
        close();
    }

    // Open files in the splitter.
    void open();

    // Reset all files in the splitter, ready for reading.
    void reset();

    // Close splitter.
    void close();
};

// ----------------------------------------------------------------------------
// Class FastxJoiner
// ----------------------------------------------------------------------------

// Allows joining by id name from FASTA data stored in a IdSplitter.
//
// Construct with IdSplitter after reset() call.

// TODO(holtgrew): Could use a heap/tournament tree.

template <typename TTag>
class FastxJoiner
{
public:
    // The type of the input iterator to use.
    typedef typename seqan::DirectionIterator<std::fstream, seqan::Input>::Type TInputIterator;

    // The IdSplitter to use.
    IdSplitter * splitter;
    // Number of active files.
    unsigned numActive;
    // Buffer for id and sequence for each input file.
    seqan::StringSet<seqan::CharString> ids, seqs, quals;
    // Maps files for activeness.
    std::vector<bool> active;
    // Input iterators, one for each input file.
    std::vector<TInputIterator> inputIterators;

    FastxJoiner() : splitter(), numActive(0)
    {}

    FastxJoiner(IdSplitter & splitter) : splitter(&splitter), numActive(0)
    {
        _init();
    }

    void _init();

    template <typename TSeq>
    bool _loadNext(TSeq & id, TSeq & seq, TSeq & qual, unsigned idx);

    bool atEnd() const
    {
        return (numActive == 0);
    }

    int get(seqan::CharString & id, seqan::CharString & seq, seqan::CharString & qual);
};

// ----------------------------------------------------------------------------
// Class SamJoiner
// ----------------------------------------------------------------------------

// Allows joining by id name from FASTA data stored in a IdSplitter.
//
// Construct with IdSplitter after reset() call.

// TODO(holtgrew): Could use a heap/tournament tree.

// Compare two BAM alignment records by query name, tie is broken by first/last flag, first < last.

class SamJoiner
{
public:
    // The IdSplitter to use.
    IdSplitter * splitter;
    // Number of active files.
    unsigned numActive;
    // Buffer for id and sequence for each input file.
    seqan::String<seqan::BamAlignmentRecord> records;
    // Maps files for activeness.
    std::vector<bool> active;
    // Input BAM files, one for each input file.
    std::vector<seqan::BamFileIn *> bamFileIns;

    // One of the identical BAM headers.
    seqan::BamHeader header;

    SamJoiner() : splitter(), numActive(0)
    {}

    SamJoiner(IdSplitter & splitter, seqan::BamFileOut * outPtr) :
            splitter(&splitter), numActive(0)
    {
        init(outPtr);
    }

    ~SamJoiner()
    {
        for (unsigned i = 0; i < bamFileIns.size(); ++i)
            delete bamFileIns[i];
    }

    void init(seqan::BamFileOut * outPtr);

    bool _loadNext(seqan::BamAlignmentRecord & record, unsigned idx);

    bool atEnd() const
    {
        return (numActive == 0);
    }

    // Get next BAM alignment record to lhs.  If it is paired-end, load the second mate as well.
    int get(seqan::BamAlignmentRecord & record);
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function FastxJoiner::_init()
// ----------------------------------------------------------------------------

template <typename TTag>
void FastxJoiner<TTag>::_init()
{
    resize(ids, splitter->files.size());
    resize(seqs, splitter->files.size());
    resize(quals, splitter->files.size());
    active.resize(splitter->files.size());

    for (unsigned i = 0; i < splitter->files.size(); ++i)
    {
        inputIterators.push_back(directionIterator(*splitter->files[i], seqan::Input()));
        active[i] = _loadNext(ids[i], seqs[i], quals[i], i);
        numActive += (active[i] != false);
    }
}

// ----------------------------------------------------------------------------
// Function FastxJoiner::_loadNext()
// ----------------------------------------------------------------------------

template <typename TTag>
template <typename TSeq>
bool FastxJoiner<TTag>::_loadNext(TSeq & id, TSeq & seq, TSeq & qual, unsigned idx)
{
    if (seqan::atEnd(inputIterators[idx]))
        return false;
    readRecord(id, seq, qual, inputIterators[idx], TTag());
    return true;
}

// ----------------------------------------------------------------------------
// Function FastxJoiner::get()
// ----------------------------------------------------------------------------

template <typename TTag>
int FastxJoiner<TTag>::get(seqan::CharString & id, seqan::CharString & seq, seqan::CharString & qual)
{
    unsigned idx = std::numeric_limits<unsigned>::max();
    for (unsigned i = 0; i < length(ids); ++i)
    {
        if (!active[i])
            continue;
        if (idx == std::numeric_limits<unsigned>::max() || strnum_cmp(toCString(ids[i]), toCString(ids[idx])) < 0)
            idx = i;
    }
    if (idx == std::numeric_limits<unsigned>::max())
        return 1;

    // We use double-buffering and the input parameters as buffers.
    active[idx] = _loadNext(id, seq, qual, idx);
    swap(id, ids[idx]);
    swap(seq, seqs[idx]);
    swap(qual, quals[idx]);
    numActive -= !active[idx];

    return 0;
}

// ----------------------------------------------------------------------------
// Function ltBamAlignmentRecord()
// ----------------------------------------------------------------------------

// Original comparison function for strings by Heng Li from samtools.

/* The MIT License

   Copyright (c) 2008-2018 Genome Research Ltd (GRL).

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

inline int strnum_cmp(const char *a, const char *b)
{
    char *pa, *pb;
    pa = (char*)a; pb = (char*)b;
    while (*pa && *pb) {
        if (isdigit(*pa) && isdigit(*pb)) {
            long ai, bi;
            ai = strtol(pa, &pa, 10);
            bi = strtol(pb, &pb, 10);
            if (ai != bi) return ai<bi? -1 : ai>bi? 1 : 0;
        } else {
            if (*pa != *pb) break;
            ++pa; ++pb;
        }
    }
    if (*pa == *pb)
        return (pa-a) < (pb-b)? -1 : (pa-a) > (pb-b)? 1 : 0;
    return *pa<*pb? -1 : *pa>*pb? 1 : 0;
}

inline bool ltBamAlignmentRecord(seqan::BamAlignmentRecord const & lhs,
                                 seqan::BamAlignmentRecord const & rhs)
{
    int res = strnum_cmp(toCString(lhs.qName), toCString(rhs.qName));
    return (res < 0) || (res == 0 && hasFlagFirst(lhs));
}

#endif  // #ifndef APPS_MASON2_EXTERNAL_SPLIT_MERGE_H_
