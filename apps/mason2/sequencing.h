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
// Simulator for the sequencing process.
// ==========================================================================

#ifndef APPS_MASON2_SEQUENCING_H_
#define APPS_MASON2_SEQUENCING_H_

#include <stdexcept>
#include <random>

#include <seqan/bam_io.h>
#include <seqan/seq_io.h>
#include <seqan/modifier.h>
#include <seqan/sequence.h>

#include "mason_options.h"
#include "methylation_levels.h"

// ============================================================================
// Forwards
// ============================================================================

class IlluminaModel;
class Roche454Model;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

typedef seqan::Dna5String TRead;
typedef seqan::CharString TQualities;
typedef std::mt19937 TRng;
typedef seqan::Infix<seqan::Dna5String const>::Type TFragment;
typedef seqan::String<seqan::CigarElement<> > TCigarString;

// ----------------------------------------------------------------------------
// Class SequencingSimulationInfo
// ----------------------------------------------------------------------------

// Composition of verbose information for sequencing simulation.
//
// Objects of this type store information such as the CIGAR string and the original sampled sequence for debug and
// evaluation purposes.

struct SequencingSimulationInfo
{
    // Originally sampled sequence, that together with the errors introduced by sequencing gives the read sequence.
    TRead sampleSequence;

    // The CIGAR string, MXID for matches, mismatches, insertions, deletions with respect to the reference.
    TCigarString cigar;

    // Whether or not this comes from the forward strand.
    bool isForward;

    // The contig and haplotype ID.
    int rID, hID;

    // The begin position of the sequence.
    int beginPos;

    // Number of bases covering SNPs and overlapping with indels.
    int snpCount, indelCount;

    SequencingSimulationInfo() : isForward(false), rID(-1), hID(-1), beginPos(-1), snpCount(0), indelCount(0)
    {}

    // Returns number in reference covered by this read.
    int lengthInRef() const
    {
        int result = 0;
        for (unsigned i = 0; i < length(cigar); ++i)
            if (cigar[i].operation == 'M' || cigar[i].operation == 'M' || cigar[i].operation == 'D')
                result += cigar[i].count;
        return result;
    }

    template <typename TStream>
    void serialize(TStream & stream) const
    {
        stream << "SEQUENCE=" << rID << " HAPLOTYPE=" << hID << " BEGIN_POS=" << beginPos
               << " SAMPLE_SEQUENCE=" << sampleSequence << " CIGAR=";
        for (unsigned i = 0; i < length(cigar); ++i)
            stream << cigar[i].count << cigar[i].operation;
        stream << " STRAND=" << (isForward ? 'F' : 'R');
        stream << " NUM_SNPS=" << snpCount << " NUM_INDELS=" << indelCount;
    }
};

// ----------------------------------------------------------------------------
// Class SequencingSimulator
// ----------------------------------------------------------------------------

// Responsible for the simulation of sequencing.
//
// We pick the read length independently of the error since the sequencing simulator might have context dependent
// errors.

class SequencingSimulator
{
public:
    // Sequencing direction from left or right side.
    enum Direction
    {
        LEFT,
        RIGHT
    };

    static Direction direction(Direction s)
    {
        return (s == LEFT) ? RIGHT : LEFT;
    }

    // Strand for sequencing.
    enum Strand
    {
        FORWARD,
        REVERSE
    };

    static Strand toggle(Strand s)
    {
        return (s == FORWARD) ? REVERSE : FORWARD;
    }

    // The random number generators for sequences and methylation/BS-treatment simulation.
    TRng & rng;
    TRng & methRng;
    // Overall sequencing options.
    SequencingOptions const * seqOptions;

    // Buffer for the materialization of BS-seq treated fragments.
    seqan::Dna5String methFrag;

    SequencingSimulator(TRng & rng, TRng & methRng, SequencingOptions const & _options) :
            rng(rng), methRng(methRng), seqOptions(&_options)
    {}

    virtual ~SequencingSimulator() = 0;

    // Pick read length for the fragment to be simulated.
    virtual unsigned readLength() = 0;

    // Simulate paired-end sequencing from a fragment.
    //
    // If BS-seq is enabled in seqOptions->bsSeqOptions then levels must be != 0.
    void simulatePairedEnd(TRead & seqL, TQualities & qualsL, SequencingSimulationInfo & infoL,
                           TRead & seqR, TQualities & qualsR, SequencingSimulationInfo & infoR,
                           TFragment const & frag,
                           MethylationLevels const * levels = 0);

    // Simulate single-end sequencing from a fragment.
    //
    // If BS-seq is enabled in seqOptions->bsSeqOptions then levels must be != 0.
    void simulateSingleEnd(TRead & seq, TQualities & quals, SequencingSimulationInfo & info,
                           TFragment const & frag,
                           MethylationLevels const * levels = 0);

    // Actually simulate read and qualities from fragment and direction forward/reverse strand.
    //
    // seq -- target sequence of the read to simulate
    // quals -- target qualities of the read to simulate
    // frag -- source fragment
    // strand -- the strand of the fragment, coordinates are relative to forward and will be affected by this flag
    // dir -- whether this is the left or right read
    virtual void simulateRead(TRead & seq, TQualities & quals, SequencingSimulationInfo & info,
                              TFragment const & frag, Direction dir, Strand strand) = 0;

private:
    // Simulate BS-seq treatment on forward/reverse strand of frag with the given methylation levels.
    //
    // The result is a DNA string with the translations.
    void _simulateBSTreatment(seqan::Dna5String & methFragment,
                              TFragment const & frag,
                              MethylationLevels const & levels,
                              bool reverse);

    // Implementation of the single-end and paired-end sequencing.  The functions without an underscore forward here.
    void _simulatePairedEnd(TRead & seqL, TQualities & qualsL, SequencingSimulationInfo & infoL,
                            TRead & seqR, TQualities & qualsR, SequencingSimulationInfo & infoR,
                            TFragment const & frag, bool isForward);

    // Simulate single-end sequencing from a fragment.
    //
    // If BS-seq is enabled in seqOptions->bsSeqOptions then levels must be != 0.
    void _simulateSingleEnd(TRead & seq, TQualities & quals, SequencingSimulationInfo & info,
                            TFragment const & frag, bool isForward);
};

inline SequencingSimulator::~SequencingSimulator() {}

// ----------------------------------------------------------------------------
// Class IlluminaSequencingSimulator
// ----------------------------------------------------------------------------

// Illumina read simulation.

class IlluminaSequencingSimulator : public SequencingSimulator
{
public:
    // Configuration for Illumina sequencing.
    IlluminaSequencingOptions illuminaOptions;

    // Storage for the Illumina simulation.
    std::unique_ptr<IlluminaModel> model;

    IlluminaSequencingSimulator(TRng & rng, TRng & methRng, SequencingOptions const & seqOptions,
                                IlluminaSequencingOptions const & illuminaOptions);

    // Pick read length for the fragment to be simulated.
    virtual unsigned readLength()
    {
        return illuminaOptions.readLength;
    }

    // Actually simulate read and qualities from fragment and direction forward/reverse strand.
    virtual void simulateRead(TRead & seq, TQualities & quals, SequencingSimulationInfo & info,
                              TFragment const & frag, Direction dir, Strand strand);

private:
    // Initialize the model.
    void _initModel();

    // Simulate PHRED qualities from the CIGAR string.
    void _simulateQualities(TQualities & quals, TCigarString const & cigar);

    // Simulate CIGAR string.  We can do this with position specific parameters only and thus independent of any
    // context.
    void _simulateCigar(TCigarString & cigar);
};

// ----------------------------------------------------------------------------
// Class Roche454SequencingSimulator
// ----------------------------------------------------------------------------

// 454 read simulation.

class Roche454SequencingSimulator : public SequencingSimulator
{
public:
    // Configuration for Roche454 sequencing.
    Roche454SequencingOptions roche454Options;

    // Precomputed model data for 454 Sequencing.
    std::unique_ptr<Roche454Model> model;

    Roche454SequencingSimulator(TRng & rng, TRng & methRng,
                                SequencingOptions const & seqOptions,
                                Roche454SequencingOptions const & roche454Options);

    // Initialize the threshold matrix.
    void _initModel();

    // Pick read length for the sequence to be sampled from fragments.
    virtual unsigned readLength();

    // Actually simulate read and qualities from fragment and direction forward/reverse strand.
    virtual void simulateRead(TRead & seq, TQualities & quals, SequencingSimulationInfo & info,
                              TFragment const & frag, Direction dir, Strand strand);
};

// ----------------------------------------------------------------------------
// Class SangerSequencingSimulator
// ----------------------------------------------------------------------------

// Sanger read simulation.

class SangerSequencingSimulator : public SequencingSimulator
{
public:
    // Configuration for Sanger sequencing.
    SangerSequencingOptions sangerOptions;

    SangerSequencingSimulator(TRng & rng, TRng & methRng,
                              SequencingOptions const & seqOptions,
                              SangerSequencingOptions const & sangerOptions) :
            SequencingSimulator(rng, methRng, seqOptions), sangerOptions(sangerOptions)
    {}

    // Pick read length for the sequence to be sampled from fragments.
    virtual unsigned readLength();

    // Actually simulate read and qualities from fragment and direction forward/reverse strand.
    virtual void simulateRead(TRead & seq, TQualities & quals, SequencingSimulationInfo & info,
                              TFragment const & frag, Direction dir, Strand strand);

    // Simulate CIGAR string.  We can do this with position specific parameters only and thus independent of any
    // context.
    void _simulateCigar(TCigarString & cigar, unsigned sampleLength);

    void _simulateQualities(TQualities & quals, TCigarString const & cigar, unsigned sampleLength);
};

// ----------------------------------------------------------------------------
// Class SequencingSimulatorFactory
// ----------------------------------------------------------------------------

// Factory for SequencingSimulator objects.op

class SequencingSimulatorFactory
{
public:
    TRng & rng;
    TRng & methRng;

    SequencingOptions const & seqOptions;
    IlluminaSequencingOptions const & illuminaOptions;
    Roche454SequencingOptions const & roche454Options;
    SangerSequencingOptions const & sangerOptions;

    SequencingSimulatorFactory(TRng & rng, TRng & methRng,
                               SequencingOptions const & seqOptions,
                               IlluminaSequencingOptions const & illuminaOptions,
                               Roche454SequencingOptions const & roche454Options,
                               SangerSequencingOptions const & sangerOptions) :
            rng(rng), methRng(methRng), seqOptions(seqOptions), illuminaOptions(illuminaOptions),
            roche454Options(roche454Options), sangerOptions(sangerOptions)
    {}

    std::unique_ptr<SequencingSimulator> make();
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function appendOrientation()
// ----------------------------------------------------------------------------

// Returns (a, b) where a is the difference in resulting read length and b is the difference in used up input/reference
// sequence length.

inline std::pair<int, int> appendOperation(TCigarString & cigar, char op)
{
    // Canceling out of events.  This can only happen if the CIGAR string is not empty.
    if (!empty(cigar) && ((back(cigar).operation == 'I' && op == 'D') ||
                          (back(cigar).operation == 'D' && op == 'I')))
    {
        if (back(cigar).count > 1)
            back(cigar).count -= 1;
        else
            eraseBack(cigar);
        return (op == 'D') ? std::make_pair(-1, 0) : std::make_pair(0, -1);
    }

    // No canceling out of events.  The read length increases by one if the operation is no deletion and one base of
    // input sequence is used up if the operation is not an insertion.
    if (!empty(cigar) && back(cigar).operation == op)
        back(cigar).count += 1;
    else
        appendValue(cigar, seqan::CigarElement<>(op, 1));
    return std::make_pair((op != 'D'), (op != 'I'));
}

#endif  // #ifndef APPS_MASON2_SEQUENCING_H_
