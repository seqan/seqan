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

#include "sequencing.h"

// ===========================================================================
// Class IlluminaSequencingOptions
// ===========================================================================

class IlluminaModel
{
public:
    // Probabilities for a mismatch at a given position.
    seqan::String<double> mismatchProbabilities;

    // Standard deviations for the normal distributions of base qualities for the mismatch case.
    seqan::String<double> mismatchQualityMeans;
    // Standard deviations for the normal distributions of base qualities for the mismatch case.
    seqan::String<double> mismatchQualityStdDevs;

    // Standard deviations for the normal distributions of base qualities for the non-mismatch case.
    seqan::String<double> qualityMeans;
    // Standard deviations for the normal distributions of base qualities for the non-mismatch case.
    seqan::String<double> qualityStdDevs;

    IlluminaModel()
    {}
};

// ===========================================================================
// Class IlluminaSequencingSimulator
// ===========================================================================

// ---------------------------------------------------------------------------
// Constructor IlluminaSequencingSimulator::IlluminaSequencingSimulator
// ---------------------------------------------------------------------------

IlluminaSequencingSimulator::IlluminaSequencingSimulator(TRng & rng,
                                                         TRng & methRng,
                                                         SequencingOptions const & seqOptions,
                                                         IlluminaSequencingOptions const & illuminaOptions) :
        SequencingSimulator(rng, methRng, seqOptions), illuminaOptions(illuminaOptions),
        model(new IlluminaModel())
{
    this->_initModel();
}

// ---------------------------------------------------------------------------
// Function IlluminaSequencingSimulator::_initModel()
// ---------------------------------------------------------------------------

void IlluminaSequencingSimulator::_initModel()
{
    // Compute mismatch probabilities, piecewise linear function.
    resize(model->mismatchProbabilities, illuminaOptions.readLength);
    // Compute probability at raise point.
    double y_r = 2 * illuminaOptions.probabilityMismatch - illuminaOptions.positionRaise * illuminaOptions.probabilityMismatchBegin - illuminaOptions.probabilityMismatchEnd + illuminaOptions.probabilityMismatchEnd * illuminaOptions.positionRaise;
    if (illuminaOptions.verbosity >= 2)
    {
        std::cerr << "Illumina error curve:\n"
                  << "  (0, " << illuminaOptions.probabilityMismatchBegin << ") -- (" << illuminaOptions.positionRaise << ", " << y_r << ") -- (1, " << illuminaOptions.probabilityMismatchEnd << ")\n";
    }
    // std::cout << "y_r = " << y_r << std::endl;
    // Compute mismatch probability at each base.
    if (!empty(illuminaOptions.probabilityMismatchFile))
    {
        // Open file.
        std::fstream file;
        file.open(toCString(illuminaOptions.probabilityMismatchFile), std::ios_base::in);
        if (!file.is_open())
        {
            std::cerr << "Failed to load mismatch probabilities from " << illuminaOptions.probabilityMismatchFile << std::endl;
            // return 1;
        }
        // Load probabilities.
        double x;
        file >> x;
        unsigned i;
        for (i = 0; i < illuminaOptions.readLength && !file.eof(); ++i) {
            model->mismatchProbabilities[i] = x;
            file >> x;
        }
        if (i != illuminaOptions.readLength)
        {
            std::cerr << "Not enough mismatch probabilites in " << illuminaOptions.probabilityMismatchFile << " (" << i << " < " << illuminaOptions.readLength << ")!" << std::endl;
            // return 1;
        }
    } else {
        // Use piecewise linear function for mismatch probability simulation.
        for (unsigned i = 0; i < illuminaOptions.readLength; ++i) {
            double x = static_cast<double>(i) / (illuminaOptions.readLength - 1);
            if (x < illuminaOptions.positionRaise) {
                double b = illuminaOptions.probabilityMismatchBegin;
                double m = (y_r - illuminaOptions.probabilityMismatchBegin) / illuminaOptions.positionRaise;
                model->mismatchProbabilities[i] = m * x + b;
                // std::cout << "model->mismatchProbabilities[" << i << "] = " << model->mismatchProbabilities[i] << std::endl;
            } else {
                double b = y_r;
                double m = (illuminaOptions.probabilityMismatchEnd - y_r) / (1 - illuminaOptions.positionRaise);
                x -= illuminaOptions.positionRaise;
                model->mismatchProbabilities[i] = m * x + b;
                // std::cout << "model->mismatchProbabilities[" << i << "] = " << model->mismatchProbabilities[i] << std::endl;
            }
        }
    }
    if (illuminaOptions.probabilityMismatchScale != 1.0) {
        for (unsigned i = 0; i < illuminaOptions.readLength; ++i)
            model->mismatchProbabilities[i] *= illuminaOptions.probabilityMismatchScale;
    }

    // Compute match/mismatch means and standard deviations.
    resize(model->mismatchQualityMeans, illuminaOptions.readLength);
    for (unsigned i = 0; i < illuminaOptions.readLength; ++i) {
        double b = illuminaOptions.meanMismatchQualityBegin;
        double x = static_cast<double>(i) / (illuminaOptions.readLength - 1);
        double m = (illuminaOptions.meanMismatchQualityEnd - illuminaOptions.meanMismatchQualityBegin);
        model->mismatchQualityMeans[i] = m * x + b;
        // std::cout << "model->mismatchQualityMeans[" << i << "] = " << model->mismatchQualityMeans[i] << std::endl;
    }
    resize(model->mismatchQualityStdDevs, illuminaOptions.readLength);
    for (unsigned i = 0; i < illuminaOptions.readLength; ++i) {
        double b = illuminaOptions.stdDevMismatchQualityBegin;
        double x = static_cast<double>(i) / (illuminaOptions.readLength - 1);
        double m = (illuminaOptions.stdDevMismatchQualityEnd - illuminaOptions.stdDevMismatchQualityBegin);
        model->mismatchQualityStdDevs[i] = m * x + b;
        // std::cout << "model->mismatchQualityStdDevs[" << i << "] = " << model->mismatchQualityStdDevs[i] << std::endl;
    }
    resize(model->qualityMeans, illuminaOptions.readLength);
    for (unsigned i = 0; i < illuminaOptions.readLength; ++i) {
        double b = illuminaOptions.meanQualityBegin;
        double x = static_cast<double>(i) / (illuminaOptions.readLength - 1);
        double m = (illuminaOptions.meanQualityEnd - illuminaOptions.meanQualityBegin);
        model->qualityMeans[i] = m * x + b;
        // std::cout << "model->qualityMeans[" << i << "] = " << model->qualityMeans[i] << std::endl;
    }
    resize(model->qualityStdDevs, illuminaOptions.readLength);
    for (unsigned i = 0; i < illuminaOptions.readLength; ++i) {
        double b = illuminaOptions.stdDevQualityBegin;
        double x = static_cast<double>(i) / (illuminaOptions.readLength - 1);
        double m = (illuminaOptions.stdDevQualityEnd - illuminaOptions.stdDevQualityBegin);
        model->qualityStdDevs[i] = m * x + b;
        // std::cout << "model->qualityStdDevs[" << i << "] = " << model->qualityStdDevs[i] << std::endl;
    }
}

// ---------------------------------------------------------------------------
// Function _simulateRead()
// ---------------------------------------------------------------------------

namespace {

// Simulate the characters that polymorphisms turn into and inserted characters.
//
// Through the usage of ModifiedString, we will always go from the left to the right end.
template <typename TFrag>
void _simulateSequence(TRead & read, TRng & rng, TFrag const & frag,
                       TCigarString const & cigar)
{
    clear(read);

    typedef typename seqan::Iterator<TFrag>::Type TFragIter;
    TFragIter it = begin(frag, seqan::Standard());

    for (unsigned i = 0; i < length(cigar); ++i)
    {
        //unsigned numSimulate = 0;
        if (cigar[i].operation == 'M')
        {
            for (unsigned j = 0; j < cigar[i].count; ++j, ++it)
                appendValue(read, *it);
            continue;
        }
        else if (cigar[i].operation == 'D')
        {
            it += cigar[i].count;
            continue;
        }

        // Otherwise, we have insertions or mismatches.
        for (unsigned j = 0; j < cigar[i].count; ++j)
        {
            // Pick a value between 0 and 1.
            std::uniform_real_distribution<double> dist(0, 1);
            double x = 1.0;
            while (x == 1.0)
                x = dist(rng);
            int num = static_cast<int>(x / 0.25);

            // NOTE: We can only insert CGAT, but we can have a polymorphism to N.

            if (cigar[i].operation == 'I')
                appendValue(read, seqan::Dna5(num));
            else
                appendValue(read, seqan::Dna5(num + (num == ordValue(*it))));
        }

        if (cigar[i].operation == 'X')
            it += cigar[i].count;
    }
}

}  // namespace (anonymous)

// ---------------------------------------------------------------------------
// Function IlluminaSequencingSimulator::simulateRead()
// ---------------------------------------------------------------------------

// Actually simulate read and qualities from fragment and direction forward/reverse strand.
void IlluminaSequencingSimulator::simulateRead(TRead & seq, TQualities & quals, SequencingSimulationInfo & info,
                                               TFragment const & frag, Direction dir, Strand strand)
{
    // std::cerr << "simulateRead(" << (char const *)(dir == LEFT ? "L" : "R") << ", " << (char const *)(strand == FORWARD ? "-->" : "<--") << ")\n";
    // Simulate sequencing operations.
    TCigarString cigar;
    _simulateCigar(cigar);
    unsigned lenInRef = 0;
    _getLengthInRef(lenInRef, cigar);

    if (lenInRef > length(frag))
    {
        throw std::runtime_error("Illumina read is too long, increase fragment length");
    }

    // Simulate sequence (materialize mismatches and insertions).
    typedef seqan::ModifiedString<seqan::ModifiedString<TFragment, seqan::ModView<seqan::FunctorComplement<seqan::Dna5> > >, seqan::ModReverse> TRevCompFrag;
    if ((dir == LEFT) && (strand == FORWARD))
    {
        _simulateSequence(seq, rng, prefix(frag, lenInRef), cigar);
    }
    else if ((dir == LEFT) && (strand == REVERSE))
    {
        seqan::Prefix<TFragment>::Type holder(prefix(frag, lenInRef));
        _simulateSequence(seq, rng, TRevCompFrag(holder), cigar);
    }
    else if ((dir == RIGHT) && (strand == FORWARD))
    {
        _simulateSequence(seq, rng, suffix(frag, length(frag) - lenInRef), cigar);
    }
    else  // ((dir == RIGHT) && (strand == REVERSE))
    {
        seqan::Suffix<TFragment>::Type holder(suffix(frag, length(frag) - lenInRef));
        _simulateSequence(seq, rng, TRevCompFrag(holder), cigar);
    }

    // Simulate qualities.
    _simulateQualities(quals, cigar);
    SEQAN_ASSERT_EQ(length(seq), length(quals));

    // // Reverse qualities if necessary.
    // if (strand == REVERSE)
    //     reverse(quals);

    // Write out extended sequencing information info if configured to do so.  We always write out the sample position
    // and alignment information.
    info.cigar = cigar;
    unsigned len = 0;
    _getLengthInRef(len, cigar);
    info.beginPos = (dir == LEFT) ? beginPosition(frag) : (beginPosition(frag) + length(frag) - len);
    info.isForward = (strand == FORWARD);
    // std::cerr << "  beginPos=" << info.beginPos - beginPosition(frag) << "\n";

    if (seqOptions->embedReadInfo)
    {
        if (dir == LEFT)
            info.sampleSequence = prefix(frag, len);
        else
            info.sampleSequence = suffix(frag, length(frag) - len);
        if (strand == REVERSE)
            reverseComplement(info.sampleSequence);
    }
    // std::cerr << "  frag  =" << frag << "\n";
    // std::cerr << "  fragRC=" << TRevCompFrag(frag) << "\n";
    // std::cerr << "  seq=" << seq << "\tquals=" << quals << "\n";
}

// ---------------------------------------------------------------------------
// Function IlluminaSequencingSimulator::_simulateQualities()
// ---------------------------------------------------------------------------

// Simulate PHRED qualities from the CIGAR string.
void IlluminaSequencingSimulator::_simulateQualities(TQualities & quals, TCigarString const & cigar)
{
    clear(quals);

    unsigned pos = 0;
    for (unsigned i = 0; i < length(cigar); ++i)
    {
        for (unsigned j = 0; j < cigar[i].count; ++j)
        {
            int q = 0;
            if (cigar[i].operation == 'M')
            {
                std::normal_distribution<double> dist(model->qualityMeans[pos], model->qualityStdDevs[pos]);
                q = static_cast<int>(dist(rng));
                ++pos;
            }
            else if (cigar[i].operation == 'I' || cigar[i].operation == 'X')
            {
                std::normal_distribution<double> dist(model->mismatchQualityMeans[pos], model->mismatchQualityStdDevs[pos]);
                q = static_cast<int>(dist(rng));
                ++pos;
            }
            else
            {
                // Deletion/padding, no quality required.
                continue;
            }
            q = std::max(0, std::min(q, 40));  // limit quality to 0..40
            appendValue(quals, (char)('!' + q));
        }
    }
}

// ---------------------------------------------------------------------------
// Function IlluminaSequencingSimulator::_simulateCigar()
// ---------------------------------------------------------------------------

// Simulate CIGAR string.  We can do this with position specific parameters only and thus independent of any
// context.
void IlluminaSequencingSimulator::_simulateCigar(TCigarString & cigar)
{
    clear(cigar);
    unsigned len = this->readLength();
    std::uniform_real_distribution<double> dist(0, 1);

    for (int i = 0; i < (int)len;)
    {
        double x = dist(rng);
        double pMismatch = model->mismatchProbabilities[i];
        double pInsert   = illuminaOptions.probabilityInsert;
        double pDelete   = illuminaOptions.probabilityDelete;
        double pMatch    = 1.0 - pMismatch - pInsert - pDelete;

        // Simulate mutation/insertion/deletion events.  If possible we reuse the last CIGAR entry.  Adjacent
        // insertion/deletion pairs cancel each other out.

        // TODO(holtgrew): No indel at beginning or ending! Same for other simulators!

        if (x < pMatch)  // match
            i += appendOperation(cigar, 'M').first;
        else if (x < pMatch + pMismatch)  // point polymorphism
            i += appendOperation(cigar, 'X').first;
        else if (x < pMatch + pMismatch + pInsert) // insertion
            i += appendOperation(cigar, 'I').first;
        else  // deletion
            i += appendOperation(cigar, 'D').first;
    }
}

// ============================================================================
// Class SequencingSimulatorFactory
// ============================================================================

// ----------------------------------------------------------------------------
// Function SequencingSimulatorFactory::make()
// ----------------------------------------------------------------------------

std::unique_ptr<SequencingSimulator> SequencingSimulatorFactory::make()
{
    std::unique_ptr<SequencingSimulator> res;

    switch (seqOptions.sequencingTechnology)
    {
        case SequencingOptions::ILLUMINA:
            res.reset(new IlluminaSequencingSimulator(rng, methRng, seqOptions, illuminaOptions));
            break;
        case SequencingOptions::SANGER:
            res.reset(new SangerSequencingSimulator(rng, methRng, seqOptions, sangerOptions));
            break;
        case SequencingOptions::ROCHE_454:
            res.reset(new Roche454SequencingSimulator(rng, methRng, seqOptions, roche454Options));
            break;
    }

    return res;
}
