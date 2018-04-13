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
// Class SangerSequencingSimulator
// ===========================================================================

// ---------------------------------------------------------------------------
// Helper Function _simulateSequence()
// ---------------------------------------------------------------------------

namespace {

// TODO(holtgrew): Copy-and-paste from IlluminaSequencingSimulator.
//
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
    }
}

}  // namespace (anonymous)

// ---------------------------------------------------------------------------
// Function SangerSequencingSimulator::readLength()
// ---------------------------------------------------------------------------

unsigned SangerSequencingSimulator::readLength()
{
    if (sangerOptions.readLengthIsUniform)
    {
        // Pick uniformly.
        double minLen = sangerOptions.readLengthMin;
        double maxLen = sangerOptions.readLengthMax;
        std::uniform_real_distribution<double> dist(minLen, maxLen);
        double len = dist(rng);
        return static_cast<unsigned>(round(len));
    }
    else
    {
        // Pick normally distributed.
        std::normal_distribution<double> dist(sangerOptions.readLengthMean, sangerOptions.readLengthError);
        double len = dist(rng);
        return static_cast<unsigned>(round(len));
    }
}

// ---------------------------------------------------------------------------
// Function SangerSequencingSimulator::simulateRead()
// ---------------------------------------------------------------------------

// Actually simulate read and qualities from fragment and direction forward/reverse strand.
void SangerSequencingSimulator::simulateRead(
        TRead & seq, TQualities & quals, SequencingSimulationInfo & info,
        TFragment const & frag, Direction dir, Strand strand)
{
    // TODO(holtgrew): Pick better names for read length here.
    // Pick sampled length.
    unsigned sampleLength = this->readLength();

    if (sampleLength > length(frag))
    {
        throw std::runtime_error("Sanger read is too long, increase fragment length");
    }

    // Simulate CIGAR string.
    TCigarString cigar;
    this->_simulateCigar(cigar, sampleLength);

    // Simulate sequence (materialize mismatches and insertions).
    typedef seqan::ModifiedString<seqan::ModifiedString<TFragment, seqan::ModView<seqan::FunctorComplement<seqan::Dna5> > >, seqan::ModReverse> TRevCompFrag;
    if ((dir == LEFT) && (strand == FORWARD))
    {
        _simulateSequence(seq, rng, prefix(frag, sampleLength), cigar);
    }
    else if ((dir == LEFT) && (strand == REVERSE))
    {
        seqan::Prefix<TFragment>::Type holder(prefix(frag, sampleLength));
        _simulateSequence(seq, rng, TRevCompFrag(holder), cigar);
    }
    else if ((dir == RIGHT) && (strand == FORWARD))
    {
        _simulateSequence(seq, rng, suffix(frag, length(frag) - sampleLength), cigar);
    }
    else  // ((dir == RIGHT) && (strand == REVERSE))
    {
        seqan::Suffix<TFragment>::Type holder(suffix(frag, length(frag) - sampleLength));
        _simulateSequence(seq, rng, TRevCompFrag(holder), cigar);
    }

    // Simulate Qualities.
    this->_simulateQualities(quals, cigar, sampleLength);

    // Reverse qualities if necessary.
    if (strand == REVERSE)
        reverse(quals);

    // Write out sequencing information info if configured to do so.
    if (seqOptions->embedReadInfo)
    {
        info.cigar = cigar;
        unsigned len = 0;
        _getLengthInRef(len, cigar);
        if (dir == LEFT)
            info.sampleSequence = prefix(frag, len);
        else
            info.sampleSequence = suffix(frag, length(frag) - len);
        info.isForward = (strand == FORWARD);
        if (strand == REVERSE)
            reverseComplement(info.sampleSequence);
    }
}

// ---------------------------------------------------------------------------
// Function SangerSequencingSimulator::_simulateCigar()
// ---------------------------------------------------------------------------

// Simulate CIGAR string.  We can do this with position specific parameters only and thus independent of any
// context.
void SangerSequencingSimulator::_simulateCigar(TCigarString & cigar, unsigned sampleLength)
{
    clear(cigar);
    std::uniform_real_distribution<double> dist(0, 1);

    for (unsigned i = 0; i < sampleLength;)
    {
        double x = dist(rng);
        double pos = 1.0 * i / (sampleLength - 1);
        double pMismatch = sangerOptions.probabilityMismatchBegin + pos * (sangerOptions.probabilityMismatchEnd - sangerOptions.probabilityMismatchBegin);
        double pInsert   = sangerOptions.probabilityInsertBegin + pos * (sangerOptions.probabilityInsertEnd - sangerOptions.probabilityInsertBegin);
        double pDelete   = sangerOptions.probabilityDeleteBegin + pos * (sangerOptions.probabilityDeleteEnd - sangerOptions.probabilityDeleteBegin);
        double pMatch    = 1.0 - pMismatch - pInsert - pDelete;

        // Simulate mutation/insertion/deletion events.  If possible we reuse the last CIGAR entry.  Adjacent
        // insertion/deletion pairs cancel each other out.  We count i up to the input read length, thus using the
        // "second" member of appendOperation()'s result.

        if (x < pMatch)  // match
        {
            i += appendOperation(cigar, 'M').second;
        }
        else if (x < pMatch + pMismatch)  // point polymorphism
        {
            i += appendOperation(cigar, 'X').second;
        }
        else if (x < pMatch + pMismatch + pInsert) // insertion
        {
            if (i == 0 || i + 1 == sampleLength)
                continue;  // no indel at beginning/end
            i += appendOperation(cigar, 'I').second;
        }
        else  // deletion
        {
            if (i == 0 || i + 1 == sampleLength)
                continue;  // no indel at beginning/end
            i += appendOperation(cigar, 'D').second;
        }
    }
}

// ---------------------------------------------------------------------------
// Function SangerSequencingSimulator::_simulateQualities()
// ---------------------------------------------------------------------------

void SangerSequencingSimulator::_simulateQualities(
        TQualities & quals, TCigarString const & cigar, unsigned sampleLength)
{
    clear(quals);

    unsigned pos = 0;   // Position in result.
    unsigned rPos = 0;  // Position in fragment.
    for (unsigned i = 0; i < length(cigar); ++i)
    {
        for (unsigned j = 0; j < cigar[i].count; ++j, ++pos)
        {
            double mean = 0.0, stdDev = 0.0;
            double relPos = 1.0 * rPos / (1.0 * sampleLength);

            rPos += (cigar[i].operation != 'D');

            if (cigar[i].operation == 'D')
            {
                continue;  // No quality to give out.
            }
            else if (cigar[i].operation == 'I' || cigar[i].operation == 'X')
            {
                mean = sangerOptions.qualityMatchStartMean + relPos * (sangerOptions.qualityMatchEndMean - sangerOptions.qualityMatchStartMean);
                stdDev = sangerOptions.qualityMatchStartStdDev + relPos * (sangerOptions.qualityMatchEndStdDev - sangerOptions.qualityMatchStartStdDev);
            }
            else  // cigar[i].operation == 'M'
            {
                mean = sangerOptions.qualityErrorStartMean + relPos * (sangerOptions.qualityErrorEndMean - sangerOptions.qualityErrorStartMean);
                stdDev = sangerOptions.qualityErrorStartStdDev + relPos * (sangerOptions.qualityErrorEndStdDev - sangerOptions.qualityErrorStartStdDev);
            }

            std::normal_distribution<double> dist(mean, stdDev);
            int q = static_cast<int>(dist(rng));
            q = std::max(0, std::min(40, q));
            appendValue(quals, (char)('!' + q));
        }
    }
}
