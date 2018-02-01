// ==========================================================================
//                                 SeqCons
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

#include "seqcons_app.h"

#include <cctype>
#include <iostream>

#include <seqan/realign.h>
#include <seqan/sequence.h>
#include <seqan/store.h>
#include <seqan/seq_io.h>

#include "seqcons_options.h"

namespace {

// ---------------------------------------------------------------------------
// Function endsWithIgnoreCase()
// ---------------------------------------------------------------------------

bool endsWithIgnoreCase(std::string str, std::string suffix)
{
    std::transform(str.begin(), str.end(), str.begin(), ::tolower);
    std::transform(suffix.begin(), suffix.end(), suffix.begin(), ::tolower);
    return seqan::endsWith(str, suffix);
}

// ---------------------------------------------------------------------------
// Function trimAfterSpace()
// ---------------------------------------------------------------------------

void trimAfterSpace(seqan::CharString & s)
{
    unsigned i = 0;
    for (; i < length(s); ++i)
        if (isspace(s[i]))
            break;
    resize(s, i);
}

}

// ----------------------------------------------------------------------------
// Class SeqConsAppImpl
// ----------------------------------------------------------------------------

class SeqConsAppImpl
{
public:
    SeqConsAppImpl(SeqConsOptions const & options) : options(options)
    {}

    void run();

private:
    // Load reads without alignment information from FASTA file.  Will create one pseudo-contig and put all reads at the
    // first position.
    void loadReads(char const * fileName);
    // Load alignments from SAM file, will set reference to string of N of sufficient length.
    void loadAlignments(char const * fileName);

    // Perform the consensus alignment, optionally interpreting coordinates.
    void performConsensusAlignment(bool useContigID, bool usePositions, bool useGlobalAlignment);
    // Perform realignment on store only.
    void performRealignment();

    // Write out consensus sequence to file.
    void writeConsensus();
    // Write out alignments to file.
    void writeAlignments();

    // The fragment store for the data.
    seqan::FragmentStore<> store;

    // Configuration.
    SeqConsOptions options;
};

void SeqConsAppImpl::writeConsensus()
{
    if (options.verbosity >= 1)
        std::cerr << "Writing consensus to " << options.outputFileConsensus << " ...";
    seqan::SeqFileOut seqFileOut;
    if (!open(seqFileOut, options.outputFileConsensus.c_str()))
        throw std::runtime_error("Could not open consensus output file for writing.");
    for (unsigned contigID = 0; contigID < length(store.contigStore); ++contigID)
    {
        std::stringstream ss;
        ss << "consensus_" << contigID;
        writeRecord(seqFileOut, ss.str(), store.contigStore[contigID].seq);
    }
    if (options.verbosity >= 1)
        std::cerr << " OK\n";
}

void SeqConsAppImpl::writeAlignments()
{
    if (options.verbosity >= 1)
        std::cerr << "Writing alignments to " << options.outputFileAlignment << " ...";

    if (endsWithIgnoreCase(options.outputFileAlignment, ".txt"))
    {
        std::fstream out(options.outputFileAlignment.c_str(), std::ios::binary | std::ios::out);
        seqan::AlignedReadLayout layout;
        layoutAlignment(layout, store);
        for (unsigned contigID = 0; contigID < length(store.contigStore); ++contigID)
        {
            int endPos = 0;
            for (unsigned i = 0; i < length(store.alignedReadStore); ++i)
                if (store.alignedReadStore[i].contigId == contigID)
                    endPos = std::max(endPos, (int)store.alignedReadStore[i].endPos);
            out << ">consensus_" << contigID << "\n";
            printAlignment(out, layout, store, /*contigID=*/contigID, /*beginPos=*/0, /*endPos=*/endPos, 0, 100);
        }
    }
    else  // ends in .sam
    {
        seqan::BamFileOut out;
        if (!open(out, options.outputFileAlignment.c_str()))
            throw std::runtime_error("Could not open output file.");
        writeRecords(out, store);
    }
    if (options.verbosity >= 1)
        std::cerr << " OK\n";
}

void SeqConsAppImpl::run()
{
    // Load read or alignment data.
    if (options.verbosity >= 1)
        std::cerr << "\n__LOADING DATA_______________________________________________________________\n"
                  << '\n';
    if (endsWithIgnoreCase(options.inputFile, ".sam"))
        loadAlignments(options.inputFile.c_str());
    else  // sequence file
        loadReads(options.inputFile.c_str());

    if (options.verbosity >= 1)
        std::cerr << "\n__COMPUTATION________________________________________________________________\n"
                  << '\n';
    // Perform the consensus or realignment computation.
    switch (options.operation)
    {
        case SeqConsOptions::ALN_CONSENSUS:
            performConsensusAlignment(false, false, true);
            break;

        case SeqConsOptions::OVL_CONSENSUS:
            performConsensusAlignment(false, false, false);
            break;

        case SeqConsOptions::CTG_CONSENSUS:
            performConsensusAlignment(true, false, false);
            break;

        case SeqConsOptions::POS_CONSENSUS:
            performConsensusAlignment(true, true, false);
            break;

        case SeqConsOptions::REALIGN:
            performRealignment();
            break;

        case SeqConsOptions::NOP:
        default:
            // do nothing, will just write out store
            break;
    }

    // Write the consensus and/or the alignments.
    if (options.verbosity >= 1)
        std::cerr << "\n__WRITING RESULT_____________________________________________________________\n"
                  << '\n';
    if (!options.outputFileConsensus.empty())
        writeConsensus();
    if (!options.outputFileAlignment.empty())
        writeAlignments();
}

void SeqConsAppImpl::loadReads(char const * fileName)
{
    // Allocate space for one contig in the store.
    resize(store.contigStore, 1);
    std::stringstream ns;
    ns << "consensus_1";
    appendValue(store.contigNameStore, ns.str());
    // Load reads from sequence file.
    seqan::SeqFileIn seqFileIn;
    if (options.verbosity >= 1)
        std::cerr << "Loading reads from " << fileName << "...";
    if (!open(seqFileIn, fileName))
        throw std::runtime_error("Problem opening input sequence file for reading.");

    int maxPos = 0;
    seqan::CharString id, seq;
    while (!atEnd(seqFileIn))
    {
        try
        {
            readRecord(id, seq, seqFileIn);
        }
        catch (seqan::ParseError const & e)
        {
            throw std::runtime_error(std::string("Problem reading from input sequence file: ") +
                                     e.what());
        }
        trimAfterSpace(id);
        int64_t readID = appendRead(store, seq, id);
        appendAlignedRead(store, readID, /*contigID=*/0, /*beginPos=*/0, (int)length(seq));
        maxPos = std::max(maxPos, (int)length(seq));
    }

    // Create sequence of Ns with sufficient size.
    resize(store.contigStore[0].seq, 'N');

    if (options.verbosity >= 1)
        std::cerr << "OK\n";
}

void SeqConsAppImpl::loadAlignments(char const * fileName)
{
    if (options.verbosity >= 1)
        std::cerr << "  Loading alignments from " << fileName << "...";
    seqan::BamFileIn bamFileIn;
    if (!open(bamFileIn, fileName))
        throw std::runtime_error("Problem opening input alignment file for reading.");
    try
    {
        readRecords(store, bamFileIn);
    }
    catch (seqan::ParseError const & e)
    {
        throw std::runtime_error(e.what());
    }
    if (options.verbosity >= 1)
        std::cerr << "OK\n";
}

void SeqConsAppImpl::performConsensusAlignment(bool useContigID, bool usePositions, bool useGlobalAlignment)
{
    // Setup the consensus alignment options.
    seqan::ConsensusAlignmentOptions caOptions;
    caOptions.useContigID = useContigID;
    caOptions.usePositions = usePositions;
    caOptions.useGlobalAlignment = useGlobalAlignment;
    caOptions.runRealignment = false;  // will run manually.
    if (options.verbosity >= 3)
        caOptions.verbosity = 3;

    caOptions.overlapMaxErrorRate = (int)options.overlapMaxErrorRate;
    caOptions.overlapMinLength = options.overlapMinLength;
    caOptions.overlapMinCount = options.overlapMinCount;
    caOptions.posDelta = options.overlapWindowSize;

    caOptions.kMerSize = options.kMerSize;
    caOptions.kMerMaxOcc = options.kMerMaxOcc;

    // Perform the consensus alignment.
    double startTime = seqan::sysTime();
    if (options.verbosity >= 1)
        std::cerr << "Performing consensus computation...";
    if (options.verbosity >= 3)
        std::cerr << "\n";
    consensusAlignment(store, caOptions);
    if (options.verbosity >= 1)
        std::cerr << " OK\n";
    if (options.verbosity >= 2)
        std::cerr << "\t=> consensus step took " << seqan::sysTime() - startTime << "s\n";

    // Finally, compute realignment of the resulting multi-read alignment.
    performRealignment();
}

void SeqConsAppImpl::performRealignment()
{
    double startTime = seqan::sysTime();
    if (options.verbosity >= 1)
        std::cerr << "Performing realignment...";
    if (options.verbosity >= 3)
        std::cerr << "\n";
    for (unsigned contigID = 0; contigID < length(store.contigStore); ++contigID)
        reAlignment(store, contigID, /*method=*/2, options.reAlignmentBandwidth, /*includeReference=*/false);
    if (options.verbosity >= 1)
        std::cerr << " OK\n";
    if (options.verbosity >= 2)
        std::cerr << "\t=> realignment took " << seqan::sysTime() - startTime << "s\n";
}

// ----------------------------------------------------------------------------
// Class SeqConsApp
// ----------------------------------------------------------------------------

SeqConsApp::SeqConsApp(SeqConsOptions const & options) : impl(new SeqConsAppImpl(options))
{}

SeqConsApp::~SeqConsApp()
{}

void SeqConsApp::run()
{
    impl->run();
}
