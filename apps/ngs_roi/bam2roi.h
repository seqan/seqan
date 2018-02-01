// ==========================================================================
//                   NGS: Regions of Interest Analysis
// ==========================================================================
// Copyright (c) 2012-2018, Bernd Jagla, Institut Pasteur
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
// Author: Bernd Jagla <bernd.jagla@pasteur.fr>
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================


#ifndef SANDBOX_MY_SANDBOX_APPS_SAMBAMSTATS_BAM2ROI_H_
#define SANDBOX_MY_SANDBOX_APPS_SAMBAMSTATS_BAM2ROI_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

#define VERSION "0.1"


// ----------------------------------------------------------------------------
// Class Options
// ----------------------------------------------------------------------------

// Options for the bam2roi Application.

struct Options
{
    // Verbosity of output: 0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    int verbosity;

    // Paths to input read files.
    seqan::CharString inputFileName;

    // Paths to output read files of accepted reads.
    seqan::CharString outputFileName;

	// if true the experiment is strandspecific.
	bool strandspecific;
    // Whether or not to expect paired-end data, shortcut to !empty(inputFileName2).
	// TODO: Do we need (see TODO file)
    //bool pairedEnd() const
    //{
    //    return !empty(inputFileName);
    //}

    Options() :
        verbosity(0), strandspecific(false)
    {}
};
#endif  // #ifndef SANDBOX_MY_SANDBOX_APPS_SAMBAMSTATS_BAM2ROI_H_
