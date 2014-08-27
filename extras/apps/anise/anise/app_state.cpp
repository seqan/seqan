// ==========================================================================
//                                   ANISE
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
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

#include "app_state.h"

#include <seqan/basic.h>

#include "file_name_tokens.h"
#include "streaming_exception.h"

// ----------------------------------------------------------------------------
// Class AppState
// ----------------------------------------------------------------------------

void AppState::loadOrInit(TemporaryFileManager & tmpMgr)
{
    try
    {
        load(tmpMgr);
    }
    catch (std::runtime_error const & e)
    {
        // In case of error, (re-)initialize.
        numSites = -1;
        stepNo = 0;
        superStep = SuperStep::INITIAL;
        libraryInfos.resize(1);
        save(tmpMgr);
    }
}

void AppState::load(TemporaryFileManager & tmpMgr)
{
    std::fstream f;
    tmpMgr.open(f, std::ios::binary | std::ios::in, GLOBAL_STATE_TOKEN, GLOBAL_STATE_EXT);
    if (!f.good())
        throw AniseIOException() << "Could not open global state file for reading.";

    std::string buffer;
    buffer.resize(1024);
    while (f.good() && f.peek() == '#')  // Skip comments.
        f.getline(&buffer[0], buffer.size());

    std::vector<AppState::SuperStep> vSteps = { AppState::SuperStep::INITIAL,
                                                AppState::SuperStep::ASSEMBLING,
                                                AppState::SuperStep::FINISHING,
                                                AppState::SuperStep::DONE };

    std::vector<BamLibraryInfo::Orientation> vOrientation = { BamLibraryInfo::F_PLUS,
                                                              BamLibraryInfo::F_MINUS,
                                                              BamLibraryInfo::R_PLUS,
                                                              BamLibraryInfo::R_MINUS };
    std::string tmp;

    f >> buffer;  // "NUM_SITES"
    f >> numSites;

    f >> buffer;  // "SUPER_STEP"
    f >> tmp;
    superStep = vSteps[0];
    for (unsigned i = 0; i < vSteps.size(); ++i)
        if (tmp == AppState::superStepStr(vSteps[i]))
            superStep = vSteps[i];

    f >> buffer;  // "STEP_NO"
    f >> stepNo;

    f >> buffer;  // "NUM_ORPHANS"
    f >> numOrphans;

    BamLibraryInfo info;
    f >> buffer;  // "LIBRARY_INFO";
    f >> info.median;
    f >> info.stdDev;
    f >> info.maxNormalISize;
    f >> tmp;
    info.defaultOrient = vOrientation[0];
    for (unsigned i = 0; i < vOrientation.size(); ++i)
        if (tmp == BamLibraryInfo::orientationStr(vOrientation[0]))
            info.defaultOrient = vOrientation[i];
    f >> info.avgReadLen;

    libraryInfos.resize(1);
    libraryInfos[0] = info;
}

void AppState::save(TemporaryFileManager & tmpMgr)
{
    std::fstream f;
    tmpMgr.open(f, std::ios::binary | std::ios::out, GLOBAL_STATE_TOKEN, GLOBAL_STATE_EXT);
    if (!f.good())
        throw AniseIOException() << "Could not open global state file for writing.";

    SEQAN_CHECK(libraryInfos.size() > 0u, "Must have at least one library info.");

    BamLibraryInfo const & info = libraryInfos[0];

    f << "#ANISE GLOBAL APP STATE\n"
      << "NUM_SITES\t" << numSites << "\n"
      << "SUPER_STEP\t" << AppState::superStepStr(superStep) << "\n"
      << "STEP_NO\t" << stepNo << "\n"
      << "NUM_ORPHANS\t" << numOrphans << "\n"
      << "LIBRARY_INFO\t" << info.median << "\t" << info.stdDev << "\t" << info.maxNormalISize
      << "\t" << BamLibraryInfo::orientationStr(info.defaultOrient) << "\t" << info.avgReadLen << "\n";
}
