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

#ifndef SANDBOX_HERBARIUM_APPS_ANISE_APP_STATE_H_
#define SANDBOX_HERBARIUM_APPS_ANISE_APP_STATE_H_

#include <vector>

#include "anise/library_info.h"
#include "anise/temporary_file_manager.h"

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class AppState
// ----------------------------------------------------------------------------

class AppState
{
public:

    enum class SuperStep
    {
        INITIAL,
        ASSEMBLING,
        FINISHING,
        DONE
    };

    AppState() : numSites(-1), superStep(SuperStep::INITIAL), stepNo(-1), numOrphans(-1)
    {}

    // Load the global app state from temporary directory if it exists.  Otherwise, create it.
    void loadOrInit(TemporaryFileManager & tmpMgr);
    // Load the global app state from temporary directory.
    void load(TemporaryFileManager & tmpMgr);
    // Save the global state to temporary directory.
    void save(TemporaryFileManager & tmpMgr);

    // Return true if we still need to perform the given step.
    bool inSuperStep(SuperStep ss) const
    {
        return (superStep == ss);
    }

    bool isInitial() const    { return inSuperStep(SuperStep::INITIAL); }
    bool isAssembling() const { return inSuperStep(SuperStep::ASSEMBLING); }
    bool isFinishing() const  { return inSuperStep(SuperStep::FINISHING); }
    bool isDone() const       { return inSuperStep(SuperStep::DONE); }

    static char const * superStepStr(SuperStep ss)
    {
        switch (ss)
        {
            case SuperStep::INITIAL:
                return "INITIAL";
            case SuperStep::ASSEMBLING:
                return "ASSEMBLING";
            case SuperStep::FINISHING:
                return "FINISHING";
            case SuperStep::DONE:
                return "DONE";
        }

        return "<INVALID>";
    }

    // The number of available sites.
    int numSites;
    // The current super step.
    SuperStep superStep;
    // If superStep is ASSEMBLING then the stepNo, -1 otherwise.
    int stepNo;
    // Number of orphans.
    int numOrphans;
    // Library infos.  Currently, only one entry.
    std::vector<BamLibraryInfo> libraryInfos;
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

#endif  // #ifndef SANDBOX_HERBARIUM_APPS_ANISE_APP_STATE_H_
