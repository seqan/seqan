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

#include "assembly_substep.h"

#include <atomic>
#include <iostream>
#include <mutex>
#include <vector>
#include <thread>

#include "shared/progress_indicator.h"

#include "anise/anise_options.h"
#include "anise/app_state.h"
#include "anise/assembler.h"
#include "anise/parallel_utils.h"
#include "anise/realignment.h"
#include "anise/site_data.h"
#include "anise/site_state.h"
#include "anise/store_utils.h"
#include "anise/temporary_file_manager.h"
#include "anise/time_log.h"

// --------------------------------------------------------------------------
// Class AssemblySubstepImpl
// --------------------------------------------------------------------------

class AssemblySubstepImpl
{
public:

    AssemblySubstepImpl(TemporaryFileManager & tmpMgr, AniseOptions const & options) :
            tmpMgr(tmpMgr), options(options)/*, overlapper(options), msa(tmpMgr, options), scaffolder(appState, options)*/
    {
        // Load application state.
        appState.load(tmpMgr);
    }

    void run()
    {
        // Progress bar.
        ProgressBar pb(std::cerr, 0, numActiveSites(), (options.verbosity == AniseOptions::NORMAL));
        pb.setLabel("  assembling sites");
        pb.updateDisplay();

        // Job generator to atomically generate jobs.
        JobGenerator gen(appState.numSites);

        // This lambda contains the code for one thread.  It grabs the next site and performs the assembly for it.
        auto assembleLoop = [&](){
            int siteID = -1;
            // Perform assembly for each site.
            while ((siteID = gen()) != -1)
            {
                if (!siteActive(siteID))
                    continue;  // Skip sites marked as inactive already.
                if (options.debugSiteID != -1 && options.debugSiteID != siteID)
                    continue;  // Skip, not selected.

                // Load previous site data.
                SiteData prevSiteData;
                prevSiteData.load(tmpMgr, siteID, appState.stepNo - 1);

                // Stop assembly if collected too many reads.
                std::string message;
                if (doStop(message, prevSiteData))
                {
                    deactivateSite(siteID, message.c_str());
                    TimeLog::instance().log("STOPPING SITE", appState.stepNo, siteID);
                    pb.advanceBy(1);
                    continue;
                }

                // Construct current site data using the assembler.
                SiteData siteData;
                withTimeLog("ASSEMBLY", appState.stepNo, siteID, [&]() {
                        performAssemblyStep(siteData, prevSiteData, options, appState, tmpMgr);
                        siteData.state.stepNo += 1;
                    });

                // Save current step's data and update the progress display.
                siteData.save(tmpMgr);
                pb.advanceBy(1);
            }
        };

        forkJoin(options.numThreads, assembleLoop);

        pb.finish();
    }

private:

    // Whether or not to stop given the SiteData.
    bool doStop(std::string & message, SiteData const & siteData) const
    {
        std::stringstream ss;

        // Check whether too many in absolute numbers.
        if (siteData.readSet.bamRecords.size() > (unsigned)(options.stopReadCount))
        {
            ss << "deactivating: too many alignments " << siteData.readSet.bamRecords.size()
               << " > " << options.stopReadCount;
            message = ss.str();
            return true;
        }

        // Count number of freshly mapped reads.
        unsigned numNewReads = 0;
        int thisStepNo = appState.stepNo - 1;
        std::set<seqan::CharString> seen;
        for (auto const & record : siteData.readSet.bamRecords)
        {
            if (seen.count(record.qName))
                continue;
            seqan::BamTagsDict tagsDict(const_cast<seqan::CharString &>(record.tags));
            unsigned idx = 0;
            int mappedStep = 0;
            if (!findTagKey(idx, tagsDict, "mS") || !extractTagValue(mappedStep, tagsDict, idx))
                continue;  // no such step
            if (mappedStep == thisStepNo)
            {
                numNewReads += 2;
                seen.insert(record.qName);
            }
        }
        if ((int)numNewReads > options.stopInitialReadCount)  // TODO(holtgrew): Dedicated option.
        {
            // std::cerr << "new alignments " << numNewReads << "\n";
            ss << "deactivating: too many new alignments " << numNewReads
               << " > " << options.stopInitialReadCount;
            message = ss.str();
            return true;
        }

        return false;
    }

    // Return number of active sites.
    //
    // This is done by looking at each site state file.
    int numActiveSites()
    {
        int result = 0;

        // We need to lock for access to TemporaryFileManager.
        std::lock_guard<std::mutex> lock(fileSystemMutex);

        for (int siteID = 0; siteID < appState.numSites; ++siteID)
        {

            AssemblySiteState siteState;
            siteState.load(tmpMgr, siteID);
            result += !!siteState.active;
        }

        return result;
    }

    // Return true if the site has been marked as not active already.
    bool siteActive(int siteID)
    {
        // We need to lock for access to TemporaryFileManager.
        std::lock_guard<std::mutex> lock(fileSystemMutex);

        AssemblySiteState siteState;
        siteState.load(tmpMgr, siteID);
        return siteState.active;
    }

    // Deactivate site in AssemblySiteState on disk.
    void deactivateSite(int siteID, char const * comment = "")
    {
        // We need to lock for access to TemporaryFileManager.
        std::lock_guard<std::mutex> lock(fileSystemMutex);

        AssemblySiteState siteState;
        siteState.load(tmpMgr, siteID);
        if (strlen(comment) != 0)
            siteState.comment = comment;
        siteState.active = false;
        siteState.save(tmpMgr);
    }

    // The application state.
    AppState appState;
    // TemporaryFileManager for opening the files in the temporary directory.
    TemporaryFileManager & tmpMgr;
    // Application configuration.
    AniseOptions const & options;
    // Mutex for file access.
    std::mutex fileSystemMutex;
};

// --------------------------------------------------------------------------
// Class AssemblySubstep
// --------------------------------------------------------------------------

AssemblySubstep::AssemblySubstep(TemporaryFileManager & tmpMgr, AniseOptions const & options) :
        impl(new AssemblySubstepImpl(tmpMgr, options))
{}

AssemblySubstep::~AssemblySubstep()
{}

void AssemblySubstep::run()
{
    impl->run();
}
