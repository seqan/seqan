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

#ifndef SANDBOX_HERBARIUM_APPS_ANISE_TEMPORARY_FILE_MANAGER_H_
#define SANDBOX_HERBARIUM_APPS_ANISE_TEMPORARY_FILE_MANAGER_H_

#include <fstream>  // TODO(holtgrew): Get rid of openmode problem with <iosfwd>
#include <memory>
#include <string>

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// --------------------------------------------------------------------------
// Class TemporaryFileManager
// --------------------------------------------------------------------------

// Provides functionality to open temporary files and removes everything on destruction.

class TemporaryFileManagerImpl;

class TemporaryFileManager
{
public:
    explicit TemporaryFileManager(char const * tempDir = "");
    ~TemporaryFileManager();  // for pimpl

    // Create the temporary directory and go to state ACTIVE.
    void init(char const * tempDir = "");

    // Create a file name, similar to open() but does not open the file.  The file is assumed to be open by the caller
    // and registered in the to be removed files.
    std::string fileName(char const * token, char const * suffix, int siteID = -1, int stepNo = -1);

    // Open a file with a name token, possibly related to a site and step.  An site id or step number of -1 indicates
    // that no such connection exists.
    void open(std::fstream & file, std::ios::openmode mode, char const * token, char const * suffix,
              int siteID = -1, int stepNo = -1);

    // Reap all remaining temporary files.
    void reapAll();

    // Reap all files for the given step.
    void reapStepFiles(int stepNo);

    // Reap the files for the given step and site.
    void reapSiteFiles(int stepNo, int siteID);

    // Remove all files in temporary directory.
    void cleanup();

private:
    std::unique_ptr<TemporaryFileManagerImpl> impl;
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

#endif  // #ifndef SANDBOX_HERBARIUM_APPS_ANISE_TEMPORARY_FILE_MANAGER_H_
