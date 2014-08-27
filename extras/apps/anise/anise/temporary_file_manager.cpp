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

#include "temporary_file_manager.h"

#include <cstring>
#include <mutex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>

#include <sys/stat.h>  // mkdir()
#include <unistd.h>    // unlink()
#include <sys/types.h>
#include <dirent.h>

// ----------------------------------------------------------------------------
// Class TemporaryFileManagerImpl
// ----------------------------------------------------------------------------

class TemporaryFileManagerImpl
{
public:
    TemporaryFileManagerImpl(char const * tempDir) : tempDir(tempDir), tempDirCreated(false)
    {}

    void init(char const * tempDir = nullptr);

    std::string fileName(char const * token, char const * suffix, int siteID = -1, int stepNo = -1);

    void open(std::fstream & file, std::ios::openmode mode, char const * token, char const * suffix,
              int siteID = -1, int stepNo = -1);

    void reapAll();

    void reapStepFiles(int stepNo);

    void reapSiteFiles(int stepNo, int siteID);

    // Removes temporary directory with all contents.
    void cleanup();

private:

    struct FileRecord
    {
        std::string path;
        std::string token;
        std::string suffix;
        int siteID;
        int stepNo;

        FileRecord(std::string const & path, std::string const & token, std::string const & suffix,
                   int siteID, int stepNo) :
                path(path), token(token), suffix(suffix), siteID(siteID), stepNo(stepNo)
        {}

        bool operator<(FileRecord const & other) const
        {
            auto tpl = [](FileRecord const r) {
                return std::make_tuple(r.path, r.token, r.suffix, r.siteID, r.stepNo);
            };
            return tpl(*this) < tpl(other);
        }
    };

    // Mutex to use for locking file structure.
    std::recursive_mutex mutex;
    // Temporary directory.
    std::string tempDir;
    // Whether or not we already created the temporary directory.
    bool tempDirCreated;
    // List of files that were created.
    std::set<FileRecord> fileRecords;
};

void TemporaryFileManagerImpl::init(char const * _tempDir)
{
    std::lock_guard<std::recursive_mutex> lock(mutex);

    if (tempDirCreated)
        return;  // already done.
    if (tempDir != _tempDir)
        tempDir = _tempDir;
    if (tempDir.empty())
        throw std::runtime_error("Cannot open temporary directory if not given!");

    // Create temporary directory if it does not exist already.
    if (mkdir(tempDir.c_str(), 0700) != 0 && errno != EEXIST)
        throw std::runtime_error("Could not create temporary directory!");
    tempDirCreated = true;
}

std::string TemporaryFileManagerImpl::fileName(char const * token, char const * suffix, int siteID, int stepNo)
{
    std::lock_guard<std::recursive_mutex> lock(mutex);

    std::stringstream ss;
    ss << tempDir << "/";
    if (siteID != -1 && stepNo != -1)
        ss << "site_" << siteID << "_step_" << stepNo << "_" << token << suffix;
    else if (siteID != -1)
        ss << "site_" << siteID << "_" << token << suffix;
    else if (stepNo != -1)
        ss << "step_" << stepNo << "_" << token << suffix;
    else
        ss << token << suffix;

    return ss.str();
}

void TemporaryFileManagerImpl::open(std::fstream & file, std::ios::openmode mode, char const * token, char const * suffix,
                                    int siteID, int stepNo)
{
    std::lock_guard<std::recursive_mutex> lock(mutex);

    if (!tempDirCreated)
        throw std::runtime_error("Cannot open file before creating temporary directory!");

    std::string path = fileName(token, suffix, siteID, stepNo);

    file.open(path.c_str(), mode);
    if (!file.good())
        throw std::runtime_error("Problem opening temporary file.");

    fileRecords.insert(FileRecord(path, token, suffix, siteID, stepNo));
}

void TemporaryFileManagerImpl::reapAll()
{
    std::lock_guard<std::recursive_mutex> lock(mutex);

    for (auto const & record : fileRecords)
        if (remove(record.path.c_str()) != 0)
            throw std::runtime_error("Problem removing temporary file.");
    if (unlink(tempDir.c_str()) != 0)
        throw std::runtime_error("Problem removing temporary directory.");
}

void TemporaryFileManagerImpl::reapStepFiles(int stepNo)
{
    std::lock_guard<std::recursive_mutex> lock(mutex);

    decltype(fileRecords) tmp;
    for (auto & record : fileRecords)
        if (record.stepNo == stepNo)
        {
            if (remove(record.path.c_str()) != 0)
                throw std::runtime_error("Problem removing temporary file.");
        }
        else
        {
            tmp.insert(record);
        }
    swap(tmp, fileRecords);
}

void TemporaryFileManagerImpl::reapSiteFiles(int stepNo, int siteID)
{
    std::lock_guard<std::recursive_mutex> lock(mutex);

    decltype(fileRecords) tmp;
    for (auto & record : fileRecords)
        if (record.stepNo == stepNo && record.siteID == siteID)
        {
            if (remove(record.path.c_str()) != 0)
                throw std::runtime_error("Problem removing temporary file.");
        }
        else
        {
            tmp.insert(record);
        }
    swap(tmp, fileRecords);
}

void TemporaryFileManagerImpl::cleanup()
{
    // List temporary directory and remove contents.
    DIR * dp;
    struct dirent * ep;     
    dp = opendir(tempDir.c_str());
    std::string path;

    if (dp != NULL)
    {
        while ((ep = readdir(dp)))
        {
            if (strcmp(ep->d_name, ".") == 0 || strcmp(ep->d_name, "..") == 0)
                continue;
            path = tempDir + "/" + ep->d_name;
            if (remove(path.c_str()) != 0)
                throw std::runtime_error("Problem removing temporary file!");
        }

        closedir(dp);
    }
    else
    {
        perror ("Couldn't open the directory");
    }

    // Remove temporary directory.
    if (remove(tempDir.c_str()) != 0)
        throw std::runtime_error("Could not remove temporary directory!");
}

// ----------------------------------------------------------------------------
// Class TemporaryFileManager
// ----------------------------------------------------------------------------

TemporaryFileManager::TemporaryFileManager(char const * tempDir) :
        impl(new TemporaryFileManagerImpl(tempDir))
{}

TemporaryFileManager::~TemporaryFileManager()
{}

void TemporaryFileManager::init(char const * tempDir)
{
    impl->init(tempDir);
}

std::string TemporaryFileManager::fileName(char const * token, char const * suffix,
                                           int siteID, int stepNo)
{
    return impl->fileName(token, suffix, siteID, stepNo);
}

void TemporaryFileManager::open(std::fstream & file, std::ios::openmode mode, char const * token, char const * suffix,
                                int siteID, int stepNo)
{
    impl->open(file, mode, token, suffix, siteID, stepNo);
}

void TemporaryFileManager::reapAll()
{
    impl->reapAll();
}

void TemporaryFileManager::reapStepFiles(int stepNo)
{
    impl->reapStepFiles(stepNo);
}

void TemporaryFileManager::reapSiteFiles(int stepNo, int siteID)
{
    impl->reapSiteFiles(stepNo, siteID);
}

void TemporaryFileManager::cleanup()
{
    impl->cleanup();
}
