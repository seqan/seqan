// ==========================================================================
//                         Mason - A Read Simulator
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

#include "external_split_merge.h"

// ---------------------------------------------------------------------------
// Function IdSplitter::open()
// ---------------------------------------------------------------------------

void IdSplitter::open()
{
    close();
    for (unsigned i = 0; i < numContigs; ++i)
    {
#if defined(PLATFORM_WINDOWS_VS)
	    char fileNameBuffer[1000];
		char filePathBuffer[1000];
		//  Gets the temp path env string (no guarantee it's a valid path).
		DWORD dwRetVal = 0;
		dwRetVal = GetTempPath(1000,            // length of the buffer
							   filePathBuffer); // buffer for path
		if (dwRetVal > 1000 || (dwRetVal == 0))
		{
			std::cerr << "GetTempPath failed" << std::endl;
			exit(1);
		}

		UINT uRetVal   = 0;
		uRetVal = GetTempFileName(filePathBuffer,   // directory for tmp files
								  TEXT("MASON_"),   // temp file name prefix
								  0,                // create unique name
								  fileNameBuffer);  // buffer for name

		if (uRetVal == 0)
		{
			std::cerr << "GetTempFileName failed" << std::endl;
			exit(1);
		}

		// The file exists now and is already closed.

		DeleteFile(fileNameBuffer);
		files.push_back(fopen(fileNameBuffer, "w+b"));
#elif defined(PLATFORM_WINDOWS_MINGW)
		char fileNameBuffer[1000];
		tmpnam(fileNameBuffer);
        files.push_back(fopen(fileNameBuffer, "w+b"));
        remove(fileNameBuffer);
#else  // POSIX (Linux/Mac)
		// Create temporary file using POSIX API, open with <cstdio>, delete, close POSIX.
        char const * tmpdir = 0;
        if ((tmpdir = getenv("TMPDIR")) == NULL)
            tmpdir = "/tmp";
        std::string pathTpl = tmpdir;
        pathTpl += "/MASON_XXXXXX";

        int fd = mkstemp(&pathTpl[0]);
        files.push_back(fopen(pathTpl.c_str(), "w+b"));
        remove(pathTpl.c_str());
        ::close(fd);
#endif

        if (!files.back())
        {
            std::cerr << "ERROR: Could not open temporary file!\n";
            exit(1);
        }
    }
}

// ---------------------------------------------------------------------------
// Function IdSplitter::reset()
// ---------------------------------------------------------------------------

void IdSplitter::reset()
{
    for (unsigned i = 0; i < files.size(); ++i)
        if (files[i] != 0)
        {
            SEQAN_ASSERT(!ferror(files[i]));
            fflush(files[i]);
            int res = fseek(files[i], 0L, SEEK_SET);
            (void)res;
            SEQAN_ASSERT_EQ(res, 0);
            SEQAN_ASSERT(!ferror(files[i]));
        }
}

// ---------------------------------------------------------------------------
// Function IdSplitter::close()
// ---------------------------------------------------------------------------

void IdSplitter::close()
{
    for (unsigned i = 0; i < files.size(); ++i)
        if (files[i])
        {
            fclose(files[i]);
            files[i] = 0;
        }
    files.clear();
}

// ---------------------------------------------------------------------------
// Function SamJoiner::init()
// ---------------------------------------------------------------------------

void SamJoiner::init()
{
    resize(records, splitter->files.size());
    active.resize(splitter->files.size());

    for (unsigned i = 0; i < splitter->files.size(); ++i)
    {
        readers.push_back(new TReader(splitter->files[i]));

        // We use a separate header structure and name stores and caches.  Since the headers of all files are equal, we
        // will write out the first one only.
        seqan::BamHeader tmpHeader;
        seqan::StringSet<seqan::CharString> tmpNameStore;
        seqan::NameStoreCache<seqan::StringSet<seqan::CharString> > tmpNameStoreCache(tmpNameStore);
        seqan::BamIOContext<seqan::StringSet<seqan::CharString> > tmpContext(tmpNameStore, tmpNameStoreCache);
        if (readRecord(tmpHeader, tmpContext, *readers[i], seqan::Sam()) != 0)
            throw MasonIOException("Could not load SAM header from temporary file.");
        if (i == 0u)
        {
            header = tmpHeader;
            nameStore = seqan::nameStore(tmpContext);
        }

        active[i] = _loadNext(records[i], i);
        numActive += (active[i] != false);
    }

    // Refresh name store cache.
    refresh(nameStoreCache);
}

// ---------------------------------------------------------------------------
// Function SamJoiner::_loadNext()
// ---------------------------------------------------------------------------

bool SamJoiner::_loadNext(seqan::BamAlignmentRecord & record, unsigned idx)
{
    if (seqan::atEnd(*readers[idx]))
        return false;
    if (readRecord(record, context, *readers[idx], seqan::Sam()) != 0)
    {
        std::cerr << "ERROR: Problem reading temporary data.\n";
        exit(1);
    }
    return true;
}

// ---------------------------------------------------------------------------
// Function SamJoiner::get()
// ---------------------------------------------------------------------------

int SamJoiner::get(seqan::BamAlignmentRecord & record)
{
    unsigned idx = seqan::maxValue<unsigned>();
    for (unsigned i = 0; i < length(records); ++i)
    {
        if (!active[i])
            continue;
        if (idx == seqan::maxValue<unsigned>() || ltBamAlignmentRecord(records[i], records[idx]))
            idx = i;
    }
    if (idx == seqan::maxValue<unsigned>())
        return 1;

    // We use double-buffering and the input parameters as buffers.
    using std::swap;
    active[idx] = _loadNext(record, idx);
    swap(record, records[idx]);
    numActive -= !active[idx];

    return 0;
}

// ---------------------------------------------------------------------------
// Function ContigPicker::pick()
// ---------------------------------------------------------------------------

std::pair<int, int> ContigPicker::pick()
{
    // Pick reference id.
    int rID = 0;
    if (lengthSums.size() > 1u)
    {
        seqan::Pdf<seqan::Uniform<__int64> > pdf(0, lengthSums.back() - 1);
        __int64 x = pickRandomNumber(rng, pdf);
        for (unsigned i = 0; i < lengthSums.size(); ++i)
        {
            if (x >= lengthSums[i])
                rID = i + 1;
            if (x < lengthSums[i])
                break;
        }
    }

    // Pick haplotype id.
    int hID = pickRandomNumber(rng, seqan::Pdf<seqan::Uniform<int> >(0, numHaplotypes - 1));

    return std::make_pair(rID, hID);
}
