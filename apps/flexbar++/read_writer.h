// ==========================================================================
//                             helper_functions.h
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
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
// Author: Benjamin Menkuec <benjamin@menkuec.de>
// ==========================================================================

#pragma once

#include <string>
#include <future>

class OutputStreams
{
    using TSeqStream = std::unique_ptr<seqan::SeqFileOut>;
    using TStreamPair = std::pair<TSeqStream, TSeqStream>;
    std::map<int, TStreamPair> fileStreams;
    const std::string basePath;
    std::string extension;

    template < typename TStream, template<typename> class TRead, typename TSeq,
        typename = std::enable_if_t < std::is_same<TRead<TSeq>, Read<TSeq>>::value || std::is_same<TRead<TSeq>, ReadMultiplex<TSeq>>::value  > >
        inline void writeRecord(TStream& stream, TRead<TSeq>&& read, bool = false)
    {
        seqan::writeRecord(*(stream.first), std::move(read.id), std::move(read.seq));
    }

    template <typename TStream, template<typename> class TRead, typename TSeq,
        typename = std::enable_if_t < std::is_same<TRead<TSeq>, ReadPairedEnd<TSeq>>::value || std::is_same<TRead<TSeq>, ReadMultiplexPairedEnd<TSeq>>::value  > >
        inline void writeRecord(TStream& stream, TRead<TSeq>&& read)
    {
        seqan::writeRecord(*(stream.first), std::move(read.id), std::move(read.seq));
        seqan::writeRecord(*(stream.second), std::move(read.idRev), std::move(read.seqRev));
    }

    //Adds a new output streams to the collection of streams.
    void _addStream(TSeqStream& stream, const std::string fileName, int id, bool useDefault)
    {
        std::string path = getBaseFilename();
        if (fileName != "")
            path += "_";
        if (useDefault)
            path += "_result";

        path += fileName + extension;
        stream = std::make_unique<seqan::SeqFileOut>(path.c_str());
    }


public:
    // The correct file extension is determined from the base path, according to the available
    // file extensions of the SeqFileOut and used for all stored files.
    OutputStreams(const std::string& base, bool /*noQuality*/) : basePath(base)
    {
        std::vector<std::string> tmpExtensions = seqan::SeqFileOut::getFileExtensions();
        for(const auto& tmpExtension : tmpExtensions)
        {
            if (seqan::endsWith(basePath, tmpExtension))
            {
                extension = tmpExtension;
                break;
            }
        }
    }

    inline std::string getBaseFilename(void) const
    {
        return prefix(basePath, length(basePath) - length(extension));
    }

    void addStream(const std::string fileName, const int streamIndex, const bool useDefault)
    {
        _addStream(fileStreams[streamIndex].first, fileName, streamIndex, useDefault);
    }
    
    void addStreams(const std::string fileName1, const std::string fileName2, const int streamIndex, const bool useDefault)
    {
        _addStream(fileStreams[streamIndex].first, fileName1, streamIndex, useDefault);
        _addStream(fileStreams[streamIndex].second, fileName2, streamIndex, useDefault);
    }

    void addStream(const std::string fileName, const int id)
    {
        addStream(fileName, id, false);
    }

    void addStreams(const std::string fileName1, const std::string fileName2, const int id)
    {
        addStreams(fileName1, fileName2, id, false);
    }

    //This method takes a String of integers and checks if these integers are
    //already associated with a stream. If not, a new stream is added and the opened
    //file is named according to the list of names. One or two files are created.
    template <typename TNames>
    void updateStreams(const TNames& names, const bool pair)
    {

        for (unsigned i = 0; i < length(names) + 1; ++i)
        {
            const unsigned streamIndex = i;
            // If no stream for this id exists, create one.
            if (fileStreams.find(streamIndex) == fileStreams.end())
            {
                // If the index is 0 (unidentified) create special stream.
                // Otherwise use index to get appropriate name for output file.
                std::string file;
                if (streamIndex > 0)
                    file = names[streamIndex - 1];
                else
                    file = "unidentified";
                // Add file extension to stream and create it.
                if (pair)
                {
                    std::string file2 = file;
                    file += "_1";
                    file2 += "_2";
                    addStreams(file, file2, streamIndex);
                }
                else
                {
                    addStream(file, streamIndex);
                }
            }
        }
    }

    template <template<typename> class TRead, typename TSeq, typename TNames>
    void writeSeqs(std::vector<TRead<TSeq>>&& reads, const TNames& names)
    {
        updateStreams(names, std::is_same<TRead<TSeq>, ReadPairedEnd<TSeq>>::value || std::is_same<TRead<TSeq>, ReadMultiplexPairedEnd<TSeq>>::value);
        for(auto& read : reads)
        {
            const unsigned streamIndex = read.demuxResult;
            writeRecord(fileStreams[streamIndex], std::move(read));
        }
    }

    ~OutputStreams(){}

};


// writes read sets to hd and collects statistics, blocks if all slots are full
template<template<typename> class TRead, typename TSeq, typename _TWriteItem, typename TProgramParams, typename TOutputStreams, bool useCV>
struct ReadWriter
{
public:
    using TWriteItem = _TWriteItem;

private:
    const TProgramParams& _programParams;
    TOutputStreams& _outputStreams;
    std::vector<std::atomic<TWriteItem*>> _tlsReadSets;
    std::thread _thread;
    std::atomic_bool _run;
    unsigned int _sleepMS;
    std::chrono::time_point<std::chrono::steady_clock> _startTime;
    std::chrono::time_point<std::chrono::steady_clock> _lastScreenUpdate;
    std::tuple_element_t<2, TWriteItem> _stats;
    std::condition_variable slotEmptyCV;
    std::condition_variable readAvailableCV;

public:
    ReadWriter(const TProgramParams& programParams, TOutputStreams& outputStreams, unsigned int sleepMS)
        : _programParams(programParams), _outputStreams(outputStreams), _tlsReadSets(_programParams.num_threads + 1),
        _run(false), _sleepMS(sleepMS), _startTime(std::chrono::steady_clock::now()), _lastScreenUpdate()
    {
        for (auto& readSet : _tlsReadSets)
            readSet.store(nullptr);  // fill initialization does not work for atomics
    }
    ~ReadWriter()
    {
        _run = false;
        if (_thread.joinable())
            _thread.join();
    }
    void start()
    {
        _run = true;
        _thread = std::thread([this]()
        {
            TWriteItem* currentWriteItem;
            std::mutex waitForDataDummyMutex;
            while (_run)
            {
                bool nothingToDo = true;
                for (auto& readSet : _tlsReadSets)
                {
                    if (readSet.load() != nullptr)
                    {
                        currentWriteItem = readSet.load();
                        readSet.store(nullptr); // make the slot free again
                        if (useCV)
                            slotEmptyCV.notify_one();
                        nothingToDo = false;

                        //std::this_thread::sleep_for(std::chrono::milliseconds(1));  // used for debuggin slow hd case

                        const auto t1 = std::chrono::steady_clock::now();
                        _outputStreams.writeSeqs(std::move(*std::get<0>(*currentWriteItem)), std::get<1>(*currentWriteItem));
                        std::get<0>(*currentWriteItem)->clear();
                        delete std::get<0>(*currentWriteItem); // delete written data
                        _stats += std::get<2>(*currentWriteItem);
                        delete currentWriteItem;

                        // terminal output
                        const auto ioTime = std::chrono::duration_cast<std::chrono::duration<float>>(std::chrono::steady_clock::now() - t1).count();
                        _stats.ioTime += ioTime;
                        const auto deltaLastScreenUpdate = std::chrono::duration_cast<std::chrono::duration<float>>(std::chrono::steady_clock::now() - _lastScreenUpdate).count();
                        if (deltaLastScreenUpdate > 1)
                        {
                            const auto deltaTime = std::chrono::duration_cast<std::chrono::duration<float>>(std::chrono::steady_clock::now() - _startTime).count();
                            if (_programParams.showSpeed)
                                std::cout << "\rReads processed: " << _stats.readCount << "   (" << static_cast<int>(_stats.readCount / deltaTime) << " Reads/s)";
                            else
                                std::cout << "\rReads processed: " << _stats.readCount;
                            _lastScreenUpdate = std::chrono::steady_clock::now();
                        }
                    }
                }
                if (nothingToDo)
                {
                    //std::cout << std::this_thread::get_id() << "-4" << std::endl;
                    if (useCV)
                    {
                        std::unique_lock<std::mutex> lk(waitForDataDummyMutex);
                        readAvailableCV.wait(lk);
                    }
                    else
                        std::this_thread::sleep_for(std::chrono::milliseconds(_sleepMS));
                }
            }
        });
    }
    void writeReads(TWriteItem* writeItem)     // blocks until item could be added
    {
        while (true)
        {
            for (auto& reads : _tlsReadSets)
            {
                if (reads.load() == nullptr)
                {
                    TWriteItem* temp = nullptr;
                    if (reads.compare_exchange_strong(temp, writeItem))
                    {
                        if (useCV)
                            readAvailableCV.notify_one();
                        return;
                    }
                }
            }
            //std::cout << std::this_thread::get_id() << "-3" << std::endl;
            if (useCV)
            {
                std::mutex waitForEmptySlotMutex;
                std::unique_lock<std::mutex> lk(waitForEmptySlotMutex);
                slotEmptyCV.wait(lk);
            }
            else
                std::this_thread::sleep_for(std::chrono::milliseconds(_sleepMS));
        }
    }
    void getStats(std::tuple_element_t<2, TWriteItem>& stats)
    {
        _run = false;
        if (_thread.joinable())
            _thread.join();
        stats = _stats;
    }
    bool idle() noexcept
    {
        for (auto& readSet : _tlsReadSets)
            if (readSet.load() != nullptr)
                return false;
        return true;
    }
};