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

//TODO: clean up this object
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

public:
    // Constructor saves a base directory path used for all outputs.
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

    inline seqan::CharString getBaseFilename(void) const
    {
        return prefix(basePath, length(basePath) - length(extension));
    }

    //Adds a new output streams to the collection of streams.
    void addStream(seqan::CharString fileName, int id, bool useDefault)
    {
        //Prepend basePath and append file extension to the filename.
        seqan::CharString path = getBaseFilename();
        if (fileName != "")
            seqan::append(path, "_");
        if (useDefault)
            append(path, "_result");

        seqan::append(path, fileName);
        seqan::append(path, extension);
        char* file = seqan::toCString(path);
        fileStreams[id].first = std::make_unique<seqan::SeqFileOut>(file);
    }

    void addStream(seqan::CharString fileName, int id)
    {
        addStream(fileName, id, false);
    }
    //Adds two new output streams to the collection of streams.
    void addStreams(seqan::CharString fileName1, seqan::CharString fileName2, int id, bool useDefault)
    {
        //Prepend basePath and append file extension to the filename.
        seqan::CharString path1 = prefix(basePath, length(basePath) - length(extension));
        seqan::CharString path2 = prefix(basePath, length(basePath) - length(extension));
        if (fileName1 != "")
            seqan::append(path1, "_");
        if (fileName2 != "")
            seqan::append(path2, "_");

        if (useDefault)
        {
            append(path1, "_result");
            append(path2, "_result");
        }

        seqan::append(path1, fileName1);
        seqan::append(path1, extension);
        seqan::append(path2, fileName2);
        seqan::append(path2, extension);
        char* file1 = seqan::toCString(path1);
        char* file2 = seqan::toCString(path2);
        fileStreams[id] = std::make_pair(std::make_unique<seqan::SeqFileOut>(file1), std::make_unique<seqan::SeqFileOut>(file2));
    }

    void addStreams(seqan::CharString fileName1, seqan::CharString fileName2, int id)
    {
        addStreams(fileName1, fileName2, id, false);
    }

    //This method takes a String of integers and checks if these integers are
    //already associated with a stream. If not, a new stream is added and the opened
    //file is named according to the list of names. One or two files are created.
    template <typename TNames>
    void updateStreams(TNames& names, bool pair)
    {

        for (unsigned i = 0; i < length(names) + 1; ++i)
        {
            const unsigned streamIndex = i;
            // If no stream for this id exists, create one.
            if (fileStreams.find(streamIndex) == fileStreams.end())
            {
                // If the index is 0 (unidentified) create special stream.
                // Otherwise use index to get appropriate name for output file.
                seqan::CharString file;
                if (streamIndex > 0)
                {
                    file = names[streamIndex - 1];
                }
                else
                {
                    file = seqan::CharString("unidentified");
                }
                // Add file extension to stream and create it.
                if (pair)
                {
                    // Create a new subfolder at basePath/[barcodeID, unidentified].
                    seqan::CharString folderPath(basePath);
                    seqan::append(folderPath, file);
                    // Turn file target from [barcodeID,unidentified]
                    // to subfolder [barcodeID, unidentified]/[barcodeID, unidentified]
                    seqan::CharString filePath(file);
                    seqan::append(file, "/");
                    seqan::append(file, filePath);
                    // To differentiate between the paired reads, add index to the file name.
                    seqan::CharString file2 = file;
                    seqan::append(file, seqan::CharString("_1"));
                    seqan::append(file2, seqan::CharString("_2"));
                    this->addStreams(file, file2, streamIndex);
                }
                else
                {
                    //std::cerr << file << " " << streamIndex << std::endl;
                    this->addStream(file, streamIndex);
                }
            }
        }
    }

    template <template<typename> class TRead, typename TSeq, typename TNames>
    void writeSeqs(std::vector<TRead<TSeq>>&& reads, TNames& names)
    {
        updateStreams(names, false);
        for(auto& read : reads)
        {
            const unsigned streamIndex = read.demuxResult;
            writeRecord(fileStreams[streamIndex], std::move(read));
        }
    }

    ~OutputStreams(){}

};


// writes read sets to hd and collects statistics, blocks if all slots are full
template<template<typename> class TRead, typename TSeq, typename _TWriteItem, typename TProgramParams, typename TOutputStreams>
struct ReadWriter
{
    using TWriteItem = _TWriteItem;

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
            while (_run)
            {
                bool nothingToDo = true;
                for (auto& readSet : _tlsReadSets)
                {
                    if (readSet.load() != nullptr)
                    {
                        currentWriteItem = readSet.load();
                        readSet.store(nullptr); // make the slot free again
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
                        return;
                    }
                }
            }
            //std::cout << std::this_thread::get_id() << "-3" << std::endl;
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
};