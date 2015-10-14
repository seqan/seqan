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