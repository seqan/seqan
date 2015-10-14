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

// reads read sets from hd and puts them into slots, waits if no free slots are available
template<template<typename> class TRead, typename TSeq, typename TProgramParams, typename TInputFileStreams, bool useCV>
struct ReadReader
{
    using TReadSet = std::vector<TRead<TSeq>>;
    ReadReader(const TProgramParams& programParams, TInputFileStreams& inputFileStreams, unsigned int sleepMS)
        :_programParams(programParams), _inputFileStreams(inputFileStreams), _tlsReadSets(_programParams.num_threads + 1),
        _eof(false), _sleepMS(sleepMS), _numReads(0)
    {
        for (auto& readSet : _tlsReadSets)
            readSet.store(nullptr);  // fill initialization does not work for atomics
    }
    ~ReadReader()
    {
        if (_thread.joinable())
            _thread.join();
    }
    void start()
    {
        _thread = std::thread([this]()
        {
            std::unique_ptr <TReadSet> currentReadSet;
            std::mutex waitForFreeSlotDummyMutex;
            while (!_eof || currentReadSet)
            {
                if (!currentReadSet)  // load new reads from hd
                {
                    currentReadSet = std::make_unique<TReadSet>(_programParams.records);

                    readReads(*currentReadSet, _programParams.records, _inputFileStreams);
                    loadMultiplex(*currentReadSet, _programParams.records, _inputFileStreams.fileStreamMultiplex);
                    _numReads += currentReadSet->size();
                    if (currentReadSet->empty() || _numReads >= _programParams.firstReads)    // no more reads available or maximum read number reached -> dont do further reads
                    {
                        _eof = true;
                    }
                }
                for (auto& readSet : _tlsReadSets)  // now insert reads into a slot
                {
                    if (currentReadSet && readSet.load() == nullptr)
                    {
                        readSet.store(currentReadSet.release());
                        if (useCV)
                            readAvailableCV.notify_one();
                        break;
                    }
                }
                if (currentReadSet)  // no empty slow was found, wait a bit
                {
                    if (useCV)
                    {
                        std::unique_lock<std::mutex> lk(waitForFreeSlotDummyMutex);
                        slotEmptyCV.wait(lk);
                    }
                    else
                        std::this_thread::sleep_for(std::chrono::milliseconds(_sleepMS));
                }
            }
        });
    }
    inline bool eof() const noexcept
    {
        return _eof;
    }
    bool idle() noexcept
    {
        if (!_eof)
            return false;
        for (auto& readSet : _tlsReadSets)
            if (readSet.load())
                return false;
        return true;
    }
    bool getReads(TReadSet** reads) noexcept
    {
        TReadSet* temp = nullptr;
        while (true)
        {
            for (auto& readSet : _tlsReadSets)
            {
                if ((temp = readSet.load()) != nullptr)
                {
                    if (readSet.compare_exchange_strong(temp, nullptr))
                    {
                        *reads = temp;
                        if (useCV)
                            slotEmptyCV.notify_one();
                        return true;
                    }
                }
            }
            if (_eof)
                return false;
            //std::cout << std::this_thread::get_id() << "-2" << std::endl;
            if (useCV)
            {
                std::mutex waitForReadDummyMutex;
                std::unique_lock<std::mutex> lk(waitForReadDummyMutex);
                readAvailableCV.wait(lk);
            }
            else
                std::this_thread::sleep_for(std::chrono::milliseconds(_sleepMS));
        }
        return false;
    };
private:
    const TProgramParams& _programParams;
    TInputFileStreams& _inputFileStreams;
    std::vector<std::atomic<TReadSet*>> _tlsReadSets;
    std::thread _thread;
    std::atomic_bool _eof;
    unsigned int _sleepMS;
    unsigned int _numReads;
    std::condition_variable slotEmptyCV;
    std::condition_variable readAvailableCV;
};

