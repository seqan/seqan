// ==========================================================================
// Author: Benjamin Menkuec <benjamin@menkuec.de>
// ==========================================================================


#pragma once

#include <string>
#include <future>

#include "semaphore.h"

// reads read sets from hd and puts them into slots, waits if no free slots are available
template<template<typename> class TRead, typename TSeq, typename TProgramParams, typename TInputFileStreams, bool useSemaphore>
struct ReadReader
{
public:
    using TReadSet = std::vector<TRead<TSeq>>;
private:
    const TProgramParams& _programParams;
    TInputFileStreams& _inputFileStreams;
    std::vector<std::atomic<TReadSet*>> _tlsReadSets;
    std::thread _thread;
    std::atomic_bool _eof;
    unsigned int _sleepMS;
    unsigned int _numReads;
    LightweightSemaphore slotEmptySemaphore;
    LightweightSemaphore readAvailableSemaphore;

public:
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
        /*
        - check if empty slot is available, if yes, read data into it
            - set eof=true if no more data available or max_number_of_reads reached
        - return from thread if eof == true
        - go to sleep if no free slot is available
        */
        _thread = std::thread([this]()
        {
            std::unique_ptr <TReadSet> currentReadSet;
            bool noEmptySlot = true;
            while (true)
            {
                noEmptySlot = true;
                for (auto& readSet : _tlsReadSets) 
                {
                    if (readSet.load() == nullptr)
                    {
                        noEmptySlot = false;
                        currentReadSet = std::make_unique<TReadSet>(_programParams.records);
                        try {
                            readReads(*currentReadSet, _programParams.records, _inputFileStreams);
                        }
                        catch (std::exception& e){
                            std::cout << "exception while reading :" << e.what() << " after read " << _numReads << std::endl;
                            throw(e);
                        }
                        loadMultiplex(*currentReadSet, _programParams.records, _inputFileStreams.fileStreamMultiplex);
                        _numReads += currentReadSet->size();
                        if (currentReadSet->empty() || _numReads >= _programParams.firstReads)    // no more reads available or maximum read number reached -> dont do further reads
                            _eof = true;
                        readSet.store(currentReadSet.release());

                        if (useSemaphore && _eof)  // wakeup all potentially waiting threads so that they can be joined
                            readAvailableSemaphore.signal(_programParams.num_threads);
                        else if (useSemaphore)
                            readAvailableSemaphore.signal();
                    }
                    if (_eof)
                        return;
                }
                if (noEmptySlot)  // no empty slow was found, wait a bit
                {
                    if (useSemaphore)
                        slotEmptySemaphore.wait();
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
    /*
    do the following steps sequentially
    - check if any slot contains data, if yes take data out and return true
    - check if eof is reached, if yes return false
    - go to sleep until data is available
    */
    bool getReads(std::unique_ptr<TReadSet>& reads) noexcept
    {
        TReadSet* temp = nullptr;
        while (true)
        {
            bool eof = _eof;
            for (auto& readSet : _tlsReadSets)
            {
                if ((temp = readSet.load()) != nullptr)
                {
                    if (readSet.compare_exchange_strong(temp, nullptr))
                    {
                        reads.reset(temp);
                        if (useSemaphore)
                            slotEmptySemaphore.signal();
                        return true;
                    }
                }
            }
            //std::cout << std::this_thread::get_id() << "-2" << std::endl;
            if (eof) // only return if _eof == true AND all the slots are empty -> therefore read eof BEFORE checking the slots
                return false;
            else if (useSemaphore)
                readAvailableSemaphore.wait();
            else
                std::this_thread::sleep_for(std::chrono::milliseconds(_sleepMS));
        }
        return false;
    };
};

