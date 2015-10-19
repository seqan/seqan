// ==========================================================================
// Author: Benjamin Menkuec <benjamin@menkuec.de>
// ==========================================================================


#ifndef HELPERFUNCTIONS_H_
#define HELPERFUNCTIONS_H_

#include <string>
#include <future>

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include "read.h"


// seqan->std interface functions

std::string prefix(const std::string& str, unsigned int len) noexcept
{
    return str.substr(0, len);
}

inline std::string seqanToStd(const seqan::Dna5QString& rhs) noexcept
{
    std::string ret;
    ret.resize(length(rhs));
    char c;
    std::transform(begin(rhs), end(rhs), ret.begin(), [&c](const auto& element){
        seqan::assign(c, element);
        return c;
    });
    return ret;
}

void append(std::string& str1, const std::string& str2)
{
    str1 += str2;
}

unsigned int length(const std::string& str) noexcept
{
    return str.size();
}

template <typename T>
unsigned int length(const std::vector<T>& vec) noexcept
{
    return vec.size();
}

template <typename T>
void resize(std::vector<T>& vec, unsigned int len)
{
    vec.resize(len);
}


void insert(std::string& dest, unsigned int k, const std::string& token)
{
    dest.insert(k, token);
}

//

class SemaphoreTest
{
private:
    std::condition_variable condition_;
    std::atomic<unsigned int> count_;
    std::atomic_bool allowAll_;
    std::mutex mutex_;
public:
    SemaphoreTest()
    {
        count_.store(0);
        allowAll_ = false;
    }

    void notify()
    {
        ++count_;
        condition_.notify_one();
    }

    void allowAll()
    {
        std::unique_lock<std::mutex> lock(mutex_);
        allowAll_ = true;
        condition_.notify_all();
    }

    void wait()
    {
        while (!allowAll_)
        {
            auto currentVal = count_.load();
            if (currentVal != 0)
            {
                unsigned int expected = currentVal;
                if (count_.compare_exchange_strong(expected, currentVal -1))
                    return;
            }
            else
            {
                std::unique_lock<std::mutex> lock(mutex_);
                if(!allowAll_)
                    condition_.wait(lock);
            }
        }
    }
};

// always use the forward read for barcode detection
template <template <typename> class TRead, typename TSeq, typename = std::enable_if_t<std::is_same<TRead<TSeq>, Read<TSeq>>::value || std::is_same<TRead<TSeq>, ReadPairedEnd<TSeq>>::value>>
std::string getPrefix(const TRead<TSeq>& read, unsigned len) noexcept
{
    return static_cast<const std::string>(prefix(read.seq, len));
}

template <template <typename> class TRead, typename TSeq, typename = std::enable_if_t<std::is_same<TRead<TSeq>, ReadMultiplex<TSeq>>::value || std::is_same<TRead<TSeq>, ReadMultiplexPairedEnd<TSeq>>::value>>
std::string getPrefix(const TRead<TSeq>& read, unsigned len, bool = false) noexcept
{
    (void)len;
    return seqanToStd(read.demultiplex);
}

template <typename T>
std::vector<T> operator+(const std::vector<T>& a, const std::vector<T>& b) noexcept
{
    assert(a.size() == b.size());

    std::vector<T> result;
    result.reserve(a.size());

    std::transform(a.begin(), a.end(), b.begin(),
        std::back_inserter(result), std::plus<T>());
    return result;
}


template <typename TDest, typename TSource>
void insertAfterFirstToken(TDest& dest, TSource&& source)
{
    unsigned const len = length(dest);
    for (unsigned k = 0;k < len;++k)
        if (dest[k] == ' ' || k == (len - 1))
        {
            insert(dest, k, std::forward<TSource>(source));
            break;
        }
}

template <class F, class... Ts>
void for_each_argument(F f, Ts&&... a) {
    // destructor of temps blocks until all threads are finished
    std::make_tuple(f(std::forward<Ts>(a))...); 
}

template<typename Trem, typename TremVal, typename TRead>
auto _eraseSeqs(const Trem& rem, const TremVal remVal, std::vector<TRead>& reads) noexcept
{
    const auto oldSize = reads.size();
    auto it = rem.cbegin();
    reads.erase(std::remove_if(reads.begin(), reads.end(),
        [&rem, &it, remVal](const auto& read) {(void)read;return *(it++) == remVal;}), reads.end());
    return oldSize - reads.size();
}

template<typename Trem, typename TremVal, typename... TContainer>
auto _eraseSeqs(const Trem& rem, const TremVal remVal, TContainer&&... container) noexcept
{
    const auto numRemoveElements = std::count(begin(rem), end(rem), remVal);
    auto eraseElements = [&rem, numRemoveElements, remVal](auto& seq)  // erase Elements using the remove erase idiom
    {
        const auto beginAddr = &*seqan::begin(seq);
        std::remove_if(seqan::begin(seq), seqan::end(seq),
            [&rem, &beginAddr, remVal](const auto& element) {
            return rem[&element - beginAddr] == remVal;});
        resize(seq, length(seq) - numRemoveElements);
        return 0;
    };
    for_each_argument(eraseElements, std::forward<TContainer>(container)...);
    return numRemoveElements;
}


#endif // HELPERFUNCTIONS_H_
