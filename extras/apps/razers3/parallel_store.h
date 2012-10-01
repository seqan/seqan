#ifndef APPS_RAZERS_PARALLEL_STORE_H
#define APPS_RAZERS_PARALLEL_STORE_H

#if defined(PLATFORM_GCC) && __GNUC__ >= 4 && __GNUC_MINOR__ >= 3
#include <parallel/algorithm>
#include <parallel/numeric>
#else
#include <algorithm>
#include <numeric>
#endif  // #ifdef PLATFORM_GCC

using namespace seqan;

struct Parallel_;
typedef Tag<Parallel_> Parallel;

#if defined(PLATFORM_GCC) && __GNUC__ >= 4 && __GNUC_MINOR__ >= 3

// use MCSTL which is part of the GCC since version 4.3

template <typename TAlign, typename TFunctorLess>
inline void
sortAlignedReads(TAlign & alignStore, TFunctorLess const & less, Parallel const &)
{
    __gnu_parallel::sort(
        begin(alignStore, Standard()),
        end(alignStore, Standard()),
        less);
}

template <typename TAlign, typename TFunctorLess>
inline void
sortAlignedReads(TAlign const & alignStore, TFunctorLess const & less, Parallel const &)
{
    __gnu_parallel::sort(
        begin(const_cast<TAlign &>(alignStore), Standard()),
        end(const_cast<TAlign &>(alignStore), Standard()),
        less);
}

template <typename TIntString>
inline void
partialSum(TIntString & intString)
{
    __gnu_parallel::partial_sum(
        begin(intString, Standard()),
        end(intString, Standard()),
        begin(intString, Standard()));
}

#else  // #ifdef PLATFORM_GCC

// sequential fallback

template <typename TAlign, typename TFunctorLess>
inline void
sortAlignedReads(TAlign & alignStore, TFunctorLess const & less, Parallel const &)
{
    std::sort(begin(alignStore, Standard()), end(alignStore, Standard()), less);
}

template <typename TAlign, typename TFunctorLess>
inline void
sortAlignedReads(TAlign const & alignStore, TFunctorLess const & less, Parallel const &)
{
    std::sort(begin(alignStore, Standard()), end(alignStore, Standard()), less);
}

template <typename TIntString>
inline void
partialSum(TIntString & intString)
{
    std::partial_sum(
        begin(intString, Standard()),
        end(intString, Standard()),
        begin(intString, Standard()));
}

#endif  // #ifdef PLATFORM_GCC

#endif  // ifndef APPS_RAZERS_PARALLEL_STORE_H
