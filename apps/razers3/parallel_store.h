#ifndef APPS_RAZERS_PARALLEL_STORE_H
#define APPS_RAZERS_PARALLEL_STORE_H

#if defined(PLATFORM_GNU) && defined(_OPENMP)
#include <parallel/algorithm>
#include <parallel/numeric>
#else
#include <algorithm>
#include <numeric>
#endif  // #ifdef PLATFORM_GNU

#include <seqan/parallel.h>

namespace seqan {

#if defined(PLATFORM_GNU) && defined(_OPENMP)

// use MCSTL which is part of the GCC since version 4.3

template <typename TAlign, typename TFunctorLess>
inline void
sortAlignedReads(TAlign & alignStore, TFunctorLess const & less, Parallel)
{
    __gnu_parallel::sort(
        begin(alignStore, Standard()),
        end(alignStore, Standard()),
        less);
}

template <typename TAlign, typename TFunctorLess>
inline void
sortAlignedReads(TAlign const & alignStore, TFunctorLess const & less, Parallel)
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

#else  // #ifdef PLATFORM_GNU

// sequential fallback

template <typename TAlign, typename TFunctorLess>
inline void
sortAlignedReads(TAlign & alignStore, TFunctorLess const & less, Parallel)
{
    std::sort(begin(alignStore, Standard()), end(alignStore, Standard()), less);
}

template <typename TAlign, typename TFunctorLess>
inline void
sortAlignedReads(TAlign const & alignStore, TFunctorLess const & less, Parallel)
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

#endif  // #ifdef PLATFORM_GNU

}

#endif  // ifndef APPS_RAZERS_PARALLEL_STORE_H
