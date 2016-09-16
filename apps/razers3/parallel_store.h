#ifndef APPS_RAZERS_PARALLEL_STORE_H
#define APPS_RAZERS_PARALLEL_STORE_H

#if defined(COMPILER_GCC) && defined(_OPENMP)
#include <parallel/algorithm>
#include <parallel/numeric>
#else
#include <algorithm>
#include <numeric>
#endif  // #ifdef COMPILER_GCC

#include <seqan/parallel.h>

namespace seqan {

#if defined(COMPILER_GCC) && defined(_OPENMP)

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

#else  // #ifdef COMPILER_GCC

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

#endif  // #ifdef COMPILER_GCC

}

#endif  // ifndef APPS_RAZERS_PARALLEL_STORE_H
