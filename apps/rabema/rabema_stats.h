// ==========================================================================
//                      RABEMA Read Alignment Benchmark
// ==========================================================================
// Copyright (C) 2010-1012 Manuel Holtgrewe, FU Berlin
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License along
// with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// ==========================================================================
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_APPS_RABEMA_RABEMA_STATS_H_
#define SEQAN_APPS_RABEMA_RABEMA_STATS_H_

#include <seqan/sequence.h>

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

namespace seqan {
struct Tsv_;
typedef Tag<Tsv_> Tsv;
}

// ----------------------------------------------------------------------------
// Class RabemaStats
// ----------------------------------------------------------------------------

// Container for RABEMA results.

struct RabemaStats
{
    // Number of intervals that were to find.
    uint64_t intervalsToFind;
    // Number of intervals that were actually found.
    uint64_t intervalsFound;
    // SAM records that did not correspond to alignments below the configured maximal error rate.
    uint64_t invalidAlignments;

    // Total number of reads.
    uint64_t totalReads;
    // Number of reads with alignments below maximal error rate.
    uint64_t mappedReads;
    // Total number of reads with GSI records equals number of normalized intervals to find.
    uint64_t readsInGsi;
    // Normalized number of found intervals.
    double normalizedIntervals;
    // Number of additional alignments in SAM file with low enough error rate but no GSI record.
    unsigned additionalHits;

    // The following arrays are indexed by the integer value of the error rate.

    // Number of intervals that were to find for each error rate.
    String<unsigned> intervalsToFindForErrorRate;
    // Number of found intervals for each error rate.
    String<unsigned> intervalsFoundForErrorRate;

    // The following values are normalized towards all intervals.

    // Normalized number of intervals to find for each error rate.
    String<double> normalizedIntervalsToFindForErrorRate;
    // Normalized number of intervals found for each error rate.
    String<double> normalizedIntervalsFoundForErrorRate;

    RabemaStats() :
        intervalsToFind(0), intervalsFound(0), invalidAlignments(0), totalReads(0), mappedReads(0),
        readsInGsi(0), normalizedIntervals(0), additionalHits(0)
    {}

    RabemaStats(unsigned maxErrorRate) :
        intervalsToFind(0), intervalsFound(0), invalidAlignments(0), totalReads(0), mappedReads(0),
        readsInGsi(0), normalizedIntervals(0),
        additionalHits(0)
    {
        resize(intervalsToFindForErrorRate, maxErrorRate + 1, 0);
        resize(intervalsFoundForErrorRate, maxErrorRate + 1, 0);
        resize(normalizedIntervalsToFindForErrorRate, maxErrorRate + 1, 0.0);
        resize(normalizedIntervalsFoundForErrorRate, maxErrorRate + 1, 0.0);
    }

};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function updateMaximalErrorRate()
// ----------------------------------------------------------------------------

void updateMaximalErrorRate(RabemaStats & stats, unsigned maxErrorRate)
{
    if (length(stats.intervalsToFindForErrorRate) <= maxErrorRate + 1)
        resize(stats.intervalsToFindForErrorRate, maxErrorRate + 1, 0);
    if (length(stats.intervalsFoundForErrorRate) <= maxErrorRate + 1)
        resize(stats.intervalsFoundForErrorRate, maxErrorRate + 1, 0);
    if (length(stats.normalizedIntervalsToFindForErrorRate) <= maxErrorRate + 1)
        resize(stats.normalizedIntervalsToFindForErrorRate, maxErrorRate + 1, 0.0);
    if (length(stats.normalizedIntervalsFoundForErrorRate) <= maxErrorRate + 1)
        resize(stats.normalizedIntervalsFoundForErrorRate, maxErrorRate + 1, 0.0);
}

// ----------------------------------------------------------------------------
// Function write()                                                       [Raw]
// ----------------------------------------------------------------------------

template <typename TStream>
void write(TStream & stream, RabemaStats const & stats, /*int maxError,*/ seqan::Raw const & /*tag*/)
{
    stream << "Intervals to find:              " << stats.intervalsToFind << '\n'
           << "Intervals found:                " << stats.intervalsFound << '\n';
    if (stats.intervalsToFind > 0u)
        stream << "Intervals found [%]             " << (100.0 * stats.intervalsFound / stats.intervalsToFind) << '\n';
    else
        stream << "Intervals found [%]             " << 0 << '\n';
    stream << "Invalid alignments:             " << stats.invalidAlignments << '\n'
           << "Additional Hits:                " << stats.additionalHits << '\n'
           << '\n'
           << "Number of reads:                " << stats.totalReads << '\n'
           << "Number of reads with intervals: " << stats.readsInGsi << '\n'
           << '\n'
           << "Mapped reads:                   " << stats.mappedReads << '\n'
           << "Mapped reads [% of total]:      " << 100.0 * stats.mappedReads / stats.totalReads << '\n'
           << "Mapped reads [% of mappable]:   " << 100.0 * stats.mappedReads / stats.readsInGsi  << '\n'
           << '\n'
           << "Normalized intervals found:     " << stats.normalizedIntervals << '\n';
    if (stats.readsInGsi > 0u)
        stream << "Normalized intervals found [%]: " << (100.0 * stats.normalizedIntervals / stats.readsInGsi) << "\n\n";
    else
        stream << "Normalized intervals found [%]: " << 0 << "\n\n";
    stream << "Found Intervals By Error Rate\n"
           << "\n";
    char buffer[1000];
    snprintf(buffer, 1000, "  ERR\t%8s\t%8s\t%8s\t%8s\t%10s\t%10s\n", "#max", "#found", "%found", "norm max", "norm found", "norm found [%]");
    stream << buffer
           << "------------------------------------------------------------------------------------------------------\n";
    for (unsigned i = 0; i < length(stats.intervalsToFindForErrorRate); ++i)
    {
        // if (maxError != -1  && (int)i != maxError)
        //     continue;
        double percFoundIntervals = 100.0 * stats.intervalsFoundForErrorRate[i] / stats.intervalsToFindForErrorRate[i];
        if (stats.intervalsToFindForErrorRate[i] == 0u)
            percFoundIntervals = 0;
        double percFoundNormalizedIntervals = 100.0 * stats.normalizedIntervalsFoundForErrorRate[i] / stats.normalizedIntervalsToFindForErrorRate[i];
        if (stats.normalizedIntervalsToFindForErrorRate[i] == 0)
            percFoundNormalizedIntervals = 0;
        snprintf(buffer, 1000,"%5u\t%8d\t%8d\t%8.2f\t%8.2f\t%10.2f\t%10.2f\n", i, stats.intervalsToFindForErrorRate[i], stats.intervalsFoundForErrorRate[i],
                percFoundIntervals, stats.normalizedIntervalsToFindForErrorRate[i], stats.normalizedIntervalsFoundForErrorRate[i],
                percFoundNormalizedIntervals);
        stream << buffer;
    }
    stream << '\n';
}

// ----------------------------------------------------------------------------
// Function write()                                                       [Tsv]
// ----------------------------------------------------------------------------

template <typename TStream>
int write(TStream & stream, RabemaStats const & stats, int maxError, CharString const & benchmarkCategory,
          bool oracleMode, CharString const & distanceMetric, seqan::Tsv const & /*tag*/)
{
    stream << "##Rabema Results\n"
           << "##\n"
           << "##category\t" << benchmarkCategory << '\n'
           << "##max_distance\t" << maxError<< '\n'
           << "##oracle_mode\t" << (oracleMode ? (char const *) "yes" : (char const *) "no") << "\n"
           << "##distance_metric\t" << distanceMetric << '\n'
           << "##\n"
           << "##intervals_to_find\t" << stats.intervalsToFind << '\n'
           << "##intervals_found\t" << stats.intervalsFound << '\n'
           << "##intervals_found_percent\t" << (100.0 * stats.intervalsFound / stats.intervalsToFind) << '\n'
           << "##invalid_alignments\t" << stats.invalidAlignments << '\n'
           << "##additional_hits\t" << stats.additionalHits << '\n'
           << "##\n"
           << "##number_of_reads\t" << stats.totalReads << '\n'
           << "##number_of_reads_with_intervals\t" << stats.readsInGsi << '\n'
           << "##\n"
           << "##mapped_reads\t" << stats.mappedReads << '\n'
           << "##mapped_reads_percent_of_total\t" << 100.0 * stats.mappedReads / stats.totalReads << '\n'
           << "##mapped_reads_percent_of_mappable\t" << 100.0 * stats.mappedReads / stats.readsInGsi << '\n'
           << "##\n"
           << "##normalized_intervals_found\t" << stats.normalizedIntervals << '\n'
           << "##normalized_intervals_found_percent\t" << (100.0 * stats.normalizedIntervals / stats.readsInGsi) << '\n'
           << "##\n"
           << "##Found Intervals By Distance\n"
           << "##\n"
           << "#error_rate\tnum_max\tnum_found\tpercent_found\tnorm_max\tnorm_found\tpercent_norm_found\n";
    char buffer[1000];
    for (unsigned i = 0; i < length(stats.intervalsToFindForErrorRate); ++i)
    {
        // if (maxError != -1  && (int)i != maxError)
        //     continue;
        double percFoundIntervals = 100.0 * stats.intervalsFoundForErrorRate[i] / stats.intervalsToFindForErrorRate[i];
        if (stats.intervalsToFindForErrorRate[i] == 0u)
            percFoundIntervals = 0;
        double percFoundNormalizedIntervals = 100.0 * stats.normalizedIntervalsFoundForErrorRate[i] / stats.normalizedIntervalsToFindForErrorRate[i];
        if (stats.normalizedIntervalsToFindForErrorRate[i] == 0)
            percFoundNormalizedIntervals = 0;
        snprintf(buffer, 1000, "%u\t%d\t%d\t%.2f\t%.2f\t%1.2f\t%.2f\n", i, stats.intervalsToFindForErrorRate[i], stats.intervalsFoundForErrorRate[i],
                percFoundIntervals, stats.normalizedIntervalsToFindForErrorRate[i], stats.normalizedIntervalsFoundForErrorRate[i],
                percFoundNormalizedIntervals);
        stream << buffer;
    }

    return 0;
}

#endif  // #ifndef SEQAN_APPS_RABEMA_RABEMA_STATS_H_
