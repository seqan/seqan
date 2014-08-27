// ==========================================================================
//                                   ANISE
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

#ifndef SANDBOX_HERBARIUM_APPS_ANISE_REP_SEP_CLUSTER_LINKING_H_
#define SANDBOX_HERBARIUM_APPS_ANISE_REP_SEP_CLUSTER_LINKING_H_

#include <map>
#include <vector>
#include <set>

#include <seqan/basic.h>

#include "asm/frag_store.h"

namespace rep_sep {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class ConflictStore
// ----------------------------------------------------------------------------

// Manage a list of unordered conflict pairs between classes.

class ConflictStore
{
public:
    typedef std::pair<unsigned, unsigned> TConflict;
    typedef std::set<TConflict>           TConflictSet;
    typedef TConflictSet::const_iterator  TConflictIter;
    typedef std::map<unsigned, std::set<unsigned> > TConflictMap;

    // Stores conflicting classes, in each pair (i, j) the following holds: "i < j" since we need unordered pairs.
    TConflictSet classConflicts;
    // Stores conflicting sets for each set; required for recording splits.
    TConflictMap conflictMap;

    // Begin iterator to class conflict pair set.
    TConflictIter begin() const { return classConflicts.begin(); }
    TConflictIter end() const { return classConflicts.end(); }

    // Reset ConflictStore.
    void clear()
    {
        classConflicts.clear();
        conflictMap.clear();
    }

    // Record a class split.  The new class inherits all conflicts of the old class.
    void recordClassSplit(unsigned newClass, unsigned oldClass)
    {
        std::set<unsigned> const & oldSet = conflictMap[oldClass];
        for (std::set<unsigned>::const_iterator it = oldSet.begin(); it != oldSet.end(); ++it)
            addConflict(newClass, *it);
    }

    // Insert conflict pair into conflicts.
    void addConflict(unsigned i, unsigned j)
    {
        if (i != j)
        {
            if (i > j)
                std::swap(i, j);
            conflictMap[i].insert(j);
            conflictMap[j].insert(i);
            classConflicts.insert(std::make_pair(i, j));
        }
    }

    // Remove conflict pair into conflicts.
    void removeConflict(unsigned i, unsigned j)
    {
        if (i > j)
            std::swap(i, j);
        SEQAN_ASSERT(classConflicts.count(std::make_pair(i, j)));
        conflictMap[i].erase(j);
        conflictMap[j].erase(i);
    }

    // Return true if there is a conflict between i and j.
    bool inConflict(unsigned i, unsigned j) const
    {
        if (i > j)
            std::swap(i, j);
        return classConflicts.count(std::make_pair(i, j));
    }

    // Remove a class.
    void removeClass(unsigned i)
    {
        std::set<unsigned> const & others = conflictMap[i];
        for (std::set<unsigned>::const_iterator it = others.begin(); it != others.end(); ++it)
        {
            removeConflict(i, *it);
            conflictMap[*it].erase(i);
        }
        conflictMap.erase(i);
    }

    void print(std::ostream & out) const
    {
        out << "Conflicts\nSet:";
        for (auto pair : classConflicts)
            out << " (" << pair.first << ", " << pair.second << ")";
        out << "\nMap\n";
        for (auto const & pair : conflictMap)
        {
            out << pair.first << " -> ";
            for (auto x : pair.second)
                out << " " << x;
            out << "\n";
        }
    }

    void swap(ConflictStore & other)
    {
        using std::swap;
        swap(classConflicts, other.classConflicts);
        swap(conflictMap, other.conflictMap);
    }
};

// ----------------------------------------------------------------------------
// Class ClusterGlobalizer
// ----------------------------------------------------------------------------

class ClusterGlobalizer
{
public:
    // Global partition.
    std::vector<std::set<unsigned> > globalPartition;
    // Mapping from read id to class id, built in stop().
    std::map<unsigned, unsigned> globalMap;
    // Mapping read ids to their mate or INVALID_ID.
    std::vector<unsigned> mateIDs;
    // The value for INVALID_ID.
    static const unsigned INVALID_ID = (unsigned)-1;
    // The conflict store.
    ConflictStore conflictStore;
    // The minimal overlap from local to global class to consider it for joining.
    unsigned minOverlap;
    // Whether or not logging is enabled.
    bool logging;

    explicit ClusterGlobalizer(unsigned minOverlap = 2, bool logging = false) : minOverlap(minOverlap), logging(logging)
    {}

    void reset();

    void init(TFragmentStore const & store);

    // TODO(holtgrew): The following can go away since we do not store clusters any more.
    // Initialize globalizer with global partition and a set of class conflict.
    void init(std::vector<std::set<unsigned> > const & partition,
              std::set<std::pair<unsigned, unsigned> > const & conflicts);

    // Perform next step of the cluster globalization.
    //
    // k -- Number of local clusters.
    // localMap -- Local cluster map, assigning reads to values.
    void run(unsigned k, std::map<unsigned, unsigned> const & localMap);

    // Removes multi-class assignments.
    void stop();

    // Print global partition.
    void printGlobalPartition()
    { printPartition(globalPartition, "GLOBAL PARTITION"); }

    size_t numClusters() const { return globalPartition.size(); }
    std::vector<std::set<unsigned>> const & partition() { return globalPartition; }

private:

    // Remove empty classes.  Note that this invalidates all maps into globalPartition and conflicts.
    void purgeEmptyClasses();

    // Build partitions from local map.
    void buildPartition(std::vector<std::set<unsigned> > & partition,
                         unsigned k, std::map<unsigned, unsigned> const & localMap);

    // Removal of multi-class assignments, Algorithm 5 in (Kuchenbecker, 2011).
    void removeMultiClassAssignments(std::vector<std::set<unsigned> > const & partition);

    // Refinement of global partition given a local partition, Algorithm 6 in (Kuchenbecker, 2011).
    void refineGlobalPartition(std::vector<std::set<unsigned> > const & partition);

    // Print the given partition.
    void printPartition(std::vector<std::set<unsigned> > const & partition, char const * title) const;

    // Perform the global class assignment, Algorithm 7 in (Kuchenbecker, 2011).
    void performGlobalClassAssignment(std::vector<std::set<unsigned> > const & partition);

    // Merge global classes; modifies gp according to Algorithm 8 in (Kuchenbecker, 2011).
    void mergeGlobalClasses(std::set<unsigned> & gp);
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

}  // namespace seqan

#endif  // #ifndef SANDBOX_HERBARIUM_APPS_ANISE_REP_SEP_CLUSTER_LINKING_H_
