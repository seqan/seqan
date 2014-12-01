// ==========================================================================
//                                   ANISE
// ==========================================================================
// Copyright (c) 2006-2014, Knut Reinert, FU Berlin
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

#include "cluster_linking.h"

#include <seqan/sequence.h>

#include <iterator>

namespace rep_sep {

namespace {  // anonymous namespace

// ----------------------------------------------------------------------------
// Function intersectionSize()
// ----------------------------------------------------------------------------

// Utility function that returns the size of the intersection of two sets.

// TODO(holtgrew): Could be rewritten so no allocation is performed.

unsigned intersectionSize(std::set<unsigned> const & lhs,
                          std::set<unsigned> const & rhs)
{
    std::set<unsigned> buf;
    std::set_intersection(lhs.begin(), lhs.end(), rhs.begin(), rhs.end(),
                          std::inserter(buf, buf.end()));
    return buf.size();
}

}  // anonymous namespace

// ----------------------------------------------------------------------------
// Class ClusterGlobalizer
// ----------------------------------------------------------------------------

void ClusterGlobalizer::reset()
{
    globalPartition.clear();
    globalMap.clear();
    mateIDs.clear();
    conflictStore.clear();
}

void ClusterGlobalizer::init(TFragmentStore const & store)
{
    seqan::String<unsigned> mateIdxs;  // referring
    calculateMateIndices(mateIdxs, store);

    mateIDs.clear();
    mateIDs.resize(length(mateIdxs), (unsigned)-1);
    for (unsigned readID = 0; readID < mateIDs.size(); ++readID)
        mateIDs[readID] = store.alignedReadStore[mateIdxs[readID]].readId;
}

void ClusterGlobalizer::init(std::vector<std::set<unsigned> > const & partition,
                             std::set<std::pair<unsigned, unsigned> > const & conflicts)
{
    globalPartition = partition;
    for (std::set<std::pair<unsigned, unsigned> >::const_iterator it = conflicts.begin();
         it != conflicts.end(); ++it)
        conflictStore.addConflict(it->first, it->second);
}

void ClusterGlobalizer::run(unsigned k, std::map<unsigned, unsigned> const & localMap)
{
    if (logging)
    {
        std::cerr << "Incorporating new local mapping (k == " << k << ")\n";
        for (auto const & pair : localMap)
            std::cerr << pair.first << "\t->\t" << pair.second << "\n";
    }

    std::vector<std::set<unsigned> > partition;
    buildPartition(partition, k, localMap);
    if (logging)
    {
        printPartition(partition, "partition");
        printPartition(globalPartition, "global partition");
        conflictStore.print(std::cerr);
        std::cerr << "_removeMultiClassAssignments(partition)\n";
    }
    removeMultiClassAssignments(partition);
    if (logging)
    {
        printPartition(partition, "partition");
        std::cerr << "globalPartition.size() == " << globalPartition.size() << "\n";
        printPartition(globalPartition, "global partition");
        conflictStore.print(std::cerr);
        std::cerr << "_refineGlobalPartition(partition)\n";
    }
    refineGlobalPartition(partition);
    if (logging)
    {
        printPartition(partition, "partition");
        std::cerr << "globalPartition.size() == " << globalPartition.size() << "\n";
        printPartition(globalPartition, "global partition");
        conflictStore.print(std::cerr);
        std::cerr << "_performGlobalClassAssignment(partition)\n";
    }
    performGlobalClassAssignment(partition);
    if (logging)
    {
        printPartition(partition, "partition");
        std::cerr << "globalPartition.size() == " << globalPartition.size() << "\n";
        printPartition(globalPartition, "global partition");
        conflictStore.print(std::cerr);
    }
}

// Removes multi-class assignments.
void ClusterGlobalizer::stop()
{
    // std::cerr << "Finishing globalization...\n";
    // printPartition(globalPartition, "global partition");
    // conflictStore.print(std::cerr);

    // Build multi-mapping from read id to all containing clusters.  We'll use this below for building the
    // multi-class sets.
    std::map<unsigned, std::set<unsigned> > assigned;
    typedef std::set<unsigned>::const_iterator TConstIterator;
    for (unsigned i = 0; i < globalPartition.size(); ++i)
        for (TConstIterator it = globalPartition[i].begin(); it != globalPartition[i].end(); ++it)
            assigned[*it].insert(i);

    // Build inverse map from "set of containing clusters" to "reads that are assigned to this set of clusters."
    // TODO(holtgrew): Can this be done faster than using such inverse maps?
    std::map<std::set<unsigned>, std::set<unsigned> > inverseMap;
    typedef std::map<unsigned, std::set<unsigned> >::const_iterator TMapIterator;
    for (TMapIterator it = assigned.begin(); it != assigned.end(); ++it)
        inverseMap[it->second].insert(it->first);

    // Create new classes and remove from old global ones.
    typedef std::map<std::set<unsigned>, std::set<unsigned> >::const_iterator TIt;
    for (TIt it = inverseMap.begin(); it != inverseMap.end(); ++it)
    {
        globalPartition.push_back(it->second);
        for (TConstIterator it2 = it->first.begin(); it2 != it->first.end(); ++it2)
            for (TConstIterator it3 = it->second.begin(); it3 != it->second.end(); ++it3)
                globalPartition[*it2].erase(*it3);
        // Inherit conflicts from containing class.
        unsigned newID = globalPartition.size() - 1;
        for (TConstIterator it2 = it->first.begin(); it2 != it->first.end(); ++it2)
        {
            for (std::set<unsigned>::const_iterator it3 = conflictStore.conflictMap[*it2].begin();
                 it3 != conflictStore.conflictMap[*it2].end(); ++it3)
                conflictStore.addConflict(*it3, newID);
        }
    }

    // Remove empty global classes.
    purgeEmptyClasses();

    // std::cerr << "Done with globalization...\n";
    // printPartition(globalPartition, "globalPartition");
    // conflictStore.print(std::cerr);

    // Build global map.
    for (unsigned i = 0; i < globalPartition.size(); ++i)
        for (TConstIterator it = globalPartition[i].begin(); it != globalPartition[i].end(); ++it)
            globalMap[*it] = i;

    if (logging)
    {
        std::cerr << "AFTER STOPPING\n";
        std::cerr << "globalPartition.size() == " << globalPartition.size() << "\n";
        printPartition(globalPartition, "global partition");
        conflictStore.print(std::cerr);
    }
}

void ClusterGlobalizer::purgeEmptyClasses()
{
    std::map<unsigned, unsigned> old2new;  // built during purging, used for updating conflicts

    // Remove empty classes.
    {
        std::vector<std::set<unsigned> > tmp;
        tmp.resize(globalPartition.size());
        unsigned pos = 0;
        for (unsigned i = 0; i < globalPartition.size(); ++i)
            if (!globalPartition[i].empty())
            {
                old2new[i] = pos;
                tmp[pos++].swap(globalPartition[i]);
            }
        tmp.resize(pos);
        tmp.swap(globalPartition);
    }

    // Update ConflictStore.
    {
        ConflictStore tmp;
        for (std::set<ConflictStore::TConflict>::iterator it = conflictStore.classConflicts.begin();
             it != conflictStore.classConflicts.end(); ++it)
            if (old2new.count(it->first) && old2new.count(it->second))
            {
                std::pair<unsigned, unsigned> p(old2new[it->first], old2new[it->second]);
                if (p.first > p.second)
                    std::swap(p.first, p.second);
                tmp.classConflicts.insert(p);
            }
        for (std::map<unsigned, std::set<unsigned> >::iterator it = conflictStore.conflictMap.begin();
             it != conflictStore.conflictMap.end(); ++it)
            if (old2new.count(it->first))
            {
                unsigned newID = old2new[it->first];
                for (std::set<unsigned>::const_iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2)
                    if (old2new.count(*it2))
                        tmp.conflictMap[newID].insert(old2new[*it2]);
            }
        conflictStore.swap(tmp);
    }
}

void ClusterGlobalizer::buildPartition(std::vector<std::set<unsigned> > & partition,
                                        unsigned k, std::map<unsigned, unsigned> const & localMap)
{
    partition.clear();
    partition.resize(k);
    
    typedef std::map<unsigned, unsigned>::const_iterator TIter;
    for (TIter it = localMap.begin(); it != localMap.end(); ++it)
    {
        // std::cerr << it->first << "\t" << it->second << "\n";
        SEQAN_ASSERT_LEQ(it->second, k);
        // std::cerr << it->first << " -> " << it->second << "\n" << std::flush;
        partition.at(it->second).insert(it->first);
        // std::cerr << mateIDs.at(it->first) << " -> " << it->second << "\n" << std::flush;
        if (it->first < mateIDs.size() && mateIDs.at(it->first) != INVALID_ID)
            partition.at(it->second).insert(mateIDs.at(it->first));
    }
}

void ClusterGlobalizer::removeMultiClassAssignments(std::vector<std::set<unsigned> > const & partition)
{
    // u, set of uniquely assigned partition entries
    std::set<unsigned> uniquelyAssigned;
    // G_l, linkedGlobal[i] are the linked global clusters for read i.
    std::map<unsigned, std::set<unsigned> > linkedGlobal;

    typedef std::vector<std::set<unsigned> >::const_iterator TIterator;
    for (TIterator it = partition.begin(); it != partition.end(); ++it)
    {
        std::set<unsigned> const & p = *it;
        uniquelyAssigned.clear();

        // Build set of uniquelyAssigned partition entries
        for (std::set<unsigned>::const_iterator lIt = p.begin(); lIt != p.end(); ++lIt)
        {
            unsigned l = *lIt;

            // Build G_l.
            unsigned g = 0;  // current partition id.
            for (TIterator it2 = globalPartition.begin(); it2 != globalPartition.end(); ++it2, ++g)
                if (it2->count(l))
                    linkedGlobal[l].insert(g);
            // std::cerr << "l == " << l << ", adding g == " << g << "\n";
            // If G_l is a singleton set then we add its entry to the uniquely assigned partition entries
            if (linkedGlobal[l].size() == 1u)
                uniquelyAssigned.insert(*linkedGlobal[l].begin());
        }

        // std::cerr << "partition entry #" << (it - partition.begin()) << ", unique: ";
        // for (std::set<unsigned>::const_iterator itU = uniquelyAssigned.begin(); itU != uniquelyAssigned.end(); ++itU)
        //     std::cerr << " " << *itU;
        // std::cerr << "\n";

        // Remove multi-assignments.
        for (std::set<unsigned>::const_iterator lIt = p.begin(); lIt != p.end(); ++lIt)
        {
            unsigned l = *lIt;

            std::set<unsigned> & Gl = linkedGlobal[l];  // G_l in Algorithm 5
            for (std::set<unsigned>::iterator gIt = Gl.begin(); gIt != Gl.end(); ++gIt)
                if (!uniquelyAssigned.count(*gIt))
                {
                    globalPartition[*gIt].erase(l);
                    if (l < mateIDs.size() && mateIDs.at(l) != INVALID_ID)
                        globalPartition[*gIt].erase(mateIDs.at(l));
                }
        }
    }
}

void ClusterGlobalizer::refineGlobalPartition(std::vector<std::set<unsigned> > const & partition)
{
    // G_l, ids in globalPartition
    std::map<unsigned, std::set<unsigned> > linkedGlobal;

    // Build P_g mapping.  The ids of the local partitions that the global partition g/gID shares reads with.
    std::map<unsigned, std::set<unsigned> > part;  // part[g] is P[g] from Algorithm 6
    typedef std::vector<std::set<unsigned> >::const_iterator TConstIterator;
    unsigned pID = 0;  // id of p
    for (TConstIterator it = partition.begin(); it != partition.end(); ++it, ++pID)
    {
        std::set<unsigned> const & p = *it;
        for (std::set<unsigned>::const_iterator lIt = p.begin(); lIt != p.end(); ++lIt)
        {
            unsigned l = *lIt;

            // Build G_l.
            unsigned g = 0;  // current global partition id.
            for (TConstIterator it2 = globalPartition.begin(); it2 != globalPartition.end(); ++it2, ++g)
                if (it2->count(l) && (intersectionSize(*it2, *it) >= minOverlap))
                    linkedGlobal[l].insert(g);

            std::set<unsigned> const & Gl = linkedGlobal[l];  // G_l in Algorithm 6
            for (std::set<unsigned>::const_iterator gIt = Gl.begin(); gIt != Gl.end(); ++gIt)
                part[*gIt].insert(pID);
        }
    }

    /*
      std::cerr << "G_l\n";
      for (std::map<unsigned, std::set<unsigned> >::const_iterator it = linkedGlobal.begin(); it != linkedGlobal.end(); ++it)
      {
      std::cerr << it->first << "\t";
      for (std::set<unsigned>::const_iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2)
      std::cerr << " " << *it2;
      std::cerr << "\n";
      }
      std::cerr << "\n";

      std::cerr << "P_g\n";
      for (std::map<unsigned, std::set<unsigned> >::const_iterator it = part.begin(); it != part.end(); ++it)
      {
      std::cerr << it->first << "\t";
      for (std::set<unsigned>::const_iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2)
      std::cerr << " " << *it2;
      std::cerr << "\n";
      }
      std::cerr << "\n";
    */

    // Refine partitions.
    // typedef std::vector<std::set<unsigned> >::iterator TIterator;
    for (unsigned gID = 0; gID < globalPartition.size(); ++gID)  // for g \in G
    {
        // P_g in Algorithm 6;
        std::set<unsigned> const & Pg = part[gID];
        if (Pg.size() <= 1u)
            continue;  // Skip if uniquely assigned.

        // for p \in P_g
        for (std::set<unsigned>::const_iterator itP = Pg.begin(); itP != Pg.end(); ++itP)
        {
            SEQAN_ASSERT_LT(*itP, partition.size());
            std::set<unsigned> const & p = partition[*itP];  // p in Algorithm 6.
            // Get shortcut to current global partition.
            std::set<unsigned> & g = globalPartition[gID];
            // Ignore if intersection of p and g < \Delta_\min.
            if (intersectionSize(p, g) < minOverlap)
                continue;

            std::set<unsigned> newPartitionEntry;  // g' in Algorithm 6
            for (std::set<unsigned>::const_iterator itL = p.begin(); itL != p.end(); ++itL)  // for l \in p
            {
                unsigned l = *itL;
                if (g.count(l))  // if l \in g
                {
                    g.erase(l);
                    newPartitionEntry.insert(l);
                    if (l < mateIDs.size() && mateIDs.at(l) != INVALID_ID)
                    {
                        g.erase(mateIDs.at(l));
                        newPartitionEntry.insert(mateIDs.at(l));
                    }
                }
            }
            // Add new class g.
            globalPartition.push_back(newPartitionEntry);

            // Register splitting of class.
            conflictStore.recordClassSplit(globalPartition.size() - 1, gID);
        }
    }
}

void ClusterGlobalizer::printPartition(std::vector<std::set<unsigned> > const & partition,
                                        char const * title) const
{
    std::cerr << title << "\n";
    for (unsigned i = 0; i < partition.size(); ++i)
    {
        std::cerr << i << "\t";
        for (std::set<unsigned>::const_iterator it = partition[i].begin(); it != partition[i].end(); ++it)
            std::cerr << " " << *it;
        std::cerr << "\n";
    }
    std::cerr << "\n";
}

void ClusterGlobalizer::performGlobalClassAssignment(std::vector<std::set<unsigned> > const & partition)
{
    // Map from idx in partition to set of idx in globalPartition.  For each entry X in partition, we record all
    // global partition entries Y where a read from X was assigned to an entry Y.  This is used for projecting
    // conflicts by different class membership in partition to the global classes.
    std::map<unsigned, std::set<unsigned> > assignments;

    // for p \in P do
    typedef std::vector<std::set<unsigned> >::const_iterator TConstIterator;
    unsigned pID = 0;  // id of p
    for (TConstIterator it = partition.begin(); it != partition.end(); ++it, ++pID)
    {
        std::set<unsigned> const & p = *it;

        // Build G_p: The list of global classes that share a read with the class p.
        std::set<unsigned> gp;  // G_p, ids in globalPartition
        for (std::set<unsigned>::const_iterator itL = p.begin(); itL != p.end(); ++itL)
        {
            for (unsigned i = 0; i < globalPartition.size(); ++i)
                if (globalPartition[i].count(*itL) && intersectionSize(p, globalPartition[i]) >= minOverlap)
                    gp.insert(i);
        }

        // Create new global class if p is not linked to any global class via shared reads.
        if (gp.empty())
        {
            gp.insert(globalPartition.size());
            globalPartition.resize(globalPartition.size() + 1);
        }

        // Merging of global classes according to conflicts.
        if (logging)
        {
            std::cerr << "mergeGlobalClasses(gp)\n"
                      << "\tgp\t";
            std::copy(gp.begin(), gp.end(), std::ostream_iterator<unsigned>(std::cerr, " "));
            std::cerr << "\n";
        }

        if (gp.size() > 1u)
            mergeGlobalClasses(gp);
        if (logging)
        {
            printPartition(partition, "partition");
            std::cerr << "globalPartition.size() == " << globalPartition.size() << "\n";
            printPartition(globalPartition, "global partition");
            conflictStore.print(std::cerr);
            std::cerr << "_refineGlobalPartition(partition)\n";
        }

        // Assign each read from p to all clusters that p shares a read with.
        for (std::set<unsigned>::const_iterator itL = p.begin(); itL != p.end(); ++itL)
        {
            bool foundAny = false;
            for (TConstIterator itG = globalPartition.begin(); itG != globalPartition.end(); ++itG)
            {
                if (intersectionSize(p, *itG) < minOverlap)
                    continue;
                foundAny = foundAny || (itG->count(*itL));
            }

            // typedef std::vector<std::set<unsigned> >::iterator TIterator;
            if (!foundAny)
                for (std::set<unsigned>::const_iterator itG = gp.begin(); itG != gp.end(); ++itG)
                {
                    globalPartition[*itG].insert(*itL);
                    assignments[pID].insert(*itG);
                    if (*itL < mateIDs.size() && mateIDs.at(*itL) != INVALID_ID)
                        globalPartition[*itG].insert(mateIDs.at(*itL));
                }
        }
    }

    // Project conflicts from partition to globalPartition.
    typedef std::map<unsigned, std::set<unsigned> >::const_iterator TAssIt;
    typedef std::set<unsigned>::const_iterator                      TSetIt;
    for (TAssIt it1 = assignments.begin(); it1 != assignments.end(); ++it1)
        for (TAssIt it2 = std::next(it1); it2 != assignments.end(); ++it2)
            for (TSetIt it3 = it1->second.begin(); it3 != it1->second.end(); ++it3)
                for (TSetIt it4 = it2->second.begin(); it4 != it2->second.end(); ++it4)
                    if (*it3 != *it4)
                        conflictStore.addConflict(*it3, *it4);
}

void ClusterGlobalizer::mergeGlobalClasses(std::set<unsigned> & gp)
{
    // Collect conflicting members of gp.
    std::set<unsigned> gpConflicts;
    for (std::set<unsigned>::const_iterator it1 = gp.begin(); it1 != gp.end(); ++it1)
        for (std::set<unsigned>::const_iterator it2 = it1; it2 != gp.end(); ++it2)
            if (conflictStore.inConflict(*it1, *it2))
            {
                gpConflicts.insert(*it1);
                gpConflicts.insert(*it2);
            }
    // Build copy of gp with non-conflicting members.
    std::set<unsigned> tmp = gp;  // G'_{p}
    for (std::set<unsigned>::const_iterator it = gpConflicts.begin(); it != gpConflicts.end(); ++it)
        tmp.erase(*it);

    // We can only merge anything if tmp is not empty.
    if (tmp.empty())
        return;

    // Create new class g.
    unsigned gID = globalPartition.size();
    globalPartition.resize(globalPartition.size() + 1);
    std::set<unsigned> & g = globalPartition.back();
    for (std::set<unsigned>::const_iterator it = gp.begin(); it != gp.end(); ++it)
        g.insert(globalPartition[*it].begin(), globalPartition[*it].end());

    // TODO(holtgrew): Should this go into ConflictStore?
    // Build list of pairs to remove and to add to conflicts.
    std::vector<std::pair<unsigned, unsigned> > toRemove, toAdd;
    for (ConflictStore::TConflictIter it = conflictStore.begin(); it != conflictStore.end(); ++it)
    {
        bool matchL = tmp.count(it->first), matchR = tmp.count(it->second);
        SEQAN_ASSERT_NOT((matchL && matchR));
        if (matchL || matchR)
        {
            toRemove.push_back(*it);
            std::pair<unsigned, unsigned> p(gID, matchL ? it->second : it->first);
            if (p.first > p.second)
                std::swap(p.first, p.second);
            toAdd.push_back(p);
        }
    }

    // Update conflicts.
    for (std::vector<std::pair<unsigned, unsigned> >::const_iterator it = toRemove.begin(); it != toRemove.end(); ++it)
        conflictStore.removeConflict(it->first, it->second);
    for (std::vector<std::pair<unsigned, unsigned> >::const_iterator it = toAdd.begin(); it != toAdd.end(); ++it)
        conflictStore.addConflict(it->first, it->second);

    // Update gp.
    for (std::set<unsigned>::const_iterator it = tmp.begin(); it != tmp.end(); ++it)
    {
        globalPartition[*it].clear();
        gp.erase(*it);
    }
    gp.insert(gID);
}

}  // namespace seqan
