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
// A simple dense mapping between two reference sequence id spaces.  This is
// required, for example, when loading both a BAM file and a reference FASTA
// file.  The order in the BAM file might not be the same as in the FASTA
// file.
//
// The function rebuildMapping() allows one to rebuild the mapping from two
// name store caches.
// ==========================================================================

// TODO(holtgrew): Consider for inclusion in SeqAn library.

#ifndef SEQAN_APPS_RABEMA_REF_ID_MAPPING_H_
#define SEQAN_APPS_RABEMA_REF_ID_MAPPING_H_

#include <seqan/store.h>
#include <seqan/sequence.h>

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class RefIdMapping
// ----------------------------------------------------------------------------

class RefIdMapping
{
public:
    String<unsigned> map;

    RefIdMapping() {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function length()
// ----------------------------------------------------------------------------

inline unsigned length(RefIdMapping const & mapping)
{
    return length(mapping.map);
}

// ----------------------------------------------------------------------------
// Function rebuildMapping()
// ----------------------------------------------------------------------------

template <typename TTargetNameStore, typename TTargetNameStoreCache, typename TSourceNameStore>
void rebuildMapping(RefIdMapping & mapping,
                    TTargetNameStore const & /*targetNameStore*/,
                    TTargetNameStoreCache const & targetNameStoreCache,
                    TSourceNameStore const & sourceNameStore)
{
    clear(mapping.map);
    resize(mapping.map, length(sourceNameStore), std::numeric_limits<unsigned>::max());

    for (unsigned i = 0; i < length(sourceNameStore); ++i)
    {
        unsigned idx = 0;
        if (getIdByName(idx, targetNameStoreCache, sourceNameStore[i]))
            mapping.map[i] = idx;
    }
}

#endif  // #ifndef SEQAN_APPS_RABEMA_REF_ID_MAPPING_H_
