// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
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

#include <iostream>
#include <cstdio>
#include <vector>

#define SEQAN_TEST

#include <seqan/seeds.h>

#include "test_seeds.h"
#include "test_seeds_global_seed_chain.h"
#include "test_seeds_seed_set.h"
#include "test_seeds_banded_align.h"
#include "test_seeds_memory_manager.h"

SEQAN_DEFINE_TEST(test_seed_banded_align) {
	SEQAN_SKIP_TEST;
    Main_BandedAlign();
}


SEQAN_DEFINE_TEST(test_seed_global_seed_chain) {
    Main_GlobalSeedChain();
}


SEQAN_DEFINE_TEST(test_seed_memory_manager) {
	SEQAN_SKIP_TEST;
    Main_MemoryManager();
}


SEQAN_DEFINE_TEST(test_seed_seeds) {
	SEQAN_SKIP_TEST;
    Main_Seeds();
}


SEQAN_DEFINE_TEST(test_seed_seed_set) {
	SEQAN_SKIP_TEST;
    Main_SeedSet();
}

//////////////////////////////////////////////////////////////////////////////

SEQAN_BEGIN_TESTSUITE(test_seed) {
    SEQAN_CALL_TEST(test_seed_banded_align);
    SEQAN_CALL_TEST(test_seed_global_seed_chain);
    SEQAN_CALL_TEST(test_seed_memory_manager);
    SEQAN_CALL_TEST(test_seed_seeds);
    SEQAN_CALL_TEST(test_seed_seed_set);
	SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/seeds/banded_align.h");
	SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/seeds/banded_chain_align.h");
	SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/seeds/banded_chain_align_affine.h");
	SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/seeds/seed_base.h");
	SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/seeds/seed_multi.h");
	SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/seeds/global_seed_chain.h");
	SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/seeds/memoryManager_base.h");
	SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/seeds/memoryManager_int.h");
	SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/seeds/seedSet_base.h");
	SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/seeds/seedSet_score.h");
}
SEQAN_END_TESTSUITE
