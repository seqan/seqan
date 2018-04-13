// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
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
// Author: Hannes Hauswedell <hannes.hauswedell@fu-berlin.de>
// ==========================================================================
// Tests for the blast module
// ==========================================================================

#ifndef SEQAN_TESTS_TEST_BLAST_MISC_H_
#define SEQAN_TESTS_TEST_BLAST_MISC_H_

using namespace seqan;

SEQAN_DEFINE_TEST(test_blast_program)
{
    constexpr BlastProgram  n = BlastProgram::BLASTN;
    constexpr BlastProgram  p = BlastProgram::BLASTP;
    constexpr BlastProgram  x = BlastProgram::BLASTX;
    constexpr BlastProgram tn = BlastProgram::TBLASTN;
    constexpr BlastProgram tx = BlastProgram::TBLASTX;

    BlastProgram n2 = BlastProgram::BLASTN;
    BlastProgram p2 = BlastProgram::BLASTP;

    /*** TEST frame information functions ***/
    static_assert(qHasRevComp(n),                           "static assertion failed!");
    static_assert(!qHasRevComp(p),                          "static assertion failed!");
    SEQAN_ASSERT(qHasRevComp(n2));
    SEQAN_ASSERT(!qHasRevComp(p2));

    static_assert(qIsTranslated(x),                         "static assertion failed!");
    static_assert(!qIsTranslated(p),                        "static assertion failed!");

    static_assert(qNumFrames(x) == 6u,                      "static assertion failed!");
    static_assert(qNumFrames(tx) == 6u,                     "static assertion failed!");
    static_assert(qNumFrames(p) == 1u,                      "static assertion failed!");
    static_assert(qNumFrames(n) == 2u,                      "static assertion failed!");

    static_assert(sHasRevComp(tn),                          "static assertion failed!");
    static_assert(!sHasRevComp(n),                          "static assertion failed!");
    static_assert(!sHasRevComp(x),                          "static assertion failed!");
    SEQAN_ASSERT(!sHasRevComp(n2));

    static_assert(!sIsTranslated(x),                        "static assertion failed!");

    static_assert(sNumFrames(x) == 1u,                      "static assertion failed!");
    static_assert(sNumFrames(tx) == 6u,                     "static assertion failed!");

    /*** TEST string functions ***/
    // gcc can do the following, because it provides a constexpr strcmp implementation (not required by the standard)
    // it is not really important, though, just for demonstration purposes:
#if defined(COMPILER_GCC) || defined(COMPILER_LINTEL)
    static_assert(strcmp(_programTagToString(n), "BLASTN") == 0,       "static assertion failed!");
#endif
    SEQAN_ASSERT_EQ(_programTagToString(n2), "BLASTN");
}

SEQAN_DEFINE_TEST(test_blast_context_targs)
{
    {
        BlastIOContext<Blosum62/*, DYNAMIC, DYNAMIC*/> context;

        SEQAN_ASSERT(context.tabularSpec == BlastTabularSpec::UNKNOWN);
        SEQAN_ASSERT(context.blastProgram == BlastProgram::UNKNOWN);

        // these would cause build failures, because for ::DYNAMIC the selectors are not fixed at compile time
//         static_assert(context.blastProgram == BlastProgram::UNKNOWN,    "static assertion failed!");
//         static_assert(context.tabularSpec == BlastTabularSpec::UNKNOWN, "static assertion failed!");

        // we can set them to something else
        context.blastProgram = BlastProgram::BLASTX;
        context.tabularSpec = BlastTabularSpec::COMMENTS;

        SEQAN_ASSERT(context.blastProgram == BlastProgram::BLASTX);
        SEQAN_ASSERT(context.tabularSpec == BlastTabularSpec::COMMENTS);
    }

    {
        BlastIOContext<Blosum62, BlastProgram::BLASTX, BlastTabularSpec::COMMENTS> context;

        SEQAN_ASSERT(context.blastProgram == BlastProgram::BLASTX);
        SEQAN_ASSERT(context.tabularSpec == BlastTabularSpec::COMMENTS);

        // now they work at compile time, too
        static_assert(context.blastProgram == BlastProgram::BLASTX,    "static assertion failed!");
        static_assert(context.tabularSpec == BlastTabularSpec::COMMENTS, "static assertion failed!");

        // setting the same programs at run-time is ok
        context.blastProgram = BlastProgram::BLASTX;
        context.tabularSpec = BlastTabularSpec::COMMENTS;

        // but setting to something else would fail
//         context.blastProgram = BlastProgram::BLASTP;
//         context.tabularSpec = BlastTabularSpec::NO_COMMENTS;
    }

}

#endif  // SEQAN_TESTS_TEST_BLAST_MISC_H_
