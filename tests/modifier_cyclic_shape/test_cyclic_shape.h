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
// Author: Sascha Meiers <meiers@inf.fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_TESTS_MODIFIER_CYCLIC_SHAPE_TEST_CYCLIC_SHAPE_H_
#define SEQAN_TESTS_MODIFIER_CYCLIC_SHAPE_TEST_CYCLIC_SHAPE_H_

#include <seqan/basic.h>
#include <seqan/stream.h>
#include <seqan/sequence.h>
#include <seqan/modifier.h>

using namespace seqan;


SEQAN_DEFINE_TEST(test_modifier_cyclic_shape_cyclic_shape)
{
    CyclicShape<GenericShape> shape;
    CharString tmp;

    cyclicShapeToString(tmp, shape);
    SEQAN_ASSERT_EQ(tmp, "1");
    SEQAN_ASSERT_EQ(weight(shape), 1u);
    SEQAN_ASSERT_EQ(shape.span, 1u);

    unsigned DIFF[] = {1, 1, 1, 4, 1, 5};
    typedef CyclicShape<FixedShape<1, GappedShape<HardwiredShape<1, 1, 1, 4, 1> >, 3> > TShape;

    TShape s;
    cyclicShapeToString(tmp, s);
    SEQAN_ASSERT_EQ(tmp, "0111100011000");
    SEQAN_ASSERT_EQ(weight(s), 6u);
    SEQAN_ASSERT_EQ(static_cast<unsigned>(s.span), 13u);
    for (unsigned i = 0; i < 6; ++i)
        SEQAN_ASSERT_EQ(static_cast<unsigned>(s.diffs[i]), DIFF[i]);

    stringToCyclicShape(shape, tmp);
    SEQAN_ASSERT_EQ(tmp, "0111100011000");
    SEQAN_ASSERT_EQ(weight(shape), 6u);
    SEQAN_ASSERT_EQ(shape.span, 13u);
    for (unsigned i = 0; i < 6; ++i)
        SEQAN_ASSERT_EQ(static_cast<unsigned>(shape.diffs[i]), DIFF[i]);


}


#endif  // SEQAN_TESTS_MODIFIER_CYCLIC_SHAPE_TEST_MODIFIER_CYCLIC_SHAPE_H_
