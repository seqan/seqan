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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================

#ifndef TESTS_INDEX_TEST_SHAPES_H
#define TESTS_INDEX_TEST_SHAPES_H


//////////////////////////////////////////////////////////////////////////////

namespace seqan
{

template <typename TShape1, typename TShape2>
void testShape(TShape1 shape1, TShape2 shape2, bool dump)
{
    if (dump) std::cout << std::endl;
                  // 012345678901234 len=15
    DnaString dna = "CGGTACGTAAGTTAG";
    DnaString dna1, dna2;

    TShape1 shape1b(shape1);
    TShape2 shape2b(shape2);

    SEQAN_ASSERT_EQ(length(shape1), length(shape2));
    SEQAN_ASSERT_EQ(weight(shape1), weight(shape2));

    Iterator<DnaString>::Type it = begin(dna);
    unsigned H1b, H2b;
    for (int i = length(dna); i >= 0; --i)
    {
        {
            unsigned H1 = hash(shape1, it, i);
            unsigned H2 = hash(shape2, it, i);
            if (i >= (int)length(shape1) && i != (int)length(dna))
            {
                H1b = hashNext(shape1b, it);
                H2b = hashNext(shape2b, it);
            } else {
                H1b = hash(shape1b, it, i);
                H2b = hash(shape2b, it, i);
            }

            unhash(dna1, H1, weight(shape1));
            unhash(dna2, H2, weight(shape2));
            if (dump) std::cout << std::dec << i << "\t" << std::hex << H1 << " " << H2 << ' ' << H1b << ' ' << H2b <<"\t" << dna1 << " " << dna2 << "\t" << (H1 == H2 && H1 == H1b && H1b == H2b);

            SEQAN_ASSERT_EQ(H1, H2);
            SEQAN_ASSERT_EQ(H1, H1b);
            SEQAN_ASSERT_EQ(H1b, H2b);

            if (i >= (int)length(shape1))
            {
                unsigned H1c = hash(shape1, it);
                unsigned H2c = hash(shape2, it);
                if (dump) std::cout << " " << (H1c == H2c && H1 == H1c);
            }
        }
        if (dump) std::cout << "\t";
        {
            unsigned H1 = hashUpper(shape1, it, i);
            unsigned H2 = hashUpper(shape2, it, i);

            unhash(dna1, H1, weight(shape1));
            unhash(dna2, H2, weight(shape2));
            if (dump) std::cout << std::dec << i << "\t" << std::hex << H1 << " " << H2 << "\t" << dna1 << " " << dna2 << "\t" << (H1 == H2);

            SEQAN_ASSERT_EQ(H1, H2);
        }
        if (dump) std::cout << std::endl;
        ++it;
    }
}

template <typename TShape>
void testHashInit(TShape shape)
{
    DnaString dna = "CGGTACGTAAGTTAG";
    Iterator<DnaString>::Type it = begin(dna);

    TShape shape2(shape);

    hash(shape, it);
    hashInit(shape2, it);
    hashNext(shape2, it);

    SEQAN_ASSERT_EQ(value(shape), value(shape2));
}

SEQAN_DEFINE_TEST(testShapes)
{
    Shape<Dna, SimpleShape > shapeA(6);
    testShape(shapeA, Shape<Dna, UngappedShape<6> >(), false);
    testHashInit(shapeA);
    testHashInit(Shape<Dna, UngappedShape<6> >());

                       // 012345678  len=9
    CharString pattern = "11100110100";
    Shape<Dna, GenericShape> shapeB(pattern);
    testShape(shapeB, Shape<Dna, GappedShape<HardwiredShape<1,1,3,1,2> > >(), false);
    testHashInit(shapeB);
    testHashInit(Shape<Dna, GappedShape<HardwiredShape<1,1,3,1,2> > >());

    pattern = "11110011";
    Shape<Dna, OneGappedShape> shapeC(pattern);
    testShape(shapeC, Shape<Dna, GappedShape<HardwiredShape<1,1,1,3,1> > >(), false);
    testHashInit(shapeC);
}

//////////////////////////////////////////////////////////////////////////////


} //namespace seqan

#endif //#ifndef SEQAN_HEADER_...
