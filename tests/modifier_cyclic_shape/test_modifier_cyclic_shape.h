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

#ifndef SEQAN_TESTS_MODIFIER_CYCLIC_SHAPE_TEST_MODIFIER_CYCLIC_SHAPE_H_
#define SEQAN_TESTS_MODIFIER_CYCLIC_SHAPE_TEST_MODIFIER_CYCLIC_SHAPE_H_

#include <seqan/basic.h>
#include <seqan/stream.h>
#include <seqan/sequence.h>
#include <seqan/modifier.h>

using namespace seqan;


// test iteration over a modified string
template<typename TString, typename TSpec>
String<typename Value<TString>::Type>
test_iterator(TString & str, CyclicShape<TSpec> const & shape)
{
    typedef ModifiedString<TString, ModCyclicShape<CyclicShape<TSpec> > > TModStr;
    TModStr modStr(str, shape);
    String<typename Value<TString>::Type> returnString;

    typedef typename Iterator<TModStr, Standard>::Type TIter;
    TIter it = begin(modStr, Standard());
    TIter itEnd = end(modStr, Standard());
    for (; it != itEnd; ++it)
        appendValue(returnString, *it);
    return returnString;
}

// test constructors and assignments
SEQAN_DEFINE_TEST(test_modifier_cyclic_shape_modified_string_construct)
{
    typedef CharString                                          TString;
    typedef CyclicShape<FixedShape<1,
                                   GappedShape<HardwiredShape<1> >, 0> >                   TShape; // 011
    typedef ModifiedString<TString, ModCyclicShape<TShape> >    TModString;
    typedef ModifiedString<TString const,
                           ModCyclicShape<TShape> >                                TConstModString;

    TString s = "this is a string";
    TShape shape;

    TModString mod;
    TModString mod2(s);
    TModString mod3(s, shape);
    mod = mod2;

    TConstModString con;
    TConstModString con2(s);
    con = con2;
// TODO(meiers):
// compile error due to missing assignment operator in the general
// ModifiedString, see Ticket #1118
//    con = mod2;

    SEQAN_ASSERT_EQ(mod, "hi i astin");
    SEQAN_ASSERT_EQ(mod2, "hi i astin");
    SEQAN_ASSERT_EQ(mod3, "hi i astin");
    SEQAN_ASSERT_EQ(con, "hi i astin");
    SEQAN_ASSERT_EQ(con2, "hi i astin");
}


// test some string functions
SEQAN_DEFINE_TEST(test_modifier_cyclic_shape_modified_string_functions)
{
    CharString s = "012345678901234567890123456789012";
    //                       => 01100110110011011001101100110110011

    CyclicShape<GenericShape> shape;
    stringToCyclicShape(shape, "0110011");
    ModifiedString<CharString, ModCyclicShape<CyclicShape<GenericShape> > > modStr(s, shape);

    SEQAN_ASSERT_EQ(18u, length(modStr));
    SEQAN_ASSERT_EQ('1', modStr[0]);
    SEQAN_ASSERT_EQ('5', value(modStr, 2));
    SEQAN_ASSERT_EQ('0', back(modStr));
}


SEQAN_DEFINE_TEST(test_modifier_cyclic_shape_modified_iterator)
{
    CharString s = "012345678901234567890123456789012";

    CyclicShape<GenericShape> shape;
    stringToCyclicShape(shape, "0110011");

    typedef ModifiedString<CharString,
                           ModCyclicShape<CyclicShape<GenericShape> > > TModString;
    TModString modStr(s, shape);

    typedef Iterator<TModString, Standard>::Type TStandardModIter;
    typedef Iterator<TModString, Rooted>::Type   TRootedModIter;
    TStandardModIter a1 = begin(modStr);
    TStandardModIter a2 = begin(modStr, Standard());
    TStandardModIter a3 = begin(modStr, Rooted());
    TRootedModIter   a4 = begin(modStr);
    TRootedModIter   a5 = begin(modStr, Rooted());

    TStandardModIter e1 = end(modStr);
    TStandardModIter e2 = end(modStr, Standard());
    TStandardModIter e3 = end(modStr, Rooted());
    TStandardModIter e4 = end(modStr);
    TStandardModIter e5 = end(modStr, Rooted());

    // TODO(meiers): Some of the following functions throw compile errors
    //              so they are commented out. These problems are
    //              related to the design of ModifiedIterator... before
    //              spending a lot of work on solving them for this class
    //              alone I'd rather wait for the redesign of ModStrings
    //              in general.
    // goBegin(e1);
    // goBegin(e2);
    // goBegin(e3);
    // goBegin(e4);
    // goBegin(e5);
    // SEQAN_ASSERT_EQ(a1, e1);
    // SEQAN_ASSERT_EQ(a1, e2);
    // SEQAN_ASSERT_EQ(a1, e3);
    // SEQAN_ASSERT_EQ(a1, e4);
    // SEQAN_ASSERT_EQ(a1, e5);
    // goEnd(a1);
    // goEnd(a2);
    // goEnd(a3);
    // goEnd(a4);
    // goEnd(a5);
    // SEQAN_ASSERT_EQ(a1, end(modStr));
    // SEQAN_ASSERT_EQ(a1, end(modStr));
    // SEQAN_ASSERT_EQ(a1, end(modStr));
    // SEQAN_ASSERT_EQ(a1, end(modStr));
    // SEQAN_ASSERT_EQ(a1, end(modStr));
}


// 3 Tests for Generic Shapes
SEQAN_DEFINE_TEST(test_modifier_cyclic_shape_iterator_generic_alloc_charstring)
{
    CharString const STRING1 = "This is an alloc string that I now use the modifier on";
    // 00110100011010001101000110100011010001101000110100011010
    CyclicShape<GenericShape> shape;
    stringToCyclicShape(shape, "0011010");
    CharString s = test_iterator(STRING1, shape);
    SEQAN_ASSERT_EQ(s, "isin l sr ta nwe hodf o");
}
SEQAN_DEFINE_TEST(test_modifier_cyclic_shape_iterator_generic_mod_charstring)
{
    CharString const STRING1 = "This is an alloc string that I now use the modifier on";
    // 011000101100010110001011000101100010110001011000101100 <==
    ModifiedString<CharString const, ModReverse> revStr(STRING1);
    CyclicShape<GenericShape> shape;
    stringToCyclicShape(shape, "0011010");
    CharString s = test_iterator(revStr, shape);
    SEQAN_ASSERT_EQ(s, " riomees  It gi clnasih");
}
SEQAN_DEFINE_TEST(test_modifier_cyclic_shape_iterator_generic_infix_charstring)
{
    CharString const STRING1 = "This is an alloc string that I now use the modifier on";
    //       00110100011010001101000110100011010
    Infix<CharString const>::Type inf = infix(STRING1, 5, 40);
    CyclicShape<GenericShape> shape;
    stringToCyclicShape(shape, "0011010");
    CharString s = test_iterator(inf, shape);
    SEQAN_ASSERT_EQ(s, " a ocsngt Inus ");
}


// 3 Tests for Fixed Shapes
SEQAN_DEFINE_TEST(test_modifier_cyclic_shape_iterator_fixed_alloc_charstring)
{
    CharString const STRING1 = "This is an alloc string that I now use the modifier on";
    // 00110100011010001101000110100011010001101000110100011010
    typedef GappedShape<HardwiredShape<1, 2> > TInnerShape;
    typedef CyclicShape<FixedShape<2, TInnerShape, 1> > TShape;
    TShape shape;
    CharString s = test_iterator(STRING1, shape);
    SEQAN_ASSERT_EQ(s, "isin l sr ta nwe hodf o");
}
SEQAN_DEFINE_TEST(test_modifier_cyclic_shape_iterator_fixed_mod_charstring)
{
    CharString const STRING1 = "This is an alloc string that I now use the modifier on";
    // 011000101100010110001011000101100010110001011000101100 <==
    ModifiedString<CharString const, ModReverse> revStr(STRING1);
    CyclicShape<GenericShape> shape;
    stringToCyclicShape(shape, "0011010");
    CharString s = test_iterator(revStr, shape);
    SEQAN_ASSERT_EQ(s, " riomees  It gi clnasih");
}
SEQAN_DEFINE_TEST(test_modifier_cyclic_shape_iterator_fixed_infix_charstring)
{
    CharString const STRING1 = "This is an alloc string that I now use the modifier on";
    //    => 00110100011010001101000110100011010
    Infix<CharString const>::Type inf = infix(STRING1, 5, 40);
    CyclicShape<GenericShape> shape;
    stringToCyclicShape(shape, "0011010");
    CharString s = test_iterator(inf, shape);
    SEQAN_ASSERT_EQ(s, " a ocsngt Inus ");
}


#endif  // SEQAN_TESTS_MODIFIER_CYCLIC_SHAPE_TEST_MODIFIER_CYCLIC_SHAPE_H_
