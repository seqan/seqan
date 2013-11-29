// ========================================================================== 
//                 SeqAn - The Library for Sequence Analysis
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
// Author: Sascha Meiers <meiers@inf.fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_HEADER_CYCLIC_SHAPE_H
#define SEQAN_HEADER_CYCLIC_SHAPE_H

#include <seqan/index/shape_base.h>
#include <seqan/index/shape_gapped.h>

namespace seqan
{

    // ==========================================================================
    // Forwards
    // ==========================================================================

    template <typename TSpec> // needed for GappedShape<HardwiredShape<...> >
    struct GappedShape;

    template <unsigned LeftOffset, typename THardWiredShape, unsigned RightOffset>
    struct FixedShape;

    // GenericShape defined in shape_gapped.h via typedef

    typedef GappedShape<Default> GenericShape;

    //              .-------------- CyclicShape -----------.
    //             /                                        \
    //       GenericShape              or      .------- FixedShape -----.
    //                                        /             |            \
    //                                       /              |             \
    //                                     LOff ,      GappedShape ,     ROff
    //                                                      |
    //                                          .---- HardwiredShape -.
    //                                         /      /                \
    //                                        P00,   P01,     ...      P19


    // ==========================================================================
    // Class CyclicShape
    // ==========================================================================

/*!
 * @class CyclicShape
 * @headerfile <seqan/modifier.h>
 *
 * @brief A pattern of zeros and ones to mark "don't care"-positions in a text.
 *
 * @signature CyclicShape<TShapeSpec>
 *
 * @tparam TShapeSpec The specializing type. Default: @link GenericCyclicShape GenericShape @endlink.
 *
 * Unlike @link Shape @endlink, the CyclicShape does not perform hashing of q-grams.
 * It is instead useful to modify a text in such a way that zero
 * position are ignored. The pattern is applied repeatedly on the whole
 * text, as the example shows.
 *
 * The length of one pattern is saved as the member variable `span` (here 3), the
 * number of ones in the function `weight(shape)` (here 2).
 *
 * @code{.txt}
 * "this is a normal sequence"
 *  110110110110110110110110110 // CyclicShape 110
 *
 * => "ths s  nrml eqene"
 * @endcode
 *
 * Note that CyclicShapes can start and end with zero characters, which
 * is not allowed in Shape.
 *
 * @see ModCyclicShapeModifiedString
 */
    template <typename TCyclicSpec = GenericShape>
    class CyclicShape;
    // member: diffs, loffset, span


    // --------------------------------------------------------------------------
    // Metafunctions
    // --------------------------------------------------------------------------

/*!
 * @mfn CyclicShape#Size
 * @headerfile seqan/modifier.h
 *
 * @brief Size type for parameters used in CyclicShape.
 *
 * @signature Size<CyclicShape<TSpec> >::Type
 *
 * @tparam TSpec The CyclicShape specialization.
 *
 * @return TReturn Currently the return type `unsigned char`.
 *
 * @section Remarks
 */
    template <typename TSpec>
    struct Size<CyclicShape<TSpec> >
    {
        typedef unsigned char Type;
    };

    // --------------------------------------------------------------------------
    // Function weight
    // --------------------------------------------------------------------------

/*!
 * @fn CyclicShape#weight
 * @headerfile seqan/modifier.h
 *
 * @brief Converts a 0/1 string to a CyclicShape.
 *
 * @signature weight(cyclicShape)
 *
 * @tparam TSpec Specialisation of the CyclicShape.
 * @param cyclicShape CyclicShape<TSpec> object
 * @return weight of the CyclicShape (number of 1s). Return type is
 * `Size<CyclicShape<TSpec> >::Type`.
 */
    template <typename TSpec>
    inline typename Size<CyclicShape<TSpec> >::Type
    weight(CyclicShape<TSpec> const & s)
    {
        return static_cast<Size<CyclicShape<> >::Type> (length(s.diffs));
    }

    // ==========================================================================
    // Class Generic CyclicShape
    // ==========================================================================

/*!
 * @class GenericCyclicShape Generic CyclicShape
 * @extends CyclicShape
 * @headerfile <seqan/modifier.h>
 *
 * @brief CyclicShape initiated at runtime.
 *
 * @signature CyclicShape<GenericShape>
 *
 * @section Remarks
 *
 * Generic Shapes can be easily created using `stringToCyclicShape`,
 * but their usage comes with a significant loss in efficiency.
 * For longer sequences, prefer the @link FixedCylicShape @endlink.
 */
    template <>
    class CyclicShape<GenericShape>
    {
    public:
        String<Size<CyclicShape>::Type> diffs;
        Size<CyclicShape>::Type loffset, span;

        CyclicShape(): loffset(0), span(1)
        {
            // Pattern "1"
            appendValue(diffs, 1);
        }
        CyclicShape(CyclicShape const & other) : diffs(other.diffs), loffset(other.loffset), span(other.span)
        {}
    };

    // --------------------------------------------------------------------------
    // Function stringToCyclicShape()
    // --------------------------------------------------------------------------

/*!
 * @fn GenericCyclicShape#stringToCyclicShape
 * @headerfile seqan/modifier.h
 *
 * @brief Converts a 0/1 string to a Generic CyclicShape.
 *
 * @signature stringToCyclicShape(shape, TString bitmap)
 *
 * @tparam TString A string type.
 * @param shape Generic CyclicShape
 * @param bitmap 0/1 string. CyclicShapes may start and end with zeros,
 * but must contain at least one 1.
 *
 */
    template <typename TString>
    inline bool
    stringToCyclicShape(CyclicShape<GenericShape> & shape, TString const & bitmap)
    {
        // TODO: assert that length of bitmap is less than max value of Size<CyclicShape>::Size    SEQAN_ASSERT_GEQ(..., length(bitmap));

        typename Iterator<TString const>::Type it    = begin(bitmap);
        typename Iterator<TString const>::Type itBeg = begin(bitmap);
        typename Iterator<TString const>::Type itEnd = end(bitmap);

        Size<CyclicShape<GenericShape> >::Type countOnes = 0;
        for(; it != itEnd; ++it)
        {
            if(*it == '1') ++countOnes;
        }
        SEQAN_ASSERT_GT(countOnes, 0u);

        resize(shape.diffs, countOnes);

        typename Iterator<String<Size<CyclicShape<GenericShape> >::Type> >::Type diffIter = begin(shape.diffs);

        for(it = itBeg; *it != '1'; ++it)
        {}
        shape.loffset = (Size<CyclicShape<GenericShape> >::Type)(it - itBeg);

        countOnes = 1; // now used as the distance between 1s
        for(++it; it != itEnd; ++it)
        {
            if(*it == '1')
            {
                *diffIter = countOnes;
                ++diffIter;
                countOnes = 1;
            }
            else
                ++countOnes;
        }

        *diffIter = countOnes + shape.loffset;
        shape.span = (typename Size<CyclicShape<GenericShape> >::Type) (it - itBeg);
        return true;
    }


    // --------------------------------------------------------------------------
    // Function cyclicShapeToString()
    // --------------------------------------------------------------------------

/*!
 * @fn CyclicShape#cyclicShapeToString
 *
 * @brief Converts a given cyclic shape into a sequence of '1' (relevant position) and
 *        '0' (irrelevant position).
 *
 * @signature cyclicShapeToString(bitmap, cyclicShape)
 *
 * @param cyclicShape CyclicShape object. Types: @link CyclicShape @endlink
 * @param bitmap The resulting sequence object. Types: @link String @endlink
 *
 * @see stringToCyclicShape
 * @see shapeToString Equivalent for Shapes
 */
	template <typename TShapeString, typename TSpec>
	inline void
	cyclicShapeToString(
                  TShapeString &bitmap,
                  CyclicShape<TSpec> const &me)
	{
		clear(bitmap);
		if (weight(me) == 0) return;

        typename Size<CyclicShape<TSpec> >::Type i,j;
        for(i=0; i < me.loffset; ++i)
            appendValue(bitmap, '0');
        for(i=0; i < weight(me)-1; ++i)
        {
            appendValue(bitmap, '1');
            for(j=1; j < static_cast<typename Size<CyclicShape<TSpec> >::Type>(me.diffs[i]); ++j)
                appendValue(bitmap, '0');
        }
        appendValue(bitmap, '1');
        for(i=length(bitmap); i < me.span; ++i)
            appendValue(bitmap, '0');
    }


    // --------------------------------------------------------------------------
    // Function carePositions()
    // --------------------------------------------------------------------------

    template <typename TString, typename TSpec>
    inline void
    carePositions(TString & output, CyclicShape<TSpec> const & shape)
    {
        typedef typename Size<CyclicShape<TSpec> >::Type TPos;

        resize(output, weight(shape));
        TPos val = shape.loffset;
        output[0] = shape.loffset;

        for(TPos i=1; i< weight(shape); ++i) {
            val += shape.diffs[i-1];
            output[i] = val;
        }
    }


    // --------------------------------------------------------------------------
    // Function cyclicShapeToSuffixLengths()
    // --------------------------------------------------------------------------
    
/*!
 * @fn CyclicShape#cyclicShapeToSuffixLengths
 *
 * @brief Calculates the number of real characters of the gapped modification of suffixes
 *        shorter than the shape span.
 *
 * @signature cyclicShapeToSuffixLengths(TString & suffLengths, TShape const & cyclicShape)
 *
 * @param cyclicShape CyclicShape object. Types: @link CyclicShape @endlink
 * @param suffLengths String to be filled. Value type should be an integral type. Fixed length arrays should also work.
 *
 * @tparam TString A string type. The Value type should be integral, i.e. unsigned or int.
 * @tparam TShape CyclicShape.
 *
 * Given a CyclicShape, it calculates a String as long as the shape's span
 * that contains at position i the length of a modified suffix that starts
 * at position N-i. This is useful when short modified suffixes are lexicographically
 * compared as in the construction methods for gapped suffix arrays.
 * 
 * @see stringToCyclicShape
 * @see shapeToString Equivalent for Shapes
 */

    template <typename TString,
        typename TShape>
    inline void cyclicShapeToSuffixLengths(TString & suffLengths, TShape const & shape)
    {
        typedef typename Size<TString>::Type TSize;
        typedef typename Iterator<TString, Standard>::Type TIter;

        TSize w = weight(shape);
        TSize s = shape.span;
        TSize o = shape.loffset;
        
        SEQAN_ASSERT_GEQ(s, 1u);
        //SEQAN_ASSERT_EQ(length(suffLengths), s); // Disabled to support arrays too.
        
        // build cummulative some of distances first
        String<TSize> sums;
        resize(sums, w);
        sums[0] = o;
        for(unsigned i=0; i < w-1; ++i)
            sums[i+1] = sums[i]+shape.diffs[i];
        
        // write suffixLengths
        unsigned sumPos = 0;
        for(unsigned i=0; i < s; ++i)
            suffLengths[i] = (sumPos < w && i == sums[sumPos]) ? sumPos++ : sumPos;
    }
    

    // ==========================================================================
    // Class Fixed CyclicShape
    // ==========================================================================


/*!
 * @class FixedCyclicShape Fixed CyclicShape
 * @extends CyclicShape
 * @headerfile <seqan/modifier.h>
 *
 * @brief CyclicShape defined at compile time.
 *
 * @signature template <unsigned L, typename THardwiredShape, unsigned R>
 *            CycledShape<FixedShape<L, GappedShape<THardwiredShape>, R> >
 *
 * @tparam L Left offset. Number of leading zeros of the Shape.
 * @tparam R Right offset. Number of trailing zeros of the Shape.
 * @tparam THardwiredShape A specialization of @link HardwiredShape @endlink
 *
 * @section Remarks
 *
 * Fixed CylcicShapes contain their information at compile time, so in
 * most cases no object must be created.
 * The notation is similar to the one of @link HardwiredShape @endlink: Template paramters
 * code for the distance between ones. Additionally you have to specify the number of leading
 * and trailing zeros.
 *
 * The notation is chosen in such a way that predefined Shapes like
 * PatternHunter can be plugged into a CyclicShape.
 *
 * @section Example
 *
 * @code{.cpp}
 * // CyclicShape 01101100 is defined as:
 * typedef CyclicShape<FixedShape<1, GappedShape<HardwiredShape<1,2,1> >, 2> > TMyShape;
 *
 * //The predefined ShapePatternHunter can also be used:
 * typedef CyclicShape<FixedShape<0, Patternhunter, 0> > TCyclicPatternHunter;
 * @endcode
 *
 */
    template <unsigned LeftOffset, unsigned RightOffset,
        int P00, int P01, int P02, int P03, int P04,
        int P05, int P06, int P07, int P08, int P09,
        int P10, int P11, int P12, int P13, int P14,
        int P15, int P16, int P17, int P18, int P19 >
    class CyclicShape<FixedShape<LeftOffset, GappedShape<HardwiredShape<
        P00,P01,P02,P03,P04,P05,P06,P07,P08,P09,
        P10,P11,P12,P13,P14,P15,P16,P17,P18,P19	>
        >, RightOffset> >
    {
    public:
    typedef HardwiredShape<
            P00,P01,P02,P03,P04,P05,P06,P07,P08,P09,
            P10,P11,P12,P13,P14,P15,P16,P17,P18,P19>      THardWiredShape;

        enum { loffset = LeftOffset };
        enum { span = LeftOffset + LENGTH<THardWiredShape>::VALUE + RightOffset };
        enum { roffset = RightOffset };

        static const int diffs[];

        CyclicShape()
        {}

        CyclicShape(CyclicShape const & /* static members already exist */)
        {}
    };

    // --------------------------------------------------------------------------
    // Function Fixed CyclicShape weight()
    // --------------------------------------------------------------------------

    template <unsigned L, typename THardwiredShape, unsigned R>
    struct WEIGHT<CyclicShape<FixedShape<L, GappedShape<THardwiredShape>, R> > >
    {
        enum {VALUE = WEIGHT<THardwiredShape>::VALUE };
        typedef CyclicShape<FixedShape<L, GappedShape<THardwiredShape>, R> > TCyclicShape;
        typedef typename Size<TCyclicShape>::Type Type;
    };

    template <unsigned L, typename THardwiredShape, unsigned R>
    inline typename Size<CyclicShape<FixedShape<L, GappedShape<THardwiredShape>, R> > >::Type
    weight(CyclicShape<FixedShape<L, GappedShape<THardwiredShape>, R> > const & /* s */)
    {
        typedef CyclicShape<FixedShape<L, GappedShape<THardwiredShape>, R> > TCyclicShape;
        return WEIGHT<TCyclicShape>::VALUE;
    }

    // --------------------------------------------------------------------------
    // static const int CyclicShape<FixedShape<L, TSpec, R> >::diffs[]
    //
    // ugly work around to derive the static const array diffs[] from
    // the static const array HardwiredShape::DIFFS. In diffs, the first
    // zero-entry of DIFFS had to be changed.
    // Note: This requires that the non-used entries of
    // HardwiredShape::DIFFS are zero.
    // --------------------------------------------------------------------------

    template <unsigned entry, typename TFixedS>
    struct _arrEntry;

    // normal positions in diff[]
    template <unsigned entry, unsigned L, unsigned R,
        int P00, int P01, int P02, int P03, int P04,
        int P05, int P06, int P07, int P08, int P09,
        int P10, int P11, int P12, int P13, int P14,
        int P15, int P16, int P17, int P18, int P19>
    struct _arrEntry<entry, FixedShape<L, GappedShape<HardwiredShape
        <P00,P01,P02,P03,P04,P05,P06,P07,P08,P09,
        P10,P11,P12,P13,P14,P15,P16,P17,P18,P19> >, R> >
    {
        enum {VALUE = entry };
    };

    // zero positions in diff[]
    template <unsigned L, unsigned R,
        int P00, int P01, int P02, int P03, int P04,
        int P05, int P06, int P07, int P08, int P09,
        int P10, int P11, int P12, int P13, int P14,
        int P15, int P16, int P17, int P18, int P19>
    struct _arrEntry<0, FixedShape<L, GappedShape<HardwiredShape
        <P00,P01,P02,P03,P04,P05,P06,P07,P08,P09,
        P10,P11,P12,P13,P14,P15,P16,P17,P18,P19> >, R> >
    {
        enum {VALUE = L + R + 1};
    };

#define ENTRY(X) _arrEntry<X, FixedShape<L, GappedShape<HardwiredShape\
       <P00,P01,P02,P03,P04,P05,P06,P07,P08,P09,\
        P10,P11,P12,P13,P14,P15,P16,P17,P18,P19> >, R> >::VALUE

    template <unsigned L, unsigned R,
        int P00, int P01, int P02, int P03, int P04,
        int P05, int P06, int P07, int P08, int P09,
        int P10, int P11, int P12, int P13, int P14,
        int P15, int P16, int P17, int P18, int P19 >
    const int CyclicShape<FixedShape<L, GappedShape<HardwiredShape<
        P00,P01,P02,P03,P04,P05,P06,P07,P08,P09,
        P10,P11,P12,P13,P14,P15,P16,P17,P18,P19	>
        >, R> >::diffs[] =
    { ENTRY(P00), ENTRY(P01), ENTRY(P02), ENTRY(P03), ENTRY(P04),
      ENTRY(P05), ENTRY(P06), ENTRY(P07), ENTRY(P08), ENTRY(P09),
      ENTRY(P10), ENTRY(P11), ENTRY(P12), ENTRY(P13), ENTRY(P14),
      ENTRY(P15), ENTRY(P16), ENTRY(P17), ENTRY(P18), ENTRY(P19),
      ENTRY(0) };

#undef ENTRY
    
   
} //namespace


#endif
