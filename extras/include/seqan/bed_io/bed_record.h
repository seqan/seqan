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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

#ifndef CORE_INCLUDE_SEQAN_BED_IO_BED_RECORD_H_
#define CORE_INCLUDE_SEQAN_BED_IO_BED_RECORD_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct Bed3_;
typedef Tag<Bed3_> Bed3;

struct Bed4_;
typedef Tag<Bed4_> Bed4;

struct Bed5_;
typedef Tag<Bed5_> Bed5;

struct Bed6_;
typedef Tag<Bed6_> Bed6;

struct Bed12_;
typedef Tag<Bed12_> Bed12;

// ----------------------------------------------------------------------------
// Class BedRecord
// ----------------------------------------------------------------------------

/**
.Class.BedRecord
..cat:BED I/O
..summary:Data structure for storing BED records.
..signature:BedRecord<TSpec>
..description:
BED files allow the easy representation of intervals on the genome.
Originally, they were designed for tracks in the UCSC genome browser.
The original format has 12 columns but often variants using fewer
columns with interpreted data are used and the rest is kept as application
dependent data.
..description:
The BedRecord class allows for storing BED records.
The various subclasses provide access to 3, 4, 5, 6, or 12 fields of the BED format.
For example, a $BedRecord<Bed5>$ has members variables for the first 5 columns of a BED file.
The remaining data is stored as the @Shortcut.CharString@ member variable $data$.
..remarks:
The $ref$ field is the name of the reference as loaded from the BED file.
The $rID$ field can be used to store a numeric reference id.
When loading without using a @Class.BedIOContext@, the $rID$ field remains set to @Memvar.BedRecord#INVALID_REFID@, otherwise the field is set to a different value.
..remarks:
Note that while the BED file format is 1-based, the coordinates in the BedRecord are 0-based.
..param.TSpec:The specialization to use.
...default:$Bed12$
..include:seqan/bed_io.h

.Memfunc.BedRecord#BedRecord
..class:Class.BedRecord
..summary:Constructor.
..description:Default constructor.
..signature:BedRecord::BedRecord()

.Memvar.BedRecord#INVALID_REFID
..class:Class.BedRecord
..summary:Constant for invalid references.
..signature:static const __int32 INVALID_REFID = -1

.Memvar.BedRecord#INVALID_POS
..class:Class.BedRecord
..summary:Constant for invalid positions.
..signature:static const __int32 INVALID_POS = -1

.Memvar.BedRecord#ref
..class:Class.BedRecord
..summary:Name of the interval's reference (@Shortcut.CharString@).
..signature:CharString ref

.Memvar.BedRecord#rID
..class:Class.BedRecord
..summary:Numeric id of the interval's reference ($__int32$, defaults to $INVALID_REFID$).
..signature:__int32 rID

.Memvar.BedRecord#beginPosition
..class:Class.BedRecord
..summary:Begin position on the reference.
..signature:__int32 beginPosition

.Memvar.BedRecord#endPosition
..class:Class.BedRecord
..summary:End position on the reference.
..signature:__int32 endPosition

.Memvar.BedRecord#data
..class:Class.BedRecord
..summary:Any data after the last position.
..signature:CharString data
*/

template <typename TSpec = Bed12>
class BedRecord;

// ----------------------------------------------------------------------------
// Class Bed3 BedRecord
// ----------------------------------------------------------------------------

/**
.Spec.Bed3 BedRecord
..cat:BED I/O
..general:Class.BedRecord
..summary:BedRecord with 3 fields.
..signature:class BedRecord<Bed3>
..description:
This BedRecord specialization stores the first three fields (ref, beginPos, endPos) of a BED file.
..include:seqan/bed_io.h
*/

template <>
class BedRecord<Bed3>
{
public:
    static const int INVALID_REFID = -1;
    static const int INVALID_POS = -1;

    // The chromosome name.
    CharString ref;
    // The id of the chromosome, -1 if not translated.
    __int32 rID;
    // The start position.
    __int32 beginPos;
    // The end position;
    __int32 endPos;
    // The remaining data from the file, unparsed.
    CharString data;

    BedRecord() : rID(INVALID_REFID), beginPos(INVALID_POS), endPos(INVALID_POS)
    {}

    void _clear()
    {
        rID = INVALID_REFID;
        beginPos = INVALID_POS;
        endPos = INVALID_POS;
        clear(ref);
        clear(data);
    }
};

// ----------------------------------------------------------------------------
// Class Bed4 BedRecord
// ----------------------------------------------------------------------------

/**
.Spec.Bed4 BedRecord
..cat:BED I/O
..general:Spec.Bed3 BedRecord
..summary:BedRecord with 4 fields.
..signature:class BedRecord<Bed3>
..description:
This BedRecord specialization stores the first four fields (ref, beginPos, endPos, name) of a BED file.
..include:seqan/bed_io.h

.Memvar.Bed4 BedRecord#name
..class:Spec.Bed4 BedRecord
..summary:The name of the interval (@Shortcut.CharString@).
*/

template <>
class BedRecord<Bed4> : public BedRecord<Bed3>
{
public:
    // The name of the feature.
    CharString name;

    BedRecord() : BedRecord<Bed3>()
    {}

    void _clear()
    {
        BedRecord<Bed3>::_clear();
        clear(name);
    }
};

// ----------------------------------------------------------------------------
// Class Bed5 BedRecord
// ----------------------------------------------------------------------------

/**
.Spec.Bed5 BedRecord
..cat:BED I/O
..general:Spec.Bed4 BedRecord
..summary:BedRecord with 5 fields.
..signature:class BedRecord<Bed5>
..description:
This BedRecord specialization stores the first five fields (ref, beginPos, endPos, name, score) of a BED file.
..include:seqan/bed_io.h

.Memvar.Bed5 BedRecord#score
..class:Spec.Bed5 BedRecord
..summary:The score of the interval (stored as @Shortcut.CharString@ to allow more flexible annotation).
..remarks:Storing the score as a @Shortcut.CharString@ is provided for compatibility with bedtools.
*/

template <>
class BedRecord<Bed5> : public BedRecord<Bed4>
{
public:
    // The score of the feature, stored as CharString for compatibility with bedtools.
    CharString score;

    BedRecord() : BedRecord<Bed4>()
    {}

    void _clear()
    {
        BedRecord<Bed4>::_clear();
        clear(score);
    }
};

// ----------------------------------------------------------------------------
// Class Bed6 BedRecord
// ----------------------------------------------------------------------------

/**
.Spec.Bed6 BedRecord
..cat:BED I/O
..general:Spec.Bed5 BedRecord
..summary:BedRecord with 6 fields.
..signature:class BedRecord<Bed6>
..description:
This BedRecord specialization stores the first six fields (ref, beginPos, endPos, name, score, strand) of a BED file.
..include:seqan/bed_io.h

.Memvar.Bed6 BedRecord#strand
..class:Spec.Bed6 BedRecord
..summary:The strand of the interval (stored as $char$, one of $.$, '-', and $+$).
..remarks:Defaults to '.'.
*/

template <>
class BedRecord<Bed6> : public BedRecord<Bed5>
{
public:
    // The strand of the feature.  One of '.', '-', and '+'.
    char strand;

    BedRecord() : BedRecord<Bed5>(), strand('.')
    {}

    void _clear()
    {
        BedRecord<Bed5>::_clear();
        strand = '.';
    }
};

// ----------------------------------------------------------------------------
// Class BedRgb
// ----------------------------------------------------------------------------

/**
.Class.BedRgb
..cat:BED I/O
..signature:class BedRgb
..summary:RGB color for @Spec.Bed12 BedRecord@.

.Memfunc.BedRgb#BedRgb
..class:Class.BedRgb
..signature:BedRgb::BedRgb()
..signature:BedRgb::BedRgb(red, green, blue)
..summary:Default constructor and initialization of integer RGB values.
..param.red:Integer red value $0-255$.
..param.green:Integer green value $0-255$.
..param.blue:Integer blue value $0-255$.

.Memvar.BedRgb#red
..class:Class.BedRgb
..summary:Red value of RGB color (default is $0$).

.Memvar.BedRgb#green
..class:Class.BedRgb
..summary:Green value of RGB color (default is $0$).

.Memvar.BedRgb#blue
..class:Class.BedRgb
..summary:Blue value of RGB color (default is $0$).
*/

class BedRgb
{
public:
    __int32 red, green, blue;

    BedRgb() : red(0), green(0), blue(0)
    {}

    BedRgb(int red, int green, int blue) : red(red), green(green), blue(blue)
    {}

    bool operator==(BedRgb const & other) const
    {
        return red == other.red && green == other.green && blue == other.blue;
    }

    bool operator!=(BedRgb const & other) const
    {
        return !(*this == other);
    }
};

// ----------------------------------------------------------------------------
// Class Bed12 BedRecord
// ----------------------------------------------------------------------------

/**
.Spec.Bed12 BedRecord
..cat:BED I/O
..general:Spec.Bed6 BedRecord
..summary:BedRecord with 12 fields.
..signature:class BedRecord<Bed12>
..description:
This BedRecord specialization stores all fields of a BED file.
..include:seqan/bed_io.h

.Memvar.Bed12 BedRecord#thickBegin
..class:Spec.Bed12 BedRecord
..summary:The begin position of thick drawing ($__int32$).

.Memvar.Bed12 BedRecord#thickEnd
..class:Spec.Bed12 BedRecord
..summary:The end position of thick drawing ($__int32$).

.Memvar.Bed12 BedRecord#itemRgb
..class:Spec.Bed12 BedRecord
..summary:RGB color of item (@Class.BedRgb@).

.Memvar.Bed12 BedRecord#blockCount
..class:Spec.Bed12 BedRecord
..summary:The number of blocks.

.Memvar.Bed12 BedRecord#blockSizes
..class:Spec.Bed12 BedRecord
..summary:The sizes of the blocks (@Spec.Alloc String@ of $__int32$).

.Memvar.Bed12 BedRecord#blockBegins
..class:Spec.Bed12 BedRecord
..summary:The begin positions of the blocks (@Spec.Alloc String@ of $__int32$).
*/

template <>
class BedRecord<Bed12> : public BedRecord<Bed6>
{
public:
    // The starting position of thick line for feature.
    __int32 thickBegin;
    // The end position of thick line for feature.
    __int32 thickEnd;
    // The color of the item.
    BedRgb itemRgb;
    // The number of blocks/exons for the feature.
    __int32 blockCount;
    // The list of block size.
    String<__int32> blockSizes;
    // List of block starts.
    String<__int32> blockBegins;

    BedRecord() : BedRecord<Bed6>(), thickBegin(INVALID_POS), thickEnd(INVALID_POS), blockCount(0)
    {}

    void _clear()
    {
        BedRecord<Bed6>::_clear();
        thickBegin = INVALID_POS;
        thickEnd = INVALID_POS;
        blockCount = 0;
        itemRgb = BedRgb();
        clear(blockSizes);
        clear(blockBegins);
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

/**
.Function.BedRecord#clear
..class:Class.BedRecord
..signature:void clear(record)
..summary:Reset BED record to state after default initialization.
..param.record:@Class.BedRecord@ to reset.
...type:Class.BedRecord
..include:seqan/bed_io.h
*/

template <typename TSpec>
void clear(BedRecord<TSpec> & record)
{
    record._clear();
}

}  // namespace seqan

#endif  // #ifndef CORE_INCLUDE_SEQAN_BED_IO_BED_RECORD_H_
