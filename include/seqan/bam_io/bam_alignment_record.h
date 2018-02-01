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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// The class BamAlignmentRecord, flag checking methods, flag constants.
// ==========================================================================

#ifndef INCLUDE_SEQAN_BAM_IO_BAM_RECORD_H_
#define INCLUDE_SEQAN_BAM_IO_BAM_RECORD_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

class BamAlignmentRecord;
inline void clear(BamAlignmentRecord & record);

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/*!
 * @enum BamFlags
 * @headerfile <seqan/bam_io.h>
 * @brief Shortcuts to the bitmask flags for BAM/SAM files.
 *
 * @signature enum BamFlags;
 *
 * @val BamFlags BAM_FLAG_MULTIPLE = 0x0001;
 * @brief Template has multiple fragments in sequencing.
 *
 * @val BamFlags BAM_FLAG_ALL_PROPER = 0x0002;
 * @brief All fragments in the template are properly mapped.
 *
 * @val BamFlags BAM_FLAG_UNMAPPED = 0x0004;
 * @brief This fragment is unmapped.
 *
 * @val BamFlags BAM_FLAG_NEXT_UNMAPPED = 0x0008;
 * @brief Next fragment in template is unmapped.
 *
 * @val BamFlags BAM_FLAG_RC = 0x0010;
 * @brief Fragment is reverse-complemented.
 *
 * @val BamFlags BAM_FLAG_NEXT_RC = 0x0020;
 * @brief Next fragment in template is reverse-complemented.
 *
 * @val BamFlags BAM_FLAG_FIRST = 0x0040;
 * @brief This fragment is the first one in its template.
 *
 * @val BamFlags BAM_FLAG_LAST = 0x0080;
 * @brief This fragment is the last one in its template (second in case of paired sequencing).
 *
 * @val BamFlags BAM_FLAG_SECONDARY = 0x0100;
 * @brief Secondary alignment.
 *
 * @val BamFlags BAM_FLAG_QC_NO_PASS = 0x0200;
 * @brief Does not pass quality controls.
 *
 * @val BamFlags BAM_FLAG_DUPLICATE = 0x0400;
 * @brief PCR or optical duplicate.
 *
 * @var BamFlags BAM_FLAG_SUPPLEMENTARY = 0x0800;
 * @brief Supplementary alignment.
 */

enum BamFlags
{
    BAM_FLAG_MULTIPLE      = 0x0001,
    BAM_FLAG_ALL_PROPER    = 0x0002,
    BAM_FLAG_UNMAPPED      = 0x0004,
    BAM_FLAG_NEXT_UNMAPPED = 0x0008,
    BAM_FLAG_RC            = 0x0010,
    BAM_FLAG_NEXT_RC       = 0x0020,
    BAM_FLAG_FIRST         = 0x0040,
    BAM_FLAG_LAST          = 0x0080,
    BAM_FLAG_SECONDARY     = 0x0100,
    BAM_FLAG_QC_NO_PASS    = 0x0200,
    BAM_FLAG_DUPLICATE     = 0x0400,
    BAM_FLAG_SUPPLEMENTARY = 0x0800
};

template <typename TValue>
struct BamTypeChar
{
    enum
    {
        VALUE =
            (IsSameType<TValue, char>::VALUE)?              'A':
            (IsSameType<TValue, int8_t>::VALUE)?            'c':
            (IsSameType<TValue, uint8_t>::VALUE)?           'C':
            (IsSameType<TValue, int16_t>::VALUE)?           's':
            (IsSameType<TValue, uint16_t>::VALUE)?          'S':
            (IsSameType<TValue, int32_t>::VALUE)?           'i':
            (IsSameType<TValue, uint32_t>::VALUE)?          'I':
            (IsSameType<TValue, float>::VALUE)?             'f':
//          (IsSameType<TValue, double>::VALUE)?            'd':
            (IsSequence<TValue>::VALUE)?                    'Z':
                                                            '?'
    };
};

// List of primitive BAM types (ordered by expected usage frequency)
typedef TagList<int32_t,
        TagList<uint32_t,
        TagList<float,
        TagList<int16_t,
        TagList<uint16_t,
        TagList<char,
        TagList<uint8_t,
        TagList<int8_t
//      TagList<double
        > > > > > > > > BamTagTypes;

// ----------------------------------------------------------------------------
// Class BamAlignmentRecord
// ----------------------------------------------------------------------------

/*!
 * @class BamAlignmentRecord
 * @headerfile <seqan/bam_io.h>
 * @implements FormattedFileRecordConcept
 * @signature class BamAlignmentRecord;
 * @brief Represent a record from a BAM or SAM file.
 *
 * @section Remarks
 *
 * While also used to represent SAM records, the type is called <tt>BamAlignmentRecord</tt> since the data directly
 * reflects a BAM records (0-based positions, identify references by id, and tags are stored in BAM format.
 *
 * @see BamFileIn
 * @see BamFileOut
 */

/*!
 * @var uint32_t BamAlignmentRecord::INVALID_POS
 * @brief Static member with invalid sentinel/position value (-1).
 *
 * @var uint32_t BamAlignmentRecord::INVALID_REFID
 * @brief Static member with invalid sentinel/position value (-1).
 *
 * @var uint32_t BamAlignmentRecord::INVALID_LEN
 * @brief Static member with invalid/sentinel reference ids (0 as in BAM/SAM).
 *
 * @var CharString BamAlignmentRecord::qName;
 * @brief The query/read name.
 *
 * Note that the reads of a template all of the same query name and are differentiated by their position
 * and the <tt>BAM_FLAG_FIRST</tt>/<tt>BAM_FLAG_LAST</tt> flag values.
 *
 * @var uint16_t BamAlignmentRecord::flag;
 * @brief The flag of this mapping.
 *
 * See @link BamFlags @endlink for flag constants and also see the <tt>hasFlag*()</tt> functions.
 *
 * @var int32_t BamAlignmentRecord::rID
 * @brief ID of reference for this fragment mapping (0-based, <tt>INVALID_REFID</tt> for '*' in SAM).
 *
 * @var int32_t BamAlignmentRecord::beginPos
 * @brief Begin position of the alignment (0-based, <tt>INVALID_POS</tt> for '0' in SAM).
 *
 * @var uint8_t BamAlignmentRecord::mapQ;
 * @brief Mapping quality (255 for '*').
 *
 * @var uint16_t BamAlignmentRecord::bin;
 * @brief The bin of the alignment, automatically computed when writing BAM.
 *
 * @var TCigarString BamAlignmentRecord::cigar;
 * @brief The CIGAR string for the BAM alignment (of type String<CigarElement<> >).
 *
 * @var int32_t BamAlignmentRecord::rNextId;
 * @brief The ID of the reference where the next fragment in this template aligns.
 *
 * <tt>INVALID_REFID</tt> for '*'.
 *
 * @var int32_t BamAlignmentRecord::pNext;
 * @brief Position on the reference where the next fragment in this template aligns.
 *
 * <tt>INVALID_POS</tt> for '*'.
 *
 * @var int32_t BamAlignmentRecord::tLen;
 * @brief The inferred template size.
 *
 * <tt>INVALID_LEN</tt> for '*'.
 *
 * @var CharString BamAlignmentRecord::seq;
 * @brief The fragment sequence.
 *
 * @var CharString BamAlignmentRecord::qual;
 * @brief The PHRED quality values of the sequence (as in SAM), empty for '*'.
 *
 * @var CharString BamAlignmentRecord::tags;
 * @brief Raw BAM tag string, use @link BamTagsDict @endlink for comfortable access.
 */

struct BamAlignmentRecordCore
{
            int32_t rID;
            int32_t beginPos;
    mutable uint32_t _l_qname:8;
            uint32_t mapQ:8;
    mutable uint32_t bin:16;
    mutable uint32_t _n_cigar:16;
            uint32_t flag:16;
    mutable int32_t _l_qseq;  // _l_qname, _n_cigar and _l_qseq for internal usage
            int32_t rNextId;
            int32_t pNext;
            int32_t tLen;
};

class BamAlignmentRecord : public BamAlignmentRecordCore
{
public:
    uint32_t _qID;  // TODO(holtgrew): Undocumented as of yet.
    String<CigarElement<> > cigar;
    CharString qName;
    IupacString seq;
    CharString qual;
    CharString tags;  // raw tags in BAM format
    CharString _buffer; // reusable internal buffer (used for I/O)

    static int32_t const INVALID_POS = -1;
    static int32_t const INVALID_REFID = -1;  // TODO(holtgrew): Rename to ...REF_ID.
    static int32_t const INVALID_LEN = 0;
    static uint32_t const INVALID_QID = 4294967295u;  // TODO(holtgrew): Undocumented as of yet.

    BamAlignmentRecord() : _qID(std::numeric_limits<unsigned>::max()) { clear(*this); }
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

/*!
 * @fn BamAlignmentRecord#clear
 * @brief Clear BamAlignmentRecord.
 *
 * @signature void clear(record);
 *
 * @param[in,out] record The BamAlignmentRecord to clear.
 *
 * Clears all strings and resets it to default initialization state.
 */

inline void
clear(BamAlignmentRecord & record)
{
    clear(record.qName);
    record.flag = 0;
    record._qID = std::numeric_limits<uint32_t>::max();
    record.rID = BamAlignmentRecord::INVALID_REFID;
    record.beginPos = BamAlignmentRecord::INVALID_POS;
    record.mapQ = 255;
    record.bin = 0;
    clear(record.cigar);
    record.rNextId = BamAlignmentRecord::INVALID_REFID;
    record.pNext = BamAlignmentRecord::INVALID_POS;
    record.tLen = BamAlignmentRecord::INVALID_LEN;
    clear(record.seq);
    clear(record.qual);
    clear(record.tags);
}

// ----------------------------------------------------------------------------
// Function hasFlagMultiple()
// ----------------------------------------------------------------------------

/*!
 * @fn BamAlignmentRecord#hasFlagMultiple
 * @headerfile <seqan/bam_io.h>
 * @brief Return true if a @link BamAlignmentRecord @endlink has the "multiple" flag set.
 *
 * @signature bool hasFlagMultiple(record);
 *
 * @param[in] record The BamAlignmentRecord to query.
 *
 * @return bool <tt>true</tt> if the flag is set, <tt>false</tt> otherwise.
 *
 * @see BamFlags
 */

inline bool
hasFlagMultiple(BamAlignmentRecord const & record)
{
    return (record.flag & BAM_FLAG_MULTIPLE) == BAM_FLAG_MULTIPLE;
}

// ----------------------------------------------------------------------------
// Function hasFlagAllProper()
// ----------------------------------------------------------------------------

/*!
 * @fn BamAlignmentRecord#hasFlagAllProper
 * @headerfile <seqan/bam_io.h>
 * @brief Return true if a @link BamAlignmentRecord @endlink has the "all properly aligned" flag set.
 *
 * @signature bool hasFlagAllProper(record);
 *
 * @param[in] record The BamAlignmentRecord to query.
 *
 * @return bool <tt>true</tt> if the flag is set, <tt>false</tt> otherwise.
 *
 * @see BamFlags
 */

inline bool
hasFlagAllProper(BamAlignmentRecord const & record)
{
    return (record.flag & BAM_FLAG_ALL_PROPER) == BAM_FLAG_ALL_PROPER;
}

// ----------------------------------------------------------------------------
// Function hasFlagUnmapped()
// ----------------------------------------------------------------------------

/*!
 * @fn BamAlignmentRecord#hasFlagUnmapped
 * @headerfile <seqan/bam_io.h>
 * @brief Return true if a @link BamAlignmentRecord @endlink has the "unmapped" flag set.
 *
 * @signature bool hasFlagUnmapped(record);
 *
 * @param[in] record The BamAlignmentRecord to query.
 *
 * @return bool <tt>true</tt> if the flag is set, <tt>false</tt> otherwise.
 *
 * @see BamFlags
 */

inline bool
hasFlagUnmapped(BamAlignmentRecord const & record)
{
    return (record.flag & BAM_FLAG_UNMAPPED) == BAM_FLAG_UNMAPPED;
}

// ----------------------------------------------------------------------------
// Function hasFlagNextUnmapped()
// ----------------------------------------------------------------------------

/*!
 * @fn BamAlignmentRecord#hasFlagNextUnmapped
 * @headerfile <seqan/bam_io.h>
 * @brief Return true if a @link BamAlignmentRecord @endlink has the "next unmapped" flag set.
 *
 * @signature bool hasFlagNextUnmapped(record);
 *
 * @param[in] record The BamAlignmentRecord to query.
 *
 * @return bool <tt>true</tt> if the flag is set, <tt>false</tt> otherwise.
 *
 * @see BamFlags
 */

inline bool
hasFlagNextUnmapped(BamAlignmentRecord const & record)
{
    return (record.flag & BAM_FLAG_NEXT_UNMAPPED) == BAM_FLAG_NEXT_UNMAPPED;
}

// ----------------------------------------------------------------------------
// Function hasFlagRC()
// ----------------------------------------------------------------------------

/*!
 * @fn BamAlignmentRecord#hasFlagRC
 * @headerfile <seqan/bam_io.h>
 * @brief Return true if a @link BamAlignmentRecord @endlink has the "reverse-complemented" flag set.
 *
 * @signature bool hasFlagRC(record);
 *
 * @param[in] record The BamAlignmentRecord to query.
 *
 * @return bool <tt>true</tt> if the flag is set, <tt>false</tt> otherwise.
 *
 * @see BamFlags
 */

inline bool
hasFlagRC(BamAlignmentRecord const & record)
{
    return (record.flag & BAM_FLAG_RC) == BAM_FLAG_RC;
}

// ----------------------------------------------------------------------------
// Function hasFlagNextRC()
// ----------------------------------------------------------------------------

/*!
 * @fn BamAlignmentRecord#hasFlagNextRC
 * @headerfile <seqan/bam_io.h>
 * @brief Return true if a @link BamAlignmentRecord @endlink has the "next reverse-complemented" flag set.
 *
 * @signature bool hasFlagNextRC(record);
 *
 * @param[in] record The BamAlignmentRecord to query.
 *
 * @return bool <tt>true</tt> if the flag is set, <tt>false</tt> otherwise.
 *
 * @see BamFlags
 */

inline bool
hasFlagNextRC(BamAlignmentRecord const & record)
{
    return (record.flag & BAM_FLAG_NEXT_RC) == BAM_FLAG_NEXT_RC;
}

// ----------------------------------------------------------------------------
// Function hasFlagFirst()
// ----------------------------------------------------------------------------

/*!
 * @fn BamAlignmentRecord#hasFlagFirst
 * @headerfile <seqan/bam_io.h>
 * @brief Return true if a @link BamAlignmentRecord @endlink has the "first in template" flag set.
 *
 * @signature bool hasFlagFirst(record);
 *
 * @param[in] record The BamAlignmentRecord to query.
 *
 * @return bool <tt>true</tt> if the flag is set, <tt>false</tt> otherwise.
 *
 * @see BamFlags
 */

inline bool
hasFlagFirst(BamAlignmentRecord const & record)
{
    return (record.flag & BAM_FLAG_FIRST) == BAM_FLAG_FIRST;
}

// ----------------------------------------------------------------------------
// Function hasFlagLast()
// ----------------------------------------------------------------------------

/*!
 * @fn BamAlignmentRecord#hasFlagLast
 * @headerfile <seqan/bam_io.h>
 * @brief Return true if a @link BamAlignmentRecord @endlink has the "last in template" flag set.
 *
 * @signature bool hasFlagLast(record);
 *
 * @param[in] record The BamAlignmentRecord to query.
 *
 * @return bool <tt>true</tt> if the flag is set, <tt>false</tt> otherwise.
 *
 * @see BamFlags
 */

inline bool
hasFlagLast(BamAlignmentRecord const & record)
{
    return (record.flag & BAM_FLAG_LAST) == BAM_FLAG_LAST;
}

// ----------------------------------------------------------------------------
// Function hasFlagSecondary()
// ----------------------------------------------------------------------------

/*!
 * @fn BamAlignmentRecord#hasFlagSecondary
 * @headerfile <seqan/bam_io.h>
 * @brief Return true if a @link BamAlignmentRecord @endlink has the "secondary" flag set.
 *
 * @signature bool hasFlagSecondary(record);
 *
 * @param[in] record The BamAlignmentRecord to query.
 *
 * @return bool <tt>true</tt> if the flag is set, <tt>false</tt> otherwise.
 *
 * @see BamFlags
 */

inline bool
hasFlagSecondary(BamAlignmentRecord const & record)
{
    return (record.flag & BAM_FLAG_SECONDARY) == BAM_FLAG_SECONDARY;
}

// ----------------------------------------------------------------------------
// Function hasFlagQCNoPass()
// ----------------------------------------------------------------------------

/*!
 * @fn BamAlignmentRecord#hasFlagQCNoPass
 * @headerfile <seqan/bam_io.h>
 * @brief Return true if a @link BamAlignmentRecord @endlink has the "did not pass QC" flag set.
 *
 * @signature bool hasFlagQCNoPass(record);
 *
 * @param[in] record The BamAlignmentRecord to query.
 *
 * @return bool <tt>true</tt> if the flag is set, <tt>false</tt> otherwise.
 *
 * @see BamFlags
 */

inline bool
hasFlagQCNoPass(BamAlignmentRecord const & record)
{
    return (record.flag & BAM_FLAG_QC_NO_PASS) == BAM_FLAG_QC_NO_PASS;
}

// ----------------------------------------------------------------------------
// Function hasFlagDuplicate()
// ----------------------------------------------------------------------------

/*!
 * @fn BamAlignmentRecord#hasFlagDuplicate
 * @headerfile <seqan/bam_io.h>
 * @brief Return true if a @link BamAlignmentRecord @endlink has the "duplicate" flag set.
 *
 * @signature bool hasFlagDuplicate(record);
 *
 * @param[in] record The BamAlignmentRecord to query.
 *
 * @return bool <tt>true</tt> if the flag is set, <tt>false</tt> otherwise.
 *
 * @see BamFlags
 */

inline bool
hasFlagDuplicate(BamAlignmentRecord const & record)
{
    return (record.flag & BAM_FLAG_DUPLICATE) == BAM_FLAG_DUPLICATE;
}

// ----------------------------------------------------------------------------
// Function hasFlagSupplementary()
// ----------------------------------------------------------------------------

/*!
 * @fn BamAlignmentRecord#hasFlagSupplementary
 * @headerfile <seqan/bam_io.h>
 * @brief Return true if a @link BamAlignmentRecord @endlink has the "supplementary" flag set.
 *
 * @signature bool hasFlagSupplementary(record);
 *
 * @param record The BamAlignmentRecord to query.
 *
 * @return bool <tt>true</tt> if the flag is set, <tt>false</tt> otherwise.
 */

inline bool
hasFlagSupplementary(BamAlignmentRecord const & record)
{
    return (record.flag & BAM_FLAG_SUPPLEMENTARY) == BAM_FLAG_SUPPLEMENTARY;
}

// ----------------------------------------------------------------------------
// Function getAlignmentLengthInRef()
// ----------------------------------------------------------------------------

/*!
 * @fn BamAlignmentRecord#getAlignmentLengthInRef
 * @headerfile <seqan/bam_io.h>
 * @brief Return the alignment length in the record's projection in the reference.
 *
 * @signature unsigned getAlignmentLengthInRef(record);
 *
 * @param[in] record The BamAlignmentRecord to compute length for.
 *
 * @return unsigned The alignment length.
 *
 * @see BamFlags
 */

inline unsigned
getAlignmentLengthInRef(BamAlignmentRecord const & record)
{
    unsigned l = 0;
    _getLengthInRef(l, record.cigar);
    return l;
}

// ----------------------------------------------------------------------------
// Function appendRawPod()
// ----------------------------------------------------------------------------

#if SEQAN_BIG_ENDIAN
inline void
enforceLittleEndian(BamAlignmentRecordCore & r)
{
    enforceLittleEndian(r.rID);
    enforceLittleEndian(r.beginPos);
    // _l_qname unchanged because 8bit
    // mapQ unchanged because 8bit
    r.bin       = htole16(r.bin);
    r._n_cigar  = htole16(r._n_cigar);
    r.flag      = htole16(r.flag);
    enforceLittleEndian(r._l_qseq);
    enforceLittleEndian(r.rNextId);
    enforceLittleEndian(r.pNext);
    enforceLittleEndian(r.tLen);
}

// This function is guarded so that we save the copy on little endian systems
template <typename TTarget>
inline void
appendRawPod(TTarget & target, BamAlignmentRecordCore r)
{
    enforceLittleEndian(r);
    appendRawPodImpl(target, r);
}
#endif

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_BAM_IO_BAM_RECORD_H_
