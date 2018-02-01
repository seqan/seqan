// ==========================================================================
//                         Mason - A Read Simulator
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

#include "vcf_materialization.h"

#include "mason_types.h"

// ----------------------------------------------------------------------------
// Function contains()
// ----------------------------------------------------------------------------

// Returns true if haystack contains needle.

bool contains(seqan::CharString const & haystack, char const * needle)
{
    return strstr(toCString(haystack), needle);
}

// ----------------------------------------------------------------------------
// Function getSVLen()
// ----------------------------------------------------------------------------

// Get integer value of SVLEN=<number> from VCF INFO column strings.

int getSVLen(seqan::CharString const & str)
{
    seqan::DirectionIterator<seqan::CharString const, seqan::Input>::Type inputIter =
            directionIterator(str, seqan::Input());

    // Parse out key/value pairs and interpret SVLEN.
    seqan::CharString key, val;
    enum { IN_KEY, IN_VALUE } state = IN_KEY;
    for (; !atEnd(inputIter); ++inputIter)
    {
        if (*inputIter == '=')
        {
            state = IN_VALUE;
            continue;
        }
        else if (*inputIter == ';')
        {
            if (key == "SVLEN")
                return seqan::lexicalCast<int>(val);

            clear(val);
            clear(key);
            state = IN_KEY;
            continue;
        }
        else if (state == IN_KEY)
        {
            appendValue(key, *inputIter);
        }
        else  // (state == IN_VALUE)
        {
            appendValue(val, *inputIter);
        }
    }

    if (key == "SVLEN")
        return seqan::lexicalCast<int>(val);

    SEQAN_FAIL("Missing INFO SVLEN %s", toCString(str));
    return 0;
}

// ----------------------------------------------------------------------------
// Function getTargetPos()
// ----------------------------------------------------------------------------

// Returns reference name and position from TARGETPOS=<ref>:<pos> in VCF INFO column strings.

std::pair<seqan::CharString, int> getTargetPos(seqan::CharString const & str)
{
    seqan::DirectionIterator<seqan::CharString const, seqan::Input>::Type inputIter =
            directionIterator(str, seqan::Input());

    // Parse out key/value pairs and interpret SVLEN.
    seqan::CharString key, val;
    enum { IN_KEY, IN_VALUE } state = IN_KEY;
    for (; !atEnd(inputIter); ++inputIter)
    {
        if (*inputIter == '=')
        {
            state = IN_VALUE;
            continue;
        }
        else if (*inputIter == ';')
        {
            if (key == "TARGETPOS")
            {
                seqan::StringSet<seqan::CharString> xs;
                strSplit(xs, val, seqan::EqualsChar<':'>());
                SEQAN_CHECK(length(xs) == 2u, "TARGETPOS has invalid format %s", toCString(val));
                SEQAN_CHECK(!empty(xs[0]) && !empty(xs[1]), "TARGETPOS has invalid format %s", toCString(val));
                return std::make_pair(xs[0], seqan::lexicalCast<int>(xs[1]));
            }

            clear(val);
            clear(key);
            state = IN_KEY;
            continue;
        }
        else if (state == IN_KEY)
        {
            appendValue(key, *inputIter);
        }
        else  // (state == IN_VALUE)
        {
            appendValue(val, *inputIter);
        }
    }

    if (key == "TARGETPOS")
    {
        seqan::StringSet<seqan::CharString> xs;
        strSplit(xs, val, seqan::EqualsChar<':'>());
        SEQAN_CHECK(length(xs) == 2u, "TARGETPOS has invalid format %s", toCString(val));
        SEQAN_CHECK(!empty(xs[0]) && !empty(xs[1]), "TARGETPOS has invalid format %s", toCString(val));
        return std::make_pair(xs[0], seqan::lexicalCast<int>(xs[1]));
    }

    SEQAN_FAIL("Missing INFO TARGETPOS %s", toCString(str));
    return std::make_pair("", 0);
}

// ----------------------------------------------------------------------------
// Function VcfMaterializer::init()
// ----------------------------------------------------------------------------

void VcfMaterializer::init()
{
    if (!empty(vcfFileName))
    {
        // Open VCF stream.
        if (!open(vcfFileIn, toCString(vcfFileName)))
            throw MasonIOException("Could not open VCF stream.");

        // Read header.
        readHeader(vcfHeader, vcfFileIn);

        // Read first VCF record.
        if (!atEnd(vcfFileIn))
            readRecord(vcfRecord, vcfFileIn);

        // Get number of haplotypes in VCF file.
        SEQAN_ASSERT_NOT(empty(vcfRecord.genotypeInfos));
        seqan::StringSet<seqan::CharString> xs;
        seqan::DirectionIterator<seqan::CharString const, seqan::Input>::Type inputIter =
                directionIterator(vcfRecord.genotypeInfos[0], seqan::Input());
        numHaplotypes = 1;
        for (; !atEnd(inputIter); ++inputIter)
            numHaplotypes += (*inputIter == '|' || *inputIter == '/');
    }
    else
    {
        numHaplotypes = 1;
    }

    // Open input FASTA file and FAI.
    if (!open(faiIndex, toCString(fastaFileName)))
    {
        if (!build(faiIndex, toCString(fastaFileName)))
            throw MasonIOException("Could not build FAI index.");

        seqan::CharString faiPath = fastaFileName;
        append(faiPath, ".fai");
        if (!save(faiIndex, toCString(faiPath)))
            throw MasonIOException("Could not write FAI index.");
    }

    // Open methylation FASTA FAI file if given.
    if (!empty(methFastaFileName))
    {
        if (!open(methFaiIndex, toCString(methFastaFileName)))
        {
            if (!build(methFaiIndex, toCString(methFastaFileName)))
                throw MasonIOException("Could not build methylation levels FAI index.");

            seqan::CharString faiPath = methFastaFileName;
            append(faiPath, ".fai");
            if (!save(methFaiIndex, toCString(faiPath)))
                throw MasonIOException("Could not write methylation levels FAI index.");
        }
    }
}

// ----------------------------------------------------------------------------
// Function VcfMaterializer::_loadLevels()
// ----------------------------------------------------------------------------

void VcfMaterializer::_loadLevels(int rID)
{
    currentLevels.clear();

    std::stringstream ssTop, ssBottom;
    ssTop << sequenceName(faiIndex, rID) << "/TOP";
    unsigned idx = 0;
    if (!getIdByName(idx, methFaiIndex, ssTop.str().c_str()))
        throw MasonIOException("Could not find top levels in methylation FASTA.");
    readSequence(currentLevels.forward, methFaiIndex, idx);
    ssBottom << sequenceName(faiIndex, rID) << "/BOT";
    if (!getIdByName(idx, methFaiIndex, ssBottom.str().c_str()))
        throw MasonIOException("Could not find bottom levels in methylation FASTA.");
    readSequence(currentLevels.reverse, methFaiIndex, idx);
}

// ----------------------------------------------------------------------------
// Function VcfMaterializer::materializeNext()
// ----------------------------------------------------------------------------

bool VcfMaterializer::materializeNext(seqan::Dna5String & seq,
                                      std::vector<SmallVarInfo> & varInfos,
                                      std::vector<std::pair<int, int> > & breakpoints,
                                      int & rID,
                                      int & haplotype)
{
    return _materializeNext(seq, 0, varInfos, breakpoints, rID, haplotype);
}

bool VcfMaterializer::materializeNext(seqan::Dna5String & seq,
                                      MethylationLevels & levels,
                                      std::vector<SmallVarInfo> & varInfos,
                                      std::vector<std::pair<int, int> > & breakpoints,
                                      int & rID,
                                      int & haplotype)
{
    return _materializeNext(seq, &levels, varInfos, breakpoints, rID, haplotype);
}

// ----------------------------------------------------------------------------
// Function VcfMaterializer::_materializeNext()
// ----------------------------------------------------------------------------

bool VcfMaterializer::_materializeNext(seqan::Dna5String & seq,
                                       MethylationLevels * levels,
                                       std::vector<SmallVarInfo> & varInfos,
                                       std::vector<std::pair<int, int> > & breakpoints,
                                       int & rID,
                                       int & haplotype)
{
    if (levels)
        SEQAN_CHECK(!empty(methFastaFileName), "Must initialize with methylation FASTA file for levels");

    if (empty(vcfFileName))
    {
        if (currRID >= (int)(numSeqs(faiIndex) - 1))
            return false;
        currRID += 1;
        rID = currRID;
        readSequence(seq, faiIndex, currRID);
        if (levels && !empty(methFastaFileName))
        {
            _loadLevels(currRID);
            swap(*levels, currentLevels);
        }

        // Build identity PositionMap.
        TJournalEntries journal;
        reinit(journal, length(seq));
        posMap.reinit(journal);
        GenomicInterval gi(0, length(seq), 0, length(seq));
        seqan::String<PositionMap::TInterval> intervals;
        appendValue(intervals, PositionMap::TInterval(gi.svBeginPos, gi.svEndPos, gi));
        createIntervalTree(posMap.svIntervalTree, intervals);
        createIntervalTree(posMap.svIntervalTreeSTL, intervals);

        return true;
    }

    // Number of sequences.
    int numSeqs = length(contigNames(context(vcfFileIn)));

    // Stop if there are no more haplotypes to materialize.
    if (currRID >= (numSeqs - 1) && nextHaplotype == numHaplotypes)
        return false;

    // Load variants for next contig if necessary.
    if (currRID == -1 || nextHaplotype == numHaplotypes)
    {
        currRID += 1;
        nextHaplotype = 0;

        _loadVariantsForContig(contigVariants, currRID);
        readSequence(contigSeq, faiIndex, currRID);
        if (levels && !empty(methFastaFileName))
            _loadLevels(currRID);
    }

    // Materialize variants for the current haplotype.
    VariantMaterializer varMat(rng, contigVariants, *methOptions);
    if (levels)
        varMat.run(seq, posMap, *levels, varInfos, breakpoints, contigSeq, currentLevels, nextHaplotype);
    else
        varMat.run(seq, posMap, varInfos, breakpoints, contigSeq, nextHaplotype);

    // Write out rID and haploty
    rID = currRID;
    haplotype = nextHaplotype++;
    return true;
}

// ----------------------------------------------------------------------------
// Function VcfMaterializer::_loadVariantsForContig()
// ----------------------------------------------------------------------------

int VcfMaterializer::_loadVariantsForContig(Variants & variants, int rID)
{
    variants.clear();

    // Compute number of haplotypes.
    SEQAN_ASSERT_NOT(empty(vcfRecord.genotypeInfos));
    seqan::StringSet<seqan::CharString> xs;
    seqan::DirectionIterator<seqan::CharString const, seqan::Input>::Type inputIter =
            directionIterator(vcfRecord.genotypeInfos[0], seqan::Input());
    numHaplotypes = 1;
    for (; !atEnd(inputIter); ++inputIter)
        numHaplotypes += (*inputIter == '|' || *inputIter == '/');

    std::vector<seqan::VcfRecord> chunk;
    while (vcfRecord.rID != -1 && vcfRecord.rID <= rID)
    {
        // Translocations are the only SVs that are stored as breakends (BND).  We collect the BNDs in chunks of 6
        // which then represent a translocation.  Requiring that the BNDs are stored in adjacent chunks of 6 records
        // is a strong limitation but supporting more generic variations would be a bit too much here.
        if (contains(vcfRecord.info, "SVTYPE=BND"))
        {
            chunk.push_back(vcfRecord);
            if (chunk.size() == 6u)
            {
                _appendToVariantsBnd(variants, chunk);
                chunk.clear();
            }
        }
        else
        {
            SEQAN_CHECK(chunk.size() == 0u, "Found chunk of != 6 BND records!");
            _appendToVariants(variants, vcfRecord);
        }

        if (atEnd(vcfFileIn))
        {
            vcfRecord.rID = -1;
            continue;
        }
        readRecord(vcfRecord, vcfFileIn);
    }

    return 0;
}

// ----------------------------------------------------------------------------
// Function VcfMaterializer::_appendToVariants()
// ----------------------------------------------------------------------------

void VcfMaterializer::_appendToVariants(Variants & variants, seqan::VcfRecord const & vcfRecord)
{
    // Compute maximal length of alternative.
    unsigned altLength = 0;
    seqan::StringSet<seqan::CharString> alts;
    strSplit(alts, vcfRecord.alt, seqan::EqualsChar<','>());
    for (unsigned i = 0; i < length(alts); ++i)
        altLength = std::max(altLength, (unsigned)length(alts[i]));

    if (contains(vcfRecord.info, "SVTYPE"))  // Structural Variant
    {
        StructuralVariantRecord svRecord;
        svRecord.rId = vcfRecord.rID;
        svRecord.pos = vcfRecord.beginPos + 1;  // given with shift of -1
        svRecord.haplotype = 0;

        SEQAN_ASSERT_EQ(length(alts), 1u);

        if (contains(vcfRecord.info, "SVTYPE=INS"))  // Insertion
        {
            svRecord.kind = StructuralVariantRecord::INDEL;
            svRecord.size = getSVLen(vcfRecord.info);
            svRecord.seq = suffix(vcfRecord.alt, 1);
        }
        else if (contains(vcfRecord.info, "SVTYPE=DEL"))  // Deletion
        {
            svRecord.kind = StructuralVariantRecord::INDEL;
            svRecord.size = getSVLen(vcfRecord.info);
        }
        else if (contains(vcfRecord.info, "SVTYPE=INV"))  // Inversion
        {
            svRecord.kind = StructuralVariantRecord::INVERSION;
            svRecord.size = getSVLen(vcfRecord.info);
        }
        else if (contains(vcfRecord.info, "SVTYPE=DUP"))  // Duplication
        {
            svRecord.kind = StructuralVariantRecord::DUPLICATION;
            svRecord.size = getSVLen(vcfRecord.info);
            std::pair<seqan::CharString, int> pos = getTargetPos(vcfRecord.info);
            unsigned idx = 0;
            if (!getIdByName(idx, contigNamesCache(context(vcfFileIn)), pos.first))
                SEQAN_FAIL("Unknown sequence %s", toCString(pos.first));
            svRecord.targetRId = idx;
            svRecord.targetPos = pos.second - 1;
        }
        else if (contains(vcfRecord.info, "SVTYPE=BND"))  // Breakend (Must be Translocation)
        {
            SEQAN_FAIL("Unexpected 'SVTYPE=BND' at this place!");
        }
        else
        {
            SEQAN_FAIL("ERROR: Unknown SVTYPE!\n");
        }

        // Split the target variants.
        SEQAN_ASSERT_NOT(empty(vcfRecord.genotypeInfos));
        seqan::DirectionIterator<seqan::CharString const, seqan::Input>::Type inputIter =
                directionIterator(vcfRecord.genotypeInfos[0], seqan::Input());
        seqan::CharString buffer;
        svRecord.haplotype = 0;
        for (; !atEnd(inputIter); ++inputIter)
            if ((*inputIter == '|' || *inputIter == '/'))
            {
                if (!empty(buffer))
                {
                    unsigned idx = std::min(seqan::lexicalCast<unsigned>(buffer), 1u);
                    if (idx != 0u)  // if not == ref
                        appendValue(variants.svRecords, svRecord);
                }
                svRecord.haplotype++;
                clear(buffer);
            }
            else
            {
                appendValue(buffer, *inputIter);
            }
        if (!empty(buffer))
        {
            unsigned idx = std::min(seqan::lexicalCast<unsigned>(buffer), 1u);
            if (idx != 0u)  // if not == ref
                appendValue(variants.svRecords, svRecord);
        }
    }
    else if (length(vcfRecord.ref) == 1u && altLength == 1u)  // SNP
    {
        SnpRecord snpRecord;
        snpRecord.rId = vcfRecord.rID;
        snpRecord.pos = vcfRecord.beginPos;

        // Split the alternatives.
        seqan::StringSet<seqan::CharString> alternatives;
        strSplit(alternatives, vcfRecord.alt, seqan::EqualsChar<','>());

        // Split the target variants.
        SEQAN_ASSERT_NOT(empty(vcfRecord.genotypeInfos));
        seqan::DirectionIterator<seqan::CharString const, seqan::Input>::Type inputIter =
                directionIterator(vcfRecord.genotypeInfos[0], seqan::Input());
        seqan::CharString buffer;
        snpRecord.haplotype = 0;
        for (; !atEnd(inputIter); ++inputIter)
        {
            if (*inputIter == ':') break;
            if ((*inputIter == '|' || *inputIter == '/'))
            {
                if (!empty(buffer))
                {
                    unsigned idx = std::min(seqan::lexicalCast<unsigned>(buffer),
                                            (unsigned)length(alternatives));
                    if (idx != 0u)  // if not == ref
                    {
                        SEQAN_ASSERT_NOT(empty(alternatives[idx - 1]));
                        snpRecord.to = alternatives[idx - 1][0];
                        appendValue(variants.snps, snpRecord);
                    }
                }
                snpRecord.haplotype++;
                clear(buffer);
            }
            else
            {
                appendValue(buffer, *inputIter);
            }
        }
        if (!empty(buffer))
        {
            unsigned idx = std::min(seqan::lexicalCast<unsigned>(buffer),
                                    (unsigned)length(alternatives));
            if (idx != 0u)  // if not == ref
            {
                SEQAN_ASSERT_NOT(empty(alternatives[idx - 1]));
                snpRecord.to = alternatives[idx - 1][0];
                appendValue(variants.snps, snpRecord);
            }
        }
    }
    else  // Small Indel
    {
        SmallIndelRecord smallIndel;
        smallIndel.rId = vcfRecord.rID;
        smallIndel.pos = vcfRecord.beginPos + 1;

        if (contains(vcfRecord.alt, ",")) return;  // only one alternative
        if (length(vcfRecord.alt) > 1u && length(vcfRecord.ref) > 1u) return; // skip MNP and COMPLEX variants

        smallIndel.haplotype = 0;
        if (length(vcfRecord.ref) == 1u)  // insertion
        {
            smallIndel.seq = suffix(vcfRecord.alt, 1);
            smallIndel.size = length(smallIndel.seq);
        }
        else  // deletion
        {
            smallIndel.size = -(int)(length(vcfRecord.ref) - 1);
        }

        // Split the target variants.
        SEQAN_ASSERT_NOT(empty(vcfRecord.genotypeInfos));
        seqan::DirectionIterator<seqan::CharString const, seqan::Input>::Type inputIter =
                directionIterator(vcfRecord.genotypeInfos[0], seqan::Input());
        seqan::CharString buffer;
        smallIndel.haplotype = 0;
        for (; !atEnd(inputIter); ++inputIter)
            if ((*inputIter == '|' || *inputIter == '/'))
            {
                if (!empty(buffer))
                {
                    unsigned idx = std::min(seqan::lexicalCast<unsigned>(buffer), 1u);
                    if (idx != 0u)  // if not == ref
                        appendValue(variants.smallIndels, smallIndel);
                }
                smallIndel.haplotype++;
                clear(buffer);
            }
            else
            {
                appendValue(buffer, *inputIter);
            }
        if (!empty(buffer))
        {
            unsigned idx = std::min(seqan::lexicalCast<unsigned>(buffer), 1u);
            if (idx != 0u)  // if not == ref
                appendValue(variants.smallIndels, smallIndel);
        }
    }
}

// ----------------------------------------------------------------------------
// Function VcfMaterializer::_appendToVariantsBnd()
// ----------------------------------------------------------------------------

void VcfMaterializer::_appendToVariantsBnd(Variants & variants, std::vector<seqan::VcfRecord> const & vcfRecords)
{
    // Sanity check.
    SEQAN_CHECK(length(vcfRecords) == 6u, "Must be 6-tuple.");
    SEQAN_CHECK(vcfRecords[0].beginPos + 1 == vcfRecords[1].beginPos, "Invalid BND 6-tuple.");
    SEQAN_CHECK(vcfRecords[2].beginPos + 1 == vcfRecords[3].beginPos, "Invalid BND 6-tuple.");
    SEQAN_CHECK(vcfRecords[4].beginPos + 1 == vcfRecords[5].beginPos, "Invalid BND 6-tuple.");
    SEQAN_CHECK(vcfRecords[1].beginPos < vcfRecords[2].beginPos, "Wrong ordering.");
    SEQAN_CHECK(vcfRecords[3].beginPos < vcfRecords[4].beginPos, "Wrong ordering.");

    // Add translocation.
    StructuralVariantRecord svRecord;
    svRecord.kind = StructuralVariantRecord::TRANSLOCATION;
    svRecord.rId = vcfRecords[0].rID;
    svRecord.pos = vcfRecords[1].beginPos;
    svRecord.targetRId = svRecord.rId;
    svRecord.targetPos = vcfRecords[5].beginPos;
    svRecord.size = vcfRecords[3].beginPos - vcfRecords[1].beginPos;
    svRecord.haplotype = 0;

    // Split the target variants.
    SEQAN_ASSERT_NOT(empty(vcfRecord.genotypeInfos));
    seqan::DirectionIterator<seqan::CharString const, seqan::Input>::Type inputIter =
            directionIterator(vcfRecord.genotypeInfos[0], seqan::Input());
    seqan::CharString buffer;
    svRecord.haplotype = 0;
    for (; !atEnd(inputIter); ++inputIter)
        if ((*inputIter == '|' || *inputIter == '/'))
        {
            if (!empty(buffer))
            {
                unsigned idx = std::min(seqan::lexicalCast<unsigned>(buffer), 1u);
                if (idx != 0u)  // if not == ref
                    appendValue(variants.svRecords, svRecord);
            }
            svRecord.haplotype++;
            clear(buffer);
        }
        else
        {
            appendValue(buffer, *inputIter);
        }
    if (!empty(buffer))
    {
        unsigned idx = std::min(seqan::lexicalCast<unsigned>(buffer), 1u);
        if (idx != 0u)  // if not == ref
            appendValue(variants.svRecords, svRecord);
    }
}
