// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
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
// Author: Rene Rahn <rene.rahn@fu-berlin.de>
// ==========================================================================

#ifndef INCLUDE_SEQAN_JOURNALED_STRING_TREE_JOURNALED_STRING_TREE_CONVERT_VCF_H_
#define INCLUDE_SEQAN_JOURNALED_STRING_TREE_JOURNALED_STRING_TREE_CONVERT_VCF_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

static const char * VCF_HEADER_REF_KEY = "reference";
static const char * VCF_HEADER_PLOIDY_KEY = "ploidy";

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Private Functions
// ============================================================================

namespace impl
{

inline const char *
pruneURIPrefix(CharString const & refFilePath)
{
    auto it = std::find(begin(refFilePath, Standard()), end(refFilePath, Standard()), "://");
    if (it == end(refFilePath, Standard()))
        return toCString(refFilePath);
    return &refFilePath[position(it, refFilePath) + 3];  // Convert to const char *
}

inline const char *
extractRefPath(VcfHeader const & header)
{
    forEach(header, [](VcfHeaderRecord const & headerRecord)
    {
        if (headerRecord.key == VCF_HEADER_REF_KEY)
            return pruneURIPrefix(headerRecord.value);
    });
}

inline unsigned
numOfHaploytpesPerSequence(VcfHeader const & header)
{
    forEach(header, [](VcfHeaderRecord const & headerRecord)
    {
        if (headerRecord.key == VCF_HEADER_PLOIDY_KEY)
            return lexicalCast<unsigned>(headerRecord.value);
    });
    return 1;  // Assume monoploid if no ploidy information is set.
}

template <typename TSequence, typename TConfig, typename TSpec>
inline void
recordGenotypes(JournaledStringTree<TSequence, TConfig, TSpec> & jst,
                VcfRecord const & record)
{
    // What if the genotype is
    // It would be cool to know if the Alphabet supports a WildCard character => Metafunction HasWildcard? -> DNA5 RNA5 IUPAC AMINO_ACID

//    auto genInfo = _extractGenotype(record); // Is stored as a StringSet<CharStrings> -> One entry per individual.
//
//    for (auto altId : genInfo)
//    {
//
//    }

}

template <typename TSequence, typename TConfig, typename TSpec>
inline bool
openJstFromFile(JournaledStringTree<TSequence, TConfig, TSpec> & jst,
               const char * fileName,
               FileOpenMode const mode,
               Vcf)
{
    //____________________________________________________________________
    // Read vcf header.

    VcfFileIn vcfFile;
    VcfHeader header;
    SEQAN_TRY
    {
        open(vcfFile, fileName, mode);
        readHeader(header, vcfFile);
    }
    SEQAN_CATCH(Exception e)
    {
        std::cerr << e.what() << std::endl;
        return false;
    }

    //____________________________________________________________________
    // Check header for polidy information.

    auto numSeqs = numOfHaploytpesPerSequence(header) * length(sampleNames(context(vcfFile)));

    //____________________________________________________________________
    // Read refrence sequence.

    const char* refPath = extractRefPath(header);
    SeqFileIn seqFile;
    SEQAN_TRY
    {
        open(seqFile, refPath, OPEN_RDONLY);
    }
    SEQAN_CATCH(Exception e)
    {
        std::cerr << e.what() << std::endl;
        return false;
    }

    TSequence hostSeq;
    CharString hostId;

    SEQAN_TRY
    {
        readRecord(hostId, hostSeq, seqFile);
    }
    SEQAN_CATCH(Exception e)
    {
        std::cerr << e.what() << std::endl;
        return false;
    }

    jst = JournaledStringTree<TSequence, TConfig, TSpec>(numSeqs);  // Copy assignment.
    setHost(jst, std::move(hostSeq));

    // Read configuration information.
    // Phased information?
    // IUPAC
    // Dna5

    VcfRecord record;
    while (!atEnd(vcfFile))
    {
        readRecord(record, vcfFile);
        recordGenotypes(jst, record);
    }


    return true;
}

}  // namespace impl

// ============================================================================
// Functions
// ============================================================================
}

#endif  // #ifndef INCLUDE_SEQAN_JOURNALED_STRING_TREE_JOURNALED_STRING_TREE_CONVERT_VCF_H_
