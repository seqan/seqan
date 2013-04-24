// ==========================================================================
//                                 bam_stats
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
// Read a BAM file and collect statistics on the alignments therein.
// ==========================================================================

// TODO(holtgrew): Check for valid insert size / corcondantness.

#include <iostream>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>      // For printing SeqAn Strings.
#include <seqan/seq_io.h>
#include <seqan/bam_io.h>
#include <seqan/misc/misc_cmdparser.h>
#include <seqan/find.h>

#if SEQAN_HAS_ZLIB

using namespace seqan;

enum Format
{
    FORMAT_AUTO,
    FORMAT_SAM,
    FORMAT_BAM
};

struct Options
{
    bool showHelp;
    bool showVersion;
    bool realign;
    bool useNM;
    bool pairedEnd;
    bool goldStandard;
    
    unsigned verbosity;
    CharString refFile;
    CharString inFile;
    CharString inFastqFile;
    CharString goldStandardFile;
    CharString bestMatchFile;
    Format inFormat;

    int insertSizeMin;
    int insertSizeMax;

    Options() :
            insertSizeMin(-1),
            insertSizeMax(-1)
    {
        showHelp = false;
        showVersion = false;
        realign = false;
        useNM = false;
        pairedEnd = false;
        goldStandard = false;
        
        verbosity = 1;
        inFormat = FORMAT_AUTO;
    }
};

namespace  seqan {
    struct _OurScoreMatrix {};
    
    template <>
    struct ScoringMatrixData_<int, Dna5, _OurScoreMatrix> {
        enum {
            VALUE_SIZE = ValueSize<Dna5>::VALUE,
            TAB_SIZE = VALUE_SIZE * VALUE_SIZE
        };
        static inline int const * getData() {
            // The user defined data table.  In this case, we use the data from BLOSUM-30.
            static int const _data[TAB_SIZE] = {
                0, -3000, -3000, -3000, -3000,
                -3000, 0, -3000, -3000, -3000,
                -3000, -3000, 0, -3000, -3000,
                -3000, -3000, -3000, 0, -3000,
                -3000, -3000, -3000, -3000, -3000,
            };
            return _data;
        }
    };    
}

void
setupCommandLineParser(CommandLineParser & parser, Options const & options)
{
    addVersionLine(parser, "1.0");
    
    addTitleLine(parser, "*************");
    addTitleLine(parser, "* bam_stats *");
    addTitleLine(parser, "*************");
    addTitleLine(parser, "");
    addTitleLine(parser, "BAM Statistics.");
    addTitleLine(parser, "");
    addTitleLine(parser, "Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>");

    addUsageLine(parser, "bam_stats [OPTIONS] REF.fasta ALIGN.bam");
    
	addSection(parser, "General Options");
    addOption(parser, CommandLineOption("v", "verbose", "Enable verbose mode (show steps).", OptionType::Bool));
    addOption(parser, CommandLineOption("vv", "very-verbose", "Enable very verbose mode (show SAM lines and actual aligments).", OptionType::Bool));

	addSection(parser, "Input Specification");
    addOption(parser, CommandLineOption("S", "input-sam", "Input file is SAM (default: auto).", OptionType::Bool, options.inFormat == FORMAT_SAM));
    addOption(parser, CommandLineOption("B", "input-bam", "Input file is BAM (default: auto).", OptionType::Bool, options.inFormat == FORMAT_BAM));
    addOption(parser, CommandLineOption("i", "input-file", "Path to input, '-' for stdin.", OptionType::String, options.inFile));
    addOption(parser, CommandLineOption("fi", "fastq-input-file", "Path to fastq input (prefix excluding _1.fastq).", OptionType::String, options.inFastqFile));
    addOption(parser, CommandLineOption("g", "gold-standard", "Input file with gold standard", OptionType::String));

	addSection(parser, "Stats Specification");    
    addOption(parser, CommandLineOption("pe","paired-end","Paired-end mode.", OptionType::Bool));
    addOption(parser, CommandLineOption("r", "realign",   "Don't trust NM values and realign instead", OptionType::Bool));
    addOption(parser, CommandLineOption("nm","use-nm-tag","Get errors from NM values instead of CIGAR string", OptionType::Bool));
    addOption(parser, CommandLineOption("", "insert-size-min", "Minimal allowed insert size.  Default: any.", OptionType::Integer));
    addOption(parser, CommandLineOption("", "insert-size-max", "Maximal allowed insert size.  Default: any.", OptionType::Integer));
    
    addSection(parser, "Output Specification");
    addOption(parser, CommandLineOption("bm",  "best-match-file", "Output file with minimal found errors for each read", OptionType::String));

    requiredArguments(parser, 2);
}

void trimSeqHeaderToId(seqan::CharString & header)
{
    unsigned i = 0;
    for (; i < length(header); ++i)
        if (isspace(header[i]))
            break;
    resize(header, i);
}

int parseCommandLineAndCheck(Options & options,
                             CommandLineParser & parser,
                             int argc,
                             char const ** argv)
{
    bool stop = !parse(parser, argc, argv);
    if (stop)
        return 1;
    if (isSetLong(parser, "help"))
    {
        options.showHelp = true;
        return 0;
    }
    if (isSetLong(parser, "version"))
    {
        options.showVersion = true;
        return 0;
    }

    if (isSetLong(parser, "verbose"))
        options.verbosity = 2;
    if (isSetLong(parser, "very-verbose"))
        options.verbosity = 3;

    getOptionValueLong(parser, "paired-end", options.pairedEnd);

    getOptionValueLong(parser, "input-file", options.inFile);
    getOptionValueLong(parser, "best-match-file", options.bestMatchFile);
    getOptionValueLong(parser, "use-nm-tag", options.useNM);
    if (isSetLong(parser, "input-sam"))
        options.inFormat = FORMAT_SAM;
    if (isSetLong(parser, "input-bam"))
        options.inFormat = FORMAT_BAM;
    
    getOptionValueLong(parser, "fastq-input-file", options.inFastqFile);
    
    if (isSetLong(parser, "gold-standard"))
    {
        options.goldStandard = true;
        getOptionValueLong(parser, "gold-standard", options.goldStandardFile);
    }
    
    options.refFile = getArgumentValue(parser, 0);
    options.inFile = getArgumentValue(parser, 1);
    options.realign = isSetLong(parser, "realign");

    if (isSet(parser, "insert-size-min"))
        getOptionValueLong(parser, "insert-size-min", options.insertSizeMin);
    if (isSet(parser, "insert-size-max"))
        getOptionValueLong(parser, "insert-size-max", options.insertSizeMax);

	return 0;
}

struct Stats
{
    __uint64 numRecords;
    __uint64 alignedRecords;
    String<unsigned> editDistanceHisto;
    String<unsigned> mismatchHisto;
    String<unsigned> insertHisto;
    String<unsigned> deletionHisto;
    String<double> avrgQuality;

    Stats() : numRecords(0), alignedRecords(0)
    {}
};

struct Read
{
    int     contigId;
    int     beginPos;
    bool    reverseComplemented;
    
    unsigned char errors;
    unsigned char seqErrors;
    unsigned char SNPs;
    unsigned char indels;
    
    Read() :
        contigId(-1),
        beginPos(0),
        reverseComplemented(false),
        errors(0),
        seqErrors(0),
        SNPs(0),
        indels(0)
    {}
};

template <typename TSource, typename TSpec, typename TReference>
int
realignBamRecord(Align<TSource, TSpec> & result, TReference & reference, BamAlignmentRecord & record, unsigned maxIndels)
{
    resize(rows(result), 2);

    unsigned len = getAlignmentLengthInRef(record);
    __int64 posBegin = record.beginPos;

    if (!empty(record.cigar) && record.cigar[0].operation == 'S')
        posBegin -= record.cigar[0].count;

    __int64 posEnd = posBegin + len + maxIndels;
    posBegin = _max(0, posBegin - maxIndels);
    posEnd = _min(posEnd, (__int64)length(reference));

    assignSource(row(result, 0), infix(reference, posBegin, posEnd));
    assignSource(row(result, 1), record.seq);

    AlignConfig<true, false, false, true> alignConfig;
    Score<int, ScoreMatrix<Dna5, _OurScoreMatrix> > scoringScheme(-3001);
    int errors =  -globalAlignment(result, scoringScheme, alignConfig, NeedlemanWunsch()) / 3000;
	   
    if (toViewPosition(row(result, 0), 0) < toViewPosition(row(result, 1), 0))
        setClippedBeginPosition(row(result, 0), toViewPosition(row(result, 1), 0));
    if (toViewPosition(row(result, 0), length(source(row(result, 0)))) >
        toViewPosition(row(result, 1), length(source(row(result, 1)))))
    setClippedEndPosition(row(result, 0), toViewPosition(row(result, 1), length(source(row(result, 1))) ));
    return errors;
}

template <typename TReference>
int realignBamRecord(TReference & reference, BamAlignmentRecord & record, unsigned maxIndels)
{
    typedef Myers<AlignTextBanded<FindInfix, NMatchesNone_, NMatchesNone_>, True, void> TAlgorithmSpec;
    
    typedef String<Dna5>                                TReadSeq;
	typedef PatternState_<TReadSeq, TAlgorithmSpec>     TPatternState;
    
    typedef Segment<TReference, InfixSegment>           TReferenceInfix;
    typedef Finder<TReferenceInfix>                     TFinder;
    
    unsigned len = getAlignmentLengthInRef(record);
    __int64 posBegin = record.beginPos;
    
    if (!empty(record.cigar) && record.cigar[0].operation == 'S')
        posBegin -= record.cigar[0].count;
    
    __int64 posEnd = posBegin + len + maxIndels;
    posBegin = _max(0, posBegin - maxIndels);
    posEnd = _min(posEnd, (__int64)length(reference));
    if (posBegin > posEnd)
        return maxIndels + 1;  // Not possible to align.
    
    TReferenceInfix refInfix(reference, posBegin, posEnd);
    TReadSeq readSeq = record.seq;
    
    TFinder finder(refInfix);
    TPatternState patternState;
    
    int distance = maxIndels + 1;
    while (find(finder, readSeq, patternState, -static_cast<int>(maxIndels)))
        distance = std::min(distance, -getScore(patternState));
    
    return distance;
}

template <typename TString, typename TReadNames, typename TReadNameCache, typename TReadSeqs, typename TNamesSuffix>
bool loadFastqFile(TString & inFastqFile, TReadNames & readNames, TReadNameCache & readNameCache,
                   TReadSeqs & readSeqs, TNamesSuffix const & namesSuffix)
{
    MultiFasta multiFasta;
    if (!open(multiFasta.concat, toCString(inFastqFile), OPEN_RDONLY)) return false;
        
    AutoSeqFormat format;
    guessFormat(multiFasta.concat, format);
    split(multiFasta, format);
    unsigned seqCount = length(multiFasta);

    String<Dna5Q>	seq;
	CharString		qual;
	CharString      seqId;
    
    // Scan fastq file.
    for (unsigned i = 0; i < seqCount; ++i) 
	{
        assignSeqId(seqId, multiFasta[i], format);
		assignSeq(seq, multiFasta[i], format);
		assignQual(qual, multiFasta[i], format);
        
        // Truncate fastq seqId.
        CharString seqIdTruncated(seqId);
        trimSeqHeaderToId(seqIdTruncated);

        // Append left/mate suffix
        append(seqIdTruncated, namesSuffix);
        
        appendName(readNames, seqIdTruncated, readNameCache);
        appendValue(readSeqs, seq);
    }
    
    close(multiFasta.concat);
    
    return true;
}

template <typename TStreamOrReader, typename TSeqString, typename TSpec, typename TFormat>
int doWork(TStreamOrReader & reader, TStreamOrReader & greader,
           StringSet<CharString> & seqIds, StringSet<TSeqString, TSpec> & seqs,
           Options const & options, TFormat const & tag)
{
    typedef NameStoreCache<StringSet<CharString> >  TNameCache;
    typedef BamIOContext<StringSet<CharString> >    TBamContext;
    
    Stats stats;
    String<__uint64> qualSum;
    
    String<Read> goldStandard;

    // record best-match errors
    String<unsigned> minErrors;
    String<unsigned> minErrorsOrigin;
    String<unsigned> matchesCount;
    
    TNameCache seqIdsCache(seqIds);

    StringSet<CharString> refNames;
    TNameCache refNamesCache(refNames);
    TBamContext context(refNames, refNamesCache);
    
    StringSet<Dna5String> readSeqs;
    StringSet<CharString> readNames;
    TNameCache readNameCache(readNames);

    BamHeader header;
    
    unsigned line;
    unsigned oldNumRefs;

    // =================================================================================================================
    
    // Read fastq file.
    if (!empty(options.inFastqFile))
    {
        if (options.pairedEnd)
        {
            CharString fastqFile1(options.inFastqFile);
            append(fastqFile1, "_1.fastq");
            
            if (options.verbosity >= 2) std::cerr << "Loading reads from file " << fastqFile1 << std::endl;
            
            if (!loadFastqFile(fastqFile1, readNames, readNameCache, readSeqs, "/1"))
            {
                std::cerr << "Could not open file " << fastqFile1 << std::endl;
                return 1;
            }

            CharString fastqFile2(options.inFastqFile);
            append(fastqFile2, "_2.fastq");
            
            if (options.verbosity >= 2) std::cerr << "Loading reads from file " << fastqFile2 << std::endl;
            
            if (!loadFastqFile(fastqFile2, readNames, readNameCache, readSeqs, "/2"))
            {
                std::cerr << "Could not open file " << fastqFile2 << std::endl;
                return 1;
            }
        }
        else
        {
            CharString fastqFile(options.inFastqFile);
            append(fastqFile, ".fastq");
            
            if (options.verbosity >= 2) std::cerr << "Loading reads from file " << fastqFile << std::endl;

            if (!loadFastqFile(fastqFile, readNames, readNameCache, readSeqs, ""))
            {
                std::cerr << "Could not open file " << fastqFile << std::endl;
                return 1;
            }
        }
    }
    
    // =================================================================================================================
    
    // Read gold standard.
    if (options.goldStandard)
    {
        BamAlignmentRecord record;

        // Read header.
        if (options.verbosity >= 2)
            std::cerr << "Reading gold standard header" << std::endl;
        if (readRecord(header, context, greader, tag) != 0)
        {
            std::cerr << "Could not read gold standard header!" << std::endl;
            return 1;
        }
        
        // Initialize mapping from record.rID to seq id from reference.
        String<unsigned> rIdToSeqId;
        for (unsigned i = 0; i < length(refNames); ++i)
        {
            unsigned idx = 0;
            if (!getIdByName(seqIds, refNames[i], idx, seqIdsCache))
            {
                std::cerr << "Invalid reference name " << refNames[i] << " from gold standard header!\n";
                return 1;
            }
            appendValue(rIdToSeqId, idx);
        }
        
        line = 0;
        oldNumRefs = length(refNames);
        
        while (!atEnd(greader))
        {
            if (readRecord(record, context, greader, tag) != 0)
            {
                std::cerr << "Could not read gold standard alignment! lineNo=" << line << std::endl;
                return 1;
            }
            
            // Update mapping from SAM rID to index of sequence in reference sequences.
            if (oldNumRefs != length(refNames))
            {
                for (unsigned i = oldNumRefs; i < length(refNames); ++i)
                {
                    unsigned idx = 0;
                    if (!getIdByName(seqIds, refNames[i], idx, seqIdsCache))
                    {
                        std::cerr << "Invalid reference name " << refNames[i] << " from gold standard record!\n";
                        return 1;
                    }
                    appendValue(rIdToSeqId, idx);
                }
                oldNumRefs = length(refNames);
            }
            
            // Check reference id.
            if (record.rID != -1 && record.rID < static_cast<int>(length(seqs)))
            {
                Read read;
                
                BamTagsDict tagsDict(record.tags);
                unsigned idx;
                int tagValue;
                
                // Read reference id.
                read.contigId = record.rID;
                
                // Read begin position.
                read.beginPos = record.beginPos;
                
                // Read reverse complemented tag.
                if (hasFlagRC(record)) read.reverseComplemented = true;
                                
                // Read read name.
                CharString readName = record.qName;
                if (hasFlagMultiple(record))
                {
                    // Append mate discriminator if paired-end read.
                    if (hasFlagFirst(record))
                        append(readName, "/1");
                    else if (hasFlagLast(record))
                        append(readName, "/2");
                }
                int readId = -1;
                if (!getIdByName(readNames, readName, readId, readNameCache))
                {
                    if (!empty(options.inFastqFile))
                        std::cerr << "WARNING: A read in the gold standard is not in the fasta file: " << readName << std::endl;
                    // Append read sequence and name.
                    readId = length(readNames);
                    appendName(readNames, readName, readNameCache);
                    appendValue(readSeqs, record.seq);
                    if (hasFlagRC(record)) reverseComplement(back(readSeqs));
                }

                // Read total errors from NM tag.
                idx = 0;
                tagValue =-1;
                if (findTagKey(idx, tagsDict, "NM") && extractTagValue(tagValue, tagsDict, idx) && tagValue != -1)
                    read.errors = tagValue;
                else
                    std::cerr << "WARNING: Could find NM tag in gold standard " << record.qName << std::endl;

                // Read sequencing errors from XE tag.
                idx = 0;
                tagValue =-1;
                if (findTagKey(idx, tagsDict, "XE") && extractTagValue(tagValue, tagsDict, idx) && tagValue != -1)
                    read.seqErrors = tagValue;
                else
                    std::cerr << "WARNING: Could find XE tag in gold standard " << record.qName << std::endl;

                // Read SNPs from XS tag.
                idx = 0;
                tagValue =-1;
                if (findTagKey(idx, tagsDict, "XS") && extractTagValue(tagValue, tagsDict, idx) && tagValue != -1)
                    read.SNPs = tagValue;
                else
                    std::cerr << "WARNING: Could find XS tag in gold standard " << record.qName << std::endl;

                // Read indels from XI tag.
                idx = 0;
                tagValue =-1;
                if (findTagKey(idx, tagsDict, "XI") && extractTagValue(tagValue, tagsDict, idx) && tagValue != -1)
                    read.indels = tagValue;
                else
                    std::cerr << "WARNING: Could find XI tag in gold standard " << record.qName << std::endl;

                // Check for inconsistencies in error tags.
//                if (read.indels + read.SNPs + read.seqErrors != read.errors)
//                    std::cerr << "WARNING: Inconsistencies in error tags " << record.qName << std::endl;

                // Add read to gold standard.
                if (readId >= (int)length(goldStandard))
                    resize(goldStandard, readId + 1);
                goldStandard[readId] = read;
            }

            line += 1;
        }
    }

    // =================================================================================================================

    // Initialize stats.
    resize(minErrors, length(readNames), MaxValue<unsigned>::VALUE, Exact());
    resize(minErrorsOrigin, length(readNames), MaxValue<unsigned>::VALUE, Exact());
    resize(matchesCount, length(readNames), 0, Exact());

    // =================================================================================================================
    
    // Read header.
    if (options.verbosity >= 2)
        std::cerr << "Reading header" << std::endl;
    if (readRecord(header, context, reader, tag) != 0)
    {
        std::cerr << "Could not read header!" << std::endl;
        return 1;
    }

    // Initialize mapping from record.rID to seq id from reference.
    String<unsigned> rIdToSeqId;
    for (unsigned i = 0; i < length(refNames); ++i)
    {
        unsigned idx = 0;
        if (!getIdByName(seqIds, refNames[i], idx, seqIdsCache))
        {
            std::cerr << "Invalid reference name " << refNames[i] << " from SAM/BAM header!\n";
            return 1;
        }
        appendValue(rIdToSeqId, idx);
    }

    // Read alignments.
    if (options.verbosity >= 2)
        std::cerr << "Reading alignments" << std::endl;
    Align<Dna5String> align;
    __int64 reads = 0;
    
    oldNumRefs = length(refNames);
    line = 0;

    // Chunk-based reading of BAM/SAM records.
    String<BamAlignmentRecord> chunk;
    // Pre-read first record.
    BamAlignmentRecord record;
    if (readRecord(record, context, reader, tag) != 0)
    {
        std::cerr << "Could not read alignment! lineNo=" << line << std::endl;
        return 1;
    }
    trimSeqHeaderToId(record.qName);  // Some mappers create invalid query names.
    if (options.verbosity >= 3)
        write2(std::cerr, record, context, Sam());

    bool done = false;
    while (!done)
    {
        // Read next chunk of records (with the same query name).
        clear(chunk);
        appendValue(chunk, record);
        done = true;
        while (!atEnd(reader))
        {
            // Read alignment record.
            ++line;
            if (line % 100000 == 0) std::cerr << '.' << std::flush;

            if (readRecord(record, context, reader, tag) != 0)
            {
                std::cerr << "Could not read alignment! lineNo=" << line << std::endl;
                return 1;
            }
            trimSeqHeaderToId(record.qName);  // Some mappers create invalid query names.
            if (options.verbosity >= 3)
                write2(std::cerr, record, context, Sam());

            // Update mapping from SAM rID to index of sequence in reference sequences.
            if (oldNumRefs != length(refNames))
            {
                for (unsigned i = oldNumRefs; i < length(refNames); ++i)
                {
                    unsigned idx = 0;
                    if (!getIdByName(seqIds, refNames[i], idx, seqIdsCache))
                    {
                        std::cerr << "Invalid reference name " << refNames[i] << " from SAM/BAM record!\n";
                        return 1;
                    }
                    appendValue(rIdToSeqId, idx);
                }
                oldNumRefs = length(refNames);
            }

            if (front(chunk).qName != record.qName)
            {
                // If we break out here and not at the atEnd(record) above then we are not done yet.
                done = false;
                break;
            }
            appendValue(chunk, record);
        }

        // Fill in any missing read seqs.
        stats.numRecords += length(chunk);
        for (unsigned i = 0; i < length(chunk); ++i)
        {
            BamAlignmentRecord & record = chunk[i];

            // retrieve (possibly) missing read sequence
            int readId = -1;
            CharString readName = record.qName;
            if (hasFlagMultiple(record))
            {
                if (hasFlagFirst(record))
                    append(readName, "/1");
                else
                    append(readName, "/2");
            }
            if (getIdByName(readNames, readName, readId, readNameCache))
            {
                record.seq = readSeqs[readId];
                if (hasFlagRC(record))
                    reverseComplement(record.seq);
            }
            else
            {
                // This should never happen with a gold standard.
                //
                // However it DOES with incorrect flags by Hobbes if the flags are incorrect.  In this case, we mark the read as unmapped, it is then ignored below.
                if (options.goldStandard)
                {
                    std::cerr << "WARNING: A read unknown to the gold standard/fasta file has been found: " << readName << std::endl;
                    record.flag |= BAM_FLAG_UNMAPPED;
                }

                if (readId == -1)
                {
                    // Add read name and entry.
                    readId = length(readNames);
                    appendValue(readSeqs, record.seq);
                    if (hasFlagRC(record))
                        reverseComplement(back(readSeqs));
                    appendName(readNames, record.qName, readNameCache);
                    appendValue(minErrors, maxValue<unsigned>());
                    appendValue(matchesCount, 0);
                }
            }
            record._qID = readId;
        }

        // Remove any unmapped records or records with unmapped mate.
        {
            String<BamAlignmentRecord> tmpChunk;
            for (unsigned i = 0; i < length(chunk); ++i)
            {
                if (hasFlagUnmapped(chunk[i]) || chunk[i].rID == -1 || hasFlagNextUnmapped(chunk[i]))
                    continue;
                appendValue(tmpChunk, chunk[i]);
            }
            move(chunk, tmpChunk);
        }

        // In case of PE data: Sort such that the mates in the pairs folow each other.
        if (!empty(chunk) && hasFlagMultiple(front(chunk)))
        {
            using std::swap;  // seqan::swap still is available
            for (unsigned i = 0; i < length(chunk); /*nop*/)
            {
                unsigned mate = i;
                for (unsigned j = i + 1; j < length(chunk); ++j)
                {
                    if (chunk[i].rID == chunk[j].rNextId && chunk[i].beginPos == chunk[j].pNext &&
                        chunk[j].rID == chunk[i].rNextId && chunk[j].beginPos == chunk[i].pNext)
                    {
                        mate = j;
                        break;
                    }
                }
                if (mate == i)
                {
                    std::cerr << "WARNING: Could not find any mate for record\n";
                    write2(std::cerr, chunk[i], context, Sam());
                    if (options.verbosity >= 3)
                    {
                        std::cerr << "Chunk (of same-query name records is\n,--\n";
                        for (unsigned k = 0; k < length(chunk); ++k)
                        {
                            std::cerr << "| ";
                            write2(std::cerr, chunk[k], context, Sam());
                        }
                        std::cerr << "`--\n";
                    }
                    // Purge this record.
                    std::rotate(begin(chunk, Standard()) + i, begin(chunk, Standard()) + i + 1, end(chunk, Standard()));
                    resize(chunk, length(chunk) - 1);
                }
                else
                {
                    swap(chunk[i + 1], chunk[mate]);
                    if (hasFlagMultiple(chunk[i]) && hasFlagLast(chunk[i]))
                    {
                        //SEQAN_CHECK(hasFlagFirst(chunk[i + 1]), "Must be first in pair!");
                        swap(chunk[i], chunk[i + 1]);
                    }
                    i += 2;
                }
            }
        }

        // Realign or compute alignment from CIGAR; Maybe interpret NM Tag.
        String<unsigned> editDistances;
        resize(editDistances, length(chunk), 0);
        for (unsigned chunkI = 0; chunkI < length(chunk); ++chunkI)
        {
            BamAlignmentRecord & record = chunk[chunkI];
            unsigned & editDistance = editDistances[chunkI];

            SEQAN_CHECK(record.rID != -1 && record.rID < static_cast<int>(length(seqs)),
                        "Must be aligned!");

            unsigned seqId = rIdToSeqId[record.rID];

            // realign in a window 10bp left and right of the original alignment
            if (options.realign)
            {
                editDistance = realignBamRecord(seqs[seqId], record, 10);
            }
            else if (options.useNM)
            {
                // trust NM instead of CIGAR string
                unsigned idx = 0;
                int nmValue =-1;
                BamTagsDict tagsDict(record.tags);
                if (findTagKey(idx, tagsDict, "NM") && extractTagValue(nmValue, tagsDict, idx) && nmValue != -1)
                    editDistance = nmValue;
                else
                    std::cerr << "WARNING: Could find NM tag in match of read " << record.qName << std::endl;
            }
            else
            {
                // convert soft clipping into match/mismatches
                if (empty(record.cigar))
                {
                    std::cerr << "ERROR: Realigning is disabled by CIGAR string missing in record!\n";
                    write2(std::cerr, record, context, Sam());
                    return 1;
                }
                if (record.cigar[0].operation == 'S')
                    record.beginPos -= record.cigar[0].count;

                for (unsigned i = 0; i < length(record.cigar); ++i)
                {
                    if (record.cigar[i].operation == 'S')
                        record.cigar[i].operation = 'M';
                }

                bamRecordToAlignment(align, seqs[seqId], record);

                if (options.verbosity >= 3)
                    std::cerr << align << std::endl;

                typedef Align<Dna5String>             TAlign;
                typedef typename Row<TAlign>::Type    TRow;
                typedef typename Iterator<TRow>::Type TRowIter;
                unsigned posReadFwd = 0;
                for (TRowIter it0 = iter(row(align, 0), 0), it1 = iter(row(align, 1), 0); !atEnd(it0); goNext(it0), goNext(it1))
                {
                    unsigned posRead = posReadFwd;
                    // is read aligned to reverse strand?
                    if (hasFlagRC(record))
                        posRead = (length(record.seq)) - posReadFwd;
                    SEQAN_ASSERT_LEQ(posRead, length(record.seq));

                    if (isGap(it0) && isGap(it1))
                        continue;
                    if (isGap(it0) || isGap(it1))
                    {
                        if (isGap(it0))
                        {
                            unsigned len = length(stats.insertHisto);
                            resize(stats.insertHisto, std::max(len, posRead + 1), 0);
                            stats.insertHisto[posRead] += 1;
                            posReadFwd += 1;
                        }
                        else
                        {
                            unsigned len = length(stats.deletionHisto);
                            resize(stats.deletionHisto, std::max(len, posRead + 1), 0);
                            stats.deletionHisto[posRead] += 1;
                        }
                        editDistance += 1;
                        continue;
                    }
                    if (convert<char>(*it0) == 'N' || convert<char>(*it1) == 'N' || *it0 != *it1)
                    {
                        unsigned len = length(stats.mismatchHisto);
                        resize(stats.mismatchHisto, std::max(len, posRead + 1), 0);
                        stats.mismatchHisto[posRead] += 1;
                        editDistance += 1;
                    }
                    posReadFwd += 1;
                }

                resize(qualSum, std::max(length(qualSum), length(record.qual)), 0);
                if (!hasFlagRC(record))
                {
                    // read aligns to forward strand
                    for (unsigned i = 0; i < length(record.qual); ++i)
                        qualSum[i] += record.qual[i] - '!';
                } else
                {
                    // read aligns to reverse strand
                    for (unsigned i = 0; i < length(record.qual); ++i)
                        qualSum[(length(record.qual) - 1) - i] += record.qual[i] - '!';
                }
                ++reads;
            }
        }

        // Filter out the reads that do not mapp acording to the gold standard or concordantly.  This has to be done
        // after realigning.
        if (!empty(chunk) && hasFlagMultiple(front(chunk)))
        {
            String<BamAlignmentRecord> tmpChunk;
            String<unsigned> tmpEditDistances;

            SEQAN_CHECK(length(chunk) % 2 == 0u, "Chunk length must be even.");
            for (unsigned i = 0; i < length(chunk); i += 2)
            {
                if (chunk[i].rID != chunk[i + 1].rID)
                    continue;  // Invalid: Not even on same reference.
                if (chunk[i].beginPos < chunk[i + 1].beginPos && (hasFlagRC(chunk[i]) || !hasFlagRC(chunk[i + 1])))
                    continue;  // Invalid: Not aligned --> <--, case where chunk[i] is left.
                if (chunk[i].beginPos > chunk[i + 1].beginPos && (!hasFlagRC(chunk[i]) || hasFlagRC(chunk[i + 1])))
                    continue;  // Invalid: Not aligned --> <--, case where chunk[i] is right.
                chunk[i].tLen = chunk[i + 1].beginPos - chunk[i].beginPos;
                chunk[i + 1].tLen = chunk[i].beginPos - chunk[i + 1].beginPos;
                __int32 tlen = abs(chunk[i].tLen);
                if ((options.insertSizeMax != -1 && tlen > options.insertSizeMax) ||
                    (options.insertSizeMin != -1 && tlen < options.insertSizeMin))
                    continue;  // Invalid: Insert size not allowed.
                appendValue(tmpChunk, chunk[i]);
                appendValue(tmpEditDistances, editDistances[i]);
                appendValue(tmpChunk, chunk[i + 1]);
                appendValue(tmpEditDistances, editDistances[i + 1]);
            }

            move(chunk, tmpChunk);
            move(editDistances, tmpEditDistances);
        }

        // Now, finally we have a chunk of alignments that we can compute statistics on.
        stats.alignedRecords += length(chunk);

        for (unsigned i = 0; i < length(chunk); ++i)
        {
            BamAlignmentRecord & record = chunk[i];
            int readId = record._qID;
            unsigned editDistance = editDistances[i];

            // Update the maps minErrors, matchesCount, and minErrorsOrigin.  These will be used to print the results.
            if (options.goldStandard)
            {
                // Check if match corresponds to the true origin in gold standard.
                if ((record.beginPos >= goldStandard[readId].beginPos - 2 * goldStandard[readId].errors &&
                     record.beginPos <= goldStandard[readId].beginPos + 2 * goldStandard[readId].errors) &&
                    record.rID == goldStandard[readId].contigId &&
                    hasFlagRC(record) == goldStandard[readId].reverseComplemented)
                {
                    minErrorsOrigin[readId] = _min(minErrors[readId], editDistance);
                }

                minErrors[readId] = _min(minErrors[readId], editDistance);
                matchesCount[readId]++;
            }
            else
            {
                // update best-match errors
                if (readId != -1)
                {
                    minErrors[readId] = _min(minErrors[readId], editDistance);
                    matchesCount[readId]++;
                }
                else
                {
                    // add read name and entry
                    appendValue(readSeqs, record.seq);
                    if (hasFlagRC(record)) reverseComplement(back(readSeqs));
                    CharString readName = record.qName;
                    if (hasFlagMultiple(record))
                    {
                        // Append mate discriminator if paired-end read.
                        if (hasFlagFirst(record))
                            append(readName, "/1");
                        else if (hasFlagLast(record))
                            append(readName, "/2");
                    }
                    appendName(readNames, readName, readNameCache);

                    // Update stats.
                    appendValue(minErrors, editDistance);
                    appendValue(matchesCount, 1);
                }
            }

            /*if (options.verbosity >= 2 && editDistance > 20)
            {
                std::cerr << "EDIT DISTANCE == " << editDistance << "\n"
                          << "ALIGNMENT\n" << align
                          << "READ NAME\t" << record.qName << "\n";
            }
            */

            if (options.verbosity >= 3)
                std::cerr << "edit distance: " << editDistance << std::endl;
            unsigned len = length(stats.editDistanceHisto);
            resize(stats.editDistanceHisto, std::max(len, editDistance + 1), 0);
            stats.editDistanceHisto[editDistance] += 1;
        }
    }

    /*
    while (!atEnd(reader))
    {
        // Read alignment record.
        ++line;
        if (line % 100000 == 0) std::cerr << '.' << std::flush;
        
        if (readRecord(record, context, reader, tag) != 0)
        {
            std::cerr << "Could not read alignment! lineNo=" << line << std::endl;
            return 1;
        }
        if (options.verbosity >= 3)
            write2(std::cerr, record, context, Sam());
        
        // Update mapping from SAM rID to index of sequence in reference sequences.
        if (oldNumRefs != length(refNames))
        {
            for (unsigned i = oldNumRefs; i < length(refNames); ++i)
            {
                unsigned idx = 0;
                if (!getIdByName(seqIds, refNames[i], idx, seqIdsCache))
                {
                    std::cerr << "Invalid reference name " << refNames[i] << " from SAM/BAM record!\n";
                    return 1;
                }
                appendValue(rIdToSeqId, idx);
            }
            oldNumRefs = length(refNames);
        }
        */

            /*
        // Now, do something with it ;)
        unsigned editDistance = 0;
        // stats.numRecords += 1;  // One more record.
        stats.alignedRecords += !hasFlagUnmapped(record);
        
        // Compute alignment.
        if (record.rID != -1 && record.rID < static_cast<int>(length(seqs)))
        {
            // retrieve (possibly) missing read sequence
            int readId;
            if (getIdByName(readNames, record.qName, readId, readNameCache))
            {
                record.seq = readSeqs[readId];
                if (hasFlagRC(record))
                    reverseComplement(record.seq);
            }
            else
            {
                // This should never happen with a gold standard.
                if (options.goldStandard)
                {
                    std::cerr << "WARNING: A read unknown to the gold standard has been found: " << record.qName << std::endl;
                    return 1;
                }
                
                readId = -1;
            }
 
            unsigned seqId = rIdToSeqId[record.rID];
            
            // realign in a window 10bp left and right of the original alignment
            if (options.realign)
            {
                editDistance = realignBamRecord(seqs[seqId], record, 10);
            }
            else
            {
                // convert soft clipping into match/mismatches
                if (record.cigar[0].operation == 'S')
                    record.beginPos -= record.cigar[0].count;

                for (unsigned i = 0; i < length(record.cigar); ++i)
                {
                    if (record.cigar[i].operation == 'S')
                        record.cigar[i].operation = 'M';
                }
                
                bamRecordToAlignment(align, seqs[seqId], record);
            }
                
            if (options.verbosity >= 3)
                std::cerr << align << std::endl;
            
            if (!options.realign)
            {
                typedef Align<Dna5String>             TAlign;
                typedef typename Row<TAlign>::Type    TRow;
                typedef typename Iterator<TRow>::Type TRowIter;
                unsigned posReadFwd = 0;
                for (TRowIter it0 = iter(row(align, 0), 0), it1 = iter(row(align, 1), 0); !atEnd(it0); goNext(it0), goNext(it1))
                {
                    unsigned posRead = posReadFwd;
                    // is read aligned to reverse strand?
                    if (hasFlagRC(record))
                        posRead = (length(record.seq)) - posReadFwd;
                    SEQAN_ASSERT_LEQ(posRead, length(record.seq));

                    if (isGap(it0) && isGap(it1))
                        continue;
                    if (isGap(it0) || isGap(it1))
                    {
                        if (isGap(it0))
                        {
                            unsigned len = length(stats.insertHisto);
                            resize(stats.insertHisto, std::max(len, posRead + 1), 0);
                            stats.insertHisto[posRead] += 1;
                            posReadFwd += 1;
                        }
                        else
                        {
                            unsigned len = length(stats.deletionHisto);
                            resize(stats.deletionHisto, std::max(len, posRead + 1), 0);
                            stats.deletionHisto[posRead] += 1;
                        }
                        editDistance += 1;
                        continue;
                    }
                    if (*it0 != *it1)
                    {
                        unsigned len = length(stats.mismatchHisto);
                        resize(stats.mismatchHisto, std::max(len, posRead + 1), 0);
                        stats.mismatchHisto[posRead] += 1;
                        editDistance += 1;
                    }
                    posReadFwd += 1;
                }
                
                resize(qualSum, std::max(length(qualSum), length(record.qual)), 0);
                if (hasFlagRC(record))
                {
                    // read aligns to forward strand
                    for (unsigned i = 0; i < length(record.qual); ++i)
                        qualSum[i] += record.qual[i] - '!';
                } else
                {
                    // read aligns to reverse strand
                    for (unsigned i = 0; i < length(record.qual); ++i)
                        qualSum[(length(record.qual) - 1) - i] += record.qual[i] - '!';
                }
                ++reads;
            }

            if (options.useNM)
            {
                // trust NM instead of CIGAR string
                unsigned idx = 0;
                int nmValue =-1;
                BamTagsDict tagsDict(record.tags);
                if (findTagKey(idx, tagsDict, "NM") && extractTagValue(nmValue, tagsDict, idx) && nmValue != -1)
                    editDistance = nmValue;
                else
                    std::cerr << "WARNING: Could find NM tag in match of read " << record.qName << std::endl;
            }
            
            
//            if (editDistance >= 10)
//            {
//                Align<Dna5String> realign;
//                unsigned realignEditDist = realignBamRecord(realign, seqs[seqId], record, 10);
//                if (editDistance != realignEditDist)
//                {
//                    std::cout << "Realignment difference " << record.qName << " (old=" << editDistance << ", new=" << realignEditDist << ")\n";
//                    std::cout << "improvement=" << (editDistance - realignEditDist) << '\t' << record.qName << '\t' << record.beginPos << '\n';
//                    std::cout << align;
//                    std::cout << realign;
//                    std::cout << std::endl;
//                }
//            }

            // Update the maps minErrors, matchesCount, and minErrorsOrigin.  These will be used to print the results.
            if (options.goldStandard)
            {
                // Check if match corresponds to the true origin in gold standard.
                if ((record.beginPos >= goldStandard[readId].beginPos - 2 * goldStandard[readId].errors &&
                     record.beginPos <= goldStandard[readId].beginPos + 2 * goldStandard[readId].errors) &&
                    record.rID == goldStandard[readId].contigId &&
                    hasFlagRC(record) == goldStandard[readId].reverseComplemented)
                {
                    minErrorsOrigin[readId] = _min(minErrors[readId], editDistance);
                }
                
                minErrors[readId] = _min(minErrors[readId], editDistance);
                matchesCount[readId]++;
            }
            else
            {
                // update best-match errors
                if (readId != -1)
                {
                    minErrors[readId] = _min(minErrors[readId], editDistance);
                    matchesCount[readId]++;
                }
                else
                {
                    // add read name and entry
                    appendValue(readSeqs, record.seq);
                    if (hasFlagRC(record)) reverseComplement(back(readSeqs));
                    appendName(readNames, record.qName, readNameCache);
                    
                    // Update stats.
                    appendValue(minErrors, editDistance);
                    appendValue(matchesCount, 1);
                }                
            }

            
            // if (options.verbosity >= 2 && editDistance > 20)
            // {
            //     std::cerr << "EDIT DISTANCE == " << editDistance << "\n"
            //               << "ALIGNMENT\n" << align
            //               << "READ NAME\t" << record.qName << "\n";
            // }

            if (options.verbosity >= 3)
                std::cerr << "edit distance: " << editDistance << std::endl;
            unsigned len = length(stats.editDistanceHisto);
            resize(stats.editDistanceHisto, std::max(len, editDistance + 1), 0);
            stats.editDistanceHisto[editDistance] += 1;
        }
    }
*/

    // compute error probabilities
    resize(stats.avrgQuality, length(qualSum));
    for (unsigned i = 0; i < length(qualSum); ++i)
        stats.avrgQuality[i] = (double)qualSum[i] / (double)reads;

    // Resize histograms so they are of equal size.
    unsigned len = length(stats.mismatchHisto);
    len = std::max(len, (unsigned)length(stats.insertHisto));
    len = std::max(len, (unsigned)length(stats.deletionHisto));
    len = std::max(len, (unsigned)length(stats.mismatchHisto));
    len = std::max(len, (unsigned)length(stats.avrgQuality));
    resize(stats.insertHisto, len, 0);
    resize(stats.deletionHisto, len, 0);
    resize(stats.mismatchHisto, len, 0);
    resize(stats.avrgQuality, len, 0);

    // Print results.
    std::cout << "RESULTS\n\n";
    std::cout << "num records     \t" << stats.numRecords << std::endl;
    std::cout << "aligned records \t" << stats.alignedRecords << std::endl;
    std::cout << "aligned record %\t" << 100.0 * stats.alignedRecords / stats.numRecords << std::endl;
    std::cout << "distance histogram" << std::endl;
    std::cout << "#distance\tnumber of reads" << std::endl;
    for (unsigned i = 0; i < length(stats.editDistanceHisto); ++i)
        std::cout << i << '\t' << stats.editDistanceHisto[i] << std::endl;
    std::cout << "read position histogram" << std::endl;
    std::cout << "#position\tmismatches\tinsertions\tdeletions\terror prob\tquality-based error prob\tquality" << std::endl;
    for (unsigned i = 0; i < length(stats.mismatchHisto); ++i)
    {
        std::cout << i << '\t';
        std::cout << stats.mismatchHisto[i] << '\t';
        std::cout << stats.insertHisto[i] << '\t';
        std::cout << stats.deletionHisto[i] << '\t';
        std::cout << (stats.mismatchHisto[i] + stats.insertHisto[i] + stats.deletionHisto[i]) / (double)reads << '\t';
        std::cout << pow(10.0, stats.avrgQuality[i] / -10.0) << '\t';
        std::cout << stats.avrgQuality[i] << std::endl;
    }
    
    // write best-match error file
    if (!empty(options.bestMatchFile))
    {
        std::ofstream bestMatchFile(toCString(options.bestMatchFile));
        // sort output by read name
        TNameCache::TSet::iterator it = readNameCache.nameSet.begin();
        TNameCache::TSet::iterator itEnd = readNameCache.nameSet.end();
        for (; it != itEnd; ++it)
        {
            CharString readName = readNames[*it];
            // Do not output single-end records in paired-end mode.
            if (options.pairedEnd && readName[length(readName) - 2] != '/')
                continue;
            // Strip any trailing "/*" from read name.
            if (length(readName) > 2u && readName[length(readName) - 2] == '/')
                resize(readName, length(readName) - 2);
            bestMatchFile << readName << '\t';
            
            if (options.goldStandard)
            {
                bestMatchFile << (unsigned)goldStandard[*it].errors     << '\t' <<
                                 (unsigned)goldStandard[*it].seqErrors  << '\t' <<
                                 (unsigned)goldStandard[*it].SNPs       << '\t' <<
                                 (unsigned)goldStandard[*it].indels     << '\t';
                
                if (minErrorsOrigin[*it] != MaxValue<unsigned>::VALUE)
                    bestMatchFile << minErrorsOrigin[*it] << '\t';
                else
                    bestMatchFile << "None" << '\t';
            }
            
            if (minErrors[*it] != MaxValue<unsigned>::VALUE)
                bestMatchFile << minErrors[*it] << '\t';
            else
                bestMatchFile << "None" << '\t';
            
            bestMatchFile << matchesCount[*it] << std::endl;
        }
    }

    return 0;
}

int main(int argc, char const ** argv)
{
    using namespace seqan;

    // -----------------------------------------------------------------------
    // Handle Command Line
    // -----------------------------------------------------------------------

    // Setup command line parser.
    CommandLineParser parser;
    Options options;
    setupCommandLineParser(parser, options);
    
    // Then, parse the command line and handle the cases where help display
    // is requested or erroneous parameters were given.
    int res = parseCommandLineAndCheck(options, parser, argc, argv);
    if (res != 0)
    {
        std::cerr << "Error with command line arguments.\n";
        return 1;
    }
    if (options.showHelp || options.showVersion)
        return 0;

    // -----------------------------------------------------------------------
    // Guess format.
    // -----------------------------------------------------------------------

    if (options.inFormat == FORMAT_AUTO)
    {
        std::ifstream guessIn(toCString(options.inFile), std::ios::binary | std::ios::in);
        if (!guessIn.good())
        {
            std::cerr << "Could not open " << options.inFile << std::endl;
            return 1;
        }
        CharString magic;
        resize(magic, 2);
        guessIn.read(&magic[0], 2);
        if (magic != "\x1f\x8b")  // Is not gzip-compressed.
            options.inFormat = FORMAT_SAM;
        else
            options.inFormat = FORMAT_BAM;
    }

    // -----------------------------------------------------------------------
    // Do Work.
    // -----------------------------------------------------------------------

    // Read reference.
    if (options.verbosity >= 2)
        std::cerr << "Reading references from " << options.refFile << std::endl;
    StringSet<CharString> seqIds;
    StringSet<Dna5String> seqs;
    String<char, MMap<> > seqMMapString;
    if (!open(seqMMapString, toCString(options.refFile), OPEN_RDONLY))
    {
        std::cerr << "Could not open " << options.refFile << std::endl;
        return 1;
    }
    double start = sysTime();
    RecordReader<String<char, MMap<> >, DoublePass<StringReader> > refReader(seqMMapString);
    if (read2(seqIds, seqs, refReader, Fasta()) != 0)
    {
        std::cerr << "Could not read reference from " << options.refFile << std::endl;
        return 1;
    }
    std::cerr << "Loading reference took" << sysTime() - start << std::endl;

    // Open SAM/BAM file and do work.
    if (options.inFormat == FORMAT_SAM)
    {
        if (options.verbosity >= 2)
            std::cerr << "Opening SAM file " << options.inFile << std::endl;
        String<char, MMap<> > samMMapString;
        if (!open(samMMapString, toCString(options.inFile), OPEN_RDONLY))
        {
            std::cerr << "Could not open " << options.inFile << std::endl;
            return 1;
        }
        RecordReader<String<char, MMap<> >, SinglePass<StringReader> > samReader(samMMapString);
        
        if (options.goldStandard)
        {
            String<char, MMap<> > goldMMapString;
            if (!open(goldMMapString, toCString(options.goldStandardFile), OPEN_RDONLY))
            {
                std::cerr << "Could not open " << options.goldStandardFile << std::endl;
                return 1;
            }
            RecordReader<String<char, MMap<> >, SinglePass<StringReader> > goldReader(goldMMapString);
            return doWork(samReader, goldReader, seqIds, seqs, options, Sam());
        }
        else
        {
            return doWork(samReader, samReader, seqIds, seqs, options, Sam());
        }
    }
    else  // options.inFormat == FORMAT_BAM
    {
        if (options.verbosity >= 2)
            std::cerr << "Opening BAM file " << options.inFile << std::endl;
        Stream<Bgzf> bamStream;
        if (!open(bamStream, toCString(options.inFile), "r"))
        {
            std::cerr << "Could not open " << options.inFile << std::endl;
            return 1;
        }
        return doWork(bamStream, bamStream, seqIds, seqs, options, Bam());
    }

    return 0;
}

#else

int main(int, char const **)
{
    std::cerr << "bam_stats can only be compiled correctly with zlib." << std::endl;
    return 0;
}

#endif  // #if SEQAN_HAS_ZLIB
