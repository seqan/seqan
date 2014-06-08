// USAGE: compute_gain GENOME.fa PRE.sam POST.sam
//
// Where PRE.sam are the mapped reads before mapping, POST.sam are the mapped reads after mapping, and GENOME.fa is the
// used FASTA file with the genome.

#include <cctype>
#include <iostream>
#include <map>

#include <seqan/basic.h>
#include <seqan/align.h>
#include <seqan/graph_align.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/seq_io.h>
#include <seqan/bam_io.h>
#include <seqan/arg_parse.h>
#include <seqan/parallel.h>
#include <seqan/random.h>

// Data structure for options.
struct Options
{
    // Verbosity: 0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    int verbosity;

    // Number of threads to use.
    int numThreads;

    // Number of elements in one chunk
    unsigned chunkSize;

    // Path to genome.
    seqan::CharString pathGenome;
    // Path to pre-correction SAM file.
    seqan::CharString pathSamPreCorrection;
    // Path to post-correction SAM file.
    seqan::CharString pathSamPostCorrection;
    // Path to post-correction FASTA/FASTQ file.
    seqan::CharString pathFastaFastqPostCorrection;

    // Path to list with fixed and introduced errors.
    seqan::CharString pathCorrectionLog;

    // Whether to log all alignments.
    bool logAll;

    // Minimal number of unclipped bases that a read must have such that it is not ignored.  0 to disable clipping.
    int minUnclippedBases;

    // Maximal error rate in initial alignment such that the read is not ignored.  0 to disable ignoring.
    int maxErrorRate;

    // Maximal number of errors in initial alignment.  Set to -1 to disable.  See also maxErrorRate.
    int maxErrorCount;

    // Whether or not to allow indels.
    bool indels;

    // Additional characters to take out of genome left and right of alignment, in percent of the original read length.
    int padding;

    // Bandwidth to use for alignment.
    int bandwidth;

    // Check sorting.
    bool checkSorting;

    // Read only the first maxChunks chunks (default, maxValue<__uint64>())
    __uint64 maxChunks;

    Options() :
            verbosity(1), numThreads(1), chunkSize(0), logAll(false), minUnclippedBases(0), maxErrorRate(0),
            maxErrorCount(0), indels(true), padding(5), bandwidth(2 * padding), checkSorting(true)
    {}
};

// Less-than string comparison from samtools.

static inline int strnum_cmp(const char *a, const char *b)
{
    char *pa, *pb;
    pa = (char*)a; pb = (char*)b;
    while (*pa && *pb) {
        if (isdigit(*pa) && isdigit(*pb)) {
            long ai, bi;
            ai = strtol(pa, &pa, 10);
            bi = strtol(pb, &pb, 10);
            if (ai != bi) return ai<bi? -1 : ai>bi? 1 : 0;
        } else {
            if (*pa != *pb) break;
            ++pa; ++pb;
        }
    }
    if (*pa == *pb)
        return (pa-a) < (pb-b)? -1 : (pa-a) > (pb-b)? 1 : 0;
    return *pa<*pb? -1 : *pa>*pb? 1 : 0;
}

// Data structure to collect statistics in.

struct Stats
{
    // Number of unmapped records before, after, and both before and after correction.
    __uint64 numUnmappedPre;
    __uint64 numUnmappedPost;
    __uint64 numIgnoredPre;

    // Number of alignment pairs that were evaluated.
    __uint64 numTotal;

    // Number of bases aligned pre-correction/post-correction.
    __uint64 numBasesPre;
    __uint64 numBasesPost;
    // Number of errors pre-correction/post-correction.
    __uint64 numErrorsPre;
    __uint64 numErrorsPost;
    // Number of errorneous reads pre-correction/post-correction.
    __uint64 numErrorReadsPre;
    __uint64 numErrorReadsPost;
    __uint64 numReads;

    // Histogram of differences (before/after correction).
    std::map<int, unsigned> histo;

    // Histogram of read error rate distribution in percent.
    std::map<int, unsigned> preErrorHisto;

    // Histogram of errors at a certain position.
    typedef seqan::Tuple<unsigned, 4> TPositionalCounts;
    TPositionalCounts zeroCounts;
    seqan::String<TPositionalCounts> preErrorsAtPos;  // 0..mismatch, 1..insert, 2..deletion in the read compared to reference
    seqan::String<TPositionalCounts> postErrorsAtPos; // 3..total number of reads having a base at that position (not cropped)

    // We can later compute the gain from actualErrorSum / diffErrorSum.
    //
    // Summed up edit distance errors in the pre-correction alignment.
    __int64 actualErrorSum;
    // Summed up difference between pre-correction and post-correction errors.
    __int64 diffErrorSum;

    // Number of introduced and removed errors.
    __uint64 numErrorsIntroduced;
    __uint64 numErrorsRemoved;

    // Reusable data structures
    seqan::BamAlignmentRecord _recordPre;
    seqan::BamAlignmentRecord _recordPost;
    seqan::String<seqan::BamAlignmentRecord> _chunkPre;
    seqan::String<seqan::Dna5String> _chunkPost;

    seqan::Dna5String _genomeInfix;
    seqan::Dna5String _preRead;
    seqan::Dna5String _postRead;

    typedef seqan::Align<seqan::Dna5String> TAlign;
    TAlign _preAlign;
    TAlign _postAlign;

    Stats() : numUnmappedPre(0), numUnmappedPost(0), numIgnoredPre(0), numTotal(0),
              numBasesPre(0), numBasesPost(0), numErrorsPre(0), numErrorsPost(0), numErrorReadsPre(0),
              numErrorReadsPost(0), numReads(0), actualErrorSum(0), diffErrorSum(0), numErrorsIntroduced(0),
              numErrorsRemoved(0)
    {
        zeroCounts[0] = 0;
        zeroCounts[1] = 0;
        zeroCounts[2] = 0;
        zeroCounts[3] = 0;
    }
};

typedef seqan::StringSet<seqan::CharString> TNameStore;
typedef seqan::NameStoreCache<TNameStore>   TNameStoreCache;

// Allows the trimming of strings after the first whitespace.

void trimSeqHeaderToId(seqan::CharString & header)
{
    unsigned i = 0;
    for (; i < length(header); ++i)
        if (isspace(header[i]))
          break;
    resize(header, i);
}

// Returns "yes"/"no" for a bool.

char const * getYesNo(bool b)
{
    return b ? "YES" : "NO";
}

// Computes gain between two SAM records and updates Stats object.

void updateStats(Stats & stats,
                 std::ofstream & correctionLog,
                 seqan::BamAlignmentRecord const & preRecord,
                 seqan::Dna5String & postRead,  // rc-ed to same as preRecord.
                 seqan::String<unsigned> const & idMap,
                 seqan::StringSet<seqan::Dna5String> const & seqs,
                 Options const & options)
{
    using namespace seqan;

    stats.numUnmappedPre += hasFlagUnmapped(preRecord);

    if (hasFlagUnmapped(preRecord))
        return;  // Ignore.
    stats.numTotal += 1;
    stats.numReads += 1;

    // In the following, we compute scores with globalAlignment() using edit distance and convert them into distances by
    // computing -score.

    unsigned seqIdx = idMap[preRecord.rID];

    // Get additional shifts in genomic position if there is soft clipping.
    int beginShift = 0, endShift = 0;
    if (!empty(preRecord.cigar) && front(preRecord.cigar).operation == 'S')
        beginShift = front(preRecord.cigar).count;
    if (length(preRecord.cigar) > 1u && back(preRecord.cigar).operation == 'S')
        endShift = back(preRecord.cigar).count;
    if (options.minUnclippedBases != 0 && ((int)length(preRecord.seq) - (beginShift + endShift) < options.minUnclippedBases))
    {
        stats.numIgnoredPre += 1;
        return;  // Ignore, too few unclipped bases.
    }

    // Compute begin and positions of genome infix, with padding, but not going over the end of the genome.
    unsigned beginPos = preRecord.beginPos;
    unsigned endPos = beginPos + length(preRecord.seq);

    Dna5String &preRead = stats._preRead;
    preRead = preRecord.seq;

    // speed-up heuristic
    bool skipAlignPre = (infix(seqs[seqIdx], beginPos, endPos) == preRead);
    bool skipAlignPost = (infix(seqs[seqIdx], beginPos, endPos) == postRead);

    int padding = (options.padding * length(preRecord.seq) + 99) / 100;
    if (options.indels)
    {
        if ((int)beginPos > padding + beginShift)
            beginPos -= (padding + beginShift);
        else
            beginPos = 0;
        endPos = preRecord.beginPos + getAlignmentLengthInRef(preRecord);
        endPos += (padding + endShift);
        if (endPos > length(seqs[seqIdx]))
            endPos = length(seqs[seqIdx]);
    }

    // Get genome infix, and convert read seqs (CharString) into Dna5Strings.
    Dna5String &genomeInfix = stats._genomeInfix;
    genomeInfix = infix(seqs[seqIdx], beginPos, endPos);

    if (hasFlagRC(preRecord))
    {
        // bring sequences to original read direction
        reverseComplement(genomeInfix);
        reverseComplement(preRead);
        reverseComplement(postRead);
    }
//    else
//    {
//        std::cout << preRead << std::endl;
//        std::cout << postRead << std::endl;
//    }

    // Create alignment objects and AlignConfig object.
    typedef Align<Dna5String> TAlign;
    TAlign &preAlign = stats._preAlign;
    TAlign &postAlign = stats._postAlign;
    
    resize(rows(preAlign), 2);
    setSource(row(preAlign, 0), genomeInfix);
    setSource(row(preAlign, 1), preRead);
    resize(rows(postAlign), 2);
    setSource(row(postAlign, 0), genomeInfix);
    setSource(row(postAlign, 1), postRead);

    // Compute distance before and after correction.
    int diffPre = 0;
    int diffPost = 0;

    if (options.indels)
    {
        // In the case of indels, we have to perform an edit distance alignment.
        AlignConfig<true, false, false, true> alignConfig;
        Score<int, Simple> scoringScheme(0, -1000, -1001);

        /////////
        // PRE //
        /////////

        // Align pre record to genome.
        if (!skipAlignPre)
            diffPre = -globalAlignment(preAlign, scoringScheme, alignConfig, NeedlemanWunsch()) / 1000;

        int diffPreRate = (int)ceil(100.0 * diffPre / length(preRead));
        if ((options.maxErrorRate != 0 && diffPreRate > options.maxErrorRate) ||
            (options.maxErrorCount != -1 && diffPre > options.maxErrorCount))
        {
            stats.numIgnoredPre += 1;
            return;  // Ignore alignment with too many errors.
        }
        stats.preErrorHisto[diffPreRate] += 1;

        // update positional error counts PRE
        resize(stats.preErrorsAtPos, std::max(length(stats.preErrorsAtPos), length(preRead)), stats.zeroCounts);
        if (diffPre != 0)
        {
            typedef Row<TAlign>::Type TRow;
            typedef Iterator<TRow, Standard>::Type TGapsIter;

            int viewBegin = toViewPosition(row(preAlign, 1), 0);
            int viewEnd = toViewPosition(row(preAlign, 1), length(preRead) - 1) + 1;

            TGapsIter genomeIter = begin(row(preAlign, 0), Standard()) + viewBegin;
            TGapsIter readIter = begin(row(preAlign, 1), Standard()) + viewBegin;

            unsigned readPos = 0;
            for (int i = viewBegin; i < viewEnd; ++i, goNext(genomeIter), goNext(readIter))
            {
                // errro types: 0..mismatch, 1..insert, 2..deletion in the read compared to reference
                if (isGap(genomeIter))
                {
                    ++stats.preErrorsAtPos[readPos][1];
                }
                else if (isGap(readIter))
                {
                    ++stats.preErrorsAtPos[readPos][2];
                    continue; // don't increment readPos here
                }
                else if (*genomeIter != *readIter)
                {
                    ++stats.preErrorsAtPos[readPos][0];
                }
                // increment number of total bases
                ++stats.preErrorsAtPos[readPos][3];
                ++readPos;
            }
            SEQAN_CHECK(atEnd(genomeIter) == atEnd(readIter), "Invalid pairwise alignment!");
        }
        else
        {
            // increment at least number of total bases (in case of no error)
            for (unsigned i = 0; i < length(preRead); ++i)
                ++stats.preErrorsAtPos[i][3];
        }

        //////////
        // POST //
        //////////

        // Align post record to genome.
        if (!skipAlignPost)
            diffPost = -globalAlignment(postAlign, scoringScheme, alignConfig, NeedlemanWunsch()) / 1000;

        if (options.verbosity >= 3 || (options.verbosity >= 2 && (abs(diffPre) > 10 || abs(diffPost) > 10)))
            std::cerr << "RECORD: " << preRecord.qName << "\n"
                      << "BEFORE, score == " << diffPre << "\n" << preAlign << '\n'
                      << "AFTER, score == " << diffPost << "\n" << postAlign << '\n'
                      << " --> BEFORE - AFTER == " << diffPre - diffPost << '\n';

        if (correctionLog.good() && (diffPre != diffPost || options.logAll))
        {
            correctionLog << "RECORD: " << preRecord.qName << "\n"
                          << "\n"
                          << "BEFORE score = " << diffPre << "\n"
                          << preAlign << "\n"
                          << "\n"
                          << "AFTER score = " << diffPost << "\n"
                          << postAlign << "\n"
                          << "------------------------------------------------------------------------------\n";
        }

        // update positional error counts POST
        resize(stats.postErrorsAtPos, std::max(length(stats.postErrorsAtPos), length(postRead)), stats.zeroCounts);
        if (diffPost != 0)
        {
            typedef Row<TAlign>::Type TRow;
            typedef Iterator<TRow, Standard>::Type TGapsIter;

            int viewBegin = toViewPosition(row(postAlign, 1), 0);
            int viewEnd = toViewPosition(row(postAlign, 1), length(postRead) - 1) + 1;

            TGapsIter genomeIter = begin(row(postAlign, 0), Standard()) + viewBegin;
            TGapsIter readIter = begin(row(postAlign, 1), Standard()) + viewBegin;

            unsigned readPos = 0, errors = 0;
            for (int i = viewBegin; i < viewEnd; ++i, goNext(genomeIter), goNext(readIter))
            {
                // errro types: 0..mismatch, 1..insert, 2..deletion in the read compared to reference
                if (isGap(genomeIter))
                {
                    ++stats.postErrorsAtPos[readPos][1];
                    ++errors;
                }
                else if (isGap(readIter))
                {
                    ++stats.postErrorsAtPos[readPos][2];
                    ++errors;
                    continue; // don't increment readPos here
                }
                else if (*genomeIter != *readIter)
                {
                    ++stats.postErrorsAtPos[readPos][0];
                    ++errors;
                }
                // increment number of total bases
                ++stats.postErrorsAtPos[readPos][3];
                ++readPos;
            }
            SEQAN_CHECK(atEnd(genomeIter) == atEnd(readIter), "Invalid pairwise alignment!");
        }
        else
        {
            // increment at least number of total bases (in case of no error)
            for (unsigned i = 0; i < length(postRead); ++i)
                ++stats.postErrorsAtPos[i][3];
        }
    }
    else
    {
        // In the case of Hamming distance, we can simply count matching and mismatching bases.

        if (options.verbosity >= 3)
            std::cerr << "infix == " << genomeInfix << "\tpre == " << preRead << "\tpost == " << postRead << '\n';
        //SEQAN_CHECK(length(genomeInfix) >= length(preRead), "Must have geq length with Hamming distance.");
        //SEQAN_CHECK(length(genomeInfix) >= length(postRead), "Must have geq length with Hamming distance.");
        //SEQAN_CHECK(length(preRead) >= length(postRead), "Must have the same length with Hamming distance.");

        // For Hamming distance, we compare at most m characters where m is the smaller one of the
        // lengths of postRead and preRead.
        unsigned minLen = std::min(length(postRead), length(preRead));
        std::string flags;

        resize(stats.preErrorsAtPos, std::max(length(stats.preErrorsAtPos), length(preRead)), stats.zeroCounts);
        resize(stats.postErrorsAtPos, std::max(length(stats.postErrorsAtPos), length(postRead)), stats.zeroCounts);

        unsigned numIntroduced = 0;
        unsigned numRemoved = 0;
        for (unsigned i = 0; i < minLen; ++i)
        {
            bool badPre = (genomeInfix[i] != preRead[i]);
            bool badPost = (genomeInfix[i] != postRead[i]);
            if (badPre && !badPost)
            {
                flags.push_back('-');
                numRemoved += 1;
            }
            else if (!badPre && badPost)
            {
                flags.push_back('+');
                numIntroduced += 1;
            }
            else
            {
                flags.push_back(' ');
            }
            stats.preErrorsAtPos[i][0] += badPre;
            stats.postErrorsAtPos[i][0] += badPost;
            ++stats.preErrorsAtPos[i][3];
            ++stats.postErrorsAtPos[i][3];

            diffPre += badPre;
            diffPost += badPost;
        }

        // If preRead is longer than postRead then matches between pre read and genome are counted as introduced
        // errors.  Mismatches are counted as removed errors.
        for (unsigned i = minLen; i < length(preRead); ++i)
        {
            if (genomeInfix[i] == preRead[i])
            {
                diffPost += 1;
                numIntroduced += 1;
            }
            else
            {
                diffPre += 1;
                numRemoved += 1;
                ++stats.preErrorsAtPos[i][0];
            }
            ++stats.preErrorsAtPos[i][3];
        }
        // If postRead is longer than preRead then matches between the post read are counted as removed errors
        // and mismatches are counted as removed errors.  If the read reaches over the genome then this is counted
        // as mismatches.
        for (unsigned i = minLen; i < length(postRead); ++i)
        {
            if (beginPos + i > length(seqs[seqIdx]))
            {
                diffPost += length(postRead) - i;
                numIntroduced += length(postRead) - i;
                break;  // End of genome.
            }
            if (genomeInfix[i] == postRead[i])
            {
                diffPre += 1;
                numRemoved += 1;
            }
            else
            {
                diffPost += 1;
                numIntroduced += 1;
                ++stats.postErrorsAtPos[i][0];
            }
            ++stats.postErrorsAtPos[i][3];
        }

        // Update histogram of error
        int diffPreRate = (int)ceil(100.0 * diffPre / length(preRead));
        stats.preErrorHisto[diffPreRate] += 1;

        if (options.verbosity >= 3 || (options.verbosity >= 2 && (abs(diffPre) > 10 || abs(diffPost) > 10)))
            std::cerr << "RECORD: " << preRecord.qName << "\n"
                      << "GENOME " << genomeInfix << "\n"
                      << "BEFORE " << preRead << "\tscore == " << diffPre << '\n'
                      << "AFTER  " << postRead << "\tscore == " << diffPost << '\n'
                      << "REMOVED    " << numRemoved << "\n"
                      << "INTRODUCED " << numIntroduced << "\n"
                      << " --> BEFORE - AFTER == " << diffPre - diffPost << '\n';

        if (correctionLog.good() && diffPre != diffPost)
        {
            correctionLog << "RECORD: " << preRecord.qName << "\n"
                          << "INTRODUCED " << numIntroduced << "\tREMOVED\t" << numRemoved << "\n"
                          << "\n"
                          << "BEFORE\n"
                          << genomeInfix << "\n";
            for (unsigned i = 0; i < minLen; ++i)
                correctionLog << ((genomeInfix[i] == preRead[i]) ? '|' : ' ');
            correctionLog << "\n"
                          << preRead << "\n\n"
                          << "AFTER\n"
                          << genomeInfix << "\n";
            for (unsigned i = 0; i < minLen; ++i)
                correctionLog << ((genomeInfix[i] == postRead[i]) ? '|' : ' ');
            correctionLog << "\n"
                          << postRead << "\n"
                          << flags << "\n"
                          << "==============================================================================\n";
        }

        stats.numErrorsRemoved += numRemoved;
        stats.numErrorsIntroduced += numIntroduced;
    }

    stats.numBasesPre += length(preRead);
    stats.numBasesPost += length(postRead);
    stats.numErrorsPre += diffPre;
    stats.numErrorsPost += diffPost;
    stats.numErrorReadsPre += (diffPre != 0);
    stats.numErrorReadsPost += (diffPost != 0);

    SEQAN_CHECK(diffPre >= 0, "Edit distance must be >= 0!");
    SEQAN_CHECK(diffPost >= 0, "Edit distance must be >= 0!");

    stats.actualErrorSum += diffPre;
    stats.diffErrorSum += diffPre - diffPost;
    stats.histo[diffPre - diffPost] += 1;
}

// Read header and perform basic validity checks.

template <typename TRecordReader>
int readHeader(seqan::BamHeader & header,
               seqan::BamIOContext<TNameStore> & context,
               TRecordReader & reader)
{
    using namespace seqan;

    int res = 0;
    res = readRecord(header, context, reader, Sam());
    if (res != 0)
    {
        std::cerr << "Could not read SAM header!\n";
        return 1;
    }

    // Check that the SAM file is sorted by QNAME.
    CharString sortOrder;
    for (unsigned i = 0; i < length(header.records); ++i)
    {
        if (header.records[i].type != BAM_HEADER_FIRST)
            continue;
        unsigned idx = 0;
        if (findTagKey(idx, "SO", header.records[i]))
            sortOrder = header.records[i].tags[idx].i2;
    }

    return 0;
}

seqan::ArgumentParser::ParseResult
parseCommandLine(Options & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("compute_gain");

    setShortDescription(parser, "Compute read correction metric GAIN.");
    setVersion(parser, "0.2");
    setDate(parser, "August 2012");

    addUsageLine(parser,
                 "[\\fIOPTIONS\\fP] \\fB-g\\fP GENOME.fa \\fB--pre\\fP \\fIPRE.sam\\fP \\fB--post\\fP "
                 "\\fIPOST.sam\\fP");

    addDescription(parser,
                   "This program computes the read correction tool metric GAIN.  It takes a genome FASTA file and "
                   "two SAM file with read alignments before and after correction.  It then computes various "
                   "statistics and computes the GAIN, based on the edit or Hamming distance for reach read before "
                   "and after correction.");

    addOption(parser, seqan::ArgParseOption("q", "quiet", "Disable most output."));
    addOption(parser, seqan::ArgParseOption("v", "verbose", "Enable more verbose output."));
    addOption(parser, seqan::ArgParseOption("vv", "very-verbose", "Enable even more verbose output."));

#ifdef _OPENMP
    addOption(parser, seqan::ArgParseOption("nt", "num-threads", "Number of threads to use.", seqan::ArgParseOption::INTEGER, "THREADS"));
    setDefaultValue(parser, "num-threads", options.numThreads);
    setMinValue(parser, "num-threads", "1");
#endif
    addOption(parser, seqan::ArgParseOption("", "chunk-size", "Chunk size.", seqan::ArgParseOption::INTEGER, "THREADS"));
    setDefaultValue(parser, "chunk-size", "10000");
    setMinValue(parser, "chunk-size", "100");

    addOption(parser, seqan::ArgParseOption("", "max-chunks", "Maximal number of chunks to read (0=disabled).", seqan::ArgParseOption::INTEGER));
    setMinValue(parser, "max-chunks", "0");
    setDefaultValue(parser, "max-chunks", "0");

    // Argument Section -- Comparison
    //
    addSection(parser, "Comparison");

    addOption(parser, seqan::ArgParseOption("", "padding", "Additional genome characters to use for alignment in percent of the origina read length.", seqan::ArgParseOption::INTEGER, "PADDING"));
    setDefaultValue(parser, "padding", 5);
    setMinValue(parser, "padding", "0");
    addOption(parser, seqan::ArgParseOption("", "bandwidth", "Bandwidth to use for alignment.", seqan::ArgParseOption::INTEGER, "BAND"));
    setDefaultValue(parser, "bandwidth", 10);
    setMinValue(parser, "bandwidth", "1");
    hideOption(parser, "bandwidth");
    addOption(parser, seqan::ArgParseOption("", "metric", "The metric type to use.", seqan::ArgParseOption::STRING,
                                            "METRIC"));
    setDefaultValue(parser, "metric", "edit");
    setValidValues(parser, "metric", "hamming edit");

    // Argument Section -- Filtration

    addSection(parser, "Input / Output");

    addOption(parser, seqan::ArgParseOption("", "min-unclipped-bases", "Reads with fewer unclipped bases are ignored. "
                                            "Set to 0 to disable ignoring because of this.",
                                            seqan::ArgParseOption::INTEGER, "NUM"));
    setDefaultValue(parser, "min-unclipped-bases", "0");
    setMinValue(parser, "min-unclipped-bases", "0");

    addOption(parser, seqan::ArgParseOption("", "max-error-rate", "Reads with a higher error rate in the initial mapping "
                                            "are ignored.  Given in percent.  Set to 0 to disable ignoring because of "
                                            "this.", seqan::ArgParseOption::INTEGER, "NUM"));
    setDefaultValue(parser, "max-error-rate", "0");
    setMinValue(parser, "max-error-rate", "0");

    addOption(parser, seqan::ArgParseOption("", "max-error-count", "Reads with a higher error count in the initial mapping "
                                            "are ignored.  Given in percent.  Set to -1 to disable ignoring because of "
                                            "this.", seqan::ArgParseOption::INTEGER, "NUM"));
    setDefaultValue(parser, "max-error-count", "-1");
    setMinValue(parser, "max-error-count", "-1");

    // Argument Section -- Input / Output

    addSection(parser, "Input / Output");

    addOption(parser, seqan::ArgParseOption("g", "genome", "Genome file.", seqan::ArgParseOption::INPUTFILE,
                                            "GENOME.fa"));
    setRequired(parser, "genome");
    setValidValues(parser, "genome", "FA FASTA");

    addOption(parser, seqan::ArgParseOption("", "pre", "Pre-correction SAM file.", seqan::ArgParseOption::INPUTFILE,
                                            "PRE.sam"));
    setRequired(parser, "pre");
    setValidValues(parser, "pre", "SAM");

    addOption(parser, seqan::ArgParseOption("", "post-sam", "Post-correction SAM file.", seqan::ArgParseOption::INPUTFILE,
                                            "POST.sam"));
    setValidValues(parser, "post-sam", "SAM");

    addOption(parser, seqan::ArgParseOption("", "post", "Post-correction FASTQ or FASTA file.", seqan::ArgParseOption::INPUTFILE,
                                            "POST.fq"));
    setValidValues(parser, "post", "FASTQ FQ FASTA FA");

    addOption(parser, seqan::ArgParseOption("", "correction-log", "Write log about introduced/removed errors to this file.",
                                            seqan::ArgParseOption::OUTPUTFILE, "OUT.txt"));
    addOption(parser, seqan::ArgParseOption("", "log-all", "Log all not only introduced/removed errors."));

    addOption(parser, seqan::ArgParseOption("", "no-check-sorting", "No checking for reads being sorted."));

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    // Extract option values.
    if (isSet(parser, "quiet"))
        options.verbosity = 0;
    else if (isSet(parser, "verbose"))
        options.verbosity = 2;
    else if (isSet(parser, "very-verbose"))
        options.verbosity = 3;
    if (!isSet(parser, "post-sam") && !isSet(parser, "post"))
    {
        std::cerr << "ERROR: Neither --post-sam nor --post was set!\n";
        return seqan::ArgumentParser::PARSE_ERROR;
    }
#ifdef _OPENMP
    getOptionValue(options.numThreads, parser, "num-threads");
#endif
    getOptionValue(options.chunkSize, parser, "chunk-size");
    getOptionValue(options.pathGenome, parser, "genome");
    getOptionValue(options.pathSamPreCorrection, parser, "pre");
    getOptionValue(options.pathSamPostCorrection, parser, "post-sam");
    getOptionValue(options.pathFastaFastqPostCorrection, parser, "post");
    getOptionValue(options.pathCorrectionLog, parser, "correction-log");
    getOptionValue(options.minUnclippedBases, parser, "min-unclipped-bases");
    getOptionValue(options.maxErrorRate, parser, "max-error-rate");
    getOptionValue(options.maxErrorCount, parser, "max-error-count");
    getOptionValue(options.maxChunks, parser, "max-chunks");
    options.logAll = isSet(parser, "log-all");
    seqan::CharString metricValue;
    getOptionValue(metricValue, parser, "metric");
    options.indels = (metricValue == "edit");

    options.checkSorting = !isSet(parser, "no-check-sorting");

    return seqan::ArgumentParser::PARSE_OK;
}

int main(int argc, char const ** argv)
{
    using namespace seqan;

    // -----------------------------------------------------------------------
    // Program Initialization
    // -----------------------------------------------------------------------

    // Checking command line parameters.
    Options options;
    seqan::ArgumentParser::ParseResult argParseRes = parseCommandLine(options, argc, argv);

    // If parsing was not successful then exit with code 1 if there were errors.
    // Otherwise, exit with code 0 (e.g. help was printed).
    if (argParseRes != seqan::ArgumentParser::PARSE_OK)
        return argParseRes == seqan::ArgumentParser::PARSE_ERROR;

    // Opening file and record reader.

    std::fstream inPre(toCString(options.pathSamPreCorrection), std::ios::in | std::ios::binary);
    if (!inPre.good())
    {
        std::cerr << "ERROR: Could not open pre-correction file.\n";
        return 1;
    }
    RecordReader<std::fstream, SinglePass<> > readerPre(inPre);
    std::fstream inPost;
    bool postSam = !empty(options.pathSamPostCorrection);  // Whether or not to read post SAM.
    if (postSam)
        inPost.open(toCString(options.pathSamPostCorrection), std::ios::in | std::ios::binary);
    else
        inPost.open(toCString(options.pathFastaFastqPostCorrection), std::ios::in | std::ios::binary);
    if (!inPost.good())
    {
        std::cerr << "ERROR: Could not open post-correction file.\n";
        return 1;
    }
    RecordReader<std::fstream, SinglePass<> > readerPost(inPost);

    // Read pre-correction SAM header.

    std::cerr << "SAM headers...\n";

    TNameStore refNameStore;
    TNameStoreCache refNameStoreCache(refNameStore);
    BamIOContext<TNameStore> context(refNameStore, refNameStoreCache);
    BamHeader header;
    int res = readHeader(header, context, readerPre);
    if (res != 0)
        return 1;  // Error message printed above.

    // Read post-correction SAM header and forget it immediately. We COULD check for consistency with pre-correction
    // header but this tool is for internal use only...

    if (postSam)
    {
        TNameStore refNameStore;
        TNameStoreCache refNameStoreCache(refNameStore);
        BamIOContext<TNameStore> context(refNameStore, refNameStoreCache);
        BamHeader header;
        int res = readHeader(header, context, readerPost);
        if (res != 0)
            return 1;  // Error message printed above.
    }

    // Read genome and compute mapping from SAM record reference ids to seqs index.

    std::cerr << "Read Genome...\n";

    StringSet<CharString> ids;
    StringSet<Dna5String> seqs;
    std::fstream inGenome(toCString(options.pathGenome), std::ios::in | std::ios::binary);
    RecordReader<std::fstream, SinglePass<> > readerGenome(inGenome);
    if (read2(ids, seqs, readerGenome, Fasta()) != 0)
    {
        std::cerr << "Could not read genome!\n";
        return 1;
    }

    String<unsigned> idMap;
    resize(idMap, length(refNameStore), maxValue<unsigned>());
    for (unsigned i = 0; i < length(ids); ++i)
    {
        trimSeqHeaderToId(ids[i]);
        unsigned idx = 0;
        if (!getIdByName(refNameStore, ids[i], idx, refNameStoreCache))
        {
            std::cerr << "Reference " << ids[i] << " not in SAM references!\n";
            return 1;
        }
        idMap[idx] = i;
    }

    // -----------------------------------------------------------------------
    // Open file for correction log.
    // -----------------------------------------------------------------------
    std::ofstream correctionLog;
    if (!empty(options.pathCorrectionLog))
        correctionLog.open(toCString(options.pathCorrectionLog));

    // -----------------------------------------------------------------------
    // Statistics Computation
    // -----------------------------------------------------------------------

    std::cerr << "Compute Statistics...\n";

    seqan::String<Stats> stats;
    resize(stats, options.numThreads);

    // Whether or not to break out of the loop below with or without an error, required in this way because we return
    // from OpenMP block.
    bool stop = false;
    bool error = false;

    // Format recognition tag used for sequence file I/O.
    seqan::AutoSeqStreamFormat seqFormatTag;

    __uint64 chunksLeftToRead = options.maxChunks;
    --chunksLeftToRead;

    // to reduce the number of threads waiting in front of the critical section
    Rng<MersenneTwister> rng(42);
    Pdf<Uniform<double> > chunkSizeNoise(options.chunkSize, 2 * options.chunkSize);

    SEQAN_OMP_PRAGMA(parallel num_threads(options.numThreads))
    while (!stop && !error)
    {
        int const tid = omp_get_thread_num();
        BamAlignmentRecord &recordPre = stats[tid]._recordPre;
        BamAlignmentRecord &recordPost = stats[tid]._recordPost;
        seqan::String<BamAlignmentRecord> &chunkPre = stats[tid]._chunkPre;
        seqan::String<seqan::Dna5String> &chunkPost = stats[tid]._chunkPost;

        clear(chunkPre);
        clear(chunkPost);

        // Read next chunk atomically.
        //
        // The stop/error variables are only updated in the critical section which implies a flush on all variables
        // anyway.  Thus there is no need to flush manually again.
        SEQAN_OMP_PRAGMA(critical (read_chunk))
        {
            int const tid = omp_get_thread_num();
            unsigned myChunkSize = (unsigned)pickRandomNumber(rng, chunkSizeNoise);
            seqan::CharString prevName;
            seqan::CharString postId;
            clear(recordPre.qName);
            clear(recordPost.qName);

            while (!stop && !error && !atEnd(readerPre) && !atEnd(readerPost))
            {
                // Read next record into chunk.
                if (readRecord(recordPre, context, readerPre, Sam()) != 0)
                {
                    std::cerr << "ERROR: Problem reading from pre correction file. (" << omp_get_thread_num() << ")\n";
                    stop = error = true;  // Stop processing for all.
                    break;
                }
                if (hasFlagSecondary(recordPre))
                    continue;  // Skip, this happens for bwasw input.

                // check sorting
                if (options.checkSorting && !empty(prevName) && strnum_cmp(toCString(prevName), toCString(recordPre.qName)) >= 0)
                {
                    std::cerr << "ERROR: Expected sorted by quername but was.\n"
                              << "But was: " << recordPre.qName << " >= " << prevName << "\n";
                    stop = error = false;
                }
                prevName = recordPre.qName;

                
                // read post-records as long as they are less than the last pre-record
                while (!atEnd(readerPost) && (empty(recordPost.qName) || strnum_cmp(toCString(recordPre.qName), toCString(recordPost.qName)) > 0))
                {
                    if (postSam)
                    {
                        if (readRecord(recordPost, context, readerPost, Sam()) != 0)
                        {
                            std::cerr << "ERROR: Problem reading from post correction file. (" << omp_get_thread_num() << ")\n";
                            stop = error = true;  // Stop processing for all.
                            break;
                        }
                        stats[tid].numUnmappedPost += hasFlagUnmapped(recordPost);
                    }
                    else
                    {
                        if (readRecord(recordPost.qName, recordPost.seq, readerPost, seqFormatTag) != 0)
                        {
                            std::cerr << "ERROR: Problem reading from post correction file. (" << omp_get_thread_num() << ")\n";
                            stop = error = true;  // Stop processing for all.
                            break;
                        }
                        trimSeqHeaderToId(recordPost.qName);
                    }
                    stats[tid].numUnmappedPre++;
                }

                // found matching pair of qnames?
                if (recordPre.qName == recordPost.qName)
                {
                    appendValue(chunkPre, recordPre);
                    appendValue(chunkPost, recordPost.seq);
                    if (postSam)
                    {
                        if (hasFlagRC(recordPost) != hasFlagRC(recordPre))
                            reverseComplement(back(chunkPost)); // transform post-read to the same orientation as pre-read
                    }
                    else
                    {
                        if (hasFlagRC(recordPre))
                            reverseComplement(back(chunkPost)); // transform post-read to the same orientation as pre-read
                    }
                    stats[tid].numUnmappedPre--;    // we overcounted unmapped pre-records by the one that matches the post-record
                }

                // Break if chunk full.
                if (length(chunkPre) > myChunkSize)
                    break;
            }

            if (atEnd(readerPre) || atEnd(readerPost) || chunksLeftToRead == 0ul)
                stop = true;  // Do not read any more, this thread finishes its computation.
            else
                --chunksLeftToRead;
        }

        // TODO(holtgrew): Defer RC for post chunk until here with array of flags.

        // Process chunk.
        for (unsigned i = 0; i < length(chunkPre); ++i)
            updateStats(stats[tid], correctionLog, chunkPre[i], chunkPost[i], idMap, seqs, options);
    }

    if (!atEnd(readerPre) || !atEnd(readerPost))
        std::cerr << "WARNING: Files not read completely!\n";
    if (error)
    {
        std::cerr << "An error occured. Bailing out.\n";
        return 1;
    }

    // Count remaining records, are unmapped/nonpresent in other file.

    if (chunksLeftToRead != 0ul)
    {
        SEQAN_CHECK(atEnd(readerPre) && atEnd(readerPost), "Both readers must be at end!");
    }

    // -----------------------------------------------------------------------
    // Write Output
    // -----------------------------------------------------------------------

    // Compute global stats.
    Stats globalStats;
    for (unsigned i = 0; i < length(stats); ++i)
    {
        globalStats.numUnmappedPre += stats[i].numUnmappedPre;
        globalStats.numIgnoredPre += stats[i].numIgnoredPre;
        globalStats.numTotal += stats[i].numTotal;
        globalStats.numBasesPre += stats[i].numBasesPre;
        globalStats.numBasesPost += stats[i].numBasesPost;
        globalStats.numErrorsPre += stats[i].numErrorsPre;
        globalStats.numErrorsPost += stats[i].numErrorsPost;
        globalStats.numErrorReadsPre += stats[i].numErrorReadsPre;
        globalStats.numErrorReadsPost += stats[i].numErrorReadsPost;
        globalStats.numReads += stats[i].numReads;
        globalStats.actualErrorSum += stats[i].actualErrorSum;
        globalStats.diffErrorSum += stats[i].diffErrorSum;
        globalStats.numErrorsIntroduced += stats[i].numErrorsIntroduced;
        globalStats.numErrorsRemoved += stats[i].numErrorsRemoved;

        for (std::map<int, unsigned>::const_iterator it = stats[i].histo.begin(); it != stats[i].histo.end(); ++it)
            globalStats.histo[it->first] += it->second;
        for (std::map<int, unsigned>::const_iterator it = stats[i].preErrorHisto.begin(); it != stats[i].preErrorHisto.end(); ++it)
            globalStats.preErrorHisto[it->first] += it->second;
        
        resize(globalStats.preErrorsAtPos, length(stats[i].preErrorsAtPos), globalStats.zeroCounts);
        for (unsigned pos = 0; pos < length(stats[i].preErrorsAtPos); ++pos)
            for (int err = 0; err < 4; ++err)
                globalStats.preErrorsAtPos[pos][err] += stats[i].preErrorsAtPos[pos][err];

        resize(globalStats.postErrorsAtPos, length(stats[i].postErrorsAtPos), globalStats.zeroCounts);
        for (unsigned pos = 0; pos < length(stats[i].postErrorsAtPos); ++pos)
            for (int err = 0; err < 4; ++err)
                globalStats.postErrorsAtPos[pos][err] += stats[i].postErrorsAtPos[pos][err];
    }

    // The quick stats are what we need for the table in the paper.
    std::cout << "QUICK STATS\n"
              << "gain\tbase error rate pre\tread error rate pre\tbase error rate post\tread error rate post";
    if (!options.indels)
        std::cout << "errors removed\terrors introduced";
    std::cout << "\n";
    fprintf(stdout, "%2.5f\t%2.5f\t%2.5f\t%2.5f\t%2.5f",
            100.0 * globalStats.diffErrorSum / globalStats.actualErrorSum,
            100.0 * globalStats.numErrorsPre / globalStats.numBasesPre,
            100.0 * globalStats.numErrorReadsPre / globalStats.numReads,
            100.0 * globalStats.numErrorsPost / globalStats.numBasesPost,
            100.0 * globalStats.numErrorReadsPost / globalStats.numReads);
    if (!options.indels)
        std::cout << "\t" << globalStats.numErrorsRemoved << "\t" << globalStats.numErrorsIntroduced;

    std::cout << "\n\n";

    // Print detailed statistics and histogram.
    std::cout << "STATISTICS\n"
              << "total read count     " << globalStats.numTotal << "\t\t(excludes unmapped reads)\n"
              << "unmapped pre count   " << globalStats.numUnmappedPre << "\n"
              << "ignored pre count    " << globalStats.numIgnoredPre << "\n"
              //<< "unmapped post count  " << globalStats.numUnmappedPost << "\n"
              //<< "unmapped both count  " << globalStats.numUnmappedBoth << "\n"
              << "\n"
              << "CORRECTION HISTOGRAM\n"
              << "\n"
              << "E.g. a value of 2 means 2 errors corrected, a value of -2 means two errors introduced.\n"
              << "\n"
              << "diff\tcount\t\tpercentage of reads\n";
    for (std::map<int, unsigned>::const_iterator it = globalStats.histo.begin(); it != globalStats.histo.end(); ++it)
        fprintf(stdout, "%3d\t%12u\t%5.2f\n", it->first, it->second, 100.0 * it->second / globalStats.numTotal);

    std::cout << "\n\nPRE-CORRECTION DISTRIBUTION\n\n"
              << "Distribution of errors before correction\n\n"
              << "error rate\tcount\t\tpercentage of reads\n";
    for (std::map<int, unsigned>::const_iterator it = globalStats.preErrorHisto.begin(); it != globalStats.preErrorHisto.end(); ++it)
        fprintf(stdout, "%3d\t%12u\t%5.2f\n", it->first, it->second, 100.0 * it->second / globalStats.numTotal);

    // Print positional error distributions PRE/POST
    std::cout << "\n\nPRE/POST-CORRECTION POSITIONAL ERROR DISTRIBUTION\n\n"
              << "position\tpre-mis\tpre-del\tpre-ins\tpre-bases\tpost-mis\tpost-del\tpost-ins\tpost-bases\n";

    unsigned maxLen = std::max(length(globalStats.preErrorsAtPos), length(globalStats.postErrorsAtPos));
    resize(globalStats.preErrorsAtPos, maxLen, globalStats.zeroCounts);
    resize(globalStats.postErrorsAtPos, maxLen, globalStats.zeroCounts);
    for (unsigned pos = 0; pos < maxLen; ++pos)
    {
        fprintf(stdout, "%3d", pos);
        for (unsigned err = 0; err < 4; ++err)
            fprintf(stdout, "\t%12u", globalStats.preErrorsAtPos[pos][err]);
        for (unsigned err = 0; err < 4; ++err)
            fprintf(stdout, "\t%12u", globalStats.postErrorsAtPos[pos][err]);
        std::cout << '\n';
    }

    // Print configuration.
    std::cout << "\n\nCONFIGURATION\n\n"
              << "  GENOME             \t" << options.pathGenome << "\n"
              << "  SAM PRE            \t" << options.pathSamPreCorrection << "\n"
              << "  SAM POST           \t" << options.pathSamPostCorrection << "\n"
              << "  FA/FQ POST         \t" << options.pathFastaFastqPostCorrection << "\n"
              << "  INDELS             \t" << getYesNo(options.indels) << "\n"
              << "  PADDING            \t" << options.padding << "\n"
              << "  BANDWIDTH          \t" << options.padding << "\n"
              << "  MIN UNCLIPPED BASES\t" << options.minUnclippedBases << "\n"
              << "  MAX ERROR RATE     \t" << options.maxErrorRate << "\n"
              << "  MAX ERROR COUNT    \t" << options.maxErrorCount << "\n"
              << "\n"
              << "  NUM THREADS        \t" << options.numThreads << "\n"
              << "  CORRECTION LOG     \t" << options.pathCorrectionLog << "\n"
              << "  LOG ALL            \t" << getYesNo(options.logAll) << "\n"
              << "  BANDWIDTH          \t" << options.bandwidth << "\n"
              << "\n";

    return 0;
}
