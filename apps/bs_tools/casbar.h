#ifndef __APPS_BS_TOOLS_CASBAR_H__
#define __APPS_BS_TOOLS_CASBAR_H__

#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>

#include <seqan/misc/svg.h>
#include <seqan/stream.h>

namespace seqan {

class Times
{
public:
    double time_all;
    double time_doBsCalling;
    double time_optimization;
    double time_IO;
    double time_convertPWToGlobal;

    static Times & instance()
    {
        static Times times;
        return times;
    }

private:
    Times() :
        time_all(0),
        time_doBsCalling(0),
        time_optimization(0),
        time_IO(0),
        time_convertPWToGlobal(0)
    {}
};

//////////////////////////////////////////////////////////////////////////////
// Default options

struct SnpStoreSpec_;
struct SnpStoreGroupSpec_;

template<>
struct FragmentStoreConfig<SnpStoreSpec_> :
    public FragmentStoreConfig<>
{
    typedef Owner<>	TReadSeqStoreSpec;
    typedef Owner<>	TAlignedReadTagStoreSpec;
    typedef Owner<> TReadNameStoreSpec;

	typedef double		TMappingQuality;    // -> mapq up to 255
};

template<>
struct FragmentStoreConfig<SnpStoreGroupSpec_> :
    public FragmentStoreConfig<>
{
    typedef Dependent<>	TReadSeqStoreSpec;
    typedef Owner<>		TAlignedReadTagStoreSpec;
    typedef Dependent<> TReadNameStoreSpec;

	typedef double		TMappingQuality;    //
};
}

using namespace seqan;


template <typename TGPos_>
struct SimplePosition
{
    typedef typename MakeSigned_<TGPos_>::Type TGPos;

    TGPos           gBegin;         // begin position in the genome
    unsigned        gseqNo;
};


struct SNPCallingOptions
{
    int         _debugLevel;                // level of verbosity
    bool        printVersion;               // print version number
    std::stringstream   programCall;        // stores snpstore program call

    // input output options
    unsigned    positionFormat;             // position format of mapped read input
                                            // 1..position space

    CharString          genomeFName;        // name of genome file
    String<CharString>  readFNames;         // list of read file names
    String<CharString>  qualityFNames;      // list of quality file names

    CharString          vcfOut;
    CharString          bedOut;
    CharString          outputLog;          // name of log output file

    // general parameters/options
    unsigned    maxPile;                    // keep at most maxPile many reads mapped to exact same position
    bool        laneSpecificMaxPile;        // do pile up correction on file by file basis (instead of on merged read set)
    bool        orientationAware;           // do pile up correction orientation aware

    bool        storeReadNames;             // store read names

    int         minMapQual;                 // min. mapping quality of reads parsed from SAM file
    bool        keepCigars;                 // remember cigar string of each match
    bool        keepSuboptimalReads;        // keep suboptimal read matches
    bool        keepMultiReads;             // keep multiply mapped reads

    bool        realign;                    // do realignment
    int         realignAddBorder;           // add flanking bases to reference in realignment (0 seems to work best)
    int         minClippedLength;           // discard read if it is not at least minClippedLength bases long after clipping
    bool        clipTagsInFile;             // helper to remember whether there were clip tags in file
    bool        softClipTagsInFile;         // helper to remember whether there was soft clipping in file

    int         asciiQualOffset;            // how to get quality values from ascii, usually subtract 33
    String<char> toIupac;                   // IUPAC code for het snp calls

     // SNP calling related
    bool        useBaseQuality;             // use base quality instead of min{base quality,mapping,avg read quality}
    unsigned    minCoverage;                // min depth at variant calling positions
    unsigned    excludeBorderPos;
    unsigned    minDifferentReadPos;        // number of different read positions that need to support the variant call

    float       avgQualT;                   // min. average quality value of variant read base
    unsigned    minMutT;                    // min. count of variant read base

    unsigned    indelCountThreshold;        // min. count of indels
    double pHetSnp;
    double pHomoSnp;


    // misc/helpers
    unsigned    maxHitLength;               // helper to remember max. alignment length
    unsigned    minCoord;                   // current min. read mapping coordinate observed
    unsigned    maxCoord;                   // current max. read mapping coordinate observed

    unsigned    windowSize;                 // genomic window size for read parsing
    unsigned    windowBuff;                 // reads within windowBuff base pairs of current window are also kept (-> overlapping windows)

    SNPCallingOptions() :
        _debugLevel(0),
        printVersion(false),
        programCall(""),
        //
        positionFormat(1),
        genomeFName(""),
        readFNames(""),
        qualityFNames(""),
        vcfOut(""),
        bedOut(""),
        outputLog(""),
        //
        maxPile(0),    // bs_change (was 0)
        laneSpecificMaxPile(true),
        orientationAware(false),
        storeReadNames(true),  // TODO change to false
        minMapQual(1),
        keepCigars(false),
        keepSuboptimalReads(false),
        keepMultiReads(false),
        realign(false),
        realignAddBorder(0),
        minClippedLength(10),
        clipTagsInFile(false),
        softClipTagsInFile(false),
        asciiQualOffset(33),
        toIupac("AMRWMCSYRSGKWYKT"),
        // SNP calling related
        useBaseQuality(true),
        minCoverage(6),
        excludeBorderPos(0),
        minDifferentReadPos(0),
        avgQualT(10),
        minMutT(3),
        indelCountThreshold(3),
        //
        pHetSnp(0.005),
        pHomoSnp(0.0005),
        //
        maxHitLength(1),
        minCoord(std::numeric_limits<unsigned>::max()),
        maxCoord(0),
        windowSize(100000),  // 10000?
        windowBuff(70)
    {}


};

struct BsInterval
{
    unsigned    startPos;
    unsigned    endPos;
    CharString  contigName;
};

struct MethCallingOptions
{
    double convRate;                // bs conversation rate
    double methConvRate;            // bs conversation rate for methylated Cs

    // Output method:
    bool outputMethProbs;           // Bayesian likelihood method
    bool outputMethStates;          // Threshold method

    String<double>  genPriors;      // Genotype prior probabilities; calculate in 'computeGenotypePriors'

    bool outputCandidates;
    bool outputAllBsStateProbs;

    // min. score best genotype must have to be called
    double minScoreToCallSnp;
    double minProbToCallSnp;

    unsigned maxCoverage;

    // Call snps without taking bs conversions into account
    bool ignoreBs;

    // Call only Snps at C positions (both strands)
    bool callOnlyCSnps;
    bool betaSampling;
    bool useMapq;

    bool helpPrint;

    unsigned threads;

    String<BsInterval> intervals;     // Intervals to analyze, must be sorted; only reads mapped into these intervals are loaded

    double statsCGMethylated;
    double statsCHGMethylated;
    double statsCHHMethylated;

    unsigned countCG;
    unsigned countCHG;
    unsigned countCHH;

    unsigned counteBViolated;
    unsigned countPlanB;
    unsigned countNoPlanB;
    unsigned countCovTooLow;
    unsigned countScoreTooLow;

    // model uniform
    bool uniformGenPriors;
    bool uniformSeqErrorsCalling;

    bool nonSimpleSubstErrors;
    bool nonSimpleInsErrors;
    bool nonSimpleDelErrors;
    double scalingFactorDelErrorsNonSimple;
    double scalingFactorDelErrorsSimple;

    // Realinging
    double delRate;
    double delErrorRate;
    double insErrorRate;
    double endGapScore;
    double scoreLimit;

    MethCallingOptions() :
        convRate(0.998),
        methConvRate(0.0),        // for the beginning assume there are no converations of methylated Cs

        outputMethProbs(true),
        outputMethStates(true),

        outputCandidates(false),
        outputAllBsStateProbs(true),

        minScoreToCallSnp(9),
        minProbToCallSnp(0.7),

        maxCoverage(500),

        ignoreBs(false),

        callOnlyCSnps(false),
        betaSampling(false),
        useMapq(false),

        helpPrint(false),
        threads(100),

        statsCGMethylated(0.0),
        statsCHGMethylated(0.0),
        statsCHHMethylated(0.0),
        countCG(0),
        countCHG(0),
        countCHH(0),

        counteBViolated(0),
        countPlanB(0),
        countNoPlanB(0),
        countCovTooLow(0),
        countScoreTooLow(0),

        uniformGenPriors(true),
        uniformSeqErrorsCalling(true),

        nonSimpleSubstErrors(false),
        nonSimpleInsErrors(false),
        nonSimpleDelErrors(false),
        scalingFactorDelErrorsNonSimple(3.5),
        scalingFactorDelErrorsSimple(5.0),

        delRate(0.0025),
        delErrorRate(0.001),
        insErrorRate(0.001),
        endGapScore(4.5),
        scoreLimit(-10)
    {}
};

// TODO types
// Single base info
struct SingleBaseInfo {

    unsigned qual;
    unsigned mapq;
    bool top;

    SingleBaseInfo():
                    qual(0),
                    mapq(0)
    {}
};


// Reference context
struct RefContext {

    unsigned pos;
    int refAllele;
    unsigned contextF;  // 0, 1, 2 (CG, CHG, CHH)
    unsigned contextR;
    CharString genomeID;
};

// BS change
struct MethylVariant {

    int genotype;
    bool genotypeCalled;
    bool bsCalled;
    long double methLevel1;
    long double methLevel2;

    long double score;

    long double genotypeProb;

    int totalCov;

    String<long double> bsStateProbs;   // Following order holds for diff. genotypes:
                                        // CC, CD, DD
                                        // GG, GH, HH
                                        // CG, DG, CH, DH
                                        // CX, DX
                                        // GX, H X

    MethylVariant()
                  : genotypeCalled(false),
                    bsCalled(false),
                    methLevel1(0.0),
                    methLevel2(0.0),
                    score(0.0),
                    genotypeProb(0.0),
                    totalCov(0)
    {}
};



//////////////////////////////////////////////////////////////////////////////
// Typedefs

// definition of a Read match
template <typename TGPos_>
struct MappedReadMatch
{
    typedef typename MakeSigned_<TGPos_>::Type TGPos;

    TGPos           gEnd;           // end position of the match in the genome              --> endPos
    unsigned        rseqNo;         // read seqNo                                           --> readId

    unsigned        gseqNo:15;      // genome seqNo     <32K sequences                      --> contigId
    unsigned        hasIndel:1;     // is 1 if read match contains indels, 0 else           --> gaps

    unsigned        editDist:3;     // Levenshtein distance <8                              --> errors
    unsigned        mScore:7;       // mapping quality  <128                                --> currently not in use anyway
    unsigned        avgQuality:6;   // avg read quality <64                                 --> score

    char            orientation;        // 'F'..forward strand, 'R'..reverse comp. strand   --> endPos > beginPos ?

};

enum CALLSNPS_ERROR {
    CALLSNPS_GFF_FAILED = 1,
    CALLSNPS_GENOME_FAILED = 2,
    CALLSNPS_QUALITY_FAILED = 3,
    CALLSNPS_OUT_FAILED = 4
};


template<typename TFile>
void
_printRecord(TFile &file, BamAlignmentRecord &record)
{
    file << "QueryName = " << record.qName << std::endl;
    file << "Flag      = " << record.flag << std::endl;
    file << "RefId     = " << record.rID << std::endl;
    file << "MapQ      = " << record.mapQ << std::endl;
    file << "bin       = " << record.bin << std::endl;
    file << "rNextId   = " << record.rNextId << std::endl;
    file << "pNext     = " << record.pNext << std::endl;
    file << "tLen      = " << record.tLen << std::endl;
    file << "seq       = " << record.seq << std::endl;
    file << "qual      = " << record.qual << std::endl;
    file << "tags      = " << record.tags << std::endl;
    for(unsigned i = 0; i < length(record.cigar); ++i)
        file << record.cigar[i].count << " " << record.cigar[i].operation << std::endl;
}


template<typename TBamTags, typename TOptions>
int
interpretBamTags(TBamTags & tags, int & editDist, bool & multi,
                int & clipLeft, int & clipRight, TOptions & options)
{
    BamTagsDict bamTags(tags);
    unsigned editDistIndex = 0;
    bool res1 = findTagKey(editDistIndex, bamTags, "NM");
    if(res1)
    {
        SEQAN_ASSERT_EQ('i', getTagType(bamTags, editDistIndex));
        extractTagValue(editDist, bamTags, editDistIndex);
    }
    else editDist = 1; // we dont know whether there are errors in the alignment, we assume there are..

    int numBest = 0;
    unsigned numBestIndex = 0;
    res1 = findTagKey(numBestIndex, bamTags, "X0");
    if(res1)
    {
        SEQAN_ASSERT_EQ('i', getTagType(bamTags, numBestIndex));
        extractTagValue(numBest, bamTags, numBestIndex);
        if(numBest > 1) multi = true;
    }

    unsigned clipIndex = 0;
    res1 = findTagKey(clipIndex, bamTags, "XC");
    if(res1)
    {
//      SEQAN_ASSERT_EQ('Z', getTagType(bamTags, clipIndex));  // XC is also used by BWA, also for clipping, but different fron ours
        if('Z' == getTagType(bamTags, clipIndex))
        {
            CharString clipLeftRight;
            extractTagValue(clipLeftRight, bamTags, clipIndex);
            // Get position of splitter char.
            unsigned x = 0;
            while (x < length(clipLeftRight) && isdigit(clipLeftRight[x]))
                ++x;
            // Extract left and right clipping count.
            seqan::CharString buffer = infix(clipLeftRight, 0, x);
            lexicalCastWithException(clipLeft, buffer);
            if (x + 1 <= length(clipLeftRight))
                buffer = infix(clipLeftRight, x + 1, length(clipLeftRight));
            else
                buffer = "0";
            lexicalCastWithException(clipRight, buffer);
            options.clipTagsInFile = true;
        }
    }

    // split read?

    // clip tags?

    // count tag?
    return 0;
}

/////////////////////////////////////////////////////////////
// read sorted(!)
template <
typename TSetContigAnchorGaps,
typename TBamFileIn,
typename TFragmentStore,
typename TContigId,
typename TContigPos,
typename TOptions
>
int readMatchesFromSamBam(
                         TSetContigAnchorGaps       &setContigAnchorGaps,     // Store contig gap anchors separately to convert at end all matches to msa
                         TBamFileIn                 &bamFileIn,
                         BamAlignmentRecord         &record,                // Need at the moment to check if record has already been read
                         TFragmentStore             &fragmentStore,             // forward/reverse fragmentStore.alignedReadStore
                         TFragmentStore             &fragmentStore1,             // to check order of reads regarding to contigs
                         TContigId                  currContigId,
                         TContigPos                 currentBegin,
                         TContigPos                 currentEnd,
                         TOptions                   &options)
{
    //std::cout << "readMatchesFromSamBam..." << std::endl;
    //bool setZero = true;
#ifdef CALL_PROFILE
    double timeStamp = sysTime();
#endif
    typedef typename TFragmentStore::TAlignedReadStore      TMatches;
    typedef typename Value<TMatches>::Type                  TMatch;
    typedef typename TFragmentStore::TAlignQualityStore     TMatchQualities;
    typedef typename Value<TMatchQualities>::Type           TMatchQuality;
    typedef typename TFragmentStore::TReadSeqStore          TReads;
    typedef typename Value<TReads>::Type                    TRead;
    typedef typename TFragmentStore::TReadStore             TReadStore;
    typedef typename Value<TReadStore>::Type                TReadStoreElement;
    typedef typename TMatch::TGapAnchors                    TReadGapAnchors;
    typedef Gaps<TRead, AnchorGaps<TReadGapAnchors> >       TReadGaps;
    //typedef typename TContig::TGapAnchors                   TContigGapAnchors;
    //typedef Gaps<TContigSeq, AnchorGaps<TContigGapAnchors> >    TContigGaps;
    typedef String<typename TFragmentStore::TContigGapAnchor>                                       TContigAnchorGaps;
    typedef Gaps<Nothing, AnchorGaps<TContigAnchorGaps> >                TContigGaps;


    //typedef Gaps<Nothing, AnchorGaps<typename TSAMContext::TContigAnchorGaps> >
    typedef String<typename TFragmentStore::TContigGapAnchor>                                       TContigAnchorGaps;

    typedef typename Id<TFragmentStore>::Type               TId;

    if(length(fragmentStore.readSeqStore)!=length(fragmentStore.alignQualityStore))
    {
        ::std::cerr << "Lengths need to be equal!!\n";
        return 10;
    }

    int readCount = length(fragmentStore.readSeqStore);
    TContigPos genomeLen = length(fragmentStore.contigStore[0].seq);

    // general stuff that is needed
    unsigned rSeq = readCount;
    Dna5String gInf;
    String<Dna5Q> curr_read;
    CharString readTemplate, temp_read;
    CharString readName, temp_str;
    TId prevRefId = 0;
    while (!atEnd(bamFileIn))
    {
        // read next record unless current one has not been handled yet
        if (empty(record.qName))
            readRecord(record, bamFileIn);

        if ( hasFlagUnmapped(record) || empty(record.cigar) || (!options.keepSuboptimalReads && hasFlagSecondary(record)))
        {
            //std::cout << "Read " << record.qName << " has Flag=" << record.flag << std::endl;
            clear(record); continue;
        }

        TId contigId;
        clear(temp_str);
        clear(temp_read);

        bool topStrand = true;
        bool hasIndel = false;
        int editDist = 0;
        int mScore;

        // Get global contigId to check, if order of reads is the same as order in contig files
        if (!getIdByName(contigId, fragmentStore1.contigNameStoreCache, contigNames(context(bamFileIn))[record.rID]))
        {
            clear(record);
            continue;
        }
        if ((TId)contigId < prevRefId)
        {
            std::cerr << "Read files need to be sorted according to chromosomes in genome file.\n";
            return CALLSNPS_GFF_FAILED;
        }

        prevRefId = contigId;
        if (contigId < (TId)currContigId)    // havent reached the sequence of interest yet
        {
            clear(record);
            continue;
        }
        if (contigId > (TId)currContigId)    // have passed the seq of interest
        {
            break;
        }
        /*if (hasFlagMultiple(record) && !hasFlagFirst(record) && hasFlagLast(record)) // read only left mates
        {
            clear(record);
            continue;
        }*/
        contigId = 0; // if we only store one chromosome at a time

        // skip whitespaces and read entry in column 2

        TContigPos beginPos = record.beginPos;
        if(beginPos > currentEnd + (TContigPos)options.windowBuff)  // we have passed the relevant match positions
        {
            if(options._debugLevel > 1)
                std::cout  << "gBegin "<< beginPos<<"  of match is too large\n";//, seeking "<<lineStart<<"\n";
            break;
        }
        if(options._debugLevel > 1)
            ::std::cout << beginPos << "\t";

        // need to calculate endPos
        TContigPos endPos;
        _getLengthInRef(endPos, record.cigar);
        endPos = beginPos + endPos;

        // check if cigar string has indels
        for(unsigned j=0; j < length(record.cigar); ++j)
            if(record.cigar[j].operation == 'D' || record.cigar[j].operation == 'N' || record.cigar[j].operation == 'I')
                hasIndel = true;

        if(options._debugLevel > 1)
            ::std::cout << endPos << "\t";
        if(endPos + (TContigPos)options.windowBuff < currentBegin)  //we havent reached a relevant read yet
        {
            clear(record);
            continue;
        }
        if(endPos > genomeLen)
            break;

        // must be parsed from tag
        mScore = (int)record.mapQ;
        if(options._debugLevel > 1)
            ::std::cout << mScore << "\t";

        if ((!hasFlagMultiple(record) && hasFlagRC(record)) ||      // se
            (hasFlagMultiple(record) &&  hasFlagRC(record) && !hasFlagLast(record)) ||     // bs right mates are simply projected on reverse complement strand
            (hasFlagMultiple(record) && !hasFlagRC(record) && hasFlagLast(record)))
        {
            topStrand = false;
        }


        if(options._debugLevel > 1)
            ::std::cout << "myID = "<<record.qName << "\n";

        TRead curr_read = record.seq;
        for(unsigned j = 0; j < length(record.qual); ++j)
        {
            int tempQual = _max(0,(int)ordValue(record.qual[j])-options.asciiQualOffset);
            assignQualityValue(curr_read[j],tempQual);
        }
        if (!topStrand)
            reverseComplement(curr_read);

        if (mScore >= options.minMapQual)        {
            if(empty(curr_read))
            {   //read sequence not found
                if(options._debugLevel>1)::std::cout << "neither quality nor read sequence found editDist = " << editDist <<"\n";
                return 1;
            }
#ifdef READ_NAME_AWARE
            if(!options.storeReadNames) clear(record.qName);
            TId readId;
            if(options.storeReadNames && !getIdByName(fragmentStore.readNameStore, record.qName, readId, fragmentStore.readNameStoreCache))
            {
                readId = length(fragmentStore.readSeqStore);
                appendValue(fragmentStore.readSeqStore,curr_read,Generous());
                appendValue(fragmentStore.readNameStore, record.qName, Generous());
            }

#else
            TId readId = length(fragmentStore.readSeqStore);
            appendValue(fragmentStore.readSeqStore,curr_read,Generous());
            if(!options.storeReadNames) clear(record.qName);
            appendValue(fragmentStore.readNameStore, record.qName, Generous());
#endif

            if(options._debugLevel > 1)
                ::std::cout<<fragmentStore.readSeqStore[rSeq]<<" with edit="<<editDist<<" at position "<< beginPos <<"\n";

            if(endPos - beginPos > (TContigPos)options.maxHitLength)
                options.maxHitLength = endPos - beginPos;

            // remember min and max positions seen
            if(beginPos < (TContigPos)options.minCoord || options.minCoord == numeric_limits<unsigned>::max()) options.minCoord = (unsigned)beginPos;
            if(endPos > (TContigPos)options.maxCoord) options.maxCoord =  (unsigned)endPos;

            // alignedReadStoreElement
            if(!topStrand)
            {
                TContigPos tmp = beginPos;
                beginPos = endPos;
                endPos = tmp;
            }
            TReadGapAnchors readGapAnchors;
            TReadGaps readGaps(record.seq, readGapAnchors);
            cigarToGapAnchorRead(readGaps, record.cigar);
            appendAlignment(fragmentStore, readId, contigId, beginPos, endPos, readGapAnchors);
            // Contig gap anchors
            //TContigGapAnchors contigGapAnchors;
            //TContigGaps contigGaps(contigGapAnchors);
            TContigAnchorGaps contigGapAnchors;
            TContigGaps contigGaps(contigGapAnchors);
            cigarToGapAnchorContig(contigGaps, record.cigar);
            appendValue(setContigAnchorGaps, contigGapAnchors);


#ifdef SNPSTORE_DEBUG
            if (beginPos < 300 || endPos < 300)
            {
                Dna5String contigInf = infix(fragmentStore.contigStore[0].seq, std::min(beginPos, endPos), std::max(beginPos, endPos));
                TContigGaps2 contigGaps2(contigInf, contigGapAnchors);
                std::cout << "readMatches..." << std::endl;
                std::cout << "record.qName: " << record.qName << std::endl;
                std::cout << "  contigGaps: " << contigGaps2 << std::endl;
                std::cout << "  readGaps:   " << readGaps << std::endl;
                std::cout << " beginPos: " << beginPos << "  endPos: " << endPos  << " topStrand: " << topStrand << std::endl;
                std::cout << " cigar: " << std::endl;
                for (unsigned i = 0; i < length(record.cigar); ++i)
                {
                    switch (record.cigar[i].operation)
                    {
                        case 'D': std::cout << 'D' << record.cigar[i].count; break;
                       // case 'N': std::cout << 'N' << record.cigar[i].count; break;
                       // case 'P': std::cout << 'P' << record.cigar[i].count; break;
                        case 'I': std::cout << 'I' << record.cigar[i].count; break;
                        case 'M': std::cout << 'M' << record.cigar[i].count; break;
                       // case 'S': std::cout << 'S' << record.cigar[i].count; break;
                       case 'X': std::cout << 'X' << record.cigar[i].count; break;
                    }
                }
                std::cout << std::endl;
            }
#endif


            // alignQualityElement
            TMatchQuality q;
            q.errors = (char)editDist;
            q.score = mScore;
            if(!options.realign && length(record.cigar)<=3) hasIndel = false;
            if(hasIndel)
                q.pairScore = 1;
            else
                q.pairScore = 0;
            appendValue(fragmentStore.alignQualityStore, q);
            // readStoreElement
            typename Value<TReadStore>::Type r;
            r.matePairId = TReadStoreElement::INVALID_ID;
            if (hasFlagMultiple(record))
            {
                if (!hasFlagRC(record) && hasFlagFirst(record)) r.matePairId = 1;       // top, forward         // TODO change to proper notation
                else if (hasFlagRC(record) && hasFlagLast(record)) r.matePairId = 2;    // top, reverse complement
                else if (!hasFlagRC(record) && hasFlagFirst(record)) r.matePairId = 1;  // bottom, reverse complement (regarding direction, is original strand)
                else if (hasFlagRC(record) && hasFlagLast(record)) r.matePairId = 2;    // bottom, forward (regarding direction, is actually rev. compl. of original)
            }
            appendValue(fragmentStore.readStore, r, Generous());

            ++rSeq;
            if(options._debugLevel > 1)
            {
                ::std::cout<<"Parsed: id= " << readId<<" name="<<record.qName<<"="<<curr_read<<" with edit="<<editDist<<" at position "<< beginPos<<"\n";
                ::std::cout << "mScore=" << mScore << " beginPos=" << beginPos << "endPos="<< endPos<<std::endl;
                if(q.pairScore==1) ::std::cout << "indel! pairScore=" << q.pairScore <<std::endl;
                if(q.pairScore==0) ::std::cout << "no indel! pairScore=" << q.pairScore <<std::endl;
                //if() ::std::cout << "reversed!" << std::endl;
            }
        }
        else
        {
            if(options._debugLevel > 1 )
            {
                ::std::cout<<"Discarded: "<<curr_read<<" with edit="<<editDist<<" at position "<< beginPos<<"\n";
                ::std::cout << "mScore = " << mScore << std::endl;
            }
        }
        clear(record);
    }

    if(options._debugLevel > 0)
        ::std::cout << ::std::endl << "Parsed "<<length(fragmentStore.alignedReadStore)<<" matches of "<<length(fragmentStore.readSeqStore)<<" reads." << ::std::endl;

#ifdef CALL_PROFILE
    Times::instance().time_IO += (sysTime() - timeStamp);
#endif

    return 0;
}










///////////////////////////////////////////////////////////////////////////////////////////////////




//looks for mismatches in alignemnt and returns positions with respect to 2nd row (read sequence)
template<typename TAlign, typename TString>
void
getMismatchMutations(TAlign & align, TString & mutations)
{
    typedef typename Row<TAlign>::Type TRow;
    typedef typename Iterator<TRow, Rooted>::Type TAlignIterator;


    TAlignIterator ali_it0_stop = iter(row(align,0),endPosition(cols(align)));
    TAlignIterator ali_it1_stop = iter(row(align,1),endPosition(cols(align)));
    TAlignIterator ali_it0 = iter(row(align,0),beginPosition(cols(align)));
    TAlignIterator ali_it1 = iter(row(align,1),beginPosition(cols(align)));


    //std::cout << "getting cigar line\n";//ali0 len = " <<ali_it0_stop-ali_it0 << " \t ali1 len = "<<ali_it1_stop-ali_it1<<"\n";
    int refPos = 0;
    //int readPos = 0;

    while(ali_it0 != ali_it0_stop && ali_it1 != ali_it1_stop)
    {
        while(ali_it0!=ali_it0_stop && ali_it1!=ali_it1_stop && !isGap(ali_it0)&& !isGap(ali_it1))
        {
            if(*ali_it1 != *ali_it0)
                appendValue(mutations,refPos);
            ++refPos;
            ++ali_it0;
            ++ali_it1;
        }
        while(ali_it0!=ali_it0_stop && isGap(ali_it0))
        {
            ++refPos;
            ++ali_it0;
            ++ali_it1;
        }
        while(isGap(ali_it1)&& ali_it1!=ali_it1_stop)
        {
            ++ali_it0;
            ++ali_it1;
        }
    }

}



//looks for position in source sequence of row1 and returns aligned position in row0
template<typename TAlign, typename TPosition>
int
getReadPos(TAlign & align, TPosition pos_row1, bool extraV = false)
{
    typedef typename Iterator<typename Row<TAlign>::Type, Rooted>::Type TAlignIterator;

    TAlignIterator ali_it0_stop = iter(row(align,0),endPosition(cols(align)));
    TAlignIterator ali_it1_stop = iter(row(align,1),endPosition(cols(align)));
    TAlignIterator ali_it0 = iter(row(align,0),beginPosition(cols(align))); // walks over read
    TAlignIterator ali_it1 = iter(row(align,1),beginPosition(cols(align))); // walks over ref

    int refPos = 0;
    int readPos = 0;
    if(extraV) std::cout << align ;
    while(ali_it0 != ali_it0_stop && ali_it1 != ali_it1_stop && refPos < pos_row1)
    {
        while(ali_it0!=ali_it0_stop && ali_it1!=ali_it1_stop && !isGap(ali_it0)&& !isGap(ali_it1) &&  refPos < pos_row1)
        {
            ++refPos;
            ++readPos;
            ++ali_it0;
            ++ali_it1;
        }
        while(ali_it0!=ali_it0_stop && isGap(ali_it0) &&  refPos < pos_row1)
        {
            ++refPos;
            ++ali_it0;
            ++ali_it1;
        }
        while(isGap(ali_it1)&& ali_it1!=ali_it1_stop &&  refPos <= pos_row1)
        {
            ++readPos;
            ++ali_it0;
            ++ali_it1;
        }
    }
    if(isGap(ali_it0)) return -1;
    else return readPos;
}

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////


struct SingleBaseVariant{
    bool called;    // did this variant pass calling criteria?
    int genotype;   // called diploid genotype (allele1 << 2 | allele2)
    int count;      // number of non-ref observations (only counting most frequent mutational base)
    int quality;    // a quality value associated with the genotype call
    int snpQuality; // a quality value associated with the SNP call
    int coverage;   // totalCoverage at position
};



template<typename TFragmentStore, typename TGroupStore, typename TMatchIterator>
void
copyFragmentStore(TGroupStore &fragStoreGroup,
                  TFragmentStore            &fragmentStore,
                  TMatchIterator            matchItBatchBegin,
                  TMatchIterator            matchItBatchEnd,
                  typename TFragmentStore::TContigPos   groupStartPos,
                  typename TFragmentStore::TContigPos   groupEndPos)
{
    //TFragmentStore fragStoreGroup = fragmentStore; //clear(fragStoreGroup.alignedReadStore); resize; arrayCopy(matchItBatchBegin,matchItBatchEnd,begin(fragStoreGroup.alignedReadStore,Standard())); // reads wont be needed anymore

    // pointers are enough
    fragStoreGroup.readSeqStore = fragmentStore.readSeqStore;
    fragStoreGroup.readStore = fragmentStore.readStore;
    fragStoreGroup.readNameStore = fragmentStore.readNameStore;
    fragStoreGroup.alignQualityStore = fragmentStore.alignQualityStore;


    // need to be copied / moved
    resize(fragStoreGroup.alignedReadStore,matchItBatchEnd-matchItBatchBegin,Exact());
    arrayCopy(matchItBatchBegin,matchItBatchEnd,begin(fragStoreGroup.alignedReadStore,Standard())); // reads wont be needed anymore
    //arrayMoveForward(matchItBatchBegin,matchItBatchEnd,begin(fragStoreGroup.alignedReadStore,Standard())); // reads wont be needed anymore

    // shorten reference sequence to the current region (groupStartPos to groupEndPos)
    // has to be copied because it will be overwritten
    // fragStoreGroup.contigStore[0].seq = infix(fragmentStore.contigStore[0].seq,groupStartPos,groupEndPos);
    typedef typename TGroupStore::TContigStore      TContigStore;           // TGenomeSet
    typedef typename Value<TContigStore>::Type      TContig;

    TContig conti;
    conti.seq = infix(fragmentStore.contigStore[0].seq,groupStartPos,groupEndPos);
    appendValue(fragStoreGroup.contigStore, conti, Generous() );
    appendValue(fragStoreGroup.contigNameStore, fragmentStore.contigNameStore[0], Generous() );

}





// check for the longest adjacent run of homopolymers
template<typename TSequence, typename TPosition>
inline typename Size<TSequence>::Type
checkSequenceContext(TSequence &reference,
                     TPosition candidatePos,
                     int indelSize)
{
    typename Size<TSequence>::Type count = 0;

    TPosition extendPos1, extendPos2;

    if(indelSize > 0) // deletion
    {
#ifdef SNPSTORE_DEBUG
        std::cout << "indelSize=" << indelSize << std::endl;
        std::cout << infix(reference,_max((int)0,(int)candidatePos-6),_min((int)candidatePos+indelSize+6,(int)length(reference)));
#endif

        // left candidate position
        extendPos1 = candidatePos > 0 ? candidatePos - 1 : 0;
        // right candidate position
        extendPos2 = candidatePos + (TPosition)indelSize < (TPosition)length(reference) ? candidatePos + (TPosition)indelSize : (TPosition)length(reference)-1;
    }
    else
    {
#ifdef SNPSTORE_DEBUG
        std::cout << "indelSize=" << indelSize << std::endl;
        std::cout << infix(reference,_max((int)0,(int)candidatePos-6),_min((int)candidatePos+6,(int)length(reference)));
#endif

        // left candidate position
        extendPos1 = candidatePos > 0 ? candidatePos - 1 : 0;
        // right candidate position
        extendPos2 = candidatePos;

    }

    typedef typename MakeSigned_<TPosition>::Type TSignedPos;
    // left candidate base
    typename Value<TSequence>::Type candBase = reference[extendPos1];

    //check to the left
    TSignedPos i = extendPos1;
    while(i >= 0 && reference[i]==candBase)
        --i;
    //check to the right
    TSignedPos j = extendPos1;
    while(j < (TSignedPos)length(reference) && reference[j]==candBase)
        ++j;
    count = j - i - 1;


    // right candidate base
    candBase = reference[extendPos2];

    //check to the left
    i = extendPos2;
    while(i >= 0 && reference[i]==candBase)
        --i;
    //check to the right
    j = extendPos2;
    while(j < (TSignedPos)length(reference) && reference[j]==candBase)
        ++j;
    count = j - i - 1 > (TSignedPos)count ? j - i - 1 : (TSignedPos)count;

#ifdef SNPSTORE_DEBUG
        std::cout << "done with seqContext" << std::endl;
#endif

    return count;

}



template<typename TFragmentStore, typename TStr>
void
_dumpMatches(TFragmentStore &fragmentStore, TStr str)
{
    //typedef typename TFragmentStore::TAlignedReadStore          TMatches;
    //typedef typename Value<TMatches>::Type                      TMatch;
    //typedef typename TFragmentStore::TAlignQualityStore         TMatchQualities;
    //typedef typename Value<TMatchQualities>::Type               TMatchQuality;
    //typedef typename TFragmentStore::TReadSeqStore              TReads;
    //typedef typename Value<TReads>::Type                        TRead;
    //typedef typename Iterator<TReads,Standard>::Type            TReadIt;
    //typedef typename Iterator<TMatchQualities,Standard>::Type   TMatchQIt;
    //typedef typename Iterator<TMatches,Standard>::Type          TMatchIt;

    std::cout << "Length of matches = " << length(fragmentStore.alignedReadStore)  << "\n";
    std::cout << "Length of reads   = " << length(fragmentStore.readSeqStore)  << "\n";
    std::cout << "Length of matchqs = " << length(fragmentStore.alignQualityStore)  << "\n";

    for(unsigned i = 0 ; i < length(fragmentStore.alignedReadStore); ++i)
    {
        char ori = (fragmentStore.alignedReadStore[i].beginPos < fragmentStore.alignedReadStore[i].endPos) ? 'F' : 'R';
        std::cout << "--"<<str<<"Match number " << i << ":\n";
        std::cout << "--"<<str<<"MatchId  = " << fragmentStore.alignedReadStore[i].id << "\n";
        std::cout << "--"<<str<<"ReadId   = " << fragmentStore.alignedReadStore[i].readId << "\n";
        std::cout << "--"<<str<<"ContigId = " << fragmentStore.alignedReadStore[i].contigId << std::flush << "\n";
        std::cout << "--"<<str<<"gBegin   = " << _min(fragmentStore.alignedReadStore[i].beginPos, fragmentStore.alignedReadStore[i].endPos) << "\n";
        std::cout << "--"<<str<<"gEnd     = " << _max(fragmentStore.alignedReadStore[i].beginPos, fragmentStore.alignedReadStore[i].endPos) << "\n";
        std::cout << "--"<<str<<"orient   = " << ori << std::flush << std::endl;
        if(length(fragmentStore.alignQualityStore) > fragmentStore.alignedReadStore[i].id)
        {
            std::cout << "--"<<str<<"EditDist = " << (int) fragmentStore.alignQualityStore[fragmentStore.alignedReadStore[i].id].errors << "\n";
            std::cout << "--"<<str<<"AvgQ     = " << (int)fragmentStore.alignQualityStore[fragmentStore.alignedReadStore[i].id].score << "\n";
        }
        std::cout << "--"<<str<<"Readseq  = " << fragmentStore.readSeqStore[fragmentStore.alignedReadStore[i].readId] << std::flush << "\n";

    }
}



///////////////////////////////////////////////////////////////////////////////////////
// Build connected subsets, check for indels and realign if necessary, then doSnpAndMethCalling
template <
    typename TFragmentStore,
    typename TContigPos,
    typename TSetContigAnchorGaps,
    typename TVcfStream,
    typename TBedStream,
    typename TMethOptions,
    typename TOptions
>
void doCheckRealignCall(
    TFragmentStore          &fragmentStore,         // forward/reverse matches
    TContigPos              startCoord,         // startCoordinate + posOnGenomeInfix = real coordinate on whole chromosome
    TContigPos              currWindowBegin,
    TContigPos              currWindowEnd,
    TSetContigAnchorGaps    &setContigAnchorGaps,
    TVcfStream              &vcfStream,
    TBedStream              &bedStream,
    TMethOptions            &methOptions,
    TOptions                &options)
{
    //std::cout << " doCheckRealignCall ..." << std::endl;
    typedef typename TFragmentStore::TAlignedReadStore              TMatches;
    typedef typename Value<TMatches>::Type                          TMatch;
    typedef typename TFragmentStore::TAlignQualityStore             TMatchQualities;
    typedef typename Iterator<TMatches,Standard>::Type              TMatchIterator;

    // for test
    typedef typename TFragmentStore::TContigStore       TContigStore;
    typedef typename Value<TContigStore>::Type          TContig;
    typedef typename TFragmentStore::TContigSeq         TContigSeq;
    typedef Gaps<TContigSeq, AnchorGaps<typename TContig::TGapAnchors> >    TContigGaps;

    TMatches &matches = fragmentStore.alignedReadStore;
    TMatchQualities &matchQualities = fragmentStore.alignQualityStore;

    std::sort(begin(matches, Standard()), end(matches, Standard()), LessGPos<TMatch>());

    TMatchIterator matchIt = begin(matches,Standard());
    TMatchIterator matchItEnd = end(matches,Standard());

    SEQAN_ASSERT_EQ(length(fragmentStore.alignedReadStore), length(setContigAnchorGaps));

    // now find connected subsets, i.e. groups of reads that overlap
    // dont realign regions unworthy of realignment (no indel reads)
    while(matchIt != matchItEnd)
    {
        TMatchIterator matchItBatchBegin = matchIt;

        TContigPos groupEndPos = _max((*matchIt).endPos,(*matchIt).beginPos);
        TContigPos groupStartPos = _min((*matchIt).endPos,(*matchIt).beginPos);

        // Translate coordinates to group-local ones.
        TContigPos groupStartCoordLocal = _max(0,(int)groupStartPos-options.realignAddBorder);

        int indelReadCount = 0; // how many reads have indels in the current group
        while(matchIt != matchItEnd && _min((*matchIt).beginPos,(*matchIt).endPos) < groupEndPos)   // Find connected subset
        {
            groupEndPos = (_max((*matchIt).beginPos,(*matchIt).endPos) > groupEndPos) ? _max((*matchIt).beginPos,(*matchIt).endPos) : groupEndPos;
            // reads wont be needed anymore! (make sure this is the case!!!)
            (*matchIt).beginPos -= groupStartCoordLocal;
            (*matchIt).endPos -= groupStartCoordLocal;
            if(matchQualities[(*matchIt).id].pairScore == 1 ) ++indelReadCount;

            ++matchIt;
        }
        TMatchIterator matchItBatchEnd = matchIt;
        unsigned numMatches = matchItBatchEnd -matchItBatchBegin;

        TContigPos groupEndCoordLocal = _min(groupEndPos+(TContigPos)options.realignAddBorder,(TContigPos)length(fragmentStore.contigStore[0].seq));

        if(numMatches >= options.minCoverage)
        {
            //make temporary fragstore for group
            // shorten reference sequence to the current region (groupStartPos to groupEndPos)

            TFragmentStore fragStoreGroup = fragmentStore;
            resize(fragStoreGroup.alignedReadStore,numMatches);
            arrayMoveForward(matchItBatchBegin,matchItBatchEnd,begin(fragStoreGroup.alignedReadStore,Standard())); // reads wont be needed anymore
            fragStoreGroup.contigStore[0].seq = infix(fragmentStore.contigStore[0].seq,groupStartCoordLocal,groupEndCoordLocal);
            //TSetContigAnchorGaps groupSetContigAnchorGaps;
            //resize(groupSetContigAnchorGaps, numMatches);
            //arrayMoveForward(contigGapsBatchBegin, contigGapsBatchEnd, begin(groupSetContigAnchorGaps,Standard()));

#ifdef CALL_PROFILE
            double timeStamp = sysTime();
 #endif
            convertPairWiseToGlobalAlignment(fragStoreGroup, setContigAnchorGaps);
#ifdef CALL_PROFILE
            Times::instance().time_convertPWToGlobal += (sysTime() - timeStamp);
#endif
            if(false)
            {
                TContigGaps contigGaps(fragStoreGroup.contigStore[0].seq, fragStoreGroup.contigStore[0].gaps);
                TContigPos maxPos = positionSeqToGap(contigGaps,length(fragStoreGroup.contigStore[0].seq)-1)+1;
                maxPos = _max(maxPos,(TContigPos)length(fragStoreGroup.contigStore[0].seq));
                std::cout << "maxPos visual = " << maxPos << std::endl;
                AlignedReadLayout layout;
                layoutAlignment(layout, fragStoreGroup);
                printAlignment(std::cout, layout, fragStoreGroup, 0, (TContigPos)(maxPos-100), (TContigPos)maxPos, 0, 150);
            }


#ifdef SNPSTORE_DEBUG
            std::cout << " groupStartCoordLocal = " << groupStartCoordLocal  << " groupEndCoordLocal=" <<  groupEndCoordLocal << std::endl;
            std::cout << " groupEndPos = " <<  groupEndPos << " groupStartPos=" <<  groupStartPos << std::endl;
            std::cout << "genomeLength= " <<  length(fragmentStore.contigStore[0].seq) << std::endl;

#endif
            groupStartPos += startCoord;    // startCoord from current window + min observed group pos
            groupEndPos += startCoord;
            TContigPos groupStartCoord = startCoord + groupStartCoordLocal;             // same as above, but adjusted to realignBorder
            TContigPos groupWindowBegin = _max(groupStartPos,currWindowBegin);          // Current window to analyze begin
            TContigPos groupWindowEnd = _min(groupEndPos,currWindowEnd);                // current window end

            //std::cout << " groupStartCoord = " << groupStartCoord << std::endl;
            //std::cout << " groupStartPos = " <<  groupStartPos << " groupEndPos=" <<  groupEndPos << std::endl;

            //the current group is formed by all reads from matchItBatchBegin until matchItBatchEnd
            if(indelReadCount >= (int)options.indelCountThreshold)
            {
                //std::cout << "groupWindowBegin: " << groupWindowBegin << "  groupWindowEnd: " << groupWindowEnd << std::endl;
                if(groupWindowBegin <= 238380 && groupWindowEnd >= 238440)
                {
                    std::cout << "Before realigning " << std::endl;
                    TContigGaps contigGaps(fragStoreGroup.contigStore[0].seq, fragStoreGroup.contigStore[0].gaps);
                    AlignedReadLayout layout;
                    layoutAlignment(layout, fragStoreGroup);
                    unsigned start = positionSeqToGap(contigGaps, 238380  - groupStartCoord);
                    unsigned end = positionSeqToGap(contigGaps,  238440 - groupStartCoord);
                    std::cout << "Start: " <<238380 << " end: " <<  238440 << std::endl;
                    printAlignment(std::cout, layout, fragStoreGroup, 0, start, end, 0, 150);
                }

                //do realignment
                doRealigning(fragStoreGroup, groupWindowBegin, groupWindowBegin, methOptions, options);   // TODO: do only for reads spanning potential indel position

                doSnpAndMethCalling(fragStoreGroup, groupStartCoord, groupWindowBegin, groupWindowEnd, true, vcfStream, bedStream, methOptions, options);
            }
            else
            {
                doSnpAndMethCalling(fragStoreGroup, groupStartCoord, groupWindowBegin, groupWindowEnd, false, vcfStream, bedStream, methOptions, options);
            }
        }
    }
}


//////////////////////////////////////////////////////////////////////////////
// Output SNPs
template <
    typename TFragmentStore,
    typename TVcfStream,
    typename TBedStream,
    typename TMethOptions,
    typename TOptions
>
void doSnpAndMethCalling(
    TFragmentStore              &fragmentStore,
    typename TFragmentStore::TContigPos startCoord,         // groupStartCoord (startCoordinate + posOnGenomeInfix = real coordinate on whole chromosome)
    typename TFragmentStore::TContigPos currStart,          // Curr. window for calling
    typename TFragmentStore::TContigPos currEnd,
    bool                    didRealign,
    TVcfStream              &vcfStream,
    TBedStream              &bedStream,
    TMethOptions            &methOptions,
    TOptions                &options)
{
    // std::cout << "doSnpAndMethCalling...." << std::endl;
    typedef typename TFragmentStore::TAlignedReadStore  TMatches;
    typedef typename Value<TMatches>::Type              TMatch;
    typedef typename TFragmentStore::TAlignQualityStore TMatchQualities;
    typedef typename TFragmentStore::TReadSeqStore      TReads;
    typedef typename Value<TReads>::Type                TRead;
    typedef typename TFragmentStore::TContigPos         TContigPos;
    typedef typename TFragmentStore::TContigSeq         TContigSeq;
    typedef typename Iterator<TMatches,Standard>::Type  TMatchIterator;
    typedef typename Value<TMatches>::Type              TReadStoreElement;


    typedef typename TFragmentStore::TContigStore       TContigStore;
    typedef typename Value<TContigStore>::Type          TContig;

    typedef Gaps<TContigSeq, AnchorGaps<typename TContig::TGapAnchors> >    TContigGaps;
    typedef Gaps<TRead, AnchorGaps<typename TMatch::TGapAnchors> >          TReadGaps;
    typedef typename Iterator<TContigGaps>::Type                            TContigGapIter;
    typedef typename Iterator<TReadGaps>::Type                              TReadGapIter;

    SEQAN_PROTIMESTART(dump_time);
    // matches need to be ordered accordign to genome position
    TReads &reads                   = fragmentStore.readSeqStore;
    TMatches &matches               = fragmentStore.alignedReadStore;
    TMatchQualities &matchQualities = fragmentStore.alignQualityStore;

    // forward match qualities
    String<int> columnQualityF;             resize(columnQualityF,5);
    String<unsigned> countF;                resize(countF,5);
    String<CharString> qualityStringF;      resize(qualityStringF,5);
    String<String<int> > mapqsF;            resize(mapqsF, 5);  // bs change
    String<String<bool> > originStringF;    resize(originStringF, 5);

    // reverse match qualities
    String<int> columnQualityR;             resize(columnQualityR,5);
    String<unsigned> countR;                resize(countR,5);
    String<CharString> qualityStringR;      resize(qualityStringR,5);
    String<String<int> > mapqsR;            resize(mapqsR, 5);
    String<String<bool> > originStringR;    resize(originStringR, 5);


    // both
    String<unsigned> count;             resize(count,5);
    String<unsigned> columnQuality;     resize(columnQuality,5);

    FunctorComplement<Dna5> f;

    // sort reads according to begin position, if not already done in after realigning
    /*if (!didRealign)*/ sortAlignedReads(fragmentStore.alignedReadStore, SortBeginPos());

    TMatchIterator matchIt  = begin(matches, Standard());
    TMatchIterator matchItEnd   = end(matches, Standard());
    if (didRealign) matchItEnd--; // exclude reference sequence

    TContigSeq reference = fragmentStore.contigStore[0].seq;
    TContigGaps referenceGaps(fragmentStore.contigStore[0].seq, fragmentStore.contigStore[0].gaps);

    TContigPos      refStart;
    String<Pair<short unsigned,short unsigned> > indelConsens;

    if (didRealign)
    {
        //unsigned numReads = length(matches)-1; // exclude reference sequence
        //refStart = (TContigPos)fragmentStore.alignedReadStore[numReads].beginPos; TODO
        //insertGaps(referenceGaps, 0, refStart); TODO
        refStart = 0;

        //TContigGaps referenceGapsTemp(fragmentStore.contigStore[0].seq, fragmentStore.alignedReadStore[numReads].gaps);
        //clear(referenceGaps);
        //referenceGaps = referenceGapsTemp;
        // for indels:
        // i1 keeps track of consensus character
        // i2 keeps track of coverage (last 8 bits) and indelcount (first 8 bits)
        resize(indelConsens,refStart + length(referenceGaps));
        for(unsigned i = 0; i < refStart + length(referenceGaps) ; ++i)
        {
            indelConsens[i].i1 = 6;
            indelConsens[i].i2 = 0;
        }
    }
    else
    {
        refStart = 0;
    }
#ifdef SNPSTORE_DEBUG
        ::std::cout << "lengthrefgaps=" << length(referenceGaps)<< std::endl;
        ::std::cout << " length(ref)=" << length(reference) << std::endl;
        std::cout << "length alignedReadStore: " << length(fragmentStore.alignedReadStore) << std::endl;
#endif


    if(options._debugLevel>1) std::cout << "Start inspecting alignment..." << std::endl;
    // now walk through the reference sequence in gaps view space,
    // i.e. position may be a gap
    // example:
    // Ref      ACCGTGCACTAGCATCATT--ACTAGCATCATA
    // Reads    ACCGTACA--AGCATCAT
    //              TACA--AGCATCATT--ACT
    //                          ATTTTACTAGCATCATA
    //

    if(false) //currStart <=  238380  && currEnd >=  238440)
    {
        std::cout << "After realigning2 " << std::endl;
        TContigGaps contigGaps(fragmentStore.contigStore[0].seq, fragmentStore.contigStore[0].gaps);
        TContigPos maxPos = positionSeqToGap(contigGaps,length(fragmentStore.contigStore[0].seq)-1)+1;
        maxPos = _max(maxPos,(TContigPos)length(fragmentStore.contigStore[0].seq));
        AlignedReadLayout layout;
        layoutAlignment(layout, fragmentStore);
        //unsigned start = positionSeqToGap(contigGaps, 238380  - startCoord);
        //unsigned end = positionSeqToGap(contigGaps,   238440 - startCoord);
        //std::cout << "Start: " << 238380  << " end: " <<   238440 << std::endl;
        //printAlignment(std::cout, layout, fragmentStore, 0, start, end, 0, 150);
        printAlignment(std::cout, layout, fragmentStore, 0, 0, 200, 0, 150);

    }


    for(TContigPos candidateViewPos = refStart; candidateViewPos < refStart + (TContigPos)length(referenceGaps); ++candidateViewPos)
    {
        // first check if reference has a gap (potential insertion in reads) at this position
        bool refGap = false;

        TContigGapIter refIt = iter(referenceGaps,candidateViewPos-refStart);
        if(isGap(refIt)) refGap = true;

        //get position in sequence space
        TContigPos candidatePos;
        candidatePos = positionGapToSeq(referenceGaps, candidateViewPos-refStart);
        //else candidatePos = candidateViewPos;   // TODO do we need reference gap positions?

        //std::cout << "candidateViewPos: " << candidateViewPos << " refStart: " << refStart << std::endl;
        //std::cout << "candidatePos : " << candidatePos << " startCoord: " << startCoord << "  currStart: " << currStart << "  currEnd: " << currEnd << std::endl;
        // not in the current window yet
        if(candidatePos + startCoord < currStart) continue;
        // not in the current window anymore
        if(candidatePos + startCoord >= currEnd) break;

        Dna5 refBase = reference[candidatePos]; // what happens if refGap==true, esp. for leading gaps?
        if(refBase=='N' || refGap) continue;

#ifdef SNPSTORE_DEBUG
        if (candidatePos < 10)
        {
            std::cout << "refStart: " << refStart << std::endl;
            std::cout << "candidateViewPos = " << candidateViewPos <<  std::endl;
            std::cout << "candidatePos = " << candidatePos << std::endl;
            std::cout << "candidatePosMitStart = " << candidatePos + startCoord << " refBase = " << refBase << std::endl;
            if(refGap) std::cout << "refGap!" << std::endl;
        }
#endif

        //find range of relevant read matches
        // CHECK: remove unnecessarily walking through same matches multiple times
        while(matchIt != matchItEnd &&  _max((*matchIt).endPos,(*matchIt).beginPos) <= candidateViewPos)
            ++matchIt;
        TMatchIterator matchRangeBegin = matchIt;
        while(matchIt != matchItEnd &&  _min((*matchIt).endPos,(*matchIt).beginPos)  <= candidateViewPos)
            ++matchIt;
        TMatchIterator matchRangeEnd = matchIt; // could remember this for next round
        matchIt = matchRangeBegin;

        int coverage = matchRangeEnd-matchRangeBegin;
#ifdef SNPSTORE_DEBUG
        std::cout <<"cov=" << coverage << std::endl;
#endif
        if(coverage<(int)options.minCoverage)
            continue; // coverage too low

        // start checking reads for this position, prepare some helpers
        Dna5 candidateBase;
        int quality;
        std::set<int> readPosMap;
        std::set<int> indelReadPosMap;

        for(unsigned t=0;t<5;++t)
        {
            countF[t] = 0;
            columnQualityF[t] = 0;
            clear(qualityStringF[t]);
            clear(mapqsF[t]);
            clear(originStringF[t]);

            countR[t] = 0;
            columnQualityR[t] = 0;
            clear(qualityStringR[t]);
            clear(mapqsR[t]);
            clear(originStringR[t]);
        }

        bool observedAtLeastOneMut = false;
        int numIndelsObservedF = 0;  // if refGap then this counts the number of insertions on forward
        int numIndelsObservedR = 0;  // if refGap then this counts the number of insertions on reverse
                        // else it counts the number of deletions
        int indelQualF = 0;
        int indelQualR = 0;

        unsigned positionCoverage = 0;   // how many reads actually span the position?

        // now check reads
        while(matchIt != matchRangeEnd)
        {
            TContigPos currViewBegin = _min((*matchIt).beginPos,(*matchIt).endPos); // gap-space
            TContigPos currViewEnd = _max((*matchIt).beginPos,(*matchIt).endPos);

            // make sure this match is really spanning the position
            if(!(currViewBegin <= candidateViewPos && candidateViewPos < currViewEnd))
            {
                ++matchIt;
                continue;
            }
            ++positionCoverage;

            char orientation = ((*matchIt).beginPos > (*matchIt).endPos) ? 'R' : 'F';

            TRead readSeq;
            if (orientation == 'F') readSeq = reads[(*matchIt).readId];
            else
            {
                reverseComplement(reads[(*matchIt).readId]);
                readSeq = reads[(*matchIt).readId];
                reverseComplement(reads[(*matchIt).readId]);
            }
            TReadGaps readGaps(readSeq,(*matchIt).gaps);
            TReadGapIter rgIt = iter(readGaps,candidateViewPos - currViewBegin);

            // check out which position is hit in this read
            int readPos;
            if(isGap(rgIt)) readPos = -1; //potential deletion in reads (insertion in reference)
            else
            {
                readPos = positionGapToSeq(readGaps,candidateViewPos - currViewBegin);
                if(orientation == 'R')
                    readPos = length(reads[(*matchIt).readId]) - readPos - 1;
            }

#ifdef SNPSTORE_DEBUG
            /*
            TContigGaps contigGaps(fragmentStore.contigStore[0].seq, fragmentStore.contigStore[0].gaps);

            setBeginPosition(contigGaps, positionGapToSeq(contigGaps, currViewBegin));
            setEndPosition(contigGaps, positionGapToSeq(contigGaps, currViewEnd));
            std::cout << "   doSnpAndMethCalling..." << std::endl;
            std::cout << "readname: " << fragmentStore.readNameStore[(*matchIt).readId] << std::endl;
            std::cout << "alignedReadStoreId: " << value(matchIt) << std::endl;
            std::cout << "  contigGaps: " << contigGaps << std::endl;
            std::cout << "  readGaps:   " << readGaps << std::endl;
            std::cout << " beginPos: " << (*matchIt).beginPos << "  endPos: " << (*matchIt).endPos << std::endl;
            std::cout << "length(readGaps): " << length(readGaps) << std::endl;
            std::cout << "candidateViewPos: " << candidateViewPos << "  currViewBegin: " << currViewBegin << "  currViewEnd: " << currViewEnd  << std::endl;

            std::cout << "Pos. in read = " << readPos  << "  beginPos of read: " << currViewBegin << std::endl;
            */
#endif
            if (false) //candidateViewPos == 1000)
            {
                TContigGaps contigGaps(fragmentStore.contigStore[0].seq, fragmentStore.contigStore[0].gaps);
                setBeginPosition(contigGaps, positionGapToSeq(contigGaps, currViewBegin));
                setEndPosition(contigGaps, positionGapToSeq(contigGaps, currViewEnd));
                std::cout << "  doSnpAndMethCalling..." << std::endl;
                //std::cout << "readId: " << (*matchIt).readId  << " length(readNameStore): " << length(fragmentStore.readNameStore) << std::endl;
                //std::cout << "readname: " << fragmentStore.readNameStore[(*matchIt).readId] << std::endl;
                std::cout << "alignedReadStoreId: " << value(matchIt) << std::endl;
                std::cout << "  contigGaps: " << contigGaps << std::endl;
                std::cout << "  readGaps:   " << readGaps << std::endl;
                CharString str;
                resize(str, length(readGaps), 'A');
                int gapReadPos = candidateViewPos - currViewBegin;
                //if(orientation == 'R') gapReadPos = length(readGaps) - gapReadPos - 1;
                str[gapReadPos] = '-';
                std::cout << "  currPos :   " << str << "  : " << gapReadPos << std::endl;
                std::cout << " beginPos: " << (*matchIt).beginPos << "  endPos: " << (*matchIt).endPos << std::endl;
                std::cout << "length(readGaps): " << length(readGaps) << std::endl;
                std::cout << "candidateViewPos: " << candidateViewPos << "  currViewBegin: " << currViewBegin << "  currViewEnd: " << currViewEnd  << std::endl;

                std::cout << "Pos. in read = " << readPos  << "  beginPos of read: " << currViewBegin << std::endl;
            }


            if(readPos != -1) //-1 indicates gap in read
            {

                if(orientation == 'R') candidateBase = f((Dna5)reads[(*matchIt).readId][readPos]);
                else candidateBase = (Dna5)reads[(*matchIt).readId][readPos];

                if(refGap)
                {
                    if(orientation == 'F')
                        ++numIndelsObservedF; // count insertions
                    else ++numIndelsObservedR;
                    if(options.minDifferentReadPos > 0)
                        if((unsigned)(length(reads[(*matchIt).readId]) - readPos) > options.excludeBorderPos  &&
                            (unsigned) readPos >= options.excludeBorderPos )
                        indelReadPosMap.insert(readPos);
                }
                else if(candidateBase != refBase)
                {
                    observedAtLeastOneMut = true;
                    if(options.minDifferentReadPos > 0)
                        if((unsigned)(length(reads[(*matchIt).readId]) - readPos) > options.excludeBorderPos  &&
                            (unsigned) readPos >= options.excludeBorderPos )
                        readPosMap.insert(readPos);
                }
                quality = getQualityValue(reads[(*matchIt).readId][readPos]);

                /*
                if(!options.useBaseQuality && quality > (int)matchQualities[(*matchIt).id].score)
                {   // dont trust the quality of this position more
                    // than the average quality of this read
                    quality = (int) matchQualities[(*matchIt).id].score;
                }
                */

                if(orientation == 'F')
                {
                    columnQualityF[ordValue(candidateBase)] += quality;
                    ++countF[ordValue(candidateBase)];
                    appendValue(qualityStringF[ordValue(candidateBase)],(char)(quality+33));
                    appendValue(mapqsF[ordValue(candidateBase)], (int)matchQualities[(*matchIt).id].score);
                    if (fragmentStore.readStore[(*matchIt).readId].matePairId != TReadStoreElement::INVALID_ID || fragmentStore.readStore[(*matchIt).readId].matePairId != 1)
                        appendValue(originStringF[ordValue(candidateBase)],  true);
                    else
                        appendValue(originStringF[ordValue(candidateBase)],  false);
                }
                else
                {
                    columnQualityR[ordValue(candidateBase)] += quality;
                    ++countR[ordValue(candidateBase)];
                    appendValue(qualityStringR[ordValue(candidateBase)],(char)(quality+33),Generous());
                    appendValue(mapqsR[ordValue(candidateBase)], (int)matchQualities[(*matchIt).id].score);
                    if (fragmentStore.readStore[(*matchIt).readId].matePairId != TReadStoreElement::INVALID_ID || fragmentStore.readStore[(*matchIt).readId].matePairId != 1)
                        appendValue(originStringR[ordValue(candidateBase)], true);
                    else
                        appendValue(originStringR[ordValue(candidateBase)], false);
                }
            }
            else
            {   //potential deletions

                if(!refGap)
                {
                    readPos = positionGapToSeq(readGaps,candidateViewPos - currViewBegin);
#ifdef SNPSTORE_DEBUG
                    std::cout <<"del readPos = " << readPos  << " readlength=" << length(reads[(*matchIt).readId]) << std::endl;
#endif
                    if(orientation == 'R')
                        readPos = length(reads[(*matchIt).readId]) - readPos;
#ifdef SNPSTORE_DEBUG
                    std::cout <<"del readPos = " << readPos  << " readlength=" << length(reads[(*matchIt).readId]) << std::endl;
#endif
                    quality = (getQualityValue(reads[(*matchIt).readId][readPos-1]) + getQualityValue(reads[(*matchIt).readId][readPos])) / 2;
                    if(orientation == 'F')
                    {
                        indelQualF += quality;
                        ++numIndelsObservedF;
                    }
                    else
                    {
                        ++numIndelsObservedR;
                        indelQualR += quality;
                    }
                    if(options.minDifferentReadPos > 0)
                    {
                        if((unsigned)(length(reads[(*matchIt).readId]) - readPos) > options.excludeBorderPos  &&
                            (unsigned) readPos >= options.excludeBorderPos )
                        indelReadPosMap.insert(readPos);
                    }
                }
            }
            ++matchIt;
        }

        matchIt = matchRangeBegin; //set iterator back to where we started from, same matches might be involved in next cand pos

#ifdef SNPSTORE_DEBUG
        if (candidatePos + startCoord == 1000)
        {
            int numIndelsObserved = numIndelsObservedF + numIndelsObservedR;
            std::cout << "posCov=" << positionCoverage << "numIndels = " << numIndelsObserved << std::endl;
            if(observedAtLeastOneMut) std::cout << "observed at least one mut " << std::endl;
        }
#endif

        bool isSnp = true;

        // coverage depth
        int refAllele = ordValue(reference[candidatePos]);
        unsigned realCoverageF = countF[0] + countF[1] +countF[2] +countF[3] +countF[4];
        unsigned realCoverageR = countR[0] + countR[1] +countR[2] +countR[3] +countR[4];

        // too few reads actually cover the position
        if(positionCoverage < options.minCoverage)
            isSnp = false;

        // is the min. number of different read positions supporting the mutation met?
        if(isSnp && options.minDifferentReadPos > 0 && readPosMap.size() < options.minDifferentReadPos)
            isSnp = false;

        if (!refGap && observedAtLeastOneMut && (realCoverageF <= (options.minCoverage/2.0) || realCoverageR <= (options.minCoverage/2.0)))  ++methOptions.countCovTooLow;

        //all observed bases match the reference allele or there were too few indels
        MethylVariant meth;
        if ( (refBase == 'C' || refBase == 'G' || (observedAtLeastOneMut && positionCoverage > options.minCoverage)) && !refGap)
        {
            if (!methOptions.outputCandidates &&
                    (positionCoverage > methOptions.maxCoverage ||                                                    // Discard positions with too high coverage
                    (realCoverageF <= (options.minCoverage/2.0) || realCoverageR <= (options.minCoverage/2.0)) ) )    // Check min. coverage on both sides
                continue;

            RefContext refContext;
            refContext.pos = candidatePos + startCoord;
            refContext.refAllele = refAllele;
            refContext.genomeID = fragmentStore.contigNameStore[0];
            // get reference context if possible
            // forward
            if (candidatePos+1 < (TContigPos)length(reference) && reference[candidatePos+1] == 'G')
                refContext.contextF = 0;
            else if (candidatePos+2 < (TContigPos)length(reference) && reference[candidatePos+2] == 'G')
                refContext.contextF = 1;
            else
                refContext.contextF = 2;
            // reverse
            if (candidatePos-1 >= 0 && reference[candidatePos-1] == 'C')
                refContext.contextR = 0;
            else if (candidatePos-2 >=0 && reference[candidatePos-2] == 'C')
                refContext.contextR = 1;
            else
                refContext.contextR = 2;

            if (false) //candidatePos + startCoord == 1000 )
            {
                std::cout << "Position: " << candidatePos + startCoord << " cov: " << positionCoverage << " refAllele: " << refAllele <<  std::endl;
                std::cout << "length: " << length(mapqsF) << std::endl;
                std::cout << "F:" << std::endl;
                for (unsigned i = 0; i < length(mapqsF); ++i)
                {
                    std::cout << "length: i " << length(mapqsF[i]) << std::endl;
                    for (unsigned j = 0; j < length(mapqsF[i]); ++j)
                    {
                        std::cout << "q: " << qualityStringF[i][j] << "  mapq: " << mapqsF[i][j] << "  q: " << static_cast<double>(ordValue(qualityStringF[i][j])-33) <<  std::endl;
                    }
                }
                std::cout << "R:" << std::endl;
                for (unsigned i = 0; i < length(mapqsR); ++i)
                {
                    std::cout << "length: i " << length(mapqsR[i]) << std::endl;
                    for (unsigned j = 0; j < length(mapqsR[i]); ++j)
                        std::cout << "q: " << qualityStringR[i][j] << "  mapq: " << mapqsR[i][j] << "  q: " << static_cast<double>(ordValue(qualityStringR[i][j])-33) << std::endl;
                }
            }
#ifdef CALL_PROFILE
            double timeStamp = sysTime();
#endif

            doBsCalling(countF, countR, qualityStringF, qualityStringR, mapqsF, mapqsR, originStringF, originStringR, refContext, methOptions, options, meth);
#ifdef CALL_PROFILE
            Times::instance().time_doBsCalling += (sysTime() - timeStamp);
#endif


            // for the beginning: write into snp output file
            if ((meth.genotypeCalled && ((meth.genotype>>2) != refContext.refAllele || (meth.genotype%4) != refContext.refAllele) ) ||  // genotype different than ref was called
                 meth.bsCalled || refAllele == 'C' || refAllele == 'G') //(methOptions.outputCandidates & (refAllele == 'C' || refAllele == 'G') ))
                writeMeth(vcfStream, bedStream, meth, qualityStringF, qualityStringR, refContext, positionCoverage, methOptions, options);
        }
    }
    CharString chrPrefix = "chr"; // should check if "chr" is already part of chromosome names (usually not)
    if(options._debugLevel>1) std::cout <<"Finished scanning window.\n"<<std::flush;
}




#endif
