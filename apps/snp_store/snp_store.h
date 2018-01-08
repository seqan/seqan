 /*==========================================================================
  SNP Calling routine of RazerS - Fast Read Mapping with Controlled Loss Rate
                   http://www.seqan.de/projects/razers.html

 ============================================================================
  Copyright (C) 2008 by Anne-Katrin Emde

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your options) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ==========================================================================*/

#ifndef SEQAN_HEADER_CALLSNPS_H
#define SEQAN_HEADER_CALLSNPS_H

// TODO(holtgrew): Using a gff_io module would increase simplicity below.

#include <iostream>
#include <fstream>
#include <cmath>

#include <seqan/misc/svg.h>
#include <seqan/stream.h>
#include <boost/math/special_functions/fpclassify.hpp>

#ifdef CORRECTED_HET
#include <boost/math/distributions.hpp>
#endif

namespace seqan
{

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
};

template<>
struct FragmentStoreConfig<SnpStoreGroupSpec_> :
    public FragmentStoreConfig<>
{
    typedef Dependent<>	TReadSeqStoreSpec;
    typedef Owner<>		TAlignedReadTagStoreSpec;
    typedef Dependent<> TReadNameStoreSpec;
};


    template <typename TGPos_>
    struct SimplePosition
    {
        typedef typename MakeSigned_<TGPos_>::Type TGPos;

        TGPos           gBegin;         // begin position in the genome
        unsigned        gseqNo;


    };

    template < bool _HAMMING_ONLY = true >
    struct SNPCallingSpec
    {
        enum { HAMMING_ONLY = _HAMMING_ONLY };              // omit verifying potential matches
    };


    struct TagMaqMethod_;
    typedef Tag<TagMaqMethod_> const MaqMethod;

    struct TagThresholdMethod_;
    typedef Tag<TagThresholdMethod_> const ThresholdMethod;


    template < typename TSpec = SNPCallingSpec<> >
    struct SNPCallingOptions
    {

        TSpec       spec;
        int         _debugLevel;                // level of verbosity
        bool        printVersion;               // print version number
        std::stringstream   programCall;        // stores snpstore program call
        std::string version;                    // version string
        std::string runID;                      // runID needed for gff output

        // input output options
        unsigned    outputFormat;               // 0 (detailed output of all snp candidate positios)
                                                // or 1 (only successful candidates)
        unsigned    inputFormat;                // 0 = razers, 1=eland, 2 = maq
        unsigned    positionFormat;             // position format of mapped read input
                                                // 1..position space

        CharString          genomeFName;        // name of genome file
        String<CharString>  readFNames;         // list of read file names
        String<CharString>  qualityFNames;      // list of quality file names

        CharString          outputSNP;          // name of snp result file
        CharString          outputIndel;        // name of indel result file
        CharString          outputLog;          // name of log output file

        CharString          inputPositionFile;  // name of position analysis input file
        CharString          outputPosition;     // name of position analysis output file
        CharString          outputCNV;          // name of cnv result file

        bool        showQualityStrings;         // output ascii qualities in SNP output

        // general parameters/options
        unsigned    maxPile;                    // keep at most maxPile many reads mapped to exact same position
        bool        laneSpecificMaxPile;        // do pile up correction on file by file basis (instead of on merged read set)
        bool        orientationAware;           // do pile up correction orientation aware

        bool        storeReadNames;             // store read names

        int         minMapQual;                 // min. mapping quality of reads parsed from SAM file
        bool        keepCigars;                 // remember cigar string of each match
        bool        keepSuboptimalReads;        // keep suboptimal read matches
        bool        keepMultiReads;             // keep multiply mapped reads
        bool        dontClip;                   // dont apply clip tags

        bool        realign;                    // do realignment
        int         realignAddBorder;           // add flanking bases to reference in realignment (0 seems to work best)
        int         minClippedLength;           // discard read if it is not at least minClippedLength bases long after clipping
        bool        clipTagsInFile;             // helper to remember whether there were clip tags in file
        bool        softClipTagsInFile;         // helper to remember whether there was soft clipping in file

        int         asciiQualOffset;            // how to get quality values from ascii, usually subtract 33
        unsigned char compMask[5];              // for comparing nucleotides
        String<char> toIupac;                   // IUPAC code for het snp calls


         // SNP calling related
        unsigned    method;                     // 0 = threshold method, 1 = Maq method
        bool        useBaseQuality;             // use base quality instead of min{base quality,mapping,avg read quality}
        unsigned    minCoverage;                // min depth at variant calling positions
        int         forceCallCount;             // force variant call if there are at least this many mutationas observed at the candidate site
        double      newQualityCalibrationFactor;  // experimental.. downweighting read quality if mapped with many errors..
        double      minExplainedColumn;         // if two most frequent bases don't make up at least this fraction of the alignment column --> discard snp, most likely noise
        unsigned    minDifferentReadPos;        // number of different read positions that need to support the variant call
        unsigned    excludeBorderPos;           // exclude this many read positions at read borders when looking at minDifferentReadPos


        // threshold method related
        float       avgQualT;                   // min. average quality value of variant read base
        float       percentageT;                // min. percentage of variant read base in alignment column
        unsigned    minMutT;                    // min. count of variant read base
        float       snpHetMax;                  // SNP is called as homozygote if variant percentage in alignment column is greater than snpHetMax

        // Maq-method related
        double      hetRate;                    // heterozygote rate
        int         numHaplotypes;              // number of haplotypes (always two..)
        double      theta;                      // theta parameter modeling error dependency in Maq model
        double      eta;                        // eta parameter of Maq model
        String<long double> cnks;               // precomputed table needed for homozygote prob computation
        String<long double> fks;                // precomputed table needed for homozygote prob computation
        String<long double> hetTable;           // precomputed table needed for heterozygote prob computation
        double      priorHetQ;                  // het. rate transformed to quality value

        // branching process corrected allele distribution (Heinrich, Krawitz, 2011)
        int         amplificationCycles;        // number of cycles during DNA amplification
        double      amplificationEfficiency;    // efficiency of amplification
        int         initialN;                   // number of DNA fragments in initial pool
        double      meanAlleleFrequency;        // mean allele frequency of reference allele at het SNP positions
        String<long double> hetTable2;          // same as hetTable but with corrected variance for amplification bias
        bool        correctedHetTable;          // use corrected het table
        bool        printHetTable;              // print corrected het table


        // indel calling related
        int         indelDepthMinOverlap;       // min. overlap of a read with a candidate indel to be considered
        int         maxPolymerRun;              // max. length of homopolymer run for indel to be callable
        bool        bothIndelStrands;           // indel needs to be observed on reads of both strands
        int         indelQualityThreshold;      // min. average quality of indel-neighboring read bases
        float       indelPercentageT;           // min. percentage of indel in alignment columns
        unsigned    indelCountThreshold;        // min. count of indels
        unsigned    indelWindow;                // heuristic to merge neighboring indels in non-realigned data
        float       indelHetMax;                // max percentage of indel-supporting reads for indel to be called heterozygous

        // misc/helpers
        unsigned    maxHitLength;               // helper to remember max. alignment length
        unsigned    minCoord;                   // current min. read mapping coordinate observed
        unsigned    maxCoord;                   // current max. read mapping coordinate observed

        unsigned    windowSize;                 // genomic window size for read parsing
        unsigned    windowBuff;                 // reads within windowBuff base pairs of current window are also kept (-> overlapping windows)

        // cnv calling related // not in use
        unsigned    expectedReadsPerBin;
        unsigned    expectedReadsSD;
        unsigned    cnvWindowSize;

        SNPCallingOptions()
        {

            _debugLevel = 0;
            printVersion = false;
            programCall << "";
            runID = ""; //

            outputFormat = 0;
            inputFormat = 0;
            positionFormat = 1;
            genomeFName = "";
            readFNames = "";
            qualityFNames = "";
            showQualityStrings = true;
            outputSNP = "";
            inputPositionFile = "";
            outputPosition = "";
            outputLog = "";
            outputIndel = "";
            outputCNV = "";

            dontClip = false;
            keepCigars = false;
            keepSuboptimalReads = false;
            keepMultiReads = false;
            minMapQual = 1;

            asciiQualOffset = 33;
            storeReadNames = false;

            maxPile = 1;
            laneSpecificMaxPile = true;
            orientationAware = false;
            realign = false;
            realignAddBorder = 0;


            for (unsigned i = 0; i < 5; ++i)
                compMask[i] = 1 << i;
            //compMask[4] = 0;
            toIupac = "AMRWMCSYRSGKWYKT";

            minClippedLength = 10;
            clipTagsInFile = false;
            softClipTagsInFile = false;

            // SNP calling related
            method = 1;                 // Maq-method is default
            forceCallCount = 10;
            minCoverage = 5;
            minDifferentReadPos = 0;
            excludeBorderPos = 0;
            minExplainedColumn = 0.8;   //


            // threshold-method related
            avgQualT = 10;
            percentageT = (float)0.25;
            minMutT = 3;
            snpHetMax = (float)0.8;
            useBaseQuality = true;      //

            // Maq-method related
            hetRate = 0.001;            // Maq default values
            theta = 0.85;
            eta = 0.03;
            priorHetQ = 0;              // will be set during het table computation
            numHaplotypes = 2;          // only 2 works...

            // amplification bias distribution
            printHetTable = false;
            correctedHetTable = false;
            amplificationCycles = 18;           // realistic values (pers. communication Krawitz)
            amplificationEfficiency = 0.3;
            initialN = 10;
            meanAlleleFrequency = 0.51;
            newQualityCalibrationFactor = 0.0;  // off..

            // indel-calling related
            maxPolymerRun = 100;        // off
            bothIndelStrands = false;
            indelQualityThreshold = 1;  // basically off.. doesn't seem to be very helpful anyway (room for improvement)

            indelDepthMinOverlap = 0;
            indelPercentageT = (float) 0.25;
            indelCountThreshold = 3;
            indelWindow = 0;            // off
            indelHetMax = 0.70f;

            windowSize = 1000000;
            windowBuff = 70;
            minCoord = std::numeric_limits<unsigned>::max();
            maxCoord = 0;
            maxHitLength = 1;

            // unused cnv options
            expectedReadsPerBin = 125;  // unused anyway..
            expectedReadsSD = 25;
            cnvWindowSize = 1000;


       }
    };




//////////////////////////////////////////////////////////////////////////////
// Typedefs

    // definition of a Read match
    template <typename TGPos_>
    struct MappedReadMatch
    {
        typedef typename MakeSigned_<TGPos_>::Type TGPos;

//      TGPos       Batch   gBegin;         // begin position of the match in the genome            --> beginPos
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

//////////////////////////////////////////////////////////////////////////////
// Definitions


// sort operators
    template <typename TMatches, typename TMatchQualities>
    struct LessGStackMQ :
        public ::std::binary_function < typename Value<TMatches>::Type, typename Value<TMatchQualities>::Type, bool >
    {
        TMatchQualities &qualStore;

        LessGStackMQ(TMatchQualities &_qualStore):
            qualStore(_qualStore) {}

        inline bool operator() (
            typename Value<TMatches>::Type const &a,
            typename Value<TMatches>::Type const &b) const
        {
            typedef typename Value<TMatches>::Type TMatch;


            // contig number
            if (a.contigId < b.contigId) return true;
            if (a.contigId > b.contigId) return false;

            // begin position
            typename TMatch::TPos ba = _min(a.beginPos, a.endPos);
            typename TMatch::TPos bb = _min(b.beginPos, b.endPos);

            if (ba < bb) return true;
            if (ba > bb) return false;

            // end position
            typename TMatch::TPos ea = _max(a.beginPos, a.endPos);
            typename TMatch::TPos eb = _max(b.beginPos, b.endPos);

            if (ea < eb) return true;
            if (ea > eb) return false;

            // quality
            if (a.id == TMatch::INVALID_ID) return false;
            if (b.id == TMatch::INVALID_ID) return true;

            if (qualStore[a.id].score > qualStore[b.id].score) return true;
            if (!(qualStore[a.id].score >= qualStore[b.id].score)) return false;

            if (qualStore[a.id].errors < qualStore[b.id].errors) return true;
            if (qualStore[a.id].errors > qualStore[b.id].errors) return false;

            return a.id < b.id;
        }
    };



    template <typename TPosLen>
    struct LessPosLen : public ::std::binary_function < TPosLen, TPosLen, bool >
    {
        inline bool operator() (TPosLen const &a, TPosLen const &b) const
        {
            // read number
            if (a.i1 < b.i1) return true;
            if (a.i1 > b.i1) return false;

            return (a.i2 < b.i2);

        }
    };

    template <typename TMatches, typename TMatchQualities>
    struct LessGStackOaMQ :
        public ::std::binary_function < typename Value<TMatches>::Type, typename Value<TMatchQualities>::Type, bool >
    {
        TMatchQualities &qualStore;

        LessGStackOaMQ(TMatchQualities &_qualStore):
            qualStore(_qualStore) {}

        inline bool operator() (
            typename Value<TMatches>::Type const &a,
            typename Value<TMatches>::Type const &b) const
        {
            typedef typename Value<TMatches>::Type TMatch;

            // contig number
            if (a.contigId < b.contigId) return true;
            if (a.contigId > b.contigId) return false;

            // begin position
            typename TMatch::TPos ba = _min(a.beginPos, a.endPos);
            typename TMatch::TPos bb = _min(b.beginPos, b.endPos);
            if (ba < bb) return true;
            if (ba > bb) return false;

            // end position
            typename TMatch::TPos ea = _max(a.beginPos, a.endPos);
            typename TMatch::TPos eb = _max(b.beginPos, b.endPos);
            if (ea < eb) return true;
            if (ea > eb) return false;

            // orientation
            bool oa = a.beginPos < a.endPos;
            bool ob = b.beginPos < b.endPos;
            if (oa != ob) return oa;

            // quality
            if (a.id == TMatch::INVALID_ID) return false;
            if (b.id == TMatch::INVALID_ID) return true;
            if (qualStore[a.id].score > qualStore[b.id].score) return true;
            if (!(qualStore[a.id].score >= qualStore[b.id].score)) return false;
            if (qualStore[a.id].errors < qualStore[b.id].errors) return true;
            if (qualStore[a.id].errors > qualStore[b.id].errors) return false;
            return a.id < b.id;
        }
    };



    template <typename TReadMatch>
    struct LessId : public ::std::binary_function < TReadMatch, TReadMatch, bool >
    {
        inline bool operator() (TReadMatch const &a, TReadMatch const &b) const
        {
            // genome sequence
            return (a.readId < b.readId);

        }
    };


    // ... to sort matches according to gBegin
    template <typename TReadMatch>
    struct LessGPos : public ::std::binary_function < TReadMatch, TReadMatch, bool >
    {
        inline bool operator() (TReadMatch const &a, TReadMatch const &b) const
        {
            // genome sequence
            if (a.contigId < b.contigId) return true;
            if (a.contigId > b.contigId) return false;

            // begin position
            if (std::min(a.beginPos, a.endPos) < std::min(b.beginPos, b.endPos))
                return true;
            if (std::min(a.beginPos, a.endPos) > std::min(b.beginPos, b.endPos))
                return false;

            // Break tie by read id.
            return a.readId < b.readId;
        }
    };



    // ... to sort matches according to gEnd
    template <typename TReadMatch>
    struct LessGPosEnd : public ::std::binary_function < TReadMatch, TReadMatch, bool >
    {
        inline bool operator() (TReadMatch const &a, TReadMatch const &b) const
        {
            // genome sequence
            if (a.contigId < b.contigId) return true;
            if (a.contigId > b.contigId) return false;

            // end position
            if (std::max(a.endPos,a.beginPos) < std::max(b.endPos,b.beginPos)) return true;
            if (std::max(a.endPos,a.beginPos) > std::max(b.endPos,b.beginPos)) return false;

            return a.readId < b.readId;
        }
    };



    template <typename TReadMatch>
    struct LessGPosEndOa : public ::std::binary_function < TReadMatch, TReadMatch, bool >
    {
        inline bool operator() (TReadMatch const &a, TReadMatch const &b) const
        {
            // genome sequence
            if (a.contigId < b.contigId) return true;
            if (a.contigId > b.contigId) return false;

            // end position
            if (_max(a.endPos,a.beginPos) < _max(b.endPos,b.beginPos)) return true;
            if (_max(a.endPos,a.beginPos) > _max(b.endPos,b.beginPos)) return false;

            // orientation
            bool oa = a.beginPos < a.endPos;
            // bool ob = b.beginPos < b.endPos;
            return oa;

        }
    };

        // ... to sort quality values //
    template <typename TQual>
    struct HigherQ : public ::std::binary_function < TQual, TQual, bool >
    {
        inline bool operator() (TQual const &a, TQual const &b) const
        {
            // quality
            return ordValue(a) > ordValue(b);
        }
    };

//_____________________________________________________________________________________//
/////////////////////////////////////////////////////////////////////////////////////////


// Get reference file names

template<typename TOptions>
int getGenomeFileNameList(StringSet<CharString> & genomeFileNames, TOptions const & options)
{
    std::ifstream file;
    file.open(toCString(options.genomeFName), std::ios_base::in | std::ios_base::binary);
    if(!file.is_open())
        return CALLSNPS_GENOME_FAILED;

    typename DirectionIterator<std::ifstream, Input>::Type fileIter = directionIterator(file, Input());

    CharString nameStr;
    if (*fileIter != '>' && *fileIter != '@')
    {
        // If file does not start with a fasta header --> list of multiple reference genome files.
        if(options._debugLevel >=1)
            std::cout << std::endl << "Reading multiple genome files:" << std::endl;
        /*      //locations of genome files are relative to list file's location
        ::std::string tempGenomeFile(filename);
        size_t lastPos = tempGenomeFile.find_last_of("/\\");
        if (lastPos == tempGenomeFile.npos)
            lastPos = 0;
        else
            ++lastPos;
        ::std::string filePrefix = tempGenomeFile.substr(0,lastPos);*/
        unsigned i = 0;
        for (; !atEnd(fileIter); ++i)
        {
            clear(nameStr);
            skipUntil(fileIter, NotFunctor<IsWhitespace>());
            readUntil(nameStr, fileIter, NotFunctor<IsGraph>());
            appendValue(genomeFileNames,nameStr,Generous());
            if(options._debugLevel >= 2)
                std::cout << "Genome file #" << (i + 1) << ": " << genomeFileNames[length(genomeFileNames) - 1] << std::endl;
        }
        if(options._debugLevel >=1)
            std::cout << i << " genome files total." << std::endl;
    }
    else
    {
        // If file starts with a fasta header --> regular one-genome-file input.
        appendValue(genomeFileNames,options.genomeFName,Generous());
    }

    return 0;
}


/////////////////////////////////////////////////////////////
// read sorted(!) Gff input file containing mapped reads
template <
typename TFile,
typename TFragmentStore,
typename TReadCounts,
typename TCigarStr,
typename TGenome,
typename TGenomeIdMap,
typename TContigPos,
typename TSize,
typename TValue,
typename TOptions
>
int readMatchesFromGFF_Batch(
                             TFile                  &file,
                             TFragmentStore                 &fragmentStore,             // forward/reverse fragmentStore.alignedReadStore
                             TReadCounts                &readCounts,
                             String<Pair<int,int> >         &readClips,
                             StringSet<TCigarStr>           &readCigars,
                             TGenome                    &genome,
                             TGenomeIdMap               &gIdStringToIdNumMap,
                             TSize                  currSeqNo,
                             TContigPos             currentBegin,
                             TContigPos             currentEnd,
                             TValue                 &highestChrId,
                             TOptions               &options,
                             bool setZero = true)
{


    typedef typename TFragmentStore::TAlignedReadStore  TMatches;
    typedef typename Value<TMatches>::Type          TMatch;
    typedef typename TFragmentStore::TAlignQualityStore     TMatchQualities;
    typedef typename Value<TMatchQualities>::Type       TMatchQuality;
    //typedef typename TFragmentStore::TReadSeqStore      TReads;
    //typedef typename Value<TReads>::Type            TRead;
    typedef typename TFragmentStore::TReadStore     TReadStore;
    typedef typename Value<TReadStore>::Type        TReadStoreElement;
    typedef typename Value<TCigarStr>::Type         TCigar;
    //typedef typename Value<TReads>::Type            TRead;
    //typedef typename TFragmentStore::TContigStore       TGenomeSet;
    typedef typename Id<TFragmentStore>::Type       TId;
    //typedef typename Iterator<TMatches,Standard>::Type  TMatchIterator;


    if(length(fragmentStore.readSeqStore)!=length(fragmentStore.alignQualityStore))
    {
        ::std::cerr << "Lengths need to be equal!!\n";
        return 10;
    }
    int readCount = length(fragmentStore.readSeqStore);
    TContigPos genomeLen = length(genome);

    // general stuff that is needed
    typename TGenomeIdMap::const_iterator it;
    unsigned rSeq = readCount;
    Dna5String gInf;
    String<Dna5Q> curr_read;
    TCigarStr tmpCigarStr;
    CharString readTemplate, temp_read;
    CharString readName, temp_str;
    String<int> gAliPos;

    //  bool test = true;

    typename DirectionIterator<std::fstream, Input>::Type fileIter = directionIterator(*file, Input());

    while (!atEnd(fileIter))
    {
        // our razers gff output looks like this:
        //X       razers      read            100919085       100919120       2       +       .       ID=s_3_1_3;unique;mutations=34A;quality=I)IEIIII-7IA>IIIIII07,-%I>)&#029.2-.

        typename std::ifstream::pos_type lineStart = position(fileIter);

        TId contigId;

        // clear temporary variables
        clear(temp_str);
        clear(temp_read);

        unsigned pos= 0;
        unsigned pos2 = 0;
        int clipLeft = 0;
        int clipRight = 0;
        unsigned readCount = 0;
        bool qualityFound = false;
        bool readFound = false;
        char orientation = 'F';
        bool hasIndel = false;
        int editDist = 0;
        int mScore = 100;

        // skip whitespaces just in case (actually there shouldnt be a whitespace at the beginning of a line)
        // and read entry in column 1  --> genomeID
        skipUntil(fileIter, NotFunctor<IsWhitespace>());
        readUntil(temp_str, fileIter, NotFunctor<IsGraph>());

        //check if the genomeID is in our map of relevant genomeIDs, otherwise skip match
        it = gIdStringToIdNumMap.find(temp_str);
        if(options._debugLevel > 1)
            ::std::cout << temp_str << "\t";
        if(it != gIdStringToIdNumMap.end()) contigId = it->second;
        else
        {
            skipLine(fileIter);
            continue;
        }
        if((int)contigId < (int)highestChrId)
        {
            std::cerr << "Read files need to be sorted according to chromosomes in genome file.\n";
            return CALLSNPS_GFF_FAILED;
        }

        highestChrId = contigId;
        if(contigId < currSeqNo)    // havent reached the sequence of interest yet
        {
            skipLine(fileIter);
            continue;
        }

        if(contigId > currSeqNo)    // have passed the seq of interest
        {
            setPosition(fileIter, lineStart);
            break;
        }
        if(setZero) contigId = 0; // if we only store one chromosome at a time

        // skip whitespaces and read entry in column 2
        skipUntil(fileIter, NotFunctor<IsWhitespace>());
        clear(temp_str);
        readUntil(temp_str, fileIter, NotFunctor<IsGraph>());

        if(options._debugLevel > 1)
            ::std::cout << temp_str << "\t";

        // skip whitespaces and read entry in column 3
        skipUntil(fileIter, NotFunctor<IsWhitespace>());
        clear(temp_str);
        readUntil(temp_str, fileIter, NotFunctor<IsGraph>());
        if(options._debugLevel > 1)
            ::std::cout << temp_str << "\t";

        // skip whitespaces and read entry in column 4  --> genomic begin position
        skipUntil(fileIter, NotFunctor<IsWhitespace>());
        clear(temp_str);
        readUntil(temp_str, fileIter, NotFunctor<IsDigit>());
        TContigPos beginPos = 0;
        lexicalCastWithException(beginPos, temp_str);
        beginPos -= options.positionFormat;
        if(beginPos > currentEnd + (TContigPos)options.windowBuff)  // we have passed the relevant match positions
        {
            if(options._debugLevel > 1)
                std::cout  << "gBegin "<< beginPos<<"  of match is too large, seeking "<<lineStart<<"\n";
            setPosition(fileIter, lineStart);
            break;
        }
        if(options._debugLevel > 1)
            ::std::cout << beginPos << "\t";

        // skip whitespaces and read entry in column 5  --> genomic end position
        skipUntil(fileIter, NotFunctor<IsWhitespace>());
        clear(temp_str);
        readUntil(temp_str, fileIter, NotFunctor<IsDigit>());
        TContigPos endPos = 0;
        lexicalCastWithException(endPos, temp_str);

        if(options._debugLevel > 1)
            ::std::cout << endPos << "\t";
        if(endPos + (TContigPos)options.windowBuff < currentBegin)  //we havent reached a relevant read yet
        {
            skipLine(fileIter);
            continue;
        }

        int gMatchLen = endPos - beginPos;
        int rLen = gMatchLen;
        if(endPos > genomeLen)
        {
            setPosition(fileIter, lineStart);
            break;
        }

        // skip whitespaces and read entry in column 6  --> score (percent identity or mapping quality) or a '.'
        skipUntil(fileIter, NotFunctor<IsWhitespace>());
        if(*fileIter == '.')
        {
            mScore = 100;  //not used, but needs to be >= options.minMapQual (default 1)
            ++fileIter;
            if (atEnd(fileIter))
                return CALLSNPS_GFF_FAILED;
        }
        else
        {
            clear(temp_str);
            readUntil(temp_str, fileIter, NotFunctor<IsGraph>());
            double tmp;
            lexicalCastWithException(tmp, temp_str);
            mScore = (int)tmp;
        }
        if(options._debugLevel > 1)
            ::std::cout << mScore << "\t";

        // skip whitespaces and read entry in column 7  --> strand information: '+' or '-'
        skipUntil(fileIter, NotFunctor<IsWhitespace>());
        if (*fileIter == '+')
            orientation = 'F';
        else
            orientation = 'R';
        ++fileIter;
        if (atEnd(fileIter))
            return CALLSNPS_GFF_FAILED;

        // skip whitespaces and read entry in column 8  --> in razers output this is always a '.'
        skipUntil(fileIter, NotFunctor<IsWhitespace>());
        goNext(fileIter);
        if (atEnd(fileIter))
            return CALLSNPS_GFF_FAILED;

        // skip whitespaces and read entry in column 9  --> tags, extra information. in razers output first tag is always "ID"
        skipUntil(fileIter, NotFunctor<IsWhitespace>());
        clear(temp_str);
        readUntil(temp_str, fileIter, NotFunctor<OrFunctor<IsAlphaNum, OrFunctor<EqualsChar<'-'>, EqualsChar<'_'> > > >());
        if(options._debugLevel > 1)
            ::std::cout << temp_str << "\n";
        if(temp_str!="ID")
            ::std::cout << "first feature field should be 'ID'"<<::std::endl;

        // skip the "="
        goNext(fileIter);
        if (atEnd(fileIter))
            return CALLSNPS_GFF_FAILED;

        // read the readID
        CharString readName;
        readUntil(readName, fileIter, OrFunctor<EqualsChar<';'>, IsNewline>());
        if(options._debugLevel > 1)
            ::std::cout << "myID = "<<readName << "\n";
#ifdef SNPSTORE_DEBUG
        bool extraV = true;
#endif
        // cut out the read template from the genomic coordinates
        gInf = infix(genome, beginPos, endPos);
        if (orientation == 'R')
            reverseComplement(gInf);
        readTemplate = gInf;
        if(options._debugLevel > 1) std::cout << readTemplate << "\n";

        // process tags in a loop
        CharString current_tag;
        skipUntil(fileIter, NotFunctor<IsWhitespace>());
        bool multi = false, suboptimal = false/*, unique = true*/, splitRead = false;
        clipLeft = 0; clipRight = 0;
        bool discardRead = false;
        bool first = true;
//        bool softClipped = false;
        int softClippedLeft = 0, softClippedRight = 0;

        //int maxIndelLen = 0;
        clear(gAliPos);
        clear(tmpCigarStr);
        while (!atEnd(fileIter) && *fileIter != '\r' && *fileIter != '\n') // while in same line
        {
            // different tags are separated by ';'  ATTENTION: ascii qualities can contain ';', therefore the tag "quality" MUST be the last tag in a line!!!!!!!!
            while (*fileIter != ';')
            {
                if (*fileIter == '\r')
                    goNext(fileIter);
                if (atEnd(fileIter))
                    return CALLSNPS_GFF_FAILED;
                if (*fileIter == '\n')
                    break;
            }
            goNext(fileIter);
            if (*fileIter == '\n')  // end of line
            {
                skipLine(fileIter);
                break;
            }
            // get the current tag
            clear(current_tag);
            readUntil(current_tag, fileIter, NotFunctor<OrFunctor<IsAlphaNum, OrFunctor<EqualsChar<'-'>, EqualsChar<'_'> > > >());
#ifdef SNPSTORE_DEBUG
            if(options._debugLevel > 1)
                ::std::cout << current_tag << " in features\n";
#endif
            if(current_tag=="quality")
            {
                // add the quality to the read
                qualityFound = true;
                if(!readFound)  //read fragmentStore.alignedReadStore with 0 errors --> read == genomeInfix
                {
                    temp_read = infix(readTemplate,0,gMatchLen);  // without quality values
                    curr_read = temp_read;              // initialized with q40
                    readFound = true;
                }
                for(int i = 0; i < rLen ; ++i) //vorsicht! rLen muss hier schon bekannt sein! (i.e. quality tag is last tag!)
                {
                    goNext(fileIter);  // Skip '='.
                    if (*fileIter != '\n' && *fileIter != '\r')
                    {
                        int tempQual = _max(0, (int)ordValue(*fileIter)-options.asciiQualOffset);
                        assignQualityValue(curr_read[i],tempQual);
#ifdef SNPSTORE_DEBUG
                        if(extraV) ::std::cout << (char)(getQualityValue(curr_read[i]) );
#endif
                    }
                    else
                    {
                        // shouldnt happen
                        if(i != rLen-1 && options._debugLevel > 1) std::cout << curr_read << " gives problems\n";
                        break;
                    }
                }
            }
            else
            {
                // parse other tags
                pos = 0;
                pos2 = 0;
                if (current_tag == "unique") {/*unique = true;*/}
                else if (current_tag == "multi") {multi = true;}
                else if (current_tag == "suboptimal") {suboptimal = true;}
                else if (current_tag == "split") {splitRead = true;}
                else if (current_tag == "clip")
                {
                    options.clipTagsInFile = true;
                    if (*fileIter == '=')
                        skipOne(fileIter);
                    clear(temp_str);
                    readUntil(temp_str, fileIter, NotFunctor<IsDigit>());
                    lexicalCastWithException(clipLeft, temp_str);
                    goNext(fileIter);  // Skip ','
                    clear(temp_str);
                    readUntil(temp_str, fileIter, NotFunctor<IsDigit>());
                    lexicalCastWithException(clipRight, temp_str);
                }
                else if (current_tag == "count")
                {
                    if (*fileIter == '=')
                        skipOne(fileIter);
                    clear(temp_str);
                    readUntil(temp_str, fileIter, NotFunctor<IsDigit>());
                    lexicalCastWithException(readCount, temp_str);
                }
                else if (current_tag == "read")
                {
                    clear(curr_read);
                    readFound = true;
                    if (*fileIter == '=')
                        skipOne(fileIter);
                    readUntil(readName, fileIter, OrFunctor<EqualsChar<';'>, IsNewline>());
                    if (mScore != 100)
                        editDist = (int)((length(curr_read) * ((100.0 - mScore + 0.001)/100.0)));
                }
                else if (current_tag == "mutations")
                {
                    if (!readFound)
                    {
                        temp_read = infix(readTemplate,0,gMatchLen);
                        curr_read = temp_read;
                        readFound = true;
                    }
                    while (*fileIter == ',' || *fileIter == '=') // and add the mutated positions (misfragmentStore.alignedReadStore and insertions in read)
                    {
                        goNext(fileIter);
                        if (atEnd(fileIter))
                            return CALLSNPS_GFF_FAILED;
                        clear(temp_str);
                        readUntil(temp_str, fileIter, NotFunctor<IsDigit>());
                        lexicalCastWithException(pos, temp_str);
                        curr_read[pos - 1] = (Dna5)*fileIter;
                        skipOne(fileIter);
                        ++editDist;
                    }
                }
                else if (current_tag == "cigar")
                {
                    int alignLength = 0;
                    pos = 0; pos2 = 0;
                    int gPos = 0;
                    readFound = true;
                    while (*fileIter != ';' && *fileIter != '\r' && *fileIter != '\n')
                    {
                        if (*fileIter == '=')
                        {
                            ++fileIter;
                            if (atEnd(fileIter))
                                return CALLSNPS_GFF_FAILED;
                        }
                        clear(temp_str);
                        readUntil(temp_str, fileIter, NotFunctor<IsDigit>());
                        lexicalCastWithException(pos2, temp_str);
                        if (*fileIter == 'M')
                        {
                            unsigned k= 0;
                            while(k<pos2)
                            {
                                appendValue(gAliPos,gPos,Generous());
                                ++gPos;
                                ++k;
                            }
                            appendValue(tmpCigarStr,TCigar('M',pos2));
                            alignLength += pos2;
                            pos2 += pos;
                            append(temp_read,infix(readTemplate,pos,pos2));
                            pos = pos2;
                            skipOne(fileIter);
                            first = false;
                            continue;
                        }
                        else if (*fileIter == 'I')
                        { //insertion in the read
                            //(*mIt).editDist += pos2; will be increased in mutations loop
                            unsigned k= 0;
                            while(k<pos2)
                            {
                                appendValue(gAliPos,-gPos,Generous());//no genome positions are used up
                                ++k;
                            }
                            appendValue(tmpCigarStr,TCigar('I',pos2));
                            for(unsigned f = 0; f < pos2; ++f)
                                appendValue(temp_read, 'A');  // will be replaced with correct base in "mutations" loop
                            skipOne(fileIter);
                            rLen += pos2;
                            hasIndel = true;
                            first = false;
                //            if(maxIndelLen < pos2) maxIndelLen = pos2;
                            continue;
                        }
                        else if (*fileIter == 'D')
                        { //there is a deletion in the read
                            unsigned k= 0;
                            while(k<pos2)
                            {
                                ++gPos;
                                ++k;
                            }
                            editDist += pos2;
                            appendValue(tmpCigarStr,TCigar('D',pos2));
                            alignLength += pos2;
                            pos += pos2;
                            rLen -= pos2;
                            skipOne(fileIter);
                            hasIndel = true;
                            first = false;
               //             if(maxIndelLen < pos2) maxIndelLen = pos2;
                            continue;
                        }
                        else if (*fileIter == 'S')
                        {
                            if (first)
                                softClippedLeft = pos2;
                            else
                                softClippedRight = pos2;
                            unsigned k= 0;
                            while(k<pos2)
                            {
                                appendValue(gAliPos,gPos,Generous());
                                ++gPos;
                                ++k;
                            }
                            appendValue(tmpCigarStr,TCigar('S',pos2));
                            alignLength += pos2;
                            pos2 += pos;
//                            if(first)
                                append(temp_read,infix(readTemplate,pos,pos2));
                            pos = pos2;
                            skipOne(fileIter);
                            first = false;
                            options.softClipTagsInFile = true;
                            continue;
                        }
                    }
                    if(alignLength != endPos - beginPos)
                    {
                        std::cerr << "WARNING! Read "<<readName<<": cigar alignment length does not match genome coordinates. Discarding read.."<<std::endl;
                        //std::cout << "align length = " << alignLength << " endPos=" << endPos << " beginPos=" << beginPos << std::endl;
                        discardRead = true;
                    }
                    curr_read = temp_read;
                }
            }
            if (qualityFound)
            {
                skipLine(fileIter);
                break;
            }
        }
        if (options._debugLevel>0&&(rSeq%1000000)==0) std::cout <<rSeq<<".."<<std::flush;
        if(!discardRead && (mScore >= options.minMapQual && (!multi || options.keepMultiReads) && (!suboptimal || options.keepSuboptimalReads)))// && (!((*mIt).hasIndel==1 && options.hammingOnly)))
        {
            if(!readFound)
            {   //neither quality nor read sequence found
                if(options._debugLevel>1)::std::cout << "neither quality nor read sequence found editDist = " << editDist <<"\n";
                temp_read = infix(readTemplate,0,gMatchLen);
                curr_read = temp_read;
            }
            // make sure softClipping is taken into account in clip tags
            if(options.dontClip) // only soft clipping if other clipping is switched off
            {
                    clipRight = softClippedRight;
                    clipLeft = softClippedLeft;
            }
            else
            {
                clipRight = _max(clipRight,softClippedRight);
                clipLeft = _max(clipLeft,softClippedLeft);
            }

            if(clipLeft + clipRight > (int)length(curr_read) - (int)options.minClippedLength)
            {
                if (options._debugLevel>1) std::cout <<"Discarding read "<<readName<<", too short after clipping.."<<std::endl;
                skipUntil(fileIter, NotFunctor<IsWhitespace>());
                continue;
            }
            if(options.realign && splitRead)
            {
                if(endPos-beginPos > (int)((float)length(curr_read)*1.5))
                {
                    if (options._debugLevel>1) std::cout <<"Discarding split read "<<readName<<", deletion too large.."<<std::endl;
                    skipUntil(fileIter, NotFunctor<IsWhitespace>());
                    continue;
                }
            }
#ifdef READ_NAME_AWARE
            if(!options.storeReadNames) clear(readName);
            TId readId;
            if(options.storeReadNames && !getIdByName(fragmentStore.readNameStore, readName, readId, fragmentStore.readNameStoreCache))
            {
                readId = length(fragmentStore.readSeqStore);
                appendValue(fragmentStore.readSeqStore,curr_read,Generous());
                appendValue(fragmentStore.readNameStore, readName, Generous());
            }
            else
            {

            }
#else
            TId readId = length(fragmentStore.readSeqStore);
            appendValue(fragmentStore.readSeqStore,curr_read,Generous());
            if(!options.storeReadNames) clear(readName);
            appendValue(fragmentStore.readNameStore, readName, Generous());
#endif

#ifdef SNPSTORE_DEBUG
            if(clipLeft + clipRight > 76 )
                ::std::cerr << "clipLeft = " << clipLeft << " clipRight = "<<clipRight << "\n";
#endif

            if(options._debugLevel > 1)
                ::std::cout<<fragmentStore.readSeqStore[rSeq]<<" with edit="<<editDist<<" at position "<< beginPos <<"\n";

            if(endPos - beginPos > (TContigPos)options.maxHitLength)
                options.maxHitLength = endPos - beginPos;

            // remember min and max positions seen
            if(beginPos < (TContigPos)options.minCoord || options.minCoord == std::numeric_limits<unsigned>::max()) options.minCoord = (unsigned)beginPos;
            if(endPos > (TContigPos)options.maxCoord) options.maxCoord =  (unsigned)endPos;

            // create match m
            TMatch m;
            m.id = readId; //length(fragmentStore.alignedReadStore);
            if(orientation == 'F')
            {
                m.beginPos = beginPos;
                m.endPos = endPos;
            }
            else
            {
                m.beginPos = endPos;
                m.endPos = beginPos;
            }
            m.contigId = contigId;
            m.readId = m.id;

            // corresponding match quality attributes are stored in q
            TMatchQuality q;
            q.errors = (char)editDist;
            q.score = (char) 0;
            if(options._debugLevel > 1)
            {
                if(splitRead && hasIndel)
                    std::cout << "has indel!\n"; //TODO: how should split reads be treated in realignment?
            }
            if(!options.realign && splitRead && length(tmpCigarStr)<=3) hasIndel = false;
            if(splitRead) clipLeft = clipRight = 0;
            if(hasIndel)
                q.pairScore = 1;
            else
                q.pairScore = 0;

            typename Value<TReadStore>::Type r;
            r.matePairId = TReadStoreElement::INVALID_ID;
            if(readCount > 0) appendValue(readCounts, readCount, Generous());

            appendValue(fragmentStore.readStore, r, Generous());
            appendValue(fragmentStore.alignedReadStore, m, Generous());
            appendValue(fragmentStore.alignQualityStore, q, Generous());
            appendValue(readClips,Pair<int,int>(clipLeft,clipRight));

            if(!splitRead)
                clear(tmpCigarStr); // split reads store their cigar string explicitly

            appendValue(readCigars,tmpCigarStr);
            ++rSeq;
            if(options._debugLevel > 1)
            {
                ::std::cout<<"Parsed: id= " <<m.readId<<" name="<<readName<<"="<<curr_read<<" with edit="<<editDist<<" at position "<< beginPos<<"\n";
                ::std::cout << "mScore=" << mScore << " m.beginPos=" << m.beginPos << " m.endPos=" << m.endPos << std::endl;
                if(q.pairScore==1) ::std::cout << "indel! pairScore=" << q.pairScore <<std::endl;
                if(q.pairScore==0) ::std::cout << "no indel! pairScore=" << q.pairScore <<std::endl;

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

        skipUntil(fileIter, NotFunctor<IsWhitespace>());
    }
    if(options._debugLevel > 0)
        ::std::cout << ::std::endl << "Parsed "<<length(fragmentStore.alignedReadStore)<<" matches of "<<length(fragmentStore.readSeqStore)<<" reads." << ::std::endl;


    return 0;
}

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
// read sorted(!) Gff input file containing mapped reads
template <
typename TFragmentStore,
typename TReadCounts,
typename TCigarStr,
typename TGenome,
typename TGenomeIdMap,
typename TContigPos,
typename TSize,
typename TValue,
typename TOptions
>
int readMatchesFromSamBam_Batch(
                         BamFileIn                  & bamFileIn,
                         BamAlignmentRecord         &record,
                         TFragmentStore             &fragmentStore,             // forward/reverse fragmentStore.alignedReadStore
                         TReadCounts                &readCounts,
                         String<Pair<int,int> >     &readClips,
                         StringSet<TCigarStr>       &readCigars,
                         TGenome                    &genome,
                         TGenomeIdMap               &gIdStringToIdNumMap,
                         TSize                  currSeqNo,
                         TContigPos             currentBegin,
                         TContigPos             currentEnd,
                         TValue                 &highestChrId,
                         TOptions               &options,
                         bool                   firstCall,
                         bool                   setZero = true)
{


    typedef typename TFragmentStore::TAlignedReadStore  TMatches;
    typedef typename Value<TMatches>::Type          TMatch;
    typedef typename TFragmentStore::TAlignQualityStore     TMatchQualities;
    typedef typename Value<TMatchQualities>::Type       TMatchQuality;
    typedef typename TFragmentStore::TReadSeqStore      TReads;
    typedef typename Value<TReads>::Type            TRead;
    typedef typename TFragmentStore::TReadStore     TReadStore;
    typedef typename Value<TReadStore>::Type        TReadStoreElement;
    //typedef typename Value<TCigarStr>::Type         TCigar;
    typedef typename Value<TReads>::Type            TRead;
    //typedef typename TFragmentStore::TContigStore       TGenomeSet;
    typedef typename Id<TFragmentStore>::Type       TId;
    //typedef typename Iterator<TMatches,Standard>::Type  TMatchIterator;


    if(length(fragmentStore.readSeqStore)!=length(fragmentStore.alignQualityStore))
    {
        ::std::cerr << "Lengths need to be equal!!\n";
        return 10;
    }
    int readCount = length(fragmentStore.readSeqStore);
    TContigPos genomeLen = length(genome);

    // general stuff that is needed
    typename TGenomeIdMap::const_iterator it;
    unsigned rSeq = readCount;
    Dna5String gInf;
    String<Dna5Q> curr_read;
    TCigarStr tmpCigarStr;
    CharString readTemplate, temp_read;
    CharString readName, temp_str;
//    String<int> gAliPos;
//  BamAlignmentRecord record;
    BamHeader header;
    int res = 0;
    if (firstCall)
        readHeader(header, bamFileIn);
    //  bool test = true;
    while (!atEnd(bamFileIn))
    {

        //typename std::ifstream::pos_type lineStart = ((reader)._file).tellg();
        //lineStart = lineStart - (std::ifstream::pos_type) 1;

        // read next record unless current one has not been handled yet
        res = 0;
        if (empty(record.qName))
            readRecord(record, bamFileIn);
        //_printRecord(std::cout, record);
        if(res != 0)
        {
            ::std::cerr << "Something wrong with SAM/BAM file?" << std::endl;
            return 2;
        }

        if ( hasFlagUnmapped(record) || empty(record.cigar) || (!options.keepSuboptimalReads && hasFlagSecondary(record)))
        {
            //std::cout << "Read " << record.qName << " has Flag=" << record.flag << std::endl;
            clear(record); continue;
        }

        TId contigId;

        // clear temporary variables
        clear(temp_str);
        clear(temp_read);

//      unsigned pos= 0;
//      unsigned pos2 = 0;
        int clipLeft = 0;
        int clipRight = 0;
        unsigned readCount = 0;
//      bool qualityFound = false;
//      bool readFound = false;
        char orientation = 'F';
        bool hasIndel = false;
        int editDist = 0;
        int mScore = 100;

        // problem: contigNameStore is not filled --> what is rID set to? can i get the original string?? -> extra BamAlignmentRecord spec?
        //check if the genomeID is in our map of relevant genomeIDs, otherwise skip match
        it = gIdStringToIdNumMap.find(contigNames(context(bamFileIn))[record.rID]);
        if(options._debugLevel > 1)
            ::std::cout << record.rID << "\t";
        if(it != gIdStringToIdNumMap.end())
            contigId = it->second;
        else // check if there is "chr"
        {
            if(prefix(contigNames(context(bamFileIn))[record.rID],3) == "chr")
                it = gIdStringToIdNumMap.find(suffix(contigNames(context(bamFileIn))[record.rID],3));
            if(it != gIdStringToIdNumMap.end())
                contigId = it->second;
            else
            {
                CharString temp = "chr";
                append(temp,contigNames(context(bamFileIn))[record.rID]);
                it = gIdStringToIdNumMap.find(temp);
                if(it != gIdStringToIdNumMap.end())
                    contigId = it->second;
                else
                {
                    clear(record);
                    continue;
                }
            }
        }

        if((int)contigId < (int)highestChrId)
        {
            std::cerr << "Read files need to be sorted according to chromosomes in genome file.\n";
            return CALLSNPS_GFF_FAILED;
        }

        highestChrId = contigId;
        if(contigId < currSeqNo)    // havent reached the sequence of interest yet
        {
            clear(record);
            continue;
        }

        if(contigId > currSeqNo)    // have passed the seq of interest
        {
//          streamSeek(reader,lineStart);
//            (reader._file).seekp(lineStart);
            break;
        }
        if(setZero) contigId = 0; // if we only store one chromosome at a time

        // skip whitespaces and read entry in column 2

        TContigPos beginPos = record.beginPos;
        if(beginPos > currentEnd + (TContigPos)options.windowBuff)  // we have passed the relevant match positions
        {
            if(options._debugLevel > 1)
                std::cout  << "gBegin "<< beginPos<<"  of match is too large\n";//, seeking "<<lineStart<<"\n";
//          streamSeek(reader,lineStart);
//            (reader._file).seekp(lineStart);
            break;
        }
        if(options._debugLevel > 1)
            ::std::cout << beginPos << "\t";

        //_printRecord(std::cout, record);

        // need to calculate endPos
        TContigPos endPos = beginPos; // + length of alignment! sum of Ms and Ds in cigar string
        // check if cigar string has indels
        bool first = true;
        bool softClipped = false;
        int softClippedLeft = 0, softClippedRight = 0;
        for(unsigned j=0; j < length(record.cigar); ++j)
        {
            if(record.cigar[j].operation == 'M')
            {
                first = false;
                endPos += (TContigPos) record.cigar[j].count;
            }
            if(record.cigar[j].operation == 'D' || record.cigar[j].operation == 'N')
            {
                first = false;
                endPos += (TContigPos) record.cigar[j].count;
                hasIndel = true;
            }
            if(record.cigar[j].operation == 'I')
            {
                first = false;
                hasIndel = true;
            }
            if(record.cigar[j].operation == 'S')
            {
                softClipped = true;
                if(first)
                {
                    softClippedLeft = record.cigar[j].count;
                 // dont apply tags yet, wait until after pile up correction
                 // instead act as if clipped positions were matched characters --> change mapping cooridnates (will be adjusted again in clipReads function)
                    beginPos = _max(0,beginPos-softClippedLeft);
                }
                else
                {
                    softClippedRight = record.cigar[j].count;
                 // dont apply tags yet, wait until after pile up correction
                 // instead act as if clipped positions were matched characters --> change mapping cooridnates (will be adjusted again in clipReads function)
                    endPos += (TContigPos) record.cigar[j].count;
                }
                first = false;
                options.softClipTagsInFile = true;
            }

        }

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

        if (hasFlagRC(record)) orientation = 'R';
        else orientation = 'F';


        CharString readName = record.qName;
        if(options._debugLevel > 1)
            ::std::cout << "myID = "<<readName << "\n";
#ifdef SNPSTORE_DEBUG
        bool extraV = true;
#endif

        bool multi = false, splitRead = false;
        bool suboptimal = hasFlagSecondary(record);

        clipLeft = 0; clipRight = 0;
        TRead curr_read = record.seq;
        for(unsigned j = 0; j < length(record.qual); ++j)
        {
            int tempQual = _max(0,(int)ordValue(record.qual[j])-options.asciiQualOffset);
            assignQualityValue(curr_read[j],tempQual);
        }
        if (orientation == 'R')
            reverseComplement(curr_read);


        interpretBamTags(record.tags,editDist,multi,clipLeft,clipRight,options);

        // make sure softClipping is taken into account in clip tags
        if(options.dontClip) // only soft clipping if other clipping is switched off
        {
                clipRight = softClippedRight;
                clipLeft = softClippedLeft;
        }
        else
        {
            clipRight = _max(clipRight,softClippedRight);
            clipLeft = _max(clipLeft,softClippedLeft);
        }
        if (orientation == 'R')
        {
            int temp = clipLeft;
            clipLeft = clipRight;
            clipRight = temp;
        }



        if (options._debugLevel>0&&(rSeq%1000000)==0) std::cout <<rSeq<<".."<<std::flush;
        if( /*length(curr_read)> 30 && */ mScore >= options.minMapQual && (!multi || options.keepMultiReads) && (!suboptimal || options.keepSuboptimalReads))// && (!((*mIt).hasIndel==1 && options.hammingOnly)))
        {
            if(empty(curr_read))
            {   //read sequence not found
                if(options._debugLevel>1)::std::cout << "neither quality nor read sequence found editDist = " << editDist <<"\n";
                return 1;
            }
            if(clipLeft + clipRight > (int)length(curr_read) - (int)options.minClippedLength)
            {
                if (options._debugLevel>1) std::cout <<"Discarding read "<<readName<<", too short after clipping.."<<std::endl;
                clear(record);
                continue;
            }
            if(options.realign && splitRead)
            {
                if(endPos-beginPos > (int)((float)length(curr_read)*1.5))
                {

                    if (options._debugLevel>1) std::cout <<"Discarding split read "<<readName<<", deletion too large.."<<std::endl;
                    // should just not consider these reads for realignment instead of kicking them out entirely
                }
            }
#ifdef READ_NAME_AWARE
            if(!options.storeReadNames) clear(readName);
            TId readId;
            if(options.storeReadNames && !getIdByName(fragmentStore.readNameStore, readName, readId, fragmentStore.readNameStoreCache))
            {
                readId = length(fragmentStore.readSeqStore);
                appendValue(fragmentStore.readSeqStore,curr_read,Generous());
                appendValue(fragmentStore.readNameStore, readName, Generous());
            }
            else
            {

            }
#else
            TId readId = length(fragmentStore.readSeqStore);
            appendValue(fragmentStore.readSeqStore,curr_read,Generous());
            if(!options.storeReadNames) clear(readName);
            appendValue(fragmentStore.readNameStore, readName, Generous());
#endif

#ifdef SNPSTORE_DEBUG
            if(clipLeft + clipRight > 76 )
                ::std::cerr << "clipLeft = " << clipLeft << " clipRight = "<<clipRight << "\n";
#endif

            if(options._debugLevel > 1)
                ::std::cout<<fragmentStore.readSeqStore[rSeq]<<" with edit="<<editDist<<" at position "<< beginPos <<"\n";

            if(endPos - beginPos > (TContigPos)options.maxHitLength)
                options.maxHitLength = endPos - beginPos;

            // remember min and max positions seen
            if(beginPos < (TContigPos)options.minCoord || options.minCoord == std::numeric_limits<unsigned>::max()) options.minCoord = (unsigned)beginPos;
            if(endPos > (TContigPos)options.maxCoord) options.maxCoord =  (unsigned)endPos;

            // create match m
            TMatch m;
            m.id = readId; //length(fragmentStore.alignedReadStore);
            if(orientation == 'F')
            {
                m.beginPos = beginPos;
                m.endPos = endPos;
            }
            else
            {
                m.beginPos = endPos;
                m.endPos = beginPos;
            }
            m.contigId = contigId;
            m.readId = m.id;

            // corresponding match quality attributes are stored in q
            TMatchQuality q;
            q.errors = (char)editDist;
            q.score = (char) mScore;
            if(options._debugLevel > 1)
            {
                if(splitRead && hasIndel)
                    std::cout << "has indel!\n"; //TODO: how should split reads be treated in realignment?
            }
            if(!options.realign && splitRead && length(tmpCigarStr)<=3) hasIndel = false;
            if(splitRead) clipLeft = clipRight = 0;
            if(hasIndel)
                q.pairScore = 1;
            else
                q.pairScore = 0;

            typename Value<TReadStore>::Type r;
            r.matePairId = TReadStoreElement::INVALID_ID;
            if(readCount > 0) appendValue(readCounts, readCount, Generous());

            appendValue(fragmentStore.readStore, r, Generous());
            appendValue(fragmentStore.alignedReadStore, m, Generous());
            appendValue(fragmentStore.alignQualityStore, q, Generous());
            appendValue(readClips,Pair<int,int>(clipLeft,clipRight));

            if(!splitRead)
                clear(tmpCigarStr); // split reads store their cigar string explicitly

            appendValue(readCigars,tmpCigarStr);
            ++rSeq;
            if(options._debugLevel > 1)
            {
                ::std::cout<<"Parsed: id= " <<m.readId<<" name="<<readName<<"="<<curr_read<<" with edit="<<editDist<<" at position "<< beginPos<<"\n";
                ::std::cout << "mScore=" << mScore << " m.beginPos=" << m.beginPos << "m.endPos="<<m.endPos<<std::endl;
                if(q.pairScore==1) ::std::cout << "indel! pairScore=" << q.pairScore <<std::endl;
                if(q.pairScore==0) ::std::cout << "no indel! pairScore=" << q.pairScore <<std::endl;
                if(softClipped) ::std::cout << "SoftClipped!" << std::endl;
                if(orientation == 'R') ::std::cout << "reversed!" << std::endl;
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
//  if(length(fragmentStore.alignedReadStore) > 0)
        ::std::cout << ::std::endl << "Parsed "<<length(fragmentStore.alignedReadStore)<<" matches of "<<length(fragmentStore.readSeqStore)<<" reads." << ::std::endl;


    return 0;
}



///////////////////////////////////////////////////////////////////////////////////////////////7
// simple position stats analysis
template <typename TPositions, typename TGenomeSetSize, typename TOptions>
bool loadPositions(TPositions & positions,
                   ::std::map<CharString, TGenomeSetSize> &gIdStringToIdNumMap,
                   char const * filename,
                   TOptions & options)
{
    ::std::ifstream file;
    file.open(filename,::std::ios_base::in | ::std::ios_base::binary);
    if(!file.is_open())
        return 1;

    CharString chrId;

    int numPos = 0;
    typename ::std::map<CharString, TGenomeSetSize>::const_iterator it;
    unsigned contigId;
    typename DirectionIterator<std::ifstream, Input>::Type fileIter = directionIterator(file, Input());
    while (!atEnd(fileIter))
    {
        skipUntil(fileIter, NotFunctor<IsWhitespace>());

        // Skip whitespaces just in case (actually there shouldnt be a whitespace at the beginning of a line)
        // and read entry in column 1  --> genomeID
        clear(chrId);
        readUntil(chrId, fileIter, IsWhitespace());

        // Check if the genomeID is in our map of relevant genomeIDs, otherwise skip position.
        it = gIdStringToIdNumMap.find(chrId);
        if (options._debugLevel > 1)
            ::std::cout << chrId << "\t";
        if (it != gIdStringToIdNumMap.end())
        {
            contigId = it->second;
        }
        else
        {
            skipLine(fileIter);
            continue;
        }
        SEQAN_ASSERT_GT(length(positions), contigId);
        skipUntil(fileIter, NotFunctor<IsWhitespace>());
        seqan::CharString buffer;
        skipUntil(fileIter, NotFunctor<IsDigit>());
        unsigned pos = 0;
        lexicalCastWithException(pos, buffer);
        pos -= options.positionFormat;
        appendValue(positions[contigId],pos);
        ++numPos;
        skipLine(fileIter);
    }

    if (numPos > 0)
        return 0;
    return 1;
}





#ifdef STDLIB_VS

template<typename TVal>
double
lgamma(TVal x)
{
    // TODO: replace
    x -= 1;
    if(x < 2) return log((double)1);
    else
    {
        double f = log((double)2);
        for(int s = 3; s <= x;++s) f += log((double)s);
        return f;
    }

}

#endif


// function taken from keith b. hall, computation of probs in log-space
template<typename TValue>
inline TValue
logSum(TValue x, TValue y)
{
    // If one value is much smaller than the other, keep the larger value.
    if (x < (y - log(1e200)))
        return y;
    if (y < (x - log(1e200)))
        return x;
    double diff = x - y;
    double retVal;
    if (!std::isfinite((double)exp(diff))) // difference is too large
        return (x > y ? x : y);
    // otherwise return the sum.
    retVal = (double)(y + log((double)(1.0) + exp(diff)));
    return retVal;
}


// this is basically maq's source code translated into seqan
// see Li, H., Ruan, J. & Durbin, R. Mapping short DNA sequencing reads and calling variants
// using mapping quality scores. Genome Res. 2008.
// and http://maq.sourceforge.net
template <typename THomoTable, typename TDependencies, typename TOptions>
void computeCnks(THomoTable & cnks, TDependencies & fks, TOptions & options)
{
    typedef typename Value<THomoTable>::Type TValue;

    String<TValue> sum_a, beta, q_c, temp, lFks, lC;
    resize(sum_a,257);
    resize(beta,256);
    resize(q_c,256);
    resize(temp,256);
    resize(lFks,256);
    resize(lC, 256*256);
    resize(fks,256);
    resize(cnks, 256*256*64,0.0); // n<256, k<256, q<64



    fks[0] = lFks[0] = 1.0;
    for (int n = 1; n < 256; ++n)
    {
        fks[n] = pow(options.theta, n) * (1.0 - options.eta) + options.eta;
        lFks[n] = fks[n>>1]; //
    }

    for (int n = 1; n < 256; ++n)
        for (int k = 0; k <= n; ++k)  //
            lC[n<<8|k] = lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1); // for every combination of k errors in n reads,
                                        // (n and k share 16bit in lC)

    for (int q = 1; q < 64; ++q)
    {
        // for all quality values up to 64
        // these are the 'average' values computed from the mapped reads
        double e = pow(10.0, -q/10.0);
        double le = log(e);
        double le1 = log(1.0-e);

        for (int n = 1; n < 256; ++n)
        {
            sum_a[n+1] = 0.0;
            for (int k = n; k >= 0; --k) { // a_k = \sum_{i=k}^n C^n_k \epsilon^k (1-\epsilon)^{n-k}
                sum_a[k] = sum_a[k+1] + expl(lC[n<<8|k] + k*le + (n-k)*le1);
                beta[k] = sum_a[k+1] / sum_a[k];
                if (beta[k] > 0.99) beta[k] = 0.99;
            }

            for (int k = 0; k < n; ++k)                     // c_k
                q_c[k] = -4.343 * lFks[k] * logl(beta[k] / e);
            for (int k = 1; k < n; ++k)
                q_c[k] += q_c[k-1]; // \prod_{i=0}^k c_i

            for (int k = 0; k <= n; ++k)
            {
                temp[k] = -4.343 * logl(1.0 - expl(lFks[k] * logl(beta[k])));
                cnks[q<<16|n<<8|k] = (k > 0 ? q_c[k-1] : 0);// + temp[k];
                cnks[q<<16|n<<8|k] += (boost::math::isnan(temp[k]) ? 0 : temp[k]);
            }

        }
    }

//  std::fstream testf;
//  testf.open("cnks", std::ios_base::out);
//  if(testf.is_open())
//  {
//      for(unsigned f = 0; f < length(cnks); ++f)
//              testf << cnks[f] << " ";
//          testf <<std::endl;
//  }
//  testf.close();

}

#ifdef CORRECTED_HET
using boost::math::normal;

// amplification bias corrected heterozygote probabilities, for diploid only
// fills hetTable structure
// returns het prior in phred scale
template<typename THeteroTable, typename TOptions>
double computeHetTable(THeteroTable & hetTable, TOptions & options)
{
//    std::cout << "Compute het table!!!!" << std::endl;
    typedef typename Value<THeteroTable>::Type TValue;
    double poly_rate;
    int numHaplotypes = 2;//options.numHaplotypes;

//    options.amplificationCycles = 18;
//    options.amplificationEfficiency = 0.3;
//    options.initialN = 10;
//    options.meanAlleleFraction = 0.54;

    resize(hetTable,256*256); // for n,k < 256
    TValue sum_harmo = 0.0;
    for (int k = 1; k <= numHaplotypes - 1; ++k)
        sum_harmo += 1.0 / k;

    // prepare normal distribution with mean n*meanAlleleFraction
    String<normal> distributions;
//    resize(distributions,512);
//    std::cout << "corrections:";
    normal distribution(options.meanAlleleFrequency, 0.2); // dummy
    appendValue(distributions,distribution);
    for (int n = 1; n < 512 ; ++n)
    {
        double mean = options.meanAlleleFrequency * n;
        double standardDev = (double) sqrt((1.0-options.meanAlleleFrequency) * mean);
        long double correction = (long double) sqrt((long double) n*(2.0*(1.0/(1.0+options.amplificationEfficiency))-2.0*pow(1.0+options.amplificationEfficiency,-options.amplificationCycles-1.0)+pow(1.0+options.amplificationEfficiency,-options.amplificationCycles)-1)/(8*options.initialN));
  //      std::cout << n << ": " << correction << "  ";

        // get normal distribution
        normal distribution(mean, standardDev + correction);
        appendValue(distributions,distribution);
    }
    //std::cout << std::endl;


      //  P(1:i+1,i+1)=poisscdf(0.5:(i+1),lambda)-poisscdf(-0.5:(i),lambda);
      //  N(1:i+1,i+1) = normcdf(0.5:(i+1),mu,sigma)-normcdf(-0.5:(i),mu,sigma);
      //  V(1:i+1,i+1) =  normcdf(0.5:(i+1),mu,(sigma+correction))-normcdf(-0.5:(i),mu,(sigma+correction));

       // compute values of normal distribution with corrected variance
    for (int n1 = 0; n1 < 256; ++n1)
    {

        for (int n2 = 0; n2 < 256; ++n2)
        {
            double corrected = 0.0;
            if(n1+n2 != 0)
                corrected = cdf(distributions[n1+n2],n2+0.5) - cdf(distributions[n1+n2],n2-0.5);
            hetTable[n1<<8|n2] = logl(corrected);

        }
    }

    // tests
   // std::cout << "n1=10 , n2=20" << hetTable[10<<8|20] << std::endl;
   // std::cout << "n1=20 , n2=10" << hetTable[20<<8|10] << std::endl;


    poly_rate = options.hetRate * sum_harmo;
    double hetPriorQ = -4.343 * log(2.0 * poly_rate / (1.0 - poly_rate));

    return hetPriorQ;
}
#endif



// maqs heterozygote probabilites
// fills hetTable structure
// returns het prior in phred scale
template<typename THeteroTable, typename TOptions>
double computeHetTable(THeteroTable & hetTable, TOptions & options, MaqMethod & )
{
    typedef typename Value<THeteroTable>::Type TValue;
    double poly_rate;
    int numHaplotypes = 2;//options.numHaplotypes;

    resize(hetTable,256*256); // for n,k < 256
    TValue sum_harmo = 0.0;
    for (int k = 1; k <= numHaplotypes - 1; ++k)
        sum_harmo += 1.0 / k;
    for (int n1 = 0; n1 < 256; ++n1)
    {
        for (int n2 = 0; n2 < 256; ++n2)
        {
            long double sum = 0.0;
            double lC = lgamma(n1+n2+1) - lgamma(n1+1) - lgamma(n2+1); // \binom{n1+n2}{n1}
            for (int k = 1; k <= numHaplotypes - 1; ++k)
            {
                double pk = 1.0 / k / sum_harmo;
                double log1 = log((double)k/numHaplotypes);
                double log2 = log(1.0 - (double)k/numHaplotypes);
                sum += pk * 0.5 * (expl(log1*n2) * expl(log2*n1) + expl(log1*n1) * expl(log2*n2));
            }
            hetTable[n1<<8|n2] = lC + logl(sum);
        }
    }
    poly_rate = options.hetRate * sum_harmo;
    double hetPriorQ = -4.343 * log(2.0 * poly_rate / (1.0 - poly_rate));

    return hetPriorQ;
}



template<typename THomoTable, typename TDependencies, typename TQStrings,typename TVal>
void
getHomoProbs(THomoTable & cnks,
            TDependencies & fks,
            TQStrings & qualitiesForward,
            TQStrings & qualitiesReverse,
            int & best,
            int & secondBest,
            long double & probQ1,
            long double & probQ2, TVal
#ifdef SNPSTORE_DEBUG_CANDPOS
            candidatePos
#endif
            )
{

    //typedef typename Value<THomoTable>::Type TValue;
    typedef typename Value< typename Value<TQStrings>::Type >::Type TQuality;

#ifdef SNPSTORE_DEBUG_CANDPOS
    bool extraV = false; //|| candidatePos < 9335310
    if(candidatePos==118487871) extraV = true;
#endif

    String<double> sumE, sumF;
    resize(sumE,4,0.0);
    resize(sumF,4,0.0);

    for(unsigned i = 0; i < 4; ++i)
    {
        sort(begin(qualitiesForward[i],Standard()),end(qualitiesForward[i],Standard()),HigherQ<TQuality>());
        sort(begin(qualitiesReverse[i],Standard()),end(qualitiesReverse[i],Standard()),HigherQ<TQuality>());
        //compute average q
        double fk = 0.0;
        double qual = 0.0;
#ifdef SNPSTORE_DEBUG_CANDPOS
        if(extraV) std::cout << "F base"<<i<<": " << std::flush;
#endif
        for(unsigned j = 0; j < length(qualitiesForward[i]); ++j)
        {
            qual = static_cast<double>(ordValue(qualitiesForward[i][j])-33);
            if(qual > 30.0) qual = 30.0; /// snp rate 1 in a 1000
            //qual = rescale into regular log
            if(j>=256) fk = fks[255];
            else fk = fks[j];
            sumE[i] += fk * qual;
            sumF[i] += fk;
#ifdef SNPSTORE_DEBUG_CANDPOS
            if(extraV)
            {
                std::cout << sumE[i] << " " << std::flush;
                std::cout << sumF[i] << " " << std::flush;
            }
#endif
        }
#ifdef SNPSTORE_DEBUG_CANDPOS
        if(extraV) std::cout << std::endl;
        if(extraV) std::cout << "R base"<<i<<": " << std::flush;
#endif
        for(unsigned j = 0; j < length(qualitiesReverse[i]); ++j)
        {
            qual = static_cast<double>(ordValue(qualitiesReverse[i][j])-33);
            if(qual > 30.0) qual = 30.0;
            if(j>=256) fk = fks[255];
            else fk = fks[j];
            sumE[i] += fk * qual;
            sumF[i] += fk;
#ifdef SNPSTORE_DEBUG_CANDPOS
            if(extraV)
            {
                std::cout << sumE[i] << " " << std::flush;
                std::cout << sumF[i] << " " << std::flush;
            }
#endif
        }

    }
#ifdef SNPSTORE_DEBUG_CANDPOS
    if(extraV)
    {
        for(unsigned j = 0; j < 256; ++j)
            std::cout << fks[j] << " " << std::flush;
        std::cout << std::endl;

    }
#endif

    best = -1;
    secondBest = -1;
    double bestSum = 0.0;
    double secondBestSum = 0.0;
    for (int j = 0; j < 4; ++j)
    {
        if (sumE[j] > bestSum)
        {
            secondBestSum = bestSum;
            secondBest = best;
            bestSum = sumE[j];
            best = j;
        }
        else
        {
            if (sumE[j] > secondBestSum)
            {
                secondBestSum = sumE[j];
                secondBest = j;
            }
        }
    }
#ifdef SNPSTORE_DEBUG_CANDPOS
    if(extraV) std::cout <<"best="<<best <<" secondbest="<<secondBest << std::flush << std::endl;
#endif

    int qAvgBest = 0, qAvgSecondBest = 0;
    int countBest = 0, countSecondBest = 0;
    if(best != -1)
    {
        // normalized quality of the base with the best weighted sum of qualities
        qAvgBest = (int)(sumE[best]/sumF[best] + 0.5);
        countBest = length(qualitiesForward[best]) + length(qualitiesReverse[best]);
    }
    else countBest = 0;
    if(secondBest != -1)
    {
        qAvgSecondBest = (int)(sumE[secondBest]/sumF[secondBest] + 0.5);
        countSecondBest = length(qualitiesForward[secondBest]) + length(qualitiesReverse[secondBest]);
    }
    else countSecondBest = 0;

    //runter skalieren
    int countTotal = countBest + countSecondBest;
    if (countTotal > 255)
    {
        countBest = int(255.0 * countBest / countTotal + 0.5);
        countSecondBest = int(255.0 * countSecondBest / countTotal + 0.5);
    //  if(extraV) ::std::cout <<  "countBest = " << countBest << "\tcountSecondBest = " << countSecondBest << ::std::endl;
        countTotal = 255;
    }
#ifdef SNPSTORE_DEBUG_CANDPOS
    if(extraV)std::cout << "qAvgBest" <<  qAvgBest << " qAvgSecond"<< qAvgSecondBest << "\n";
    if(extraV)std::cout << "totalCount" <<  countTotal << " countBest"<< countBest << " countSecondBest"<< countSecondBest << "\n";
#endif
    probQ1 = ((countSecondBest > 0) ? sumE[secondBest] : 0);
    probQ1 += (boost::math::isnan(cnks[qAvgSecondBest<<16|countTotal<<8|countSecondBest])) ? 0 : cnks[qAvgSecondBest<<16|countTotal<<8|countSecondBest];
    probQ2 = ((countBest > 0) ? sumE[best] : 0);
    probQ2 += (boost::math::isnan(cnks[qAvgBest<<16|countTotal<<8|countBest])) ? 0 : cnks[qAvgBest<<16|countTotal<<8|countBest];

#ifdef SNPSTORE_DEBUG_CANDPOS
    if(extraV)std::cout << "cnkBest" <<  cnks[qAvgBest<<16|countTotal<<8|countBest] << "  bei cnkindex " <<(qAvgBest<<16|countTotal<<8|countBest)<<"\n";
    if(extraV)std::cout << "cnkSecondBest" <<  cnks[qAvgSecondBest<<16|countTotal<<8|countSecondBest] <<  "  bei cnkindex " <<(qAvgSecondBest<<16|countTotal<<8|countSecondBest)<< "\n";
    if(extraV)std::cout << "probQ1" <<  probQ1 << "\n";
    if(extraV)std::cout << "probQ2" <<  probQ2 << "\n";
#endif

//  if(extraV)
 // {
 //     for(unsigned f = 0; f < length(cnks); ++f)
 //         std::cout << cnks[f] << " ";
 //     std::cout <<std::endl;
 // }


}


///////////////////////////////////////////////////////////////////////////////////////////////////




//looks for mismatches in alignemnt and returns positions with respect to 2nd row (read sequence)
template<typename TAlign, typename TString>
void
getMismatchMutations(TAlign & align, TString & mutations)
{

    //typedef typename Source<TAlign>::Type TSource;
    //typedef typename Iterator<TSource, Rooted>::Type TStringIterator;

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


/*
struct IndelVariant{
    bool called;   // did this variant pass calling criteria?
    int indelSize; // called diploid genotype (allele1 << 2 | allele2)
    int count;     // number of supporting reads
    int quality;   // a quality value associated with the call
    int coverage;  // totalCoverage at position
    DnaString sequence;
};
*/






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


///////////////////////////////////////////////////////////////////////////////////////
// Output SNPs, do realignment if a certain number of indels is observed in the reads
template <
    typename TFragmentStore,
    typename TReadCigars,
    typename TReadCounts,
    typename TGenomeName,
    typename TFile,
    typename TOptions
>
void dumpVariantsRealignBatchWrap(
    TFragmentStore              &fragmentStore,         // forward/reverse matches
    TReadCigars             &readCigars,
    TReadCounts const           &readCounts,
    TGenomeName const           genomeID,           // genome name
    typename TFragmentStore::TContigPos startCoord,         // startCoordinate + posOnGenomeInfix = real coordinate on whole chromosome
    typename TFragmentStore::TContigPos currWindowBegin,
    typename TFragmentStore::TContigPos currWindowEnd,
    TFile                   &fileSNPs,
    TFile                   &fileIndels,
    TOptions                &options)
{

    typedef typename TFragmentStore::TAlignedReadStore  TMatches;
    typedef typename Value<TMatches>::Type          TMatch;
    typedef typename TFragmentStore::TAlignQualityStore     TMatchQualities;
    //typedef typename Value<TMatchQualities>::Type       TMatchQuality;
    //typedef typename TFragmentStore::TReadSeqStore      TReads;
    //typedef typename Value<TReads>::Type            TRead;
    //typedef typename TFragmentStore::TContigStore       TContigStore;
    typedef typename TFragmentStore::TContigPos         TContigPos;
    //typedef typename Value<TContigStore>::Type      TContig;
    //typedef typename TFragmentStore::TContigSeq         TContigSeq;
    typedef typename Iterator<TMatches,Standard>::Type  TMatchIterator;


    TMatches &matches       = fragmentStore.alignedReadStore;
    TMatchQualities &matchQualities = fragmentStore.alignQualityStore;

    ::std::sort(begin(matches, Standard()), end(matches, Standard()), LessGPos<TMatch>());

    TMatchIterator matchIt = begin(matches,Standard());
    TMatchIterator matchItEnd = end(matches,Standard());

//    unsigned minNumIndels = 0;
    unsigned minNumIndels = options.indelCountThreshold;
    //unsigned minNumIndels = 1;  // do realignment whenever there is at least one indel

//  std::fstream tmpfile;
//  tmpfile.open("Z:/seqan071010/projects/library/apps/chr4.beforeTotal.sam", ::std::ios_base::out);
//  write(tmpfile, fragmentStore, Sam());
//  tmpfile.close();

#ifdef SNPSTORE_DEBUG
    CharString strstr = "test";
    _dumpMatches(fragmentStore,strstr);
#endif

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
        while(matchIt != matchItEnd && _min((*matchIt).beginPos,(*matchIt).endPos) < groupEndPos)
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

            //FragmentStore<SnpStoreGroupSpec_> fragStoreGroup;
            //copyFragmentStore(fragStoreGroup,fragmentStore,matchItBatchBegin,matchItBatchEnd,groupStartPos,groupEndPos);

            TFragmentStore fragStoreGroup = fragmentStore;
            arrayMoveForward(matchItBatchBegin,matchItBatchEnd,begin(fragStoreGroup.alignedReadStore,Standard())); // reads wont be needed anymore
            resize(fragStoreGroup.alignedReadStore,numMatches,Exact());
            fragStoreGroup.contigStore[0].seq = infix(fragmentStore.contigStore[0].seq,groupStartCoordLocal,groupEndCoordLocal);

#ifdef SNPSTORE_DEBUG
            std::cout << "in realign wrap: groupEndPos = " <<  groupEndPos << " groupStartPos=" <<  groupStartPos << std::endl;
            std::cout << "genomeLength= " <<  length(fragmentStore.contigStore[0].seq) << std::endl;

            CharString strstre = "testgroup";
            _dumpMatches(fragStoreGroup,strstre);

#endif
            groupStartPos += startCoord;
            groupEndPos += startCoord;
            TContigPos groupStartCoord = startCoord + groupStartCoordLocal;
            groupStartPos = _max(groupStartPos,currWindowBegin);
            groupEndPos = _min(groupEndPos,currWindowEnd);

            //the current group is formed by all reads from matchItBatchBegin until matchItBatchEnd
            if(indelReadCount >= (int)minNumIndels && options.realign)
            {
                //do realignment
                dumpVariantsRealignBatch(fragStoreGroup,readCigars,
                    readCounts,genomeID,
                    groupStartCoord,groupStartPos,groupEndPos,
                    fileSNPs,fileIndels,options);
            }
            else
            {
                // todo: switch between with or without realignment in dumpVariantsRealignBatch.. problem: split reads
/*              dumpVariantsRealignBatch(fragStoreGroup,readCigars,
                    readCounts,genomeID,
                    groupStartCoord,groupStartPos,groupEndPos,
                    fileSNPs,fileIndels,options);*/
                    dumpSNPsBatch(fragStoreGroup,readCigars,
                        readCounts,genomeID,
                        groupStartCoord,groupStartPos,groupEndPos,
                        fileSNPs,options);
            }
        }

    }


}

///////////////////////////////////////////////////////////////////////
// SNP calling Maq style
template<typename TCounts, typename TQualities, typename TOptions>
inline bool
_doSnpCall(TCounts & countF,
          TCounts & countR,
          TQualities & qualF,
          TQualities & qualR,
          int &refAllele,
          TOptions & options,
          SingleBaseVariant &snp,
          MaqMethod&
#ifdef SNPSTORE_DEBUG_CANDPOS
            , int candPos
#endif
          )
{


        // the diploid reference genotype
        int genotypeRef = (refAllele<<2) | refAllele;
        int genotypeCalled = genotypeRef, qCall1 = 0;  // genotype call quality
        int qSnp = 0; // SNP call quality
        // int genotypeCalled2 = genotypeRef, qCall2 = 0;


#ifdef SNPSTORE_DEBUG_CANDPOS
        bool extraV = false;
        if(candPos == 118487871) extraV = true;
        if(extraV)
        {
            ::std::cout << "Forward qualities:\n" << std::flush;
            for(unsigned x = 0; x < length(qualF); ++x)
                ::std::cout << qualF[x] << "\t";
            ::std::cout << "\nReverse qualities:\n" << std::flush;
            for(unsigned x = 0; x < length(qualR); ++x)
                ::std::cout << qualR[x] << "\t";
            ::std::cout << "\n" << std::flush;
        }
#endif

        // do the Maq statistics
        //
        // argmax P(g|D)=P(D|g)*P(g)/P(D)
        //    g
        //

        // get pHomo for best and second best nucleotide
        int best, secondBest;
        long double pHet = 0, pHomo1 = 0, pHomo2 = 0;
        getHomoProbs(options.cnks,options.fks,qualF,qualR,best,secondBest,pHomo1,pHomo2,
#ifdef SNPSTORE_DEBUG_CANDPOS
            candPos
#else
            0
#endif
            );
        if(secondBest == -1)
        {
            if(best==refAllele) // shouldnt happen
                return false;
            secondBest = refAllele;
        }

        //get pHet
        int n = countF[best] + countR[best] + countF[secondBest] + countR[secondBest];
#ifdef SNPSTORE_DEBUG_CANDPOS
        if(extraV)
        {
            std::cout << " n = " <<n << std::endl;
            std::cout << "(countF[secondBest] + countR[secondBest]) = " << (countF[secondBest] + countR[secondBest]) << std::endl;
        }
#endif

        // should always access hetTable with n,refAllele (unless neither best nor second best is ref)
        //if(secondBest==refAllele && ) // switch best and secondBest
        //{
        //    int temp = secondBest;
        //    secondBest = best;
        //    best = temp;
        //}

        if (n > 255)
        {
            int temp2 = (int)((countF[secondBest] + countR[secondBest])*255.0/n + 0.5);
            int temp1 = (int)((countF[best] + countR[best])*255.0/n + 0.5);
#ifdef SNPSTORE_DEBUG_CANDPOS
        if(extraV)
            std::cout << "temp1 = " << temp1 << std::endl;
#endif
            pHet = options.priorHetQ - 4.343 * options.hetTable[temp2<<8|temp1];
//          pHet = options.priorHetQ - 4.343 * options.hetTable[255<<8|temp];
        }
        else
            pHet = options.priorHetQ - 4.343 * options.hetTable[(countF[secondBest] + countR[secondBest])<<8|(countF[best] + countR[best])];
//          pHet = options.priorHetQ - 4.343 * options.hetTable[n<<8|(countF[secondBest] + countR[secondBest])];

#ifdef SNPSTORE_DEBUG_CANDPOS
        if(extraV)
        {
            std::cout << "refAllele = " << refAllele << std::endl;
        	std::cout << "best = " << best << " with " << countF[best]+countR[best] << std::endl;
        	std::cout << "secondbest = " << secondBest << " with " << countF[secondBest]+countR[secondBest] << std::endl;
        	std::cout << "pHet = " << pHet << std::endl;
        	std::cout << "pHomo1 = " << pHomo1 << std::endl;
        	std::cout << "pHomo2 = " << pHomo2 << std::endl;
        }
#endif

        pHet = (pHet > 0.0) ? pHet : 0.0;
        pHomo1 = (pHomo1 > 0.0) ? pHomo1 : 0.0;
        pHomo2 = (pHomo2 > 0.0) ? pHomo2 : 0.0;

        double pRef = pHomo1;
        if(best != refAllele)
            pRef = pHomo2;
        if(best != refAllele && secondBest != refAllele)
            qSnp = 255;

        int het,/*homo1,*/homo2; //0,1,2

        //rank them and create the genotype
        if(pHet < pHomo1)
        {
            if(pHet < pHomo2)
            {
                het = 0; //het is best
                if(best==refAllele)
                    genotypeCalled = (best<<2) | secondBest;
                else
                    genotypeCalled = (secondBest<<2) | best;

                if(pHomo1<=pHomo2)    //(1)
                {
                    qCall1 = (int)(pHomo1 - pHet  + 0.5);
                    qSnp = (int)(pRef - pHet  + 0.5);
                    // homo1 = 1; //second best
                    // genotypeCalled2 = (best<<2)| best;
                    // qCall2 = (int)(pHomo2 - pHomo1 + 0.5);
                    homo2 = 2; // last
                }
                else                //(2)
                {
                    qCall1 = (int)(pHomo2 - pHet + 0.5);
                    qSnp = (int)(pRef - pHet  + 0.5);
                    homo2 = 1;
                    // genotypeCalled2 = (secondBest<<2)| secondBest;
                    // qCall2 = (int)(pHomo1 - pHomo2 + 0.5);
                    // homo1 = 2;
                }
            }
            else
            {                       //(3)
                // shouldnt happen
                homo2 = 0;
                qCall1 = (int)(pHet - pHomo2 + 0.5);
                qSnp = (int)(pRef - pHomo2  + 0.5);
                genotypeCalled = (secondBest<<2)| secondBest;
                het = 1;
                // qCall2 = (int)(pHomo1 - pHet + 0.5);
                // if(best==refAllele)
                //  genotypeCalled2 = (best<<2)| secondBest;
                // else
                //  genotypeCalled2 = (secondBest<<2) | best;
                // homo1 = 2;
            }
        }
        else
        {
            if(pHomo2 < pHomo1) //(4)
            {
                //this case shouldnt happen
                // std::cout << "Weird case 4" << std::endl;
                homo2 = 0;
                qCall1 = (int)(pHomo1 - pHomo2 + 0.5);
                qSnp = (int)(pRef - pHomo2  + 0.5);
                genotypeCalled = (secondBest<<2)| secondBest;
                // homo1 = 1;
                // qCall2 = (int)(pHet - pHomo1 + 0.5);
                // genotypeCalled2 = (best<<2)| best;
                het = 2;
            }
            else
            {
                // homo1 = 0;
                genotypeCalled = (best<<2)| best;
                if(pHet<pHomo2)     //(5)
                {
                    qCall1 = (int)(pHet - pHomo1 + 0.5);
                    qSnp = (int)(pRef - pHomo1  + 0.5);
                    het = 1;
                    // if(best==refAllele)
                    //  genotypeCalled2 = (best<<2)| secondBest;
                    // else
                    //  genotypeCalled2 = (secondBest<<2) | best;
                    // qCall2 = (int)(pHomo2 - pHet + 0.5);
                    homo2 = 2;
                }
                else        // (6)
                {
                    qCall1 = (int)(pHomo2 - pHomo1+ 0.5);
                    qSnp = (int)(pRef - pHomo1  + 0.5);
                    homo2 = 1;
                    // qCall2 = (int)(pHet - pHomo2 + 0.5);
                    // genotypeCalled2 = (secondBest<<2)| secondBest;
                    het = 2;
                }
            }
        }

        if (het != 0 && homo2 == 0) //
        {
#ifdef SNPSTORE_DEBUG_CANDPOS
            std::cout << "Second best is best homozygote?!" << std::endl;
#endif
            //return false;
            // genotypeCalled2 = genotypeCalled;
            genotypeCalled = genotypeRef;   // disable call
            homo2 = 1;
            het = 2;
            qCall1 = 0;
            // qCall2 = 0;
        }

    //}
//    if(best != refAllele && secondBest != refAllele)
//        qSnp = 255;


    unsigned totalCoverage = countF[0] + countF[1] +countF[2] +countF[3] +countF[4]
                           + countR[0] + countR[1] +countR[2] +countR[3] +countR[4];

    snp.genotype   = genotypeCalled;
    snp.count      = countF[best] + countR[best];
    snp.quality    = qCall1;
    snp.snpQuality = qSnp;
    snp.coverage   = totalCoverage;
    int consideredCount = countF[secondBest] + countR[secondBest] + countF[best] + countR[best];
    if (genotypeCalled == genotypeRef) snp.called = false;
    else snp.called = true;
    if((double)consideredCount/totalCoverage < options.minExplainedColumn)
    {
        snp.called = false;
    }

    return true;
}


template<typename TCounts, typename TQualities, typename TOptions>
inline bool
_doSnpCall(TCounts & countF,
          TCounts & countR,
          TQualities & qualF,       // columnQualityF
          TQualities & qualR,       // columnQualityR
          int &refAllele,
          TOptions & options,
          SingleBaseVariant &snp,
          ThresholdMethod&)
{

    // find potential mutation allele
    int allele1 = -1;   // most frequent allele
    int allele2 = -1;   // second most frequent allele

    unsigned maxCount=0;
    for(int k=0; k < 5; ++k)
    {
        if(countF[k]+countR[k] > maxCount)
        {
            maxCount = countF[k]+countR[k];
            allele1 = k;
        }
    }
    maxCount = 0;
    for(int k=0; k < 5; ++k)
    {
        if(k != allele1 && countF[k]+countR[k] >= maxCount)
        {
            maxCount = countF[k]+countR[k];
            allele2 = k;
        }
    }

    // No evidence of non-ref bases left... (used to happen with onthefly-pileupcorrection
    // cannot happen anymore as these positions would never be inspected)
    if(allele1==refAllele && allele2==refAllele)
    {
        ::std::cout << "No non-ref base observed. Correct??\n";
        return false;
    }

    // get the mutational allele
    int mutAllele=allele1;
    if(allele1==refAllele) mutAllele=allele2;

    unsigned mutCoverage   = countF[mutAllele] + countR[mutAllele];
    unsigned totalCoverage = countF[0] + countF[1] +countF[2] +countF[3] +countF[4]
                           + countR[0] + countR[1] +countR[2] +countR[3] +countR[4];

    // the diploid reference genotype
    int genotypeRef = (refAllele<<2) | refAllele;
    int genotypeCalled = genotypeRef;

    // threshold method
    if( mutCoverage >= options.minMutT
        && (float)mutCoverage/totalCoverage >= options.percentageT
        && (float)(qualF[mutAllele]+qualR[mutAllele])/mutCoverage >= options.avgQualT)
    {
        if((float)mutCoverage/totalCoverage <= options.snpHetMax)
            genotypeCalled = (mutAllele<<2)|refAllele; // we dont attempt real genotype calling here
        else
            genotypeCalled = (mutAllele<<2)|mutAllele; // we dont attempt real genotype calling here
    }


    snp.genotype = genotypeCalled;
    snp.count    = mutCoverage;
    snp.quality  = (qualF[mutAllele]+qualR[mutAllele])/ mutCoverage;
    snp.coverage = totalCoverage;
    if (genotypeCalled == genotypeRef) snp.called = false;
    else snp.called = true;

    return true;
}


// write to file
template<typename TFile, typename TString, typename TQualities, typename TPos, typename TOptions>
inline bool
_writeSnp(TFile &file,
       SingleBaseVariant &snp,
       TQualities &qualityStringF,
       TQualities &qualityStringR,
       int refAllele,
       TString &genomeID,
       TPos candPos,
       unsigned realCoverage,
       TOptions &options)
{
//IOREV _nodoc_ what kind of format is this?
    if (!file.is_open())
    {
        ::std::cerr << "SNP output file is not open" << ::std::endl;
        return false;
    }

    //chromosome
    file << genomeID << '\t';
    file << candPos + options.positionFormat<< '\t';
    file << (Dna5)refAllele <<'\t';
    if(options.orientationAware)
    {
        if(options.showQualityStrings)
        {
            file << "["<<qualityStringF[0] <<"]\t";
            file << "["<<qualityStringF[1] <<"]\t";
            file << "["<<qualityStringF[2] <<"]\t";
            file << "["<<qualityStringF[3] <<"]\t";
            file << "["<<qualityStringR[0] <<"]\t";
            file << "["<<qualityStringR[1] <<"]\t";
            file << "["<<qualityStringR[2] <<"]\t";
            file << "["<<qualityStringR[3] <<"]\t";
        }
        else
        {
            file << length(qualityStringF[0]) <<"\t";
            file << length(qualityStringF[1]) <<"\t";
            file << length(qualityStringF[2]) <<"\t";
            file << length(qualityStringF[3]) <<"\t";
            file << length(qualityStringR[0]) <<"\t";
            file << length(qualityStringR[1]) <<"\t";
            file << length(qualityStringR[2]) <<"\t";
            file << length(qualityStringR[3]) <<"\t";
        }
    }
    else
    {
        if(options.showQualityStrings)
        {
            file << "["<<qualityStringF[0]<<qualityStringR[0] <<"]\t";
            file << "["<<qualityStringF[1]<<qualityStringR[1] <<"]\t";
            file << "["<<qualityStringF[2]<<qualityStringR[2] <<"]\t";
            file << "["<<qualityStringF[3]<<qualityStringR[3] <<"]\t";
        }
        else
        {
            file << length(qualityStringF[0])+length(qualityStringR[0]) <<"\t";
            file << length(qualityStringF[1])+length(qualityStringR[1]) <<"\t";
            file << length(qualityStringF[2])+length(qualityStringR[2]) <<"\t";
            file << length(qualityStringF[3])+length(qualityStringR[3]) <<"\t";
        }
    }
    file << realCoverage;

    if(options.method == 1)
    {
        //genotypeCalled to string
        if(snp.called)//genotypeCalled != genotypeRef)
            file << '\t' << (char)options.toIupac[(unsigned)(snp.genotype&15)]<< '\t' << snp.quality << '\t' << snp.snpQuality;
        else
            file << "\t\t";
    }
    else
    {
        if(snp.called)
            file  << '\t' << (char)options.toIupac[(unsigned)(snp.genotype)&15] << '\t' << snp.quality;// mutAllele;// << "/" << (Dna)mutAllele;
//            file  << '\t' << (Dna)(snp.genotype & 3) << '\t' << snp.quality;// mutAllele;// << "/" << (Dna)mutAllele;
        else file << "\t\t";
    }
    file << std::endl;
    return true;
}


// write to file
template<typename TFile, typename TString, typename TQualities, typename TPos, typename TOptions>
inline bool
_writePos(TFile &file,
       TQualities &qualityStringF,
       TQualities &qualityStringR,
       unsigned delPlus,
       unsigned delMinus,
       TString &genomeID,
       TPos candPos,
       unsigned /*coverage*/,
       TOptions &options)
{
//IOREV _nodoc_ what kind of format is this?
    if (!file.is_open())
    {
        ::std::cerr << "SNP output file is not open" << ::std::endl;
        return false;
    }

    //chromosome
    file << genomeID << '\t';
    file << candPos + options.positionFormat<< '\t';
    if(options.orientationAware)
    {
        file << length(qualityStringF[0]) <<"\t";
        file << length(qualityStringF[1]) <<"\t";
        file << length(qualityStringF[2]) <<"\t";
        file << length(qualityStringF[3]) <<"\t";
        file << delPlus <<"\t";
        file << length(qualityStringR[0]) <<"\t";
        file << length(qualityStringR[1]) <<"\t";
        file << length(qualityStringR[2]) <<"\t";
        file << length(qualityStringR[3]) <<"\t";
        file << delMinus ; //<< std::endl;
    }
    else
    {
        file << length(qualityStringF[0])+length(qualityStringR[0]) <<"\t";
        file << length(qualityStringF[1])+length(qualityStringR[1]) <<"\t";
        file << length(qualityStringF[2])+length(qualityStringR[2]) <<"\t";
        file << length(qualityStringF[3])+length(qualityStringR[3]) <<"\t";
        file << delPlus+delMinus ; // << std::endl;
    }
    unsigned coverage = length(qualityStringF[0])+length(qualityStringF[1])+length(qualityStringF[2])+length(qualityStringF[3])+delPlus;
             coverage+= length(qualityStringR[0])+length(qualityStringR[1])+length(qualityStringR[2])+length(qualityStringR[3])+delMinus;
    file << '\t' << coverage << std::endl;

    //if(options.method == 1)
    //{
    //  //genotypeCalled to string
    //  if(snp.called)//genotypeCalled != genotypeRef)
    //      file << '\t' << (char)options.toIupac[(unsigned)(snp.genotype&15)]<< '\t' << snp.quality;
    //  else
    //      file << "\t\t";
    //}
    //else
    //{
    //  if(snp.called)
    //      file  << '\t' << (Dna)(snp.genotype & 3) << '\t' << snp.quality;// mutAllele;// << "/" << (Dna)mutAllele;
    //  else file << "\t\t";
    //}
    //file << std::endl;
    return true;
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


template<typename TFragmentStore, typename TId, typename TOptions>
void
realignReferenceToReadProfile(TFragmentStore & fragmentStore,
                              TId refId,
                              TOptions & options)
{
// fragment store types
    typedef typename TFragmentStore::TAlignedReadStore  TMatches;
    //typedef typename Value<TMatches>::Type              TMatch;
    typedef typename TFragmentStore::TAlignQualityStore TMatchQualities;
    //typedef typename Value<TMatchQualities>::Type       TMatchQuality;
    typedef typename Iterator<TMatches,Standard>::Type  TMatchIterator;

    typedef typename TFragmentStore::TReadSeqStore      TReads;
    typedef typename Value<TReads>::Type                TRead;
    typedef typename Iterator<TRead, Standard>::Type    TReadIter;
    typedef typename Value<TRead>::Type                 TAlphabet;
    typedef typename Size<TRead>::Type                  TSize;
    typedef typename TFragmentStore::TReadPos           TReadPos;

    typedef typename TFragmentStore::TContigStore       TContigStore;
    typedef typename Value<TContigStore>::Type          TContig;
    //typedef typename TFragmentStore::TContigPos         TContigPos;
    typedef typename TFragmentStore::TContigSeq         TContigSeq;

    typedef Gaps<TContigSeq, AnchorGaps<typename TContig::TGapAnchors> >    TContigGaps;
    //typedef Gaps<TRead, AnchorGaps<typename TMatch::TGapAnchors> >          TReadGaps;
    //typedef typename Iterator<TContigGaps>::Type                            TContigGapIter;
    typedef typename Iterator<String<typename TFragmentStore::TReadGapAnchor>, Standard>::Type TReadGapsIter;

// profile types
    typedef ProfileChar<TAlphabet> TProfile;
    typedef String<TProfile> TProfileString;
    typedef typename Iterator<TProfileString, Standard>::Type TProfIter;

    extern int SEQAN_CONSENSUS_OPEN_PENALTY_FACTOR;
    int tmp = SEQAN_CONSENSUS_OPEN_PENALTY_FACTOR;
    SEQAN_CONSENSUS_OPEN_PENALTY_FACTOR = 1;

    TReads &reads                   = fragmentStore.readSeqStore;
    TMatches &matches               = fragmentStore.alignedReadStore;
    TMatchQualities &matchQualities = fragmentStore.alignQualityStore;


//    separateAlleles(fragmentStore,options);

    // make profile of multi-read-alignment
    TSize gapChar = ValueSize<TAlphabet>::VALUE;
    TReadPos minPos = 0; // TODO: check if correct
    TContigGaps contigGaps(fragmentStore.contigStore[0].seq, fragmentStore.contigStore[0].gaps);
    TReadPos maxPos = (TReadPos)positionSeqToGap(contigGaps,length(fragmentStore.contigStore[0].seq)-1)+1;
    TProfileString multiReadProfile;
    resize(multiReadProfile, maxPos - minPos, TProfile()); // get maxPos minPos

    TProfIter it = begin(multiReadProfile, Standard() );
    TProfIter itEnd = end(multiReadProfile, Standard());
    TMatchIterator matchIt = begin(matches, Standard() );
    TMatchIterator matchItEnd = end(matches, Standard() );
    for(;matchIt != matchItEnd; ++matchIt) {
        matchIt->beginPos -= minPos;
        matchIt->endPos -= minPos;
        it = begin(multiReadProfile, Standard());
        it += _min(matchIt->beginPos,matchIt->endPos);

        if(matchIt->beginPos > matchIt->endPos)
            reverseComplement(fragmentStore.readSeqStore[matchIt->readId]);
        TRead& readSeq = fragmentStore.readSeqStore[matchIt->readId];

        TReadIter itRead = begin(readSeq, Standard() );
        TReadIter itReadEnd = end(readSeq, Standard() );

        TReadGapsIter gitRead = begin(matchIt->gaps, Standard() );
        TReadGapsIter gitReadEnd = end(matchIt->gaps, Standard() );

        TReadPos old = 0;
        int diff = 0;
        bool clippedEnd = false;
        if ((gitRead != gitReadEnd) && (gitRead->gapPos == 0)) {
      //      std::cout << "does this happen? shouldnt for semiglobal matches\n";
            old = gitRead->seqPos;
            itRead += old;
            diff -= old;
            ++gitRead;
        }
        for(;gitRead != gitReadEnd; ++gitRead) {
            TReadPos limit = gitRead->seqPos;
            int newDiff = (gitRead->gapPos - limit);
            SEQAN_ASSERT_LT(gitRead->gapPos, (int)length(multiReadProfile));
            if (diff > newDiff) {
                limit -= (diff - newDiff);
                clippedEnd = true;
            }
            // add non-gap positions to consensus
            for(;old < limit && itRead != itReadEnd && it != itEnd; ++old, ++itRead)
            {
                SEQAN_ASSERT_LT(itRead, itReadEnd);
                ++(value(it++)).count[ordValue(*itRead)];
            }
            // add gap positions to consensus
            for(;diff < newDiff; ++diff)
                ++(value(it++)).count[gapChar];
        }
        if (!clippedEnd) {
            for( ; itRead!=itReadEnd && it != itEnd;++itRead)
                ++(value(it++)).count[ordValue(*itRead)];
        }
        if(matchIt->beginPos > matchIt->endPos)
            reverseComplement(fragmentStore.readSeqStore[matchIt->readId]);

    }

    // put reference sequence into ProfileString
    TReadIter refSeqIt = begin(fragmentStore.readSeqStore[refId],Standard());
    TProfileString refProfile;
    resize(refProfile, length(fragmentStore.readSeqStore[refId]), TProfile());
    TProfIter refIt = begin(refProfile,Standard());
    for(; refIt != end(refProfile,Standard()); ++refIt, ++refSeqIt)
        (*refIt).count[0] = ordValue(*refSeqIt);


    // Realign the reference to the reads
    typedef StringSet<TProfileString, Dependent<> > TStringSet;
    TStringSet pairSet;

    appendValue(pairSet, multiReadProfile);
    appendValue(pairSet, refProfile);

    //for(TSize i = 0; i<length( pairSet[0]); ++i) {
    //  std::cout <<  pairSet[0][i] << std::endl;
    //}
    //std::cout << "_______________" << std::endl;
    //for(TSize i = 0; i<length( pairSet[1]); ++i) {
    //  std::cout <<   pairSet[1][i] << std::endl;
    //}
    //std::cout << "..............." << std::endl;

    //Score<int, WeightedConsensusScore<Score<int, FractionalScore>, Score<int, ConsensusScore> > > consScore;
    Score<int, ConsensusScore>  consScore;

    // TODO: Score<int, WeightedConsensusScore<Score<int, AffineFractionalScore>, Score<int, AffineConsensusScore> > > consScore;
    // or: quality fractional score! x/n where n=sum(q(read)) and x = sum(q(reads with base x))
    // consScore.verticalGapOpen = 10; bzw dependent on read depth? or read quality? --> higher quality --> higher penalty?
    // consScore.horizontalGapOpen = 10;
    assignProfile(consScore, multiReadProfile);

    typedef String<Fragment<> > TFragmentString;
    TFragmentString fragments;

    TReadPos bandWidth = 30;     // set somewhere else
    TReadPos leftDiag = -bandWidth;
    TReadPos rightDiag = bandWidth;
    //make sure the band contains the first and last position of both refProfile and multiReadProfile
    if(length(refProfile) < length(multiReadProfile)) rightDiag +=  length(multiReadProfile) - length(refProfile);
    if(length(refProfile) > length(multiReadProfile)) leftDiag -=  length(refProfile) - length(multiReadProfile);


    //// Debug code
//  Graph<Alignment<TStringSet, void, WithoutEdgeId> > g1(pairSet);
//  int sc1 = globalAlignment(g1, consScore, AlignConfig<false,true,true,false>(), leftDiag, rightDiag, Gotoh());
//  std::cout << sc1 << std::endl;
//  std::cout << g1 << std::endl;


 //   Score<int> scoreType = Score<int>(0, -1, -2, -10);    // (match, mismatch,gapExtend,gapOpen)
    //StringSet<TRead, Dependent<> > pairSet2;
 //   appendValue(pairSet2, fragmentStore.contigStore[0].seq);
 //   appendValue(pairSet2, fragmentStore.readSeqStore[refId]);
    //Graph<Alignment<StringSet<TRead, Dependent<> >, void, WithoutEdgeId> > g2(pairSet2);
    //int sc2 = globalAlignment(g2, scoreType, AlignConfig<false,true,true,false>(), leftDiag, rightDiag, Gotoh());
    //std::cout << sc2 << std::endl;
    //std::cout << g2 << std::endl;


    // make diploid consensus --> realign ref with that
    it = begin(multiReadProfile, Standard() );
    itEnd = end(multiReadProfile, Standard());

    String<bool> removeState;
    resize(removeState,length(multiReadProfile));
    typedef typename Iterator<String<bool>, Standard >::Type TStateIterator;
    TStateIterator sit = begin(removeState,Standard());
    TProfileString diploidConsensus;
    resize(diploidConsensus,length(multiReadProfile),TProfile());
    TProfIter dipIt = begin(diploidConsensus, Standard() );
   // std::cout << "length before removing low freq gaps: " << length(diploidConsensus) << std::endl;
    for(; it != itEnd; ++it, ++sit)
    {
        TSize getMax1 = 0;
        TSize countMax1 = 0;
        for(TSize i = 0; i < ValueSize<TProfile>::VALUE; ++i)
            if((*it).count[i] >= (*it).count[getMax1])
            {
                getMax1 = i;
                countMax1 = (*it).count[getMax1];
            }
        int getMax2 = 0;
        int countMax2 = -1;
        for(TSize i = 0; i < ValueSize<TProfile>::VALUE; ++i)
            if(i != getMax1 && (*it).count[i] >= (*it).count[getMax2])
            {
                getMax2 = i;
                countMax2 = (*it).count[getMax2];
            }
        bool remove = true;
        if(getMax1 != gapChar)
        {
            (*dipIt).count[getMax1] = (*it).count[getMax1];
            remove = false;
        }
        if(getMax1 != gapChar && countMax2 > 0 && getMax2 != gapChar)
        {
            (*dipIt).count[getMax2] = (*it).count[getMax2];
            remove = false;
        }
        if(getMax1 == gapChar && countMax2 > 0 && (float)countMax2/countMax1 > options.indelCountThreshold)
        {
            (*dipIt).count[getMax1] = (*it).count[getMax1];
            (*dipIt).count[getMax2] = (*it).count[getMax2];
            remove = false;
        }
        if(getMax2 == gapChar && countMax2 > 0 && (float)countMax2/countMax1 > options.indelCountThreshold)
        {
            (*dipIt).count[getMax2] = (*it).count[getMax2];
            remove = false;
        }
        *sit = remove;
        if(!remove)++dipIt;
    }
    resize(diploidConsensus,dipIt - begin(diploidConsensus),Exact());
    //std::cout << "length after removing low freq gaps: " << length(diploidConsensus) << std::endl;

    typedef StringSet<TProfileString, Dependent<> > TStringSet;
    TStringSet pairSet3;
    appendValue(pairSet3, diploidConsensus);
    appendValue(pairSet3, refProfile);

    //for(TSize i = 0; i<length( pairSet[0]); ++i) {
    //  std::cout <<  pairSet[0][i] << std::endl;
    //}
    //std::cout << "_______________" << std::endl;
    //for(TSize i = 0; i<length( pairSet[1]); ++i) {
    //  std::cout <<   pairSet[1][i] << std::endl;
    //}
    //std::cout << "..............." << std::endl;

    Score<int, WeightedConsensusScore<Score<int, FractionalScore>, Score<int, ConsensusScore> > > consScore3;
    assignProfile(consScore3, diploidConsensus);

    bandWidth = 30;     // set somewhere else
    leftDiag = -bandWidth;
    rightDiag = bandWidth;
    //make sure the band contains the first and last position of both refProfile and multiReadProfile
    if(length(refProfile) < length(diploidConsensus)) rightDiag +=  length(diploidConsensus) - length(refProfile);
    if(length(refProfile) > length(diploidConsensus)) leftDiag -=  length(refProfile) - length(diploidConsensus);

    TFragmentString fragments2;
    //// Debug code
    Graph<Alignment<TStringSet, void, WithoutEdgeId> > g3(pairSet3);
    //int sc3 =
    globalAlignment(g3, consScore3, AlignConfig<false,true,true,false>(), leftDiag, rightDiag, Gotoh());
    //std::cout << sc3 << std::endl;
    //std::cout << g3 << std::endl;


    SEQAN_CONSENSUS_OPEN_PENALTY_FACTOR = tmp;

}



template<typename TFragmentStore, typename TId, typename TOptions>
void
realignReferenceToDiploidConsensusProfile(TFragmentStore & fragmentStore,
                              TId refReadId,
                              TOptions & options)
{
// fragment store types
    typedef typename TFragmentStore::TAlignedReadStore  TMatches;
    typedef typename Value<TMatches>::Type              TMatch;
    typedef typename TFragmentStore::TAlignQualityStore TMatchQualities;
    //typedef typename Value<TMatchQualities>::Type       TMatchQuality;
    typedef typename Iterator<TMatches,Standard>::Type  TMatchIterator;

    typedef typename TFragmentStore::TReadSeqStore      TReads;
    typedef typename Value<TReads>::Type                TRead;
    typedef typename Iterator<TRead, Standard>::Type    TReadIter;
    typedef typename Value<TRead>::Type                 TAlphabet;
    typedef typename Size<TRead>::Type                  TSize;
    typedef typename TFragmentStore::TReadPos           TReadPos;
    typedef typename TFragmentStore::TReadGapAnchor     TGapAnchor;

    typedef typename TFragmentStore::TContigStore       TContigStore;
    typedef typename Value<TContigStore>::Type          TContig;
    //typedef typename TFragmentStore::TContigPos         TContigPos;
    typedef typename TFragmentStore::TContigSeq         TContigSeq;

    typedef Gaps<TContigSeq, AnchorGaps<typename TContig::TGapAnchors> >    TContigGaps;
    //typedef Gaps<TRead, AnchorGaps<typename TMatch::TGapAnchors> >          TReadGaps;
    //typedef typename Iterator<TContigGaps>::Type                            TContigGapIter;
    typedef typename Iterator<String<typename TFragmentStore::TReadGapAnchor>, Standard>::Type TReadGapsIter;

// profile types
    typedef ProfileChar<TAlphabet> TProfile;
    typedef String<TProfile> TProfileString;
    typedef typename Iterator<TProfileString, Standard>::Type TProfIter;

    extern int SEQAN_CONSENSUS_OPEN_PENALTY_FACTOR;
    int tmp = SEQAN_CONSENSUS_OPEN_PENALTY_FACTOR;

    TReads &reads                   = fragmentStore.readSeqStore;
    TMatches &matches               = fragmentStore.alignedReadStore;
    TMatchQualities &matchQualities = fragmentStore.alignQualityStore;


    // make profile of multi-read-alignment
    TSize gapChar = ValueSize<TAlphabet>::VALUE;
    TReadPos minPos = 0; // TODO: check if correct
    TContigGaps contigGaps(fragmentStore.contigStore[0].seq, fragmentStore.contigStore[0].gaps);
    TReadPos maxPos = (TReadPos)positionSeqToGap(contigGaps,length(fragmentStore.contigStore[0].seq)-1)+1;
    TProfileString multiReadProfile;
    resize(multiReadProfile, maxPos - minPos, TProfile()); // get maxPos minPos

    TId refMatchPosId = 0;
    bool refFound = false;
    TProfIter it = begin(multiReadProfile, Standard() );
    TProfIter itEnd = end(multiReadProfile, Standard());
    TMatchIterator matchIt = begin(matches, Standard() );
    TMatchIterator matchItEnd = end(matches, Standard() );
    int counter = 0;
    for(;matchIt != matchItEnd; ++matchIt) {
        if (matchIt->readId == refReadId)
        {
            refFound = true;
            refMatchPosId = (TId)counter;
            continue; // dont count the reference sequence in the multi-read-profile
        }
        ++counter;
        matchIt->beginPos -= minPos;
        matchIt->endPos -= minPos;
        it = begin(multiReadProfile, Standard());
        it += _min(matchIt->beginPos,matchIt->endPos);

        if(matchIt->beginPos > matchIt->endPos)
            reverseComplement(fragmentStore.readSeqStore[matchIt->readId]);
        TRead& readSeq = fragmentStore.readSeqStore[matchIt->readId];

        TReadIter itRead = begin(readSeq, Standard() );
        TReadIter itReadEnd = end(readSeq, Standard() );

        TReadGapsIter gitRead = begin(matchIt->gaps, Standard() );
        TReadGapsIter gitReadEnd = end(matchIt->gaps, Standard() );

        TReadPos old = 0;
        int diff = 0;
        bool clippedEnd = false;
        if ((gitRead != gitReadEnd) && (gitRead->gapPos == 0)) {
            std::cout << "does this happen? shouldnt for semiglobal matches\n";
            old = gitRead->seqPos;
            itRead += old;
            diff -= old;
            ++gitRead;
        }
        for(;gitRead != gitReadEnd; ++gitRead) {
            TReadPos limit = gitRead->seqPos;
            int newDiff = (gitRead->gapPos - limit);
            SEQAN_ASSERT_LT(gitRead->gapPos, (int)length(multiReadProfile));
            if (diff > newDiff) {
                limit -= (diff - newDiff);
                clippedEnd = true;
            }
            // add non-gap positions to consensus
            for(;old < limit && itRead != itReadEnd && it != itEnd; ++old, ++itRead)
            {
                SEQAN_ASSERT_LT(itRead, itReadEnd);
                ++(value(it++)).count[ordValue(*itRead)];
            }
            // add gap positions to consensus
            for(;diff < newDiff; ++diff)
                ++(value(it++)).count[gapChar];
        }
        if (!clippedEnd) {
            for( ; itRead!=itReadEnd && it != itEnd;++itRead)
                ++(value(it++)).count[ordValue(*itRead)];
        }
        if(matchIt->beginPos > matchIt->endPos)
            reverseComplement(fragmentStore.readSeqStore[matchIt->readId]);

    }
    SEQAN_ASSERT(refFound);

    // put reference sequence into ProfileString
    TReadIter refSeqIt = begin(fragmentStore.readSeqStore[refReadId],Standard());
    TProfileString refProfile;
    resize(refProfile, length(fragmentStore.readSeqStore[refReadId]), TProfile());
    TProfIter refIt = begin(refProfile,Standard());
    for(; refIt != end(refProfile,Standard()); ++refIt, ++refSeqIt)
        (*refIt).count[0] = ordValue(*refSeqIt);


    // align reference to full multi read profile

    // identify columns that are reference gaps and mainly gaps in multi read profile
    // --> check for adjacent fragments where fragmentEnd == nextFragmentBegin


    // make diploid consensus profile
    TProfileString diploidConsensus;
    resize(diploidConsensus,length(multiReadProfile),TProfile());
    TProfIter dipIt = begin(diploidConsensus, Standard() );

    // record which columns are removed from multiReadProfile to obtain diploidConsensus
    String<TReadPos> removeState;
    resize(removeState,length(multiReadProfile));
    typedef typename Iterator<String<TReadPos>, Standard >::Type TStateIterator;
    TStateIterator sit = begin(removeState,Standard());

    String<TReadPos> toFullProfile;
    resize(toFullProfile,length(multiReadProfile));
    typedef typename Iterator<String<TReadPos>, Standard >::Type TPosIterator;
    TPosIterator posIt = begin(toFullProfile,Standard());
    TReadPos pos = 0;
    // walk through and keep only the best (and sufficiently good) second best characters
    // and remove resulting gap-only columns
    it = begin(multiReadProfile, Standard() );
    itEnd = end(multiReadProfile, Standard());
    if(options._debugLevel > 1) std::cout << "length before removing low freq gaps: " << length(diploidConsensus) << std::endl;
    for(; it != itEnd; ++it)
    {
        TSize getMax1 = 0;
        TSize countMax1 = 0;
        for(TSize i = 0; i < ValueSize<TProfile>::VALUE; ++i)
            if((*it).count[i] >= (*it).count[getMax1])
            {
                getMax1 = i;
                countMax1 = (*it).count[getMax1];
            }
        int getMax2 = 0;
        int countMax2 = -1;
        for(TSize i = 0; i < ValueSize<TProfile>::VALUE; ++i)
            if(i != getMax1 && (*it).count[i] >= (*it).count[getMax2])
            {
                getMax2 = i;
                countMax2 = (*it).count[getMax2];
            }
        bool remove = true;
        if(getMax1 != gapChar)
        {
            (*dipIt).count[getMax1] = (*it).count[getMax1];
            remove = false;
        }
        if(getMax1 != gapChar && countMax2 > 0 && getMax2 != gapChar)
        {
            (*dipIt).count[getMax2] = (*it).count[getMax2];
            remove = false;
        }
//        if(getMax1 == gapChar && countMax2 > options.indelCountThreshold && (float)countMax2/countMax1 > options.indelPercentageT)
        if(getMax1 == gapChar && countMax2 > 0 && (float)countMax2/countMax1 > options.indelCountThreshold)
        {
            (*dipIt).count[getMax1] = (*it).count[getMax1];
            (*dipIt).count[getMax2] = (*it).count[getMax2];
            remove = false;
        }
//        if(getMax2 == gapChar && countMax2 > options.indelCountThreshold && (float)countMax2/countMax1 > options.indelPercentageT)
        if(getMax2 == gapChar && countMax2 > 0 && (float)countMax2/countMax1 > options.indelCountThreshold)
        {
            (*dipIt).count[getMax2] = (*it).count[getMax2];
            remove = false;
        }
//        *sit = remove;
        remove= false;
        if(!remove){
            *posIt = pos;
            ++posIt;
            ++dipIt;
        }
        else
        {
            *sit = pos;
            ++sit;
        }
        ++pos;
    }
    resize(diploidConsensus,dipIt - begin(diploidConsensus),Exact());
    resize(toFullProfile,posIt - begin(toFullProfile),Exact());
    *sit = length(multiReadProfile); ++sit;
    resize(removeState,sit - begin(removeState),Exact());
    if(options._debugLevel > 1) std::cout << "length after removing low freq gaps: " << length(diploidConsensus) << std::endl;

    typedef StringSet<TProfileString, Dependent<> > TStringSet;
    TStringSet pairSet;
    appendValue(pairSet, diploidConsensus);
    appendValue(pairSet, refProfile);

    //for(TSize i = 0; i<length( pairSet[0]); ++i) {
    //  std::cout <<  pairSet[0][i] << std::endl;
    //}
    //std::cout << "_______________" << std::endl;
    //for(TSize i = 0; i<length( pairSet[1]); ++i) {
    //  std::cout <<   pairSet[1][i] << std::endl;
    //}
    //std::cout << "..............." << std::endl;

    SEQAN_CONSENSUS_OPEN_PENALTY_FACTOR = 2;
    Score<int, WeightedConsensusScore<Score<int, FractionalScore>, Score<int, ConsensusScore> > > consScore;
    assignProfile(consScore, diploidConsensus);

    TReadPos bandWidth = 30;     // set somewhere else
    TReadPos leftDiag = -bandWidth;
    TReadPos rightDiag = bandWidth;
    //make sure the band contains the first and last position of both refProfile and multiReadProfile
    if(length(refProfile) < length(diploidConsensus)) rightDiag +=  length(diploidConsensus) - length(refProfile);
    if(length(refProfile) > length(diploidConsensus)) leftDiag -=  length(refProfile) - length(diploidConsensus);

    typedef String<Fragment<> > TFragmentString;
    TFragmentString fragments;
    //// Debug code
//  Graph<Alignment<TStringSet, void, WithoutEdgeId> > g3(pairSet3);
//  int sc3 = globalAlignment(g3, consScore, AlignConfig<false,true,true,false>(), leftDiag, rightDiag, Gotoh());
//  std::cout << sc3 << std::endl;
//  std::cout << g3 << std::endl;

//    std::cout << "leftDiag="<< leftDiag << std::endl;
//    std::cout << "rightDiag="<< rightDiag << std::endl;
    // reference can be aligned to gaps at the ends, diploidConsensus needs to be fully aligned
    globalAlignment(fragments, pairSet, consScore, AlignConfig<false,false,false,false>(), _max(leftDiag, -1 * (int) length(refProfile)), _min(rightDiag, (int) length(diploidConsensus)), Gotoh());

    // now use removeState string to retrieve reference<->multiReadProfile alignment
    // use oldPosLimits string to retrieve positions wrt multiReadProfile

    // transform coordinates back to full multi read profile


    // split fragments spanning removed columns
    TFragmentString fragments2;
    sit = begin(removeState,Standard());
    TStateIterator sitEnd = end(removeState,Standard());

    typedef typename Iterator<TFragmentString, Standard>::Type TFragIter;
    TFragIter fragIt = end(fragments, Standard());
    TFragIter fragItEnd = begin(fragments, Standard());
    if(fragIt != fragItEnd)
    {
        --fragIt;
        while (sit != sitEnd)
        {
            if(fragIt >= fragItEnd && toFullProfile[fragIt->begin1] > *sit) // removed column is before beginPos of fragment
            {
                ++sit;
            }
            else
            {
                if(toFullProfile[fragIt->begin1+fragIt->len-1]+1 <= *sit) // removed column is after endPos of fragment
                {
                    // copy fragment to new fragments with transformed coordinate
                    appendValue(fragments2,Fragment<>(0,toFullProfile[fragIt->begin1],1,fragIt->begin2,fragIt->len));
                    if(fragIt == fragItEnd)
                        break;
                    --fragIt;
                }
                else    // fragment contains removed column
                {
                    SEQAN_ASSERT_GT(*sit,toFullProfile[fragIt->begin1]);
                    int len = *sit-toFullProfile[fragIt->begin1];
                    // append new fragment: begin1, begin2 with len=removedCol-begin1 relative to whole profile coordinates
                    appendValue(fragments2,Fragment<>(0,toFullProfile[fragIt->begin1],1,fragIt->begin2,len));
                    // adapt old fragment accordingly (cut off beginning)
                    fragIt->begin2 += len;
                    fragIt->begin1 += len;
                    fragIt->len -= len;
                    SEQAN_ASSERT_GT(fragIt->len,0);
                    ++sit;
                }

            }

        }
    }
    fragIt = begin(fragments2, Standard());
    fragItEnd = end(fragments2, Standard());


    TMatch& refAli = fragmentStore.alignedReadStore[refMatchPosId];
    fragmentStore.alignedReadStore[refMatchPosId].beginPos = fragmentStore.alignedReadStore[refMatchPosId].endPos = 0;    // "disable read"
    clear(fragmentStore.alignedReadStore[refMatchPosId].gaps);    // -> reference gaps are handled separately

    // add reference back to multi-read-alignment according to new alignment
    TReadPos profilePos = 0;    // corresponds to consPos in realigner code
    TReadPos referencePos = 0;  // corresponds to readPos
    TReadPos alignPos = 0;      //

    // actually dont need these two in this case --> check and remove
    TReadPos bandOffset = 0;
    TReadPos clippedBeginPos = 0;

    TReadPos diff = 0;
    bool firstMatch = true;
    if (fragIt != fragItEnd) { // walk through segment matches (last one is leftmost one)
        do {
    //      --fragIt;
            int gapLen = fragIt->begin1 - profilePos;
            if (firstMatch) gapLen = 0; // gap between two adjacent segment matches
            // equivalent to profilePos + fraglen < nextProfilePos
            while (profilePos < (TReadPos)fragIt->begin1) { // cons stretch before newCons start
                ++profilePos;
                ++alignPos;
            }
            // equivalent to refPos + fraglen < nextRefPos
            while (referencePos < (TReadPos)fragIt->begin2) { // read stretch before matching fragment starts
                SEQAN_ASSERT_LT(referencePos, (TReadPos)length(fragmentStore.readSeqStore[refReadId]));
                // equivalent to profileDel
                if (gapLen) {
                    diff += gapLen; // add gap of length gaplen to readGaps
                    appendValue(fragmentStore.alignedReadStore[refMatchPosId].gaps, TGapAnchor(referencePos,referencePos + diff), Generous() );
                    gapLen = 0; // do this only once
                }
                //int numGaps =
                insertGap(matches, bandOffset + alignPos);
                ++referencePos;
                ++alignPos;
            }
            for (TSize i = 0; i<fragIt->len; ++i, ++profilePos, ++referencePos, ++alignPos) {
                SEQAN_ASSERT_LT(referencePos, (TReadPos)length(fragmentStore.readSeqStore[refReadId]));
                if (firstMatch) {
                    firstMatch = false;
                    fragmentStore.alignedReadStore[refMatchPosId].beginPos = bandOffset + profilePos;
                } else if (gapLen) {
                    diff += gapLen;
                    appendValue(fragmentStore.alignedReadStore[refMatchPosId].gaps, TGapAnchor(referencePos,referencePos + diff), Generous() );
                    gapLen = 0;
                }
            }
            ++fragIt;
        } while (fragIt != fragItEnd);
    }

    for (; referencePos < (TReadPos)length(fragmentStore.readSeqStore[refReadId]); ++referencePos) {
        //int numGaps =
        insertGap(matches, bandOffset + alignPos);
        ++alignPos;
    }
    fragmentStore.alignedReadStore[refMatchPosId].endPos = fragmentStore.alignedReadStore[refMatchPosId].beginPos + referencePos + diff;

    SEQAN_CONSENSUS_OPEN_PENALTY_FACTOR = tmp;

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

template<typename TMatchQuality, typename TRead, typename TOptions>
inline int
calibrateQuality(TRead & read, TMatchQuality & matchQuality, int originalQuality,TOptions & options)
{
    double epsilon = (double)(matchQuality.errors)/length(read);
    //return ( (int)((double) (epsilon*100.0) * (int) matchQuality.score / pow(2.0,(100.0*epsilon - 1.0)) + (double)originalQuality ) / (double)(100.0* epsilon + 1.0));
    return (int)( (int)((double) (epsilon*options.newQualityCalibrationFactor) * (int) matchQuality.score / pow(2.0,(20.0*options.newQualityCalibrationFactor - 1.0)) + (double)originalQuality ) / (double)(20.0* options.newQualityCalibrationFactor+ 1.0));

}


//////////////////////////////////////////////////////////////////////////////
// Output SNPs
template <
    typename TFragmentStore,
    typename TReadCounts,
    typename TReadCigars,
    typename TGenomeName,
    typename TFile,
    typename TOptions
>
void dumpVariantsRealignBatch(
    TFragmentStore              &fragmentStore,             // forward/reverse matches
    TReadCigars                 &,
    TReadCounts const           &,
    TGenomeName const           genomeID,                   // genome name
    typename TFragmentStore::TContigPos startCoord,         // startCoordinate + posOnGenomeInfix = real coordinate on whole chromosome
    typename TFragmentStore::TContigPos currStart,
    typename TFragmentStore::TContigPos currEnd,
    TFile                   &file,
    TFile                   &indelfile,
    TOptions                &options)
{

    typedef typename TFragmentStore::TAlignedReadStore  TMatches;
    typedef typename Value<TMatches>::Type              TMatch;
    typedef typename TFragmentStore::TAlignQualityStore TMatchQualities;
    //typedef typename Value<TMatchQualities>::Type       TMatchQuality;
    typedef typename Iterator<TMatches,Standard>::Type  TMatchIterator;

    typedef typename TFragmentStore::TReadSeqStore      TReads;
    typedef typename Value<TReads>::Type                TRead;

    typedef typename TFragmentStore::TContigStore       TContigStore;
    typedef typename Value<TContigStore>::Type          TContig;
    typedef typename TFragmentStore::TContigPos         TContigPos;
    typedef typename TFragmentStore::TContigSeq         TContigSeq;

    typedef Gaps<TContigSeq, AnchorGaps<typename TContig::TGapAnchors> >    TContigGaps;
    typedef Gaps<TRead, AnchorGaps<typename TMatch::TGapAnchors> >          TReadGaps;
    typedef typename Iterator<TContigGaps>::Type                            TContigGapIter;
    typedef typename Iterator<TReadGaps>::Type                              TReadGapIter;


    SEQAN_PROTIMESTART(dump_time);

    // matches need to be ordered according to genome position
    TReads &reads                   = fragmentStore.readSeqStore;
    TMatches &matches               = fragmentStore.alignedReadStore;
    TMatchQualities &matchQualities = fragmentStore.alignQualityStore;
    TContigPos genomeLen                = (TContigPos)length(fragmentStore.contigStore[0].seq);

    ::std::sort(begin(matches, Standard()), end(matches, Standard()), LessGPos<TMatch>());

    // make sure both output files are open
    if (!file.is_open() && !indelfile.is_open())
    {
        ::std::cerr << "Neither SNP nor indel output file is open" << ::std::endl;
        return;
    }

    // log file business
    ::std::ofstream logfile;
    if(options.outputLog != "")
    {
        logfile.open(toCString(options.outputLog), ::std::ios_base::out | ::std::ios_base::app);
        if (!logfile.is_open())
            ::std::cerr << "Failed to write to log file" << ::std::endl;
        logfile << "#stats for window " << currStart << " " << currEnd << " of " << genomeID << std::endl;
    }

    if(options._debugLevel > 1) ::std::cout << "Scanning chromosome " << genomeID << " window (" << currStart<<","<< currEnd << ") for SNPs..." << ::std::endl;

#ifdef SNPSTORE_DEBUG
    bool extraV = true;

    std::cout << genomeLen << " <-length genome \n";
    std::cout << length(fragmentStore.alignedReadStore) << " <-nummatches \n";
    std::cout << length(fragmentStore.readSeqStore) << " <-numreads \n";
    CharString str = "realignBatch";
    _dumpMatches(fragmentStore, str);
    std::cout << "startcoord=" << startCoord << std::endl;
#endif

//  std::fstream tmpfile;
//      tmpfile.open("Z:/seqan071010/projects/library/apps/chr4.before.sam", ::std::ios_base::out);
//      write(tmpfile, fragmentStore, Sam());
//      tmpfile.close()
;
    // convert matches in fragmentstore into a globally consistent alignment
    //Score<int> scoreType = Score<int>(0, -999, -1001, -1000); // levenshtein-score (match, mismatch, gapOpen, gapExtend)

#ifdef READS_454
    Score<int> scoreType = Score<int>(0, -3, -2, -2);   // (match, mismatch,gapExtend,gapOpen)
#else
    Score<int> scoreType = Score<int>(0, -1, -2, -10);  // (match, mismatch,gapExtend,gapOpen)
#endif

convertMatchesToGlobalAlignment(fragmentStore, scoreType, Nothing());

//  std::fstream tmpfile2;
//  tmpfile2.open("Z:/seqan071010/projects/library/apps/chr4.after.sam", ::std::ios_base::out);
//  write(tmpfile2, fragmentStore, Sam());
//  tmpfile2.close();

#ifdef SNPSTORE_DEBUG
    if(extraV)
    {
        TContigGaps contigGaps(fragmentStore.contigStore[0].seq, fragmentStore.contigStore[0].gaps);
        TContigPos maxPos = positionSeqToGap(contigGaps,length(fragmentStore.contigStore[0].seq)-1)+1;
        maxPos = _max(maxPos,(TContigPos)length(fragmentStore.contigStore[0].seq));
        std::cout << "maxPos visual = " << maxPos;
        std::cout << " genomeLen = " << genomeLen << std::endl;

        AlignedReadLayout layout;
        layoutAlignment(layout, fragmentStore);
        printAlignment(std::cout, Raw(), layout, fragmentStore, 0, (TContigPos)0, (TContigPos)maxPos, 0, 150);
    }
    ::std::cout << "done.\n" << std::flush;
    if(extraV)
    {
        CharString strstr = "befReal";
        _dumpMatches(fragmentStore,strstr);
    }
#endif

    //todo: do realignment here if options.realign, else work on global alignment
//  if(options.realign)
//  {
    // TODO: check out if a different scoring scheme makes more sense
    Score<int, WeightedConsensusScore<Score<int, FractionalScore>, Score<int, ConsensusScore> > > consScore;
    int bandWidth = 10; // ad hoc

    if(options._debugLevel > 1)
    {
        ::std::cout << "Realigning "<< length(matches)<<" reads to genome of length " <<genomeLen << std::flush;
//      ::std::cout << " StartCoord="<< startCoord << std::endl;
    }


#ifdef READS_454
    reAlign(fragmentStore,consScore,0,1,bandWidth,true);
#else
    reAlign(fragmentStore,consScore,0,1,bandWidth,true);
#endif

#ifdef SNPSTORE_DEBUG
    if(extraV)
    {
        CharString strstr = "befRefRal";
        _dumpMatches(fragmentStore,strstr);
        TContigGaps contigGaps(fragmentStore.contigStore[0].seq, fragmentStore.contigStore[0].gaps);
        TContigPos maxPos = positionSeqToGap(contigGaps,length(fragmentStore.contigStore[0].seq)-1)+1;
        maxPos = _max(maxPos,(TContigPos)length(fragmentStore.contigStore[0].seq));
        std::cout << "maxPos visual = " << maxPos;
        std::cout << " genomeLen = " << genomeLen << std::endl;
            AlignedReadLayout layout;
            layoutAlignment(layout, fragmentStore);
            printAlignment(std::cout, Raw(), layout, fragmentStore, 0, (TContigPos)0, (TContigPos)maxPos, 0, 150);
        }

#endif


    if(options._debugLevel > 1)::std::cout << "Realigning reads including reference..." << std::flush;


    unsigned numReads = length(matches)-1; // exclude reference sequence
    unsigned refId = length(matchQualities); // reference id (there may be more matchQs than matches due to pile up correction)


//
#ifndef  READS_454
    reAlign(fragmentStore,consScore,0,1,/*bandWidth*/5,false);
#else
//  reAlign(fragmentStore,consScore,0,1,/*bandWidth*/5,false);
//  realignReferenceToDiploidConsensusProfileDeleteSeqErrors(fragmentStore,refId,options);
    realignReferenceToDiploidConsensusProfile(fragmentStore,refId,options);
//    realignReferenceToDiploidConsensusProfileDeleteSeqErrors(fragmentStore,refId,options);
#endif

    if(options._debugLevel > 1) ::std::cout << "Finished realigning." << std::endl;



#ifdef SNPSTORE_DEBUG
    ::std::cout << "Realignment done.\n";
    if(extraV)
    {
        CharString strstr = "aftRefReal";
        _dumpMatches(fragmentStore,strstr);
        TContigGaps contigGaps(fragmentStore.contigStore[0].seq, fragmentStore.contigStore[0].gaps);
        TContigPos maxPos = positionSeqToGap(contigGaps,length(fragmentStore.contigStore[0].seq)-1)+1;
        maxPos = _max(maxPos,(TContigPos)length(fragmentStore.contigStore[0].seq));
        std::cout << "maxPos visual = " << maxPos;
        std::cout << " genomeLen = " << genomeLen << std::endl;
        AlignedReadLayout layout;
        layoutAlignment(layout, fragmentStore);
        printAlignment(std::cout, Raw(), layout, fragmentStore, 0, (TContigPos)0, (TContigPos)maxPos, 0, 150);
        std::fstream tmpfile3;
        tmpfile3.open("test.realigned.sam", ::std::ios_base::out);
        //write(tmpfile3, fragmentStore, Sam());
        tmpfile3.close();
    }
    ::std::cout << "done." << std::flush;

    //  std::fstream tmpfile2;
    //  tmpfile2.open("tmpfile_realigned.sam", ::std::ios_base::out);
    //  write(tmpfile2, fragmentStore, Sam());
    //  tmpfile2.close();

#endif
    //}


    // forward match qualities

    String<int> columnQualityF;         resize(columnQualityF,5);
    String<unsigned> countF;            resize(countF,5);
    String<CharString> qualityStringF;  resize(qualityStringF,5);

    // reverse match qualities
    String<int> columnQualityR;         resize(columnQualityR,5);
    String<unsigned> countR;            resize(countR,5);
    String<CharString> qualityStringR;  resize(qualityStringR,5);

    FunctorComplement<Dna5> f;

    // sort reads according to begin position
    sortAlignedReads(fragmentStore.alignedReadStore, SortBeginPos());
    TMatchIterator matchIt      = begin(matches, Standard());
    TMatchIterator matchItEnd   = end(matches, Standard());

    // look for reference sequence and move it to the end of alignedreads
    // todo: only do this when realignment was done
    bool refFound = false;
    TMatchIterator matchItKeep = matchIt;
    TMatch tempRef;
    while(matchIt != matchItEnd)
    {
        if((*matchIt).readId == refId) // this is the reference
        {
            refFound = true;
            tempRef = *matchIt;
            matchItKeep = matchIt;
            ++matchIt;
            continue;
        }
        if(refFound)
        {
            *matchItKeep = *matchIt; // matchItKeep lags behind by one match
            ++matchIt;++matchItKeep;
        }
        else ++matchIt;
    }
    *matchItKeep = tempRef;
    SEQAN_ASSERT(refFound);

#ifdef SNPSTORE_DEBUG
    if(!refFound) ::std::cout << "ref not Found!\n";
    ::std::cout << "done looking for ref." << std::flush << std::endl;
#endif

    matchIt     = begin(matches, Standard());
    matchItEnd  = end(matches, Standard());
    matchItEnd--; // exclude reference sequence

    TRead       &reference = fragmentStore.readSeqStore[fragmentStore.alignedReadStore[numReads].readId]; // last read is reference sequence
    TReadGaps   referenceGaps(reference, fragmentStore.alignedReadStore[numReads].gaps);
    TContigPos      refStart = (TContigPos)fragmentStore.alignedReadStore[numReads].beginPos;
    TContigGaps contigGaps(fragmentStore.contigStore[0].seq, fragmentStore.contigStore[0].gaps);
    SingleBaseVariant snp = {0,0,0,0,0,0};

#ifdef SNPSTORE_DEBUG
    if(extraV)
    {
        ::std::cout << "lengthrefgaps=" << length(referenceGaps)<< std::endl;
        ::std::cout << "length(genome) = " << genomeLen << " length(ref)=" << length(reference) << std::endl;
    }
#endif


    // for indels:
    // i1 keeps track of consensus character
    // i2 keeps track of coverage (last 8 bits) and indelcount (first 8 bits)
    String<Pair<short unsigned,short unsigned> > indelConsens;
    resize(indelConsens,refStart + length(referenceGaps));
    for(unsigned i = 0; i < refStart + length(referenceGaps) ; ++i)
    {
        indelConsens[i].i1 = 6;
        indelConsens[i].i2 = 0;
    }

    if(options._debugLevel>1) std::cout << "Start inspecting alignment..." << std::endl;
    // now walk through the reference sequence in gaps view space,
    // i.e. position may be a gap
    // example:
    // Ref      ACCGTGCACTAGCATCATT--ACTAGCATCATA
    // Reads    ACCGTACA--AGCATCAT
    //              TACA--AGCATCATT--ACT
    //                          ATTTTACTAGCATCATA
    for(TContigPos candidateViewPos = refStart; candidateViewPos < refStart + (TContigPos)length(referenceGaps); ++candidateViewPos)
    {
        // first check if reference has a gap (potential insertion in reads) at this position
        TContigGapIter refIt = iter(referenceGaps,candidateViewPos-refStart);
        bool refGap = false;
        if(isGap(refIt)) refGap = true;

        //get position in sequence space
        TContigPos candidatePos = positionGapToSeq(referenceGaps, candidateViewPos-refStart);

        // not in the current window yet
        if(candidatePos + startCoord < currStart) continue;
        // not in the current window anymore
        if(candidatePos + startCoord >= currEnd) break;

        Dna5 refBase = reference[candidatePos]; // what happens if refGap==true, esp. for leading gaps?
        if(refBase=='N') continue;

#ifdef SNPSTORE_DEBUG
        std::cout << "candidateViewPos = " << candidateViewPos << std::endl;
        std::cout << "candidatePos = " << candidatePos << std::endl;
        std::cout << "candidatePosMitStart = " << candidatePos + startCoord << " refBase = " << refBase << std::endl;
        if(refGap) std::cout << "refGap!" << std::endl;
        bool extraVVVV = false;
        if(candidatePos + startCoord == 19388258) extraVVVV=true;
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
        if(extraVVVV) std::cout <<"cov=" << coverage << std::endl;
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

            countR[t] = 0;
            columnQualityR[t] = 0;
            clear(qualityStringR[t]);
        }

        bool observedAtLeastOneMut = false;
        int numIndelsObservedF = 0;  // if refGap then this counts the number of insertions on forward
        int numIndelsObservedR = 0;  // if refGap then this counts the number of insertions on reverse
                        // else it counts the number of deletions
        int indelQualF = 0;
        int indelQualR = 0;

        int positionCoverage = 0;   // how many reads actually span the position?

        // now check reads
        while(matchIt != matchRangeEnd)
        {
            TContigPos currViewBegin = _min((*matchIt).beginPos,(*matchIt).endPos);
            TContigPos currViewEnd = _max((*matchIt).beginPos,(*matchIt).endPos);

            // make sure this match is really spanning the position
            if(!(currViewBegin <= candidateViewPos && candidateViewPos < currViewEnd))
            {
                ++matchIt;
                continue;
            }
            ++positionCoverage;

            char orientation = ((*matchIt).beginPos > (*matchIt).endPos) ? 'R' : 'F';

            TReadGaps readGaps(reads[(*matchIt).readId],(*matchIt).gaps);
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
            extraVVVV = true;
            std::cout << "ReadPos = " << readPos << std::endl;
#endif

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
                quality = getQualityValue(reads[(*matchIt).readId][readPos]) ;
                if(options.newQualityCalibrationFactor > 0.0001 ) quality = calibrateQuality(reads[(*matchIt).readId], matchQualities[(*matchIt).id],quality,options);

                if(!options.useBaseQuality && quality > (int)matchQualities[(*matchIt).id].score)
                {   // dont trust the quality of this position more
                    // than the average quality of this read
                    quality = (int) matchQualities[(*matchIt).id].score;
                }


                if(orientation == 'F')
                {
                    columnQualityF[ordValue(candidateBase)] += quality;
                    ++countF[ordValue(candidateBase)];
                    appendValue(qualityStringF[ordValue(candidateBase)],(char)(quality+33),Generous());
                }
                else
                {
                    columnQualityR[ordValue(candidateBase)] += quality;
                    ++countR[ordValue(candidateBase)];
                    appendValue(qualityStringR[ordValue(candidateBase)],(char)(quality+33),Generous());
                }
            }
            else
            {   //potential deletions

                if(!refGap)
                {
                    readPos = positionGapToSeq(readGaps,candidateViewPos - currViewBegin);
#ifdef SNPSTORE_DEBUG
                    if(extraVVVV) std::cout <<"del readPos = " << readPos  << "readlength=" << length(reads[(*matchIt).readId]) << std::endl;
#endif
                    if(orientation == 'R')
                        readPos = length(reads[(*matchIt).readId]) - readPos;
#ifdef SNPSTORE_DEBUG
                    if(extraVVVV) std::cout <<"del readPos = " << readPos  << "readlength=" << length(reads[(*matchIt).readId]) << std::endl;
#endif
                    quality = (getQualityValue(reads[(*matchIt).readId][readPos-1])  + getQualityValue(reads[(*matchIt).readId][readPos]) ) / 2;
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
        int numIndelsObserved = numIndelsObservedF + numIndelsObservedR;
#ifdef SNPSTORE_DEBUG
        if(extraVVVV)
        {
            std::cout << "posCov=" << positionCoverage << "numIndels = " << numIndelsObserved << std::endl;
            if(observedAtLeastOneMut) std::cout << "observed at least one mut " << std::endl;
        }
#endif

        // too few reads actually cover the position
        if(positionCoverage < (int)options.minCoverage)
            continue;

        //all observed bases match the reference allele or there were too few indels
        //if(!observedAtLeastOneMut && numIndelsObserved< options.indelCountThreshold)
        //  continue;


        // do SNP calling
        if(file.is_open() && observedAtLeastOneMut)
        {
            bool isSnp = true;

            // coverage depth
            int refAllele = ordValue(reference[candidatePos]);
            unsigned realCoverageF = countF[0] + countF[1] +countF[2] +countF[3] +countF[4];
            unsigned realCoverageR = countR[0] + countR[1] +countR[2] +countR[3] +countR[4];
            unsigned realCoverage  = realCoverageF + realCoverageR;

            // Coverage too low after discarding Ns and gaps from alignment column
            if(realCoverage<options.minCoverage) isSnp = false;

            // is the min. number of different read positions supporting the mutation met?
            if(isSnp && options.minDifferentReadPos > 0 && readPosMap.size() < options.minDifferentReadPos)
                isSnp = false;

            // do genotype calling
            if(isSnp && options.method == 1)
                isSnp = _doSnpCall(countF,countR,qualityStringF,qualityStringR,refAllele,options,snp,MaqMethod()
#ifdef SNPSTORE_DEBUG_CANDPOS
                ,(int) candidatePos + startCoord
#endif
                );
            else if(isSnp && options.method == 0)
                isSnp = _doSnpCall(countF,countR,columnQualityF,columnQualityR,refAllele,options,snp,ThresholdMethod() );

            // write SNP to file
            if(isSnp && (snp.called || options.outputFormat == 0))
                _writeSnp(file,snp,qualityStringF,qualityStringR,refAllele,genomeID,candidatePos+startCoord,realCoverage,options);
        }

        bool isIndel = true;
        if (!indelfile.is_open() || numIndelsObserved < (int)options.indelCountThreshold) // count threshold
            isIndel = false;

        if(isIndel && ((float)numIndelsObserved/(float)positionCoverage) < options.indelPercentageT)
            isIndel = false;

        if (isIndel && options.minDifferentReadPos > 0 && indelReadPosMap.size() < options.minDifferentReadPos) // minDiffReadPos criterium
            isIndel = false;

        if(isIndel && options.bothIndelStrands && ((float)numIndelsObservedF==0 || numIndelsObservedR==0))
            isIndel = false;

          // possible deletion, insertion quality will be handled later
        int avgIndelQuality = 0;
        if(isIndel && numIndelsObserved > 0) avgIndelQuality = (int)((double)indelQualF+indelQualR)/numIndelsObserved;
        if(isIndel && !refGap && (float)avgIndelQuality < options.indelQualityThreshold)
            isIndel = false;


        // do indel calling  //TODO: for 454 data look at sequence content! for >=6bp homopolymer runs needs to be more strict
        //if (indelfile.is_open() && numIndelsObserved >= (int)options.indelCountThreshold // count threshold
        //    && ((float)numIndelsObserved/(float)positionCoverage) >= options.indelPercentageT // percentage threshold
        //    && (options.minDifferentReadPos == 0 || indelReadPosMap.size() >= options.minDifferentReadPos)) // minDiffReadPos criterium
        if(isIndel)
        {
            bool indelCalled = true;
            char mostCommonBase = 5; // 5 represents gap char "-", potential deletion
            if(refGap) // potential insertion
            {
#ifdef SNPSTORE_DEBUG
                if(extraVVVV) std::cout << "potential insertion" << std::endl;
#endif
                SEQAN_ASSERT(!observedAtLeastOneMut);
                mostCommonBase = 0;
                unsigned maxCount = countF[0] + countR[0];
                for(unsigned j = 0; j < length(countF); ++j)
                    if(countF[j] + countR[j] > maxCount)
                    {
                        maxCount = countF[j] + countR[j];
                        mostCommonBase = j;
                    }
                // number of ref reads + insertion reads
                if(maxCount < options.indelCountThreshold) indelCalled = false;
               // if(!indelCalled) std::cout << "verloren1\n";
                if(((double)maxCount/positionCoverage) < options.indelPercentageT) indelCalled = false;
               // if(!indelCalled) std::cout << "verloren2\n";
                if((double)(maxCount+positionCoverage-numIndelsObserved)/positionCoverage < options.minExplainedColumn) indelCalled = false;
               // if(!indelCalled) std::cout << "maxCount = " << maxCount << " positionCoverage = "<< positionCoverage << "numIndelsObserved = " << numIndelsObserved <<" verloren3\n";

                // check quality
                //avgIndelQuality = (int) ((double)columnQualityR[mostCommonBase] +columnQualityF[mostCommonBase])/numIndelsObserved; // divide by all observed insertions?
                avgIndelQuality = (int) ((double)columnQualityR[mostCommonBase] +columnQualityF[mostCommonBase])/maxCount;  // or just candidate insertions?
                if ((float)avgIndelQuality < options.indelQualityThreshold) indelCalled = false;
               // if(!indelCalled) std::cout << "verloren4\n";
                numIndelsObserved = maxCount;
            }
            if(indelCalled)
            {
                int bothStrandsObserved = 0;
                if(numIndelsObservedF>0 && numIndelsObservedR>0) // remember that both strands were observed --> increase quality of call later
                    bothStrandsObserved = 1;
                indelConsens[candidateViewPos].i1 = avgIndelQuality << 8 | bothStrandsObserved << 4 | mostCommonBase;
#ifdef SNPSTORE_DEBUG
                if(extraVVVV) std::cout << "mosCommonBase = " << (int)mostCommonBase << std::endl;
#endif
                if(positionCoverage > 255) //downscaling if numbers get too large
                {
                    numIndelsObserved *= (int)((float)255.0/(float)positionCoverage);
                    positionCoverage = 255;
#ifdef SNPSTORE_DEBUG
                    if(extraVVVV) std::cout << "downscaled to " << numIndelsObserved << std::endl;
#endif
                }
                indelConsens[candidateViewPos].i2 = numIndelsObserved << 8 | positionCoverage;
            }
        }
    }

    //CharString chrPrefix = "";
    CharString chrPrefix = "chr"; // should check if "chr" is already part of chromosome names (usually not)


    ::std::string runID = options.runID;

#ifdef SNPSTORE_DEBUG
    // write out indels
    for(unsigned i = refStart; i < refStart + length(referenceGaps); ++i)
        std::cout << (indelConsens[i].i1 & 7);
    std::cout << std::endl;
#endif

    if(indelfile.is_open()) //indelcalling
    {
        if(options._debugLevel > 1) std::cout << "Calling indels..." << std::endl;
        TContigPos candidateViewPos = refStart;
        Dna5String insertionSeq;
        while(candidateViewPos < refStart + (TContigPos)length(referenceGaps))
        {

            bool bsi = true;
            if(candidateViewPos < refStart + (TContigPos)length(referenceGaps) &&
                (indelConsens[candidateViewPos].i1 & 7) == 6) // not a relevant position
            {
#ifdef SNPSTORE_DEBUG
                ::std::cout << candidateViewPos << "not relevant for indels" <<  std::endl;
#endif
                ++candidateViewPos;
                continue;
            }

            //get position in sequence space
            TContigPos candidatePos = positionGapToSeq(referenceGaps, candidateViewPos-refStart);
            int indelSize = 0;
            unsigned depth = 0;
            float percentage = 0.0;
            int quality = 0;
            // gap position
#ifdef SNPSTORE_DEBUG
            ::std::cout << candidateViewPos << " indel?" <<  std::endl;
#endif
            while(candidateViewPos < refStart + (TContigPos)length(referenceGaps) && // shouldnt happen actually
                ((indelConsens[candidateViewPos].i1 & 7) == 5  ||        // deletion in consens
                ((indelConsens[candidateViewPos].i1 & 7) == 6  && isGap(referenceGaps, candidateViewPos-refStart)))  // position in consensus is the same as in reference
                )                                                   // and reference is a gap (same candidatePosition as before)
            {
#ifdef SNPSTORE_DEBUG
                ::std::cout << startCoord + candidateViewPos << " del!!" <<  std::endl;
#endif
                if((indelConsens[candidateViewPos].i1 & 7) == 5)
                {
                    ++indelSize;
                    depth += (indelConsens[candidateViewPos].i2 & 255);
                    percentage += (float)((indelConsens[candidateViewPos].i2 >> 8) & 255);
                    quality += ((indelConsens[candidateViewPos].i1 >> 8) & 255);
                    bsi = bsi && (((indelConsens[candidateViewPos].i1 >> 4) & 1) == 1 );
                }
                ++candidateViewPos;

            }
            if(indelSize>0)
            {
                // check for the longest adjacent stretch of homopolymers
                int homoLength = checkSequenceContext(reference,candidatePos,indelSize);
                if(homoLength <= options.maxPolymerRun)
                {
                    percentage /= depth;                                // low coverage positions get a lower weight
                    depth = (depth + (indelSize >> 1)) / indelSize;     // coverage is spread over all positions
                    quality = (quality + (indelSize >> 1)) / indelSize; // quality is spread over all positions
                    int indelQ = (int)(quality * percentage);
                    if (!bsi) indelQ /= 2;

                    //print deletion
                    indelfile << chrPrefix << genomeID << '\t' << runID << "\tdeletion\t";
                    indelfile << candidatePos + startCoord + options.positionFormat  << '\t';
                    indelfile << candidatePos + startCoord + options.positionFormat + indelSize - 1;
                    indelfile << "\t" << percentage;
                    indelfile << "\t+\t.\tID=" << candidatePos + startCoord + options.positionFormat ;
                    indelfile << ";size=" << indelSize;
                    indelfile << ";count=" << (int)(percentage * depth + 0.00001);
                    indelfile << ";depth=" << depth;
                    indelfile << ";quality=" << indelQ;
                    indelfile << ";homorun=" << homoLength;
                    if(bsi)indelfile << ";bsi";
                    indelfile << ";seqContext=" << infix(reference,_max((int)0,(int)candidatePos-6),_min((int)candidatePos+indelSize+6,(int)length(reference)));
                    if(percentage <= options.indelHetMax) indelfile << ";geno=het";
                    else indelfile << ";geno=hom";

                    //if(splitSupport>0) indelfile << ";splitSupport=" << splitSupport;
                    indelfile << std::endl;
                }


                //reset
                candidatePos = positionGapToSeq(referenceGaps, candidateViewPos-refStart);
                indelSize = 0;
                depth = 0;
                percentage = 0.0;
                quality = 0;
                bsi = true;
            }
            clear(insertionSeq);
            while(candidateViewPos < refStart + (TContigPos)length(referenceGaps) && // shouldnt happen actually
                ((indelConsens[candidateViewPos].i1 & 7) < 5 ||          // insertion in consensus
                ((indelConsens[candidateViewPos].i1 & 7) == 6  && candidatePos == positionGapToSeq(referenceGaps, candidateViewPos-refStart)))   // position in consensus is the same as in reference
                )                                                   // and reference is a gap (same candidatePosition as before)
            {
#ifdef SNPSTORE_DEBUG
                ::std::cout << candidateViewPos << " ins!!!" <<  std::endl;
#endif
                if((indelConsens[candidateViewPos].i1 & 7) < 5)
                {
                    --indelSize;
                    depth += (indelConsens[candidateViewPos].i2 & 255);
                    percentage += (float)((indelConsens[candidateViewPos].i2 >> 8) & 255);
                    quality += ((indelConsens[candidateViewPos].i1 >> 8) & 255);
                    appendValue(insertionSeq,(Dna5)(indelConsens[candidateViewPos].i1 & 7));
                    bsi = bsi && (((indelConsens[candidateViewPos].i1 >> 4) & 1) == 1 );
                }
                ++candidateViewPos;
            }
            if(indelSize<0)
            {
                int homoLength = checkSequenceContext(reference,candidatePos,indelSize);
                if(homoLength <= options.maxPolymerRun)
                {
                    unsigned absIndelSize = -indelSize;
                    percentage /= depth;                                        // low coverage positions get a lower weight
                    depth = (depth + (absIndelSize >> 1)) / absIndelSize;       // coverage is spread over all positions
                    quality = (quality + (absIndelSize >> 1)) / absIndelSize;   // quality is spread over all positions
                    int indelQ = (int)(quality * percentage);
                    if (!bsi) indelQ /= 2;

                    //print insertion
                    indelfile << chrPrefix <<genomeID << '\t' << runID << "\tinsertion\t";
                    indelfile << candidatePos + startCoord + options.positionFormat - 1 << '\t';
                    indelfile << candidatePos + startCoord;// + options.positionFormat; //VORSICHT!!!
                    indelfile << "\t" << percentage;
                    indelfile << "\t+\t.\tID=" << candidatePos + startCoord + options.positionFormat;
                    indelfile << ";size=" << indelSize;
                    indelfile << ";count=" << (int)(percentage * depth + 0.00001);
                    indelfile << ";seq="<< insertionSeq;
                    indelfile << ";depth=" << depth;
                    indelfile << ";quality=" << indelQ;
                    indelfile << ";homorun=" << homoLength;
                    if(bsi)indelfile << ";bsi";
                    indelfile << ";seqContext=" << infix(reference,_max((int)0,(int)candidatePos-6),_min((int)candidatePos+6,(int)length(reference)));
                    if(percentage <= options.indelHetMax) indelfile << ";geno=het";
                    else indelfile << ";geno=hom";

                    //if(splitSupport>0) indelfile << ";splitSupport=" << splitSupport;
                    indelfile << std::endl;
                }
                //resetting will be done in next round
            }

        }
        if(options._debugLevel > 1) std::cout << "Finished calling indels..." << std::endl;

    }

    if(options._debugLevel>1) std::cout <<"Finished scanning window.\n"<<std::flush;

    if((options.outputLog != "") && logfile.is_open())
        logfile.close();


}




//////////////////////////////////////////////////////////////////////////////
// Output SNPs
template <
    typename TFragmentStore,
    typename TReadCigars,
    typename TReadCounts,
    typename TGenomeName,
    typename TFile,
    typename TOptions
>
void dumpSNPsBatch(
    TFragmentStore              &fragmentStore,             // forward/reverse matches
    TReadCigars             &,
    TReadCounts const           &readCounts,
    TGenomeName const           genomeID,                   // genome name
    typename TFragmentStore::TContigPos startCoord,         // startCoordinate + posOnGenomeInfix = real coordinate on whole chromosome
    typename TFragmentStore::TContigPos currStart,
    typename TFragmentStore::TContigPos currEnd,
    TFile               &file,
    TOptions            &options)
{

    typedef typename TFragmentStore::TAlignedReadStore  TMatches;
    typedef typename Value<TMatches>::Type              TMatch;
    typedef typename TFragmentStore::TAlignQualityStore TMatchQualities;
    //typedef typename Value<TMatchQualities>::Type       TMatchQuality;
    typedef typename TFragmentStore::TReadSeqStore      TReads;
    //typedef typename Value<TReads>::Type                TRead;
    typedef typename TFragmentStore::TContigPos         TContigPos;
    typedef typename TFragmentStore::TContigSeq         TContigSeq;
    typedef typename Iterator<TMatches,Standard>::Type  TMatchIterator;

    SEQAN_PROTIMESTART(dump_time);
    //options._debugLevel = 2;
    String<char> toIupac = "AMRWMCSYRSGKWYKT";
    //std::cout << "Hier\n";
    // matches need to be ordered accordign to genome position
    TReads &reads                   = fragmentStore.readSeqStore;
    TMatches &matches               = fragmentStore.alignedReadStore;
    TMatchQualities &matchQualities = fragmentStore.alignQualityStore;
    TContigSeq &genome              = fragmentStore.contigStore[0].seq;

    ::std::sort(begin(matches, Standard()), end(matches, Standard()), LessGPos<TMatch>());

    Align<String<Dna5>, ArrayGaps> align;
    Score<int> scoreType = Score<int>(0, -999, -1001, -1000);   // levenshtein-score (match, mismatch, gapOpen, gapExtend)
    resize(rows(align), 2);

    if (!file.is_open())
    {
        ::std::cerr << "Output file is not open" << ::std::endl;
        return;
    }

    ::std::ofstream logfile;
    if(options.outputLog != "")
    {
        logfile.open(toCString(options.outputLog), ::std::ios_base::out | ::std::ios_base::app);
        if (!logfile.is_open())
            ::std::cerr << "Failed to write to log file" << ::std::endl;
        logfile << "#stats for window " << currStart << " " << currEnd << " of " << genomeID << std::endl;
    }

    TMatchIterator matchIt  = begin(matches, Standard());
    TMatchIterator matchItEnd   = end(matches, Standard());
    //matchItEnd--;
    unsigned countLowerMQ = 0, countHigherMQ = 0;

    if(options._debugLevel > 1) ::std::cout << "Scanning chromosome " << genomeID << " window (" << currStart<<","<< currEnd << ") for SNPs..." << ::std::endl;

    // forward match qualities
    String<int> columnQualityF;         resize(columnQualityF,5);
    String<unsigned> countF;            resize(countF,5);
    String<CharString> qualityStringF;  resize(qualityStringF,5);

    // reverse match qualities
    String<int> columnQualityR;         resize(columnQualityR,5);
    String<unsigned> countR;            resize(countR,5);
    String<CharString> qualityStringR;  resize(qualityStringR,5);

    // both
    String<unsigned> count;             resize(count,5);
    String<unsigned> columnQuality;     resize(columnQuality,5);

#ifdef SNPSTORE_DEBUG
    bool extraV = false;
#endif
    SingleBaseVariant snp = {0,0,0,0,0,0};

    for(TContigPos candidatePos = 0; candidatePos < (TContigPos)length(genome); ++candidatePos)
    {
//      if(options._debugLevel > 1) ::std::cout << "Next pos\n";

        if(candidatePos + startCoord < currStart) continue;

        // not in the current window anymore
        if(candidatePos + startCoord >= currEnd)
            break;

        Dna5 refBase = genome[candidatePos];
        if(refBase=='N') continue;

#ifdef SNPSTORE_DEBUG
        ::std::cout << "candPos=" << candidatePos + startCoord << ::std::endl;
        if(candidatePos + startCoord == 861196)
            ::std::cout << "ab jetzt.." << ::std::flush;
#endif

        Dna5 candidateBase;
        int quality;

//      if(options._debugLevel > 1)std::cout << candidatePos+startCoord << "<-candidatePos\n";
        for(unsigned t=0;t<5;++t)
        {
            countF[t] = 0;
            columnQualityF[t] = 0;
            clear(qualityStringF[t]);

            countR[t] = 0;
            columnQualityR[t] = 0;
            clear(qualityStringR[t]);
        }

        //find range of relevant read matches
        while(matchIt != matchItEnd &&  _max((*matchIt).endPos,(*matchIt).beginPos) <= candidatePos)
            ++matchIt;
        TMatchIterator matchRangeBegin = matchIt;
        while(matchIt != matchItEnd &&  _min((*matchIt).endPos,(*matchIt).beginPos)  <= candidatePos)
            ++matchIt;
        TMatchIterator matchRangeEnd = matchIt;
        matchIt = matchRangeBegin;

        int coverage = matchRangeEnd-matchRangeBegin;
        if(coverage<(int)options.minCoverage) continue; // coverage too low

        if(options._debugLevel > 1)::std::cout << "Match range:" << matchRangeEnd - matchRangeBegin << ::std::endl;
#ifdef SNPSTORE_DEBUG
        if(extraV)
        {
            for (TMatchIterator tempIt = matchRangeBegin; tempIt != matchRangeEnd; ++tempIt)
                ::std::cout << reads[(*tempIt).readId]<<"\n";
        }
#endif
        std::set<unsigned> readPosMap;
        bool observedAtLeastOneMut = false;

        while(matchIt != matchRangeEnd)
        {
            TContigPos currentBegin = _min((*matchIt).beginPos,(*matchIt).endPos);
            TContigPos currentEnd   = _max((*matchIt).beginPos,(*matchIt).endPos);
            char orientation = ((*matchIt).beginPos > (*matchIt).endPos) ? 'R' : 'F';

#ifdef SNPSTORE_DEBUG
            if(extraV)
            {
                ::std::cout <<"currentBegin = "<<currentBegin << "\n";
                ::std::cout <<"currentEnd = "<<currentEnd << "\n";
            }
#endif
            if(!(currentBegin <= candidatePos && candidatePos < currentEnd))// this match is not really spanning the position
            {                                                               // (can happen because of indels or variable-length reads)
                ++matchIt;
                continue;
            }
            /*if(!empty(readCigars[(*matchIt).readId]))//splitRead, dont use for snp calling for now
            {
                ++matchIt;
                continue;
            }*/
            // or do use and only :
            // if(candidatePos lies in split gap) continue;

            // do edit alignment
            if((int)matchQualities[(*matchIt).id].pairScore == 1 /*!empty((*matchIt).gaps)*/) // splitReads: hamming: pairScore=0
            {
                Dna5String gInf = infix(genome, currentBegin, currentEnd);
                if (orientation == 'R')
                    reverseComplement(gInf);

                assignSource(row(align, 0), reads[(*matchIt).readId]);
                assignSource(row(align, 1), gInf);
                globalAlignment(align, scoreType);  //splitReads: get alignment from cigar string
            }

            if (orientation == 'R')
            {
                FunctorComplement<Dna5> f;

                int readPos = currentEnd - candidatePos - 1;
                if ((int)matchQualities[(*matchIt).id].pairScore == 1 /*!empty((*matchIt).gaps)*/)
                    readPos = getReadPos(align,readPos,false); //

#ifdef SNPSTORE_DEBUG
                if(extraV) std::cout << "readPosNacher = " << readPos << std::endl;
#endif
                if(readPos != -1) //-1 indicates gap
                {
                    candidateBase = f((Dna5)reads[(*matchIt).readId][readPos]);
#ifdef SNPSTORE_DEBUG
                    if(extraV) std::cout << candidateBase << "candBase\n";
#endif
                    quality = getQualityValue(reads[(*matchIt).readId][readPos]) ;
                    if(candidateBase != refBase)
                    {
                        observedAtLeastOneMut = true;
                        if(options.minDifferentReadPos > 0)
                            if((unsigned)(length(reads[(*matchIt).readId]) - readPos) > options.excludeBorderPos  && (unsigned) readPos >= options.excludeBorderPos )
                                readPosMap.insert(readPos);

                    }
                    if(options.newQualityCalibrationFactor > 0.0001) quality = calibrateQuality(reads[(*matchIt).readId], matchQualities[(*matchIt).id],quality,options);

                    if(!options.useBaseQuality && quality > (int)matchQualities[(*matchIt).id].score)   // dont trust the quality of this position more
                    {                                                                               // than the average quality of this read
                        quality = (int) matchQualities[(*matchIt).id].score;
                        ++countLowerMQ;
                    }
                    else ++countHigherMQ;
                    //if(quality < 0 || quality > 40)::std::cout << "falschQ candPos = " << candidatePos + startCoord << std::endl;

                    unsigned tmpCount = 1;
                    if(!empty(readCounts)) tmpCount = readCounts[(*matchIt).readId];
                    for (unsigned k = 0; k < tmpCount; ++k)
                    {
                        columnQualityR[ordValue(candidateBase)] += quality;
                        ++countR[ordValue(candidateBase)];
                        appendValue(qualityStringR[ordValue(candidateBase)],(char)(quality+33),Generous());
                    }
                }
            }
            else
            {
                int readPos = candidatePos - currentBegin;

                if((int)matchQualities[(*matchIt).id].pairScore == 1 /*!empty((*matchIt).gaps)*/)
                    readPos = getReadPos(align,readPos);

                if(readPos != -1) //-1 indicates gap
                {
                    candidateBase = (Dna5)reads[(*matchIt).readId][readPos];
                    quality = getQualityValue(reads[(*matchIt).readId][readPos]) ;
                    if(candidateBase != refBase)
                    {
                        observedAtLeastOneMut = true;
                        if(options.minDifferentReadPos > 0)
                            if((unsigned)(length(reads[(*matchIt).readId]) - readPos) > options.excludeBorderPos  && (unsigned) readPos >= options.excludeBorderPos )
                                readPosMap.insert(readPos);
                    }
                    if(options.newQualityCalibrationFactor > 0.0001) quality = calibrateQuality(reads[(*matchIt).readId], matchQualities[(*matchIt).id],quality,options);

                    if(!options.useBaseQuality && quality > (int) matchQualities[(*matchIt).id].score)
                    {
                        quality = (int) matchQualities[(*matchIt).id].score;
                        ++countLowerMQ;
                    }
                    else ++countHigherMQ;
                    //if(quality < 0 || quality > 40)::std::cout << "falschQ candPos = " << candidatePos + startCoord << std::endl;

                    unsigned tmpCount = 1;
                    if(!empty(readCounts)) tmpCount = readCounts[(*matchIt).readId];
                    for (unsigned k = 0; k < tmpCount; ++k)
                    {
                        ++countF[ordValue(candidateBase)];
                        columnQualityF[ordValue(candidateBase)] += quality;
                        appendValue(qualityStringF[ordValue(candidateBase)],(char)(quality+33),Generous());
                    }
                }
            }
            ++matchIt;
        }
        matchIt = matchRangeBegin; //set iterator back to where we started from, same matches might be involved in next cand pos

        if(!observedAtLeastOneMut) continue; //all observed bases match the reference allele

        // do SNP calling
        bool isSnp = true;

        // coverage depth
        int refAllele = ordValue(genome[candidatePos]);
        unsigned realCoverageF = countF[0] + countF[1] +countF[2] +countF[3] +countF[4];
        unsigned realCoverageR = countR[0] + countR[1] +countR[2] +countR[3] +countR[4];
        unsigned realCoverage  = realCoverageF + realCoverageR;

        // Coverage too low after discarding Ns and gaps from alignment column
        if(realCoverage<options.minCoverage) isSnp = false;

        // is the min. number of different read positions supporting the mutation met?
        if(isSnp && options.minDifferentReadPos > 0 && readPosMap.size() < options.minDifferentReadPos)
            isSnp = false;

        // do genotype calling
        if(isSnp && options.method == 1)
            isSnp = _doSnpCall(countF,countR,qualityStringF,qualityStringR,refAllele,options,snp,MaqMethod()
#ifdef SNPSTORE_DEBUG_CANDPOS
                ,(int) candidatePos + startCoord
#endif
            );
        else if(isSnp && options.method == 0)
            isSnp = _doSnpCall(countF,countR,columnQualityF,columnQualityR,refAllele,options,snp,ThresholdMethod() );

        // write SNP to file
        if(isSnp && (snp.called || options.outputFormat == 0))
            _writeSnp(file,snp,qualityStringF,qualityStringR,refAllele,genomeID,candidatePos+startCoord,realCoverage,options);
    }

    if(options._debugLevel>1) std::cout <<"Finished scanning window.\n"<<std::flush;

    if((options.outputLog != "") && logfile.is_open())
        logfile.close();


    return;

}




template<typename TAlign, typename TString, typename TPosition>
void
getIndels(TAlign & align,TString &insertions, TString &deletions, TPosition begin_, TPosition end_)
{

    //typedef typename Source<TAlign>::Type TSource;
    //typedef typename Iterator<TSource, Rooted>::Type TStringIterator;

    typedef typename Row<TAlign>::Type TRow;
    typedef typename Iterator<TRow, Rooted>::Type TAlignIterator;

    TAlignIterator ali_it0 = iter(row(align,0),begin_);
    TAlignIterator ali_it1 = iter(row(align,1),begin_);
    TAlignIterator ali_it0_stop = iter(row(align,0),end_);
    TAlignIterator ali_it1_stop = iter(row(align,1),end_);
    // TStringIterator readBase = begin(source(row(align,0)));
    //std::cout << "getting cigar line\n";//ali0 len = " <<ali_it0_stop-ali_it0 << " \t ali1 len = "<<ali_it1_stop-ali_it1<<"\n";
    int readPos = 0;
    int refPos = 0;
    while(ali_it0 != ali_it0_stop && ali_it1 != ali_it1_stop)
    {
        int inserted = 0;
        int deleted = 0;
        while(ali_it0!=ali_it0_stop && ali_it1!=ali_it1_stop && !isGap(ali_it0)&& !isGap(ali_it1))
        {
            ++readPos;
            ++refPos;
            ++ali_it0;
            ++ali_it1;
        }
        while(ali_it0!=ali_it0_stop && isGap(ali_it0))
        {
            ++refPos;
            ++ali_it0;
            ++ali_it1;
            ++deleted;
        }
        if(deleted>0)
        {
            appendValue(deletions,Pair<int,Pair<int,int> >(refPos-deleted,Pair<int,int>(readPos,deleted)));
//          appendValue(deletions,Pair<int,Pair<int,int> >(refPos,Pair<int,int>(readPos,deleted)));
        }
        while(isGap(ali_it1)&& ali_it1!=ali_it1_stop)
        {
            ++ali_it0;
            ++ali_it1;
            ++readPos;
            ++inserted;
        }
        if(inserted>0)
        {
            appendValue(insertions,Pair<int,Pair<int,int> >(refPos,Pair<int,int>(readPos-inserted,inserted)));
        }
    }

}


// TODO: get rid of this function, integrate indel calling into snp calling procedure as for realigned reads
template <
    typename TFragmentStore,
    typename TReadCigars,
    typename TGenome,
    typename TGenomeName,
    typename TFile,
    typename TOptions
>
void dumpShortIndelPolymorphismsBatch(
    TFragmentStore              &fragmentStore,             // forward/reverse matches
    TReadCigars             &readCigars,
    TGenome                 &genome,                // genome sequence
    TGenomeName const           genomeID,               // genome name
    typename TFragmentStore::TContigPos startCoord,         // startCoordinate + posOnGenomeInfix = real coordinate on whole chromosome
    typename TFragmentStore::TContigPos currStart,
    typename TFragmentStore::TContigPos currEnd,
    TFile                   &indelfile,
    TOptions                &options)
{

    typedef typename TFragmentStore::TAlignedReadStore      TMatches;
    typedef typename TFragmentStore::TAlignQualityStore         TMatchQualities;
    typedef typename Value<TMatches>::Type              TMatch;
    typedef typename TFragmentStore::TReadSeqStore          TReads;
    typedef typename Value<TReads>::Type                TRead;
    typedef typename Infix<TRead>::Type             TReadInf;
    typedef typename Iterator<TMatches, Standard>::Type TMatchIterator;
    typedef typename TFragmentStore::TContigPos TContigPos;
    // matches need to be ordered accordign to genome position
    TReads &reads = fragmentStore.readSeqStore;
    TMatches &matches = fragmentStore.alignedReadStore;
    TMatchQualities &matchQualities = fragmentStore.alignQualityStore;
    ::std::sort(begin(matches, Standard()), end(matches, Standard()), LessGPos<TMatch>());

    if (!indelfile.is_open())
    {
        ::std::cerr << "Failed to open indel output file" << ::std::endl;
        return;
    }

    TMatchIterator matchIt = begin(matches, Standard());
    TMatchIterator matchItEnd = end(matches, Standard());

    Align<Dna5String, ArrayGaps> align;
//  Score<int> scoreType = Score<int>(1, -3, -11, -1);  //
    Score<int> scoreType = Score<int>(0, -999, -1001, -1000);   // levenshtein-score (match, mismatch, gapOpen, gapExtend)
    resize(rows(align), 2);

    ::std::string runID = options.runID;

    typedef Pair<unsigned,int> TPosLen;
    typedef typename std::map<TPosLen,Pair<unsigned,TReadInf>, LessPosLen<TPosLen> > TIndelMap;
    typedef typename TIndelMap::iterator TIndelIt;
    typedef typename std::map<TPosLen,unsigned,LessPosLen<TPosLen> >  TSplitMap;
    typedef typename TSplitMap::iterator TSplitIt;
    typedef typename std::map<TPosLen,Pair<bool,bool>, LessPosLen<TPosLen> > THelperStrandMap;
    typedef typename THelperStrandMap::iterator TStrandIt; // for each indel candidate, store whether both strinds were seen


    // position,length and count,sequence
    TIndelMap indels; //readinf empty for deletions
    THelperStrandMap indelStrandHelper;


    //remember how many split reads supported the indel
    TSplitMap splitCounts;

    TIndelIt indelIt;
    TStrandIt strandIt;
    TSplitIt splitCountIt;
    //CharString chrPrefix = "";
    CharString chrPrefix = "chr"; // should check if "chr" is already part of chromosome names (usually not)


    TReadInf dummyInf;

#ifdef SNPSTORE_DEBUG
    bool extraV = true;
#endif

    // collect potential indels
    for(;matchIt != matchItEnd; ++matchIt)
    {
        if(matchQualities[(*matchIt).id].pairScore == 0 && empty(readCigars[(*matchIt).readId])) //if(length(cigar)>0)  dont skip!
            continue;

        String<Pair<int,Pair<int,int> > > readInserts;
        String<Pair<int,Pair<int,int> > > readDeletes;

        TRead& read = reads[(*matchIt).readId];
        int readLen = length(read);
        if(empty(readCigars[(*matchIt).readId]))// if this is not a split read --> do edit alignment
        {
#ifdef SNPSTORE_DEBUG
            if(extraV) ::std::cout << "read is edit indel mapped" << std::endl;
            if(extraV) ::std::cout << "read=" << read << " beg,end="<<(*matchIt).beginPos << ","<<(*matchIt).endPos <<::std::endl;
#endif
            assignSource(row(align, 0), reads[(*matchIt).readId]);
            assignSource(row(align, 1), infix(genome, _min((*matchIt).beginPos,(*matchIt).endPos), _max((*matchIt).beginPos,(*matchIt).endPos)));
            if ((*matchIt).beginPos > (*matchIt).endPos)
                reverseComplement(source(row(align, 0))); // check if reversing read is better for gap placement

            globalAlignment(align, scoreType, AlignConfig<false,true,true,false>(), Gotoh());
//          globalAlignment(align, scoreType, AlignConfig<false,false,false,false>(), Gotoh());
#ifdef SNPSTORE_DEBUG
            if(extraV) ::std::cout << align << std::endl;
#endif
            // transform first and last read character to genomic positions
            unsigned viewPosReadFirst  = toViewPosition(row(align, 0), 0);
            unsigned viewPosReadLast   = toViewPosition(row(align, 0), length(reads[(*matchIt).readId]) - 1);

            getIndels(align,readInserts,readDeletes, viewPosReadFirst,viewPosReadLast+1);

#ifdef SNPSTORE_DEBUG
            for (unsigned i = 0; i < length(readInserts); ++i)
                if(extraV) ::std::cout <<"ins: "<< readInserts[i].i1 << ","<<readInserts[i].i2.i1 <<","<< readInserts[i].i2.i2 << ::std::endl;
            for (unsigned i = 0; i < length(readDeletes); ++i)
                if(extraV) ::std::cout <<"del: "<<  readDeletes[i].i1 << ","<<readDeletes[i].i2.i1 <<","<< readDeletes[i].i2.i2 << ::std::endl;
#endif
        }
        else
        {

    //      if(extraV) ::std::cout << "read is split mapped" << std::endl;
            //this is where i have to get rid of adjacent insertions/deletions in edit-split-mapped reads
            typename Value<TReadCigars>::Type &cigar = readCigars[(*matchIt).readId];
            int readPos = 0;
            int refPos = 0;
            if((*matchIt).endPos > (*matchIt).beginPos)
            {
                for(unsigned i = 0; i < length(cigar); ++i)
                {
                    if(cigar[i].i1 == 'D') //deletion
                    {
                        appendValue(readDeletes,Pair<int,Pair<int,int> >(refPos,Pair<int,int>(readPos,cigar[i].i2)));
                        //::std::cout << " "<<cigar[i].i2 << " d at refPos " << refPos ;
                        refPos += cigar[i].i2;
                    }
                    if(cigar[i].i1 == 'I') //deletion
                    {
                        appendValue(readInserts,Pair<int,Pair<int,int> >(refPos,Pair<int,int>(readPos,cigar[i].i2)));
                        readPos += cigar[i].i2;
                        //::std::cout << " "<<cigar[i].i2<< " i at refPos " << refPos ;
                    }
                    if(cigar[i].i1 == 'M') //matches
                    {
                        refPos += cigar[i].i2;
                        readPos += cigar[i].i2;
                        //::std::cout <<"  "<< cigar[i].i2<< " m " ;
                    }
                }
//              ::std::cout << std::endl;::std::cout << std::endl;
            }
            else{
                for(int i = length(cigar)-1; i >= 0; --i)
                {
                    if(cigar[i].i1 == 'D') //deletion
                    {
                        appendValue(readDeletes,Pair<int,Pair<int,int> >(refPos,Pair<int,int>(readPos,cigar[i].i2)));
//                      ::std::cout << " "<<cigar[i].i2 << " d at refPos " << refPos ;
                        refPos += cigar[i].i2;
                    }
                    if(cigar[i].i1 == 'I') //deletion
                    {
                        appendValue(readInserts,Pair<int,Pair<int,int> >(refPos,Pair<int,int>(readPos,cigar[i].i2)));
//                      ::std::cout << " "<<cigar[i].i2<< " i at refPos " << refPos << " readPos" << readPos << std::endl;
                        readPos += cigar[i].i2;
                    }
                    if(cigar[i].i1 == 'M') //matches
                    {
                        refPos += cigar[i].i2;
                        readPos += cigar[i].i2;
//                      ::std::cout <<"  "<< cigar[i].i2<< " m " ;
                    }
                }
//              ::std::cout << std::endl;::std::cout << std::endl;
            }
            //cigar[read][1].i1 => I/D? cigar[read][1].i2 = indelLen, cigar[read][0].i2 = readPos
        }

        //::std::cout << align;
        //unsigned refStart = min((*matchIt).beginPos,(*matchIt).endPos);
        for (unsigned i = 0; i < length(readInserts); ++i)
        {
            // go to genomic indel position
            unsigned indelCandPos;//
            if((*matchIt).beginPos > (*matchIt).endPos)
                indelCandPos = (*matchIt).endPos + readInserts[i].i1;

            else indelCandPos = (*matchIt).beginPos + readInserts[i].i1;
#ifdef SNPSTORE_DEBUG
            if(extraV)  //62
                std::cout << "Pos=" << indelCandPos  + startCoord << " len=" <<  (readInserts[i].i2).i2 << std::endl;
#endif

            //TODO: make use of i2
            //TODO: remember strand of indel-supporting read
            indelIt = indels.find(TPosLen((unsigned)indelCandPos,-(int)(readInserts[i].i2).i2));
            if(indelIt == indels.end())
            {
                // this is the first insertion with this length found at this genomic position
                // --> remember with count 1 and also store indelsize and readInf
    //          if(extraV) ::std::cout << "new inds pos" << std::endl;

                TReadInf rInf;
                if((*matchIt).beginPos < (*matchIt).endPos)
                {
                    rInf = infix(read,
                    (readInserts[i].i2).i1,
                    (readInserts[i].i2).i1+(readInserts[i].i2).i2);
                    indelStrandHelper.insert(std::make_pair(Pair<unsigned,int>(indelCandPos,-(int)(readInserts[i].i2).i2),
                         Pair<bool,bool>(true,false)));
                }
                else
                {
                    rInf = infix(read,
                    readLen - (readInserts[i].i2).i1-(readInserts[i].i2).i2,
                    readLen - (readInserts[i].i2).i1);
                    reverseComplement(rInf);
                    indelStrandHelper.insert(std::make_pair(Pair<unsigned,int>(indelCandPos,-(int)(readInserts[i].i2).i2),
                         Pair<bool,bool>(false,true)));
                }

                indels.insert(std::make_pair(Pair<unsigned,int>(indelCandPos,-(int)(readInserts[i].i2).i2),
                                             Pair<unsigned,TReadInf>(1,rInf)));

    //          if(extraV)std::cout << rInf << " <-" << (*matchIt).id<<std::endl;
            }
            else
            {
    //          if(extraV) ::std::cout << "increase counter ins pos" << std::endl;
                ++(indelIt->second.i1);
                strandIt = indelStrandHelper.find(TPosLen((unsigned)indelCandPos,-(int)(readInserts[i].i2).i2));
                //                SEQAN_ASSERT_NEQ(strandIt, indelStrandHelper.end());

                if((*matchIt).beginPos < (*matchIt).endPos)
                    strandIt->second.i1 = true;
                else
                    strandIt->second.i2 = true;



/*              if (storeall insertion sequences)
                {

                    TReadInf rInf;
                    if((*matchIt).beginPos < (*matchIt).endPos)
                        rInf = infix(read,
                        (readInserts[i].i2).i1,
                        (readInserts[i].i2).i1+(readInserts[i].i2).i2);
                    else
                        rInf = infix(read,
                        readLen - (readInserts[i].i2).i1-(readInserts[i].i2).i2,
                        readLen - (readInserts[i].i2).i1);
                    if((*matchIt).beginPos > (*matchIt).endPos)
                        reverseComplement(rInf);
                    if(extraV)std::cout << rInf << " <-" << (*matchIt).id<<std::endl;
                    indelIt->second.i2 = rInf;
                }*/

            }
            if(!empty(readCigars[(*matchIt).readId]))// if this is a split read --> increase counter
            {
                splitCountIt = splitCounts.find(Pair<unsigned,int>(indelCandPos,-(int)(readInserts[i].i2).i2));
                if(splitCountIt == splitCounts.end())
                {
                    splitCounts.insert(std::make_pair(Pair<unsigned,int>(indelCandPos,-(int)(readInserts[i].i2).i2),1));

                }
                else
                    ++(splitCountIt->second);
            }

        }
        for (unsigned i = 0; i < length(readDeletes); ++i)
        {
            unsigned indelCandPos;// = refStart + readDeletes[i].i1;
            if((*matchIt).beginPos > (*matchIt).endPos)
            //  indelCandPos = (*matchIt).beginPos - readDeletes[i].i1 - (readDeletes[i].i2).i2;
                indelCandPos = (*matchIt).endPos + readDeletes[i].i1;
            else indelCandPos = (*matchIt).beginPos + readDeletes[i].i1;
#ifdef SNPSTORE_DEBUG
            if(extraV)
                std::cout << "Pos=" << indelCandPos  + startCoord << " len=" <<  (readDeletes[i].i2).i2;
#endif

            //TODO: make use of i2
            indelIt = indels.find(Pair<unsigned,int>(indelCandPos,(int)(readDeletes[i].i2).i2));
            if(indelIt == indels.end())
            {
                indels.insert(std::make_pair(Pair<unsigned,int>(indelCandPos,(int)(readDeletes[i].i2).i2),
                     Pair<unsigned,TReadInf>(1,dummyInf)));

                if((*matchIt).beginPos < (*matchIt).endPos)
                {
                    indelStrandHelper.insert(std::make_pair(Pair<unsigned,int>(indelCandPos,(int)(readDeletes[i].i2).i2),
                         Pair<bool,bool>(true,false)));
                }
                else
                {
                    indelStrandHelper.insert(std::make_pair(Pair<unsigned,int>(indelCandPos,(int)(readDeletes[i].i2).i2),
                         Pair<bool,bool>(false,true)));
                }


            }
            else
            {
                ++(indelIt->second.i1);
                strandIt = indelStrandHelper.find(Pair<unsigned,int>(indelCandPos,(int)(readDeletes[i].i2).i2));
                //                SEQAN_ASSERT_NEQ(strandIt,indelStrandHelper.end());

                if((*matchIt).beginPos < (*matchIt).endPos)
                   strandIt->second.i1 = true;
                else
                    strandIt->second.i2 = true;

            }
            if(!empty(readCigars[(*matchIt).readId]))// if this is a split read --> increase counter
            {
                splitCountIt = splitCounts.find(Pair<unsigned,int>(indelCandPos,(readDeletes[i].i2).i2));
                if(splitCountIt == splitCounts.end())
                {
                    splitCounts.insert(std::make_pair<Pair<unsigned,int>,unsigned>
                    (Pair<unsigned,int>(indelCandPos,(int)(readDeletes[i].i2).i2),1));

                }
                else
                    ++(splitCountIt->second);
            }

        }
    }

    // now output all indels that meet the filter criteria
    matchIt = begin(matches, Standard());
    while(matchIt != matchItEnd)
    {
        unsigned currSeqNo = (*matchIt).contigId;
        TMatchIterator currSeqMatchItBegin = matchIt;
        while(matchIt != matchItEnd)
        {
            if ((*matchIt).contigId != currSeqNo) break;
            ++matchIt;
        }
        TMatchIterator currSeqMatchItEnd = matchIt;

        matchIt = currSeqMatchItBegin;
        indelIt = indels.begin();
        strandIt = indelStrandHelper.begin();
        splitCountIt = splitCounts.begin();
        TIndelIt endIt = indels.end();
        TSplitIt splitEndIt = splitCounts.end();

        //indel-merging, possibly suboptimal
        if(options.indelWindow > 0)
        {
            while(indelIt != endIt)
            {
                unsigned currPos = indelIt->first.i1;
                unsigned oriCurrPos = currPos;
                if(indelIt->second.i1 == 0) {++indelIt;continue;}
                TIndelIt nextIt = indelIt;
                ++nextIt;
                //for all positions that are
                while(nextIt != endIt && currPos + options.indelWindow > nextIt->first.i1 && oriCurrPos + 2*options.indelWindow > nextIt->first.i1 )
                {
                    //add the number of found indel-reads to it
                    if(indelIt->first.i2 == nextIt->first.i2)
                    {
                        if(indelIt->second.i1 < nextIt->second.i1 ) //nextIT has a higher count for that position --> add counts of indelIT to nextIT
                        {
                            nextIt->second.i1 += indelIt->second.i1;
                            indelIt->second.i1 = 0;
                            currPos = indelIt->first.i1;
                        }
                        else
                        {
                            indelIt->second.i1 += nextIt->second.i1;
                            nextIt->second.i1 = 0;
                        }
                    }
                    ++nextIt;
                }
                ++indelIt;
            }
            indelIt = indels.begin();
        }


/*      for(splitCountIt = splitCounts.begin(); splitCountIt != splitEndIt; ++splitCountIt)
        {
            std::cout << splitCountIt->first.i1 << ","  << splitCountIt->first.i2 << "," << splitCountIt->second << std::endl;
        }
        splitCountIt = splitCounts.begin();*/
        for(; indelIt != endIt; ++indelIt, ++strandIt)
        {
            bool debug = false;
#ifdef SNPSTORE_DEBUG
            debug=true;
#endif
            int splitSupport = 0;
            if(splitCountIt != splitEndIt &&
                (splitCountIt->first.i1 == indelIt->first.i1) && (splitCountIt->first.i2 == indelIt->first.i2))
            {
                splitSupport = splitCountIt->second;
                ++splitCountIt;
            }

            if(indelIt->second.i1 < options.indelCountThreshold)
            {
                if(debug)::std::cout << "indel: count too low "<<indelIt->second.i1<<"\n";
                continue;
            }
            if((TContigPos)indelIt->first.i1 + startCoord < currStart || (TContigPos)indelIt->first.i1 + startCoord >= currEnd)
            {
                if(debug)::std::cout << "indel: pos outside range "<<indelIt->first.i1<<"\n";
                continue;
            }
            bool bsi = false;
            if(strandIt->second.i1 == true && strandIt->second.i2 == true)
                bsi = true;
            if(options.bothIndelStrands && !bsi)
            {
                if(debug)::std::cout << "indel: not supported by both strands \n";
                continue;
            }
            SEQAN_ASSERT_EQ(strandIt->first.i1,indelIt->first.i1);
            SEQAN_ASSERT_EQ(strandIt->first.i2,indelIt->first.i2);

            unsigned candidatePos = indelIt->first.i1;
            while(matchIt != currSeqMatchItEnd && _max((*matchIt).endPos,(*matchIt).beginPos) <= (TContigPos) candidatePos)
                ++matchIt;

            TMatchIterator matchRangeBegin = matchIt;
            while(matchIt != currSeqMatchItEnd && _min((*matchIt).endPos,(*matchIt).beginPos) <= (TContigPos) candidatePos)
                ++matchIt;
            TMatchIterator matchRangeEnd = matchIt;

            int minOverlapDepth = 0;
            if(options.indelDepthMinOverlap != 0)
            {
                for(TMatchIterator tmpIt = matchRangeBegin; tmpIt != matchRangeEnd; ++tmpIt)
                {
                    if(((TContigPos)candidatePos - options.indelDepthMinOverlap >= _min((*tmpIt).endPos,(*tmpIt).beginPos))
                        && ((TContigPos)candidatePos + options.indelDepthMinOverlap < _max((*tmpIt).endPos,(*tmpIt).beginPos)))
                        ++minOverlapDepth;
                }
            }

            int coverage = matchRangeEnd-matchRangeBegin;
            if(coverage<(int)options.minCoverage)
            {
                matchIt = matchRangeBegin;
                continue;

            }
            matchIt = matchRangeBegin;

            Dna5 refBase = genome[candidatePos];
            if(refBase=='N') continue;

            unsigned covF = 0;
            unsigned covR = 0;

            while(matchIt != matchRangeEnd)
            {

                if(!(_min((*matchIt).beginPos,(*matchIt).endPos) <= (TContigPos)candidatePos
                    && (TContigPos)candidatePos < _max((*matchIt).beginPos,(*matchIt).endPos)))
                {
                    ++matchIt;
//                  if(options._debugLevel > 1) std::cout << "How can this happen?\n";
                    continue;
                }

                if((*matchIt).beginPos < (*matchIt).endPos) ++covF;
                else ++covR;
                ++matchIt;
            }
            int depth = covF + covR;
            if(options.indelDepthMinOverlap != 0)
                depth = minOverlapDepth;
            if(depth < (int)options.minCoverage)
            {
                if(options._debugLevel > 1)
                    ::std::cout << "Coverage " << covF+covR << " after applying max pile filter and discarding Ns" << ::std::endl;
                matchIt = matchRangeBegin;
                continue;
            }

            if((float)indelIt->second.i1/depth < options.indelPercentageT)
            {
                matchIt = matchRangeBegin;
                continue;
            }

            int indelSize=indelIt->first.i2;


            if(options.outputFormat < 2) //
            {
                int homoLength = checkSequenceContext(genome,candidatePos,indelSize);
                if(indelSize > 0 ) indelfile << chrPrefix << genomeID << '\t' << runID << "\tdeletion\t";
                else indelfile << chrPrefix <<genomeID << '\t' << runID << "\tinsertion\t";
                if(indelSize > 0 ) indelfile << candidatePos + startCoord + options.positionFormat << '\t';
                else indelfile << candidatePos + startCoord + options.positionFormat - 1 << '\t';
                if(indelSize > 0 ) indelfile << candidatePos + startCoord + options.positionFormat  + indelSize - 1;
                else indelfile << candidatePos + startCoord;// + options.positionFormat; //VORSICHT!!!
                indelfile << "\t" << (float)indelIt->second.i1/depth;
                indelfile << "\t+\t.\tID=" << candidatePos + startCoord + options.positionFormat;
                indelfile << ";size=" << indelSize;
                indelfile << ";count=" << indelIt->second.i1;
                if(indelSize < 0)indelfile << ";seq="<<indelIt->second.i2;
                indelfile << ";ebiDepth=" << depth << ";depth=" << covF+covR;
                if(splitSupport>0) indelfile << ";splitSupport=" << splitSupport;
                if(homoLength > 1) indelfile << ";homorun=" << homoLength;
                if(bsi) indelfile << ";bsi";
                if( (float)indelIt->second.i1/depth <= options.indelHetMax) indelfile << ";geno=het";
                else indelfile << ";geno=hom";
                indelfile << std::endl;
            }
            else
            {

                //chromosome
                indelfile << genomeID << '\t';
                indelfile <<  candidatePos + startCoord + options.positionFormat << '\t';
                if(options.orientationAware)
                {
                    indelfile << covF  <<'\t';
                    indelfile << covR  <<'\t';
                }
                else
                {
                    indelfile << covF+covR  <<'\t';
                }
                indelfile << indelIt->second.i1 << std::endl;
            }

            matchIt = matchRangeBegin;
        }

        matchIt = currSeqMatchItEnd;
        if(options._debugLevel>1) std::cout <<"Finished scanning window for deletions.\n"<<std::flush;
    }


    return;

}





//do indel calling based on pronounced drops in coverage
template <
    typename TFragmentStore,
    typename TGenomeName,
    typename TFile,
    typename TOptions
>
void dumpCopyNumberPolymorphismsBatch(
    TFragmentStore              &fragmentStore,             // forward/reverse matches
    TGenomeName const           genomeID,               // genome name
    typename TFragmentStore::TContigPos startCoord,
    typename TFragmentStore::TContigPos currStart,
    typename TFragmentStore::TContigPos currEnd,
    TFile           &file,
    TOptions        &options)
{
    typedef typename TFragmentStore::TAlignedReadStore      TMatches;
    typedef typename Value<TMatches>::Type              TMatch;
    //typedef typename TFragmentStore::TReadSeqStore          TReads;
    //typedef typename Value<TReads>::Type                TRead;
    typedef typename Iterator<TMatches, Standard>::Type TMatchIterator;

    // matches need to be ordered accordign to genome position
    TMatches &matches = fragmentStore.alignedReadStore;
    ::std::sort(begin(matches, Standard()), end(matches, Standard()), LessGPos<TMatch>());

    if (!file.is_open())
    {
        ::std::cerr << "Failed to open cnv output file" << ::std::endl;
        return;
    }

    typedef typename Iterator<TMatches, Standard>::Type TMatchIterator;

    TMatchIterator matchIt = begin(matches, Standard());
    TMatchIterator matchItEnd = end(matches, Standard());

    ::std::string runID = options.runID;


    //VORSICHT!! windowsize must be a multiple of cnvwindowsize !!!!!
    // bin matches according to their start! positions
    for(unsigned currBinStart = currStart;  currBinStart < currEnd; currBinStart += options.cnvWindowSize)
    {

        while(matchIt != matchItEnd && _min((*matchIt).beginPos,(*matchIt).endPos) < currBinStart ) // havent reached bin begin yet
            ++matchIt;

        unsigned count = 0;
        while(matchIt != matchItEnd && _min((*matchIt).beginPos,(*matchIt).endPos) < currBinStart + options.cnvWindowSize ) //count reads
        {
            ++count;
            ++matchIt;
        }
        CharString guess = "normal";
        if (count > options.expectedReadsPerBin + 3 *options.expectedReadsSD)
            guess = "insertion";
        else if (count < options.expectedReadsPerBin - 3 *options.expectedReadsSD)
            guess = "deletion";
        //if(guess != "normal" )
        //{
        file << genomeID << "\tcoverage\t"<< guess << "\t";
        file << currBinStart+ startCoord+1  << "\t" << currBinStart + startCoord + options.cnvWindowSize << "\t";
        file << count << "\t+\t.\t.\n";
        //}
    }


    return;

}

///////////////////////////////////////////////////////////////////////////////////////
// Output SNPs, do realignment if a certain number of indels is observed in the reads
template <
    typename TFragmentStore,
    typename TPosIterator,
    typename TReadCigars,
    typename TReadCounts,
    typename TGenomeName,
    typename TFile,
    typename TOptions
>
void dumpPositionsRealignBatchWrap(
    TFragmentStore              &fragmentStore,         // forward/reverse matches
    TPosIterator & inspectPosIt,
    TPosIterator & inspectPosItEnd,
    TReadCigars             &readCigars,
    TReadCounts const           &readCounts,
    TGenomeName const           genomeID,           // genome name
    typename TFragmentStore::TContigPos startCoord,         // startCoordinate + posOnGenomeInfix = real coordinate on whole chromosome
    typename TFragmentStore::TContigPos currWindowBegin,
    typename TFragmentStore::TContigPos currWindowEnd,
    TFile                   &posFile,
    TOptions                &options)
{

    typedef typename TFragmentStore::TAlignedReadStore  TMatches;
    typedef typename Value<TMatches>::Type          TMatch;
    typedef typename TFragmentStore::TAlignQualityStore     TMatchQualities;
    //typedef typename Value<TMatchQualities>::Type       TMatchQuality;
    //typedef typename TFragmentStore::TReadSeqStore      TReads;
    //typedef typename Value<TReads>::Type            TRead;
    //typedef typename TFragmentStore::TContigStore       TContigStore;
    typedef typename TFragmentStore::TContigPos         TContigPos;
    //typedef typename Value<TContigStore>::Type      TContig;
    //typedef typename TFragmentStore::TContigSeq         TContigSeq;
    typedef typename Iterator<TMatches,Standard>::Type  TMatchIterator;


    TMatches &matches       = fragmentStore.alignedReadStore;
    TMatchQualities &matchQualities = fragmentStore.alignQualityStore;

    ::std::sort(begin(matches, Standard()), end(matches, Standard()), LessGPos<TMatch>());

    TMatchIterator matchIt = begin(matches,Standard());
    TMatchIterator matchItEnd = end(matches,Standard());

    unsigned minNumIndels = options.indelCountThreshold;

//  std::fstream tmpfile;
//  tmpfile.open("Z:/seqan071010/projects/library/apps/chr4.beforeTotal.sam", ::std::ios_base::out);
//  write(tmpfile, fragmentStore, Sam());
//  tmpfile.close();

#ifdef SNPSTORE_DEBUG
    CharString strstr = "test";
    _dumpMatches(fragmentStore,strstr);
#endif

    // now find connected subsets, i.e. groups of reads that overlap
    // dont realign regions unworthy of realignment (no indel reads)
    while(matchIt != matchItEnd)
    {
        TMatchIterator matchItBatchBegin = matchIt;
        TContigPos groupEndPos = _max((*matchIt).endPos,(*matchIt).beginPos);
        TContigPos groupStartPos = _min((*matchIt).endPos,(*matchIt).beginPos);

        TContigPos groupStartCoordLocal = _max(0,(int)groupStartPos-options.realignAddBorder);
        while(inspectPosIt != inspectPosItEnd && *inspectPosIt < groupStartPos + startCoord)
        {
            if(options.orientationAware)
                posFile << genomeID << '\t' << *inspectPosIt + options.positionFormat<< "\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0" << std::endl;
            else
                posFile << genomeID << '\t' << *inspectPosIt  + options.positionFormat<< "\t0\t0\t0\t0\t0\t0" << std::endl;
             ++inspectPosIt;
        }

        int indelReadCount = 0; // how many reads have indels in the current group
        while(matchIt != matchItEnd && _min((*matchIt).beginPos,(*matchIt).endPos) < groupEndPos)
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

        if(numMatches > 0 && inspectPosIt != inspectPosItEnd &&
            groupStartPos + startCoord <= *inspectPosIt && *inspectPosIt < groupEndPos + startCoord)
        {
            //make temporary fragstore for group
            // shorten reference sequence to the current region (groupStartPos to groupEndPos)

            //FragmentStore<SnpStoreGroupSpec_> fragStoreGroup;
            //copyFragmentStore(fragStoreGroup,fragmentStore,matchItBatchBegin,matchItBatchEnd,groupStartPos,groupEndPos);

            TFragmentStore fragStoreGroup = fragmentStore;
            arrayMoveForward(matchItBatchBegin,matchItBatchEnd,begin(fragStoreGroup.alignedReadStore,Standard())); // reads wont be needed anymore
            resize(fragStoreGroup.alignedReadStore,numMatches,Exact());
            fragStoreGroup.contigStore[0].seq = infix(fragmentStore.contigStore[0].seq,groupStartCoordLocal,groupEndCoordLocal);

#ifdef SNPSTORE_DEBUG
            std::cout << "in realign wrap: groupEndPos = " <<  groupEndPos << " groupStartPos=" <<  groupStartPos << std::endl;
            std::cout << "genomeLength= " <<  length(fragmentStore.contigStore[0].seq) << std::endl;

            CharString strstre = "testgroup";
            _dumpMatches(fragStoreGroup,strstre);

#endif
            groupStartPos += startCoord;
            groupEndPos += startCoord;
            //groupStartCoord = groupStartPos;
            TContigPos groupStartCoord = startCoord + groupStartCoordLocal;
            groupStartPos = _max(groupStartPos,currWindowBegin);
            groupEndPos = _min(groupEndPos,currWindowEnd);

            //the current group is formed by all reads from matchItBatchBegin until matchItBatchEnd
            if(indelReadCount >= (int)minNumIndels && options.realign)
            {
                //do realignment
                dumpPositionsRealignBatch(fragStoreGroup,inspectPosIt,inspectPosItEnd,
                    readCigars,readCounts,genomeID,
                    groupStartCoord,groupStartPos,groupEndPos,
                    posFile,options);
            }
            else
            {
                // todo: switch between with or without realignment in dumpSNPsBatch.. make global in any case
                dumpPosBatch(fragStoreGroup,inspectPosIt,inspectPosItEnd,
                    readCigars,readCounts,genomeID,
                    groupStartCoord,groupStartPos,groupEndPos,
                    posFile,options);
            }
        }

    }


}


//////////////////////////////////////////////////////////////////////////////
// Output SNPs
template <
    typename TFragmentStore,
    typename TPosIterator,
    typename TReadCounts,
    typename TReadCigars,
    typename TGenomeName,
    typename TFile,
    typename TOptions
>
void dumpPositionsRealignBatch(
    TFragmentStore              &fragmentStore,             // forward/reverse matches
    TPosIterator                &inspectPosIt,
    TPosIterator                &inspectPosItEnd,
    TReadCigars                 &,
    TReadCounts const           &,
    TGenomeName const           genomeID,                   // genome name
    typename TFragmentStore::TContigPos startCoord,         // startCoordinate + posOnGenomeInfix = real coordinate on whole chromosome
    typename TFragmentStore::TContigPos currStart,
    typename TFragmentStore::TContigPos currEnd,
    TFile                   &posFile,
    TOptions                &options)
{

    typedef typename TFragmentStore::TAlignedReadStore  TMatches;
    typedef typename Value<TMatches>::Type              TMatch;
    typedef typename TFragmentStore::TAlignQualityStore TMatchQualities;
    //typedef typename Value<TMatchQualities>::Type       TMatchQuality;
    typedef typename Iterator<TMatches,Standard>::Type  TMatchIterator;

    typedef typename TFragmentStore::TReadSeqStore      TReads;
    typedef typename Value<TReads>::Type                TRead;

    //typedef typename TFragmentStore::TContigStore       TContigStore;
    //typedef typename Value<TContigStore>::Type          TContig;
    typedef typename TFragmentStore::TContigPos         TContigPos;
    //typedef typename TFragmentStore::TContigSeq         TContigSeq;

    //typedef Gaps<TContigSeq, AnchorGaps<typename TContig::TGapAnchors> >    TContigGaps;
    typedef Gaps<TRead, AnchorGaps<typename TMatch::TGapAnchors> >          TReadGaps;
    //typedef typename Iterator<TContigGaps>::Type                            TContigGapIter;
    typedef typename Iterator<TReadGaps>::Type                              TReadGapIter;


    SEQAN_PROTIMESTART(dump_time);

    // matches need to be ordered according to genome position
    TReads &reads                   = fragmentStore.readSeqStore;
    TMatches &matches               = fragmentStore.alignedReadStore;
    TMatchQualities &matchQualities = fragmentStore.alignQualityStore;
    TContigPos genomeLen                = (TContigPos)length(fragmentStore.contigStore[0].seq);

    ::std::sort(begin(matches, Standard()), end(matches, Standard()), LessGPos<TMatch>());

    // make sure both output files are open
    if (!posFile.is_open())
    {
        ::std::cerr << "position output file is not open" << ::std::endl;
        return;
    }


    if(options._debugLevel > 1) ::std::cout << "Scanning chromosome " << genomeID << " window (" << currStart<<","<< currEnd << ") for SNPs..." << ::std::endl;

#ifdef READS_454
    Score<int> scoreType = Score<int>(0, -3, -2, -2);   // (match, mismatch,gapExtend,gapOpen)
#else
    Score<int> scoreType = Score<int>(0, -1, -2, -10);  // (match, mismatch,gapExtend,gapOpen)
#endif

    convertMatchesToGlobalAlignment(fragmentStore, scoreType, Nothing());

    // TODO: check out if a different scoring scheme makes more sense
    Score<int, WeightedConsensusScore<Score<int, FractionalScore>, Score<int, ConsensusScore> > > consScore;
    int bandWidth = 10; // ad hoc, but seems good, increasing is not necessarily better

    if(options._debugLevel > 1)
        ::std::cout << "Realigning "<< length(matches)<<" reads to genome of length " <<genomeLen << std::flush;



#ifdef READS_454
    reAlign(fragmentStore,consScore,0,1,bandWidth,true);
#else
    reAlign(fragmentStore,consScore,0,1,bandWidth,true);
#endif


    if(options._debugLevel > 1)::std::cout << "Realigning reads including reference..." << std::flush;


    unsigned numReads = length(matches)-1; // exclude reference sequence
    unsigned refId = length(matchQualities); // reference id (there may be more matchQs than matches due to pile up correction)


//
#ifndef  READS_454
    reAlign(fragmentStore,consScore,0,1,/*bandWidth*/5,false);
    //realignReferenceToReadProfile(fragmentStore,refId,options);
    //realignReferenceToDiploidConsensusProfile(fragmentStore,refId,options);
#else
//  reAlign(fragmentStore,consScore,0,1,/*bandWidth*/5,false);
    realignReferenceToDiploidConsensusProfile(fragmentStore,refId,options);
//    realignReferenceToDiploidConsensusProfileDeleteSeqErrors(fragmentStore,refId,options);
#endif

    if(options._debugLevel > 1) ::std::cout << "Finished realigning." << std::endl;


    // forward match qualities
    String<int> columnQualityF;         resize(columnQualityF,5);
    String<unsigned> countF;            resize(countF,5);
    String<CharString> qualityStringF;  resize(qualityStringF,5);

    // reverse match qualities
    String<int> columnQualityR;         resize(columnQualityR,5);
    String<unsigned> countR;            resize(countR,5);
    String<CharString> qualityStringR;  resize(qualityStringR,5);

    FunctorComplement<Dna5> f;

    // sort reads according to begin position
    sortAlignedReads(fragmentStore.alignedReadStore, SortBeginPos());
    TMatchIterator matchIt      = begin(matches, Standard());
    TMatchIterator matchItEnd   = end(matches, Standard());

    // look for reference sequence and move it to the end of alignedreads
    // todo: only do this when realignment was done
    bool refFound = false;
    TMatchIterator matchItKeep = matchIt;
    TMatch tempRef;
    while(matchIt != matchItEnd)
    {
        if((*matchIt).readId == refId) // this is the reference
        {
            refFound = true;
            tempRef = *matchIt;
            matchItKeep = matchIt;
            ++matchIt;
            continue;
        }
        if(refFound)
        {
            *matchItKeep = *matchIt; // matchItKeep lags behind by one match
            ++matchIt;++matchItKeep;
        }
        else ++matchIt;
    }
    *matchItKeep = tempRef;
    SEQAN_ASSERT(refFound);

    matchIt     = begin(matches, Standard());
    matchItEnd  = end(matches, Standard());
    matchItEnd--; // exclude reference sequence

    TRead       &reference = fragmentStore.readSeqStore[fragmentStore.alignedReadStore[numReads].readId]; // last read is reference sequence
    TReadGaps   referenceGaps(reference, fragmentStore.alignedReadStore[numReads].gaps);
//  TContigPos      refStart = (TContigPos)fragmentStore.alignedReadStore[numReads].beginPos;
//  TContigGaps contigGaps(fragmentStore.contigStore[0].seq, fragmentStore.contigStore[0].gaps);

    if(options._debugLevel>1) std::cout << "Start inspecting alignment..." << std::endl;
    // now walk through the reference sequence in gaps view space,
    // i.e. position may be a gap
    // example:
    // Ref      ACCGTGCACTAGCATCATT--ACTAGCATCATA
    // Reads    ACCGTACA--AGCATCAT
    //              TACA--AGCATCATT--ACT
    //                          ATTTTACTAGCATCATA
    for( ;inspectPosIt != inspectPosItEnd && *inspectPosIt < currEnd; ++inspectPosIt)
    //for(TContigPos candidateViewPos = refStart; candidateViewPos < refStart + (TContigPos)length(referenceGaps); ++candidateViewPos)
    {

        //get position in sequence space
        TContigPos candidatePos = *inspectPosIt;
        // not in the current window yet
        if(candidatePos < currStart)
        {
            if(options.orientationAware)
                posFile << genomeID << '\t' << *inspectPosIt + options.positionFormat<< "\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0" << std::endl;
            else
                posFile << genomeID << '\t' << *inspectPosIt + options.positionFormat<< "\t0\t0\t0\t0\t0\t0" << std::endl;
            continue;
        }
        candidatePos -= startCoord;

        TContigPos candidateViewPos = positionSeqToGap(referenceGaps, candidatePos);
        // not in the current window anymore
        if(candidatePos + startCoord >= currEnd) break;

        //Dna5 refBase = reference[candidatePos];   // what happens if refGap==true, esp. for leading gaps?

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

        // start checking reads for this position, prepare some helpers
        Dna5 candidateBase;
        int quality;
//        std::set<int> readPosMap;

        for(unsigned t=0;t<5;++t)
        {
            countF[t] = 0;
            columnQualityF[t] = 0;
            clear(qualityStringF[t]);

            countR[t] = 0;
            columnQualityR[t] = 0;
            clear(qualityStringR[t]);
        }

        unsigned delPlus = 0;
        unsigned delMinus = 0;

        // now check reads
        while(matchIt != matchRangeEnd)
        {
            TContigPos currViewBegin = _min((*matchIt).beginPos,(*matchIt).endPos);
            TContigPos currViewEnd = _max((*matchIt).beginPos,(*matchIt).endPos);

            // make sure this match is really spanning the position
            if(!(currViewBegin <= candidateViewPos && candidateViewPos < currViewEnd))
            {
                ++matchIt;
                continue;
            }

            char orientation = ((*matchIt).beginPos > (*matchIt).endPos) ? 'R' : 'F';

            TReadGaps readGaps(reads[(*matchIt).readId],(*matchIt).gaps);
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
            std::cout << "ReadPos = " << readPos << std::endl;
#endif

            if(readPos != -1) //-1 indicates gap in read
            {
                //if(options.minDifferentReadPos > 0)
                //    if((unsigned)(length(reads[(*matchIt).readId]) - readPos) > options.excludeBorderPos  &&
                //            (unsigned) readPos >= options.excludeBorderPos )
                //        readPosMap.insert(readPos);

                if(orientation == 'R') candidateBase = f((Dna5)reads[(*matchIt).readId][readPos]);
                else candidateBase = (Dna5)reads[(*matchIt).readId][readPos];

                quality = getQualityValue(reads[(*matchIt).readId][readPos]) ;

                if(!options.useBaseQuality && quality > (int)matchQualities[(*matchIt).id].score)
                {   // dont trust the quality of this position more
                    // than the average quality of this read
                    quality = (int) matchQualities[(*matchIt).id].score;
                }

                if(orientation == 'F')
                {
                    columnQualityF[ordValue(candidateBase)] += quality;
                    ++countF[ordValue(candidateBase)];
                    appendValue(qualityStringF[ordValue(candidateBase)],(char)(quality+33),Generous());
                }
                else
                {
                    columnQualityR[ordValue(candidateBase)] += quality;
                    ++countR[ordValue(candidateBase)];
                    appendValue(qualityStringR[ordValue(candidateBase)],(char)(quality+33),Generous());
                }
            }
            else
            {
                if(orientation == 'R') ++delMinus;
                else ++delPlus;
            }
            ++matchIt;
        }
        matchIt = matchRangeBegin; //set iterator back to where we started from, same matches might be involved in next cand pos

        //all observed bases match the reference allele or there were too few indels
        //if(!observedAtLeastOneMut && numIndelsObserved< options.indelCountThreshold)
        //  continue;


        // write out coverage info
        _writePos(posFile,qualityStringF,qualityStringR,delPlus,delMinus,genomeID,candidatePos+startCoord,coverage,options);

    }

    if(options._debugLevel>1) std::cout <<"Finished scanning window.\n"<<std::flush;


}


//////////////////////////////////////////////////////////////////////////////
// Output SNPs
template <
    typename TFragmentStore,
    typename TPosIterator,
    typename TReadCigars,
    typename TReadCounts,
    typename TGenomeName,
    typename TFile,
    typename TOptions
>
void dumpPosBatch(
    TFragmentStore              &fragmentStore,             // forward/reverse matches
    TPosIterator                &inspectPosIt,
    TPosIterator                &inspectPosItEnd,
    TReadCigars             &,
    TReadCounts const           &readCounts,
    TGenomeName const           genomeID,                   // genome name
    typename TFragmentStore::TContigPos startCoord,         // startCoordinate + posOnGenomeInfix = real coordinate on whole chromosome
    typename TFragmentStore::TContigPos currStart,
    typename TFragmentStore::TContigPos currEnd,
    TFile               &posFile,
    TOptions            &options)
{

    typedef typename TFragmentStore::TAlignedReadStore  TMatches;
    typedef typename Value<TMatches>::Type              TMatch;
    typedef typename TFragmentStore::TAlignQualityStore TMatchQualities;
    //typedef typename Value<TMatchQualities>::Type       TMatchQuality;
    typedef typename TFragmentStore::TReadSeqStore      TReads;
    //typedef typename Value<TReads>::Type                TRead;
    typedef typename TFragmentStore::TContigPos         TContigPos;
    typedef typename TFragmentStore::TContigSeq         TContigSeq;
    typedef typename Iterator<TMatches,Standard>::Type  TMatchIterator;

    SEQAN_PROTIMESTART(dump_time);
    //options._debugLevel = 2;
    String<char> toIupac = "AMRWMCSYRSGKWYKT";
    //std::cout << "Hier\n";
    // matches need to be ordered accordign to genome position
    TReads &reads                   = fragmentStore.readSeqStore;
    TMatches &matches               = fragmentStore.alignedReadStore;
    TMatchQualities &matchQualities = fragmentStore.alignQualityStore;
    TContigSeq &genome              = fragmentStore.contigStore[0].seq;

    ::std::sort(begin(matches, Standard()), end(matches, Standard()), LessGPos<TMatch>());

    Align<String<Dna5>, ArrayGaps> align;
    Score<int> scoreType = Score<int>(0, -999, -1001, -1000);   // levenshtein-score (match, mismatch, gapOpen, gapExtend)
    resize(rows(align), 2);

    if (!posFile.is_open())
    {
        ::std::cerr << "Output file is not open" << ::std::endl;
        return;
    }

    TMatchIterator matchIt  = begin(matches, Standard());
    TMatchIterator matchItEnd   = end(matches, Standard());
    //matchItEnd--;
    unsigned countLowerMQ = 0, countHigherMQ = 0;

    if(options._debugLevel > 1) ::std::cout << "Scanning chromosome " << genomeID << " window (" << currStart<<","<< currEnd << ") for SNPs..." << ::std::endl;

    // forward match qualities
    String<int> columnQualityF;         resize(columnQualityF,5);
    String<unsigned> countF;            resize(countF,5);
    String<CharString> qualityStringF;  resize(qualityStringF,5);

    // reverse match qualities
    String<int> columnQualityR;         resize(columnQualityR,5);
    String<unsigned> countR;            resize(countR,5);
    String<CharString> qualityStringR;  resize(qualityStringR,5);

    // both
    String<unsigned> count;             resize(count,5);
    String<unsigned> columnQuality;     resize(columnQuality,5);

#ifdef SNPSTORE_DEBUG
    bool extraV = false;
#endif

    for(TContigPos candidatePos = *inspectPosIt; inspectPosIt !=  inspectPosItEnd && *inspectPosIt < currEnd; ++inspectPosIt)
    {

        if(candidatePos < currStart)
        {
            if(options.orientationAware)
                posFile << genomeID << '\t' << *inspectPosIt + options.positionFormat << "\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0" << std::endl;
            else
                posFile << genomeID << '\t' << *inspectPosIt + options.positionFormat << "\t0\t0\t0\t0\t0\t0" << std::endl;
            continue;
        }
        candidatePos -= startCoord;

        //Dna5 refBase = genome[candidatePos];

#ifdef SNPSTORE_DEBUG
        ::std::cout << "candPos=" << candidatePos + startCoord << ::std::endl;
        if(candidatePos + startCoord == 861196)
            ::std::cout << "ab jetzt.." << ::std::flush;
#endif

        Dna5 candidateBase;
        int quality;

//      if(options._debugLevel > 1)std::cout << candidatePos+startCoord << "<-candidatePos\n";
        for(unsigned t=0;t<5;++t)
        {
            countF[t] = 0;
            columnQualityF[t] = 0;
            clear(qualityStringF[t]);

            countR[t] = 0;
            columnQualityR[t] = 0;
            clear(qualityStringR[t]);
        }

        //find range of relevant read matches
        while(matchIt != matchItEnd &&  _max((*matchIt).endPos,(*matchIt).beginPos) <= candidatePos)
            ++matchIt;
        TMatchIterator matchRangeBegin = matchIt;
        while(matchIt != matchItEnd &&  _min((*matchIt).endPos,(*matchIt).beginPos)  <= candidatePos)
            ++matchIt;
        TMatchIterator matchRangeEnd = matchIt;
        matchIt = matchRangeBegin;

        int coverage = matchRangeEnd-matchRangeBegin;

        if(options._debugLevel > 1)::std::cout << "Match range:" << matchRangeEnd - matchRangeBegin << ::std::endl;
#ifdef SNPSTORE_DEBUG
        if(extraV)
        {
            for (TMatchIterator tempIt = matchRangeBegin; tempIt != matchRangeEnd; ++tempIt)
                ::std::cout << reads[(*tempIt).readId]<<"\n";
        }
#endif
//        std::set<unsigned> readPosMap;

        unsigned delPlus = 0;
        unsigned delMinus = 0;

        while(matchIt != matchRangeEnd)
        {
            TContigPos currentBegin = _min((*matchIt).beginPos,(*matchIt).endPos);
            TContigPos currentEnd   = _max((*matchIt).beginPos,(*matchIt).endPos);
            char orientation = ((*matchIt).beginPos > (*matchIt).endPos) ? 'R' : 'F';

#ifdef SNPSTORE_DEBUG
            if(extraV)
            {
                ::std::cout <<"currentBegin = "<<currentBegin << "\n";
                ::std::cout <<"currentEnd = "<<currentEnd << "\n";
            }
#endif
            if(!(currentBegin <= candidatePos && candidatePos < currentEnd))// this match is not really spanning the position
            {                                                               // (can happen because of indels or variable-length reads)
                ++matchIt;
                continue;
            }
            // do edit alignment
            if((int)matchQualities[(*matchIt).id].pairScore == 1 /*!empty((*matchIt).gaps)*/) // splitReads: hamming: pairScore=0
            {
                Dna5String gInf = infix(genome, currentBegin, currentEnd);
                if (orientation == 'R')
                    reverseComplement(gInf);

                assignSource(row(align, 0), reads[(*matchIt).readId]);
                assignSource(row(align, 1), gInf);
                globalAlignment(align, scoreType);  //splitReads: get alignment from cigar string
            }

            if (orientation == 'R')
            {
                FunctorComplement<Dna5> f;

                int readPos = currentEnd - candidatePos - 1;
                if ((int)matchQualities[(*matchIt).id].pairScore == 1 /*!empty((*matchIt).gaps)*/)
                    readPos = getReadPos(align,readPos,false); //

#ifdef SNPSTORE_DEBUG
                if(extraV) std::cout << "readPosNacher = " << readPos << std::endl;
#endif
                if(readPos != -1) //-1 indicates gap
                {
                    //if(options.minDifferentReadPos > 0)
                    //    if((unsigned)(length(reads[(*matchIt).readId]) - readPos) > options.excludeBorderPos  &&
                    //        (unsigned) readPos >= options.excludeBorderPos )
                    //        readPosMap.insert(readPos);
                    candidateBase = f((Dna5)reads[(*matchIt).readId][readPos]);
#ifdef SNPSTORE_DEBUG
                    if(extraV) std::cout << candidateBase << "candBase\n";
#endif
                    quality = getQualityValue(reads[(*matchIt).readId][readPos]) ;

                    if(!options.useBaseQuality && quality > (int)matchQualities[(*matchIt).id].score)   // dont trust the quality of this position more
                    {                                                                               // than the average quality of this read
                        quality = (int) matchQualities[(*matchIt).id].score;
                        ++countLowerMQ;
                    }
                    else ++countHigherMQ;
                    //if(quality < 0 || quality > 40)::std::cout << "falschQ candPos = " << candidatePos + startCoord << std::endl;

                    unsigned tmpCount = 1;
                    if(!empty(readCounts)) tmpCount = readCounts[(*matchIt).readId];
                    for (unsigned k = 0; k < tmpCount; ++k)
                    {
                        columnQualityR[ordValue(candidateBase)] += quality;
                        ++countR[ordValue(candidateBase)];
                        appendValue(qualityStringR[ordValue(candidateBase)],(char)(quality+33),Generous());
                    }
                }
                else
                    ++delMinus;
            }
            else
            {
                int readPos = candidatePos - currentBegin;

                if((int)matchQualities[(*matchIt).id].pairScore == 1 /*!empty((*matchIt).gaps)*/)
                    readPos = getReadPos(align,readPos);

                if(readPos != -1) //-1 indicates gap
                {
                    //if(options.minDifferentReadPos > 0)
                    //    if((unsigned)(length(reads[(*matchIt).readId]) - readPos) > options.excludeBorderPos  &&
                    //        (unsigned) readPos >= options.excludeBorderPos )
                    //        readPosMap.insert(readPos);

                    candidateBase = (Dna5)reads[(*matchIt).readId][readPos];
                    quality = getQualityValue(reads[(*matchIt).readId][readPos]) ;

                    if(!options.useBaseQuality && quality > (int) matchQualities[(*matchIt).id].score)
                    {
                        quality = (int) matchQualities[(*matchIt).id].score;
                        ++countLowerMQ;
                    }
                    else ++countHigherMQ;
                    //if(quality < 0 || quality > 40)::std::cout << "falschQ candPos = " << candidatePos + startCoord << std::endl;

                    unsigned tmpCount = 1;
                    if(!empty(readCounts)) tmpCount = readCounts[(*matchIt).readId];
                    for (unsigned k = 0; k < tmpCount; ++k)
                    {
                        ++countF[ordValue(candidateBase)];
                        columnQualityF[ordValue(candidateBase)] += quality;
                        appendValue(qualityStringF[ordValue(candidateBase)],(char)(quality+33),Generous());
                    }
                }
                else
                    ++delPlus;
            }
            ++matchIt;
        }
        matchIt = matchRangeBegin; //set iterator back to where we started from, same matches might be involved in next cand pos

        // write info to file
        _writePos(posFile,qualityStringF,qualityStringR,delPlus,delMinus,genomeID,candidatePos+startCoord,coverage,options);

    }

    if(options._debugLevel>1) std::cout <<"Finished scanning window.\n"<<std::flush;

    return;

}



}

#endif
