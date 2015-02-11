#ifndef __APPS_BS_TOOLS_CASBAR_REALIGNMENT_H__
#define __APPS_BS_TOOLS_CASBAR_REALIGNMENT_H__

#include "bisar_score_data.h"
#include "casbar_score_data.h"
#include "casbar_score.h"
#include "casbar_consensus_realign.h"


using namespace std;
using namespace seqan;


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

////////////////////////////////////////////////////////////////////////////////

template <
    typename TFragmentStore,
    typename TMethOptions,
    typename TOptions
>
void doRealigning(
    TFragmentStore              &fragmentStore,             // forward/reverse matches
    typename TFragmentStore::TContigPos currStart,
    typename TFragmentStore::TContigPos currEnd,
    TMethOptions            &methOptions,
    TOptions                &options)
{
    std::cout << " doRealigning() ... " << std::endl;
    typedef typename TFragmentStore::TAlignedReadStore  TMatches;
    typedef typename Value<TMatches>::Type              TMatch;
    typedef typename Iterator<TMatches,Standard>::Type  TMatchIterator;

    SEQAN_PROTIMESTART(dump_time);

    // matches need to be ordered according to genome position
    TMatches &matches = fragmentStore.alignedReadStore;

    std::sort(begin(matches, Standard()), end(matches, Standard()), LessGPos<TMatch>());


    if(options._debugLevel > 1) ::std::cout << "Scanning chromosome " << fragmentStore.contigNameStore[0] << " window (" << currStart<<","<< currEnd << ") for SNPs..." << ::std::endl;

    int bandWidth = 10; // ad hoc
    reAlign(fragmentStore, 0, bandWidth, true, methOptions, BsSimple());


    if(options._debugLevel > 1)::std::cout << "Realigning reads including reference..." << std::flush;

    // Would diploid consensus profile would be not that clear in bs case, hence we do not do that for the beginning!
    //unsigned refId = length(matchQualities); // reference id (there may be more matchQs than matches due to pile up correction)
    //realignReferenceToDiploidConsensusProfile(fragmentStore,refId,options);

    // sort reads according to begin position
    sortAlignedReads(fragmentStore.alignedReadStore, SortBeginPos());
    TMatchIterator matchIt      = begin(matches, Standard());
    TMatchIterator matchItEnd   = end(matches, Standard());

    // look for reference sequence and move it to the end of alignedreads
    bool refFound = false;
    TMatchIterator matchItKeep = matchIt;
    TMatch tempRef;
    while(matchIt != matchItEnd)
    {
        if ((*matchIt).readId == length(fragmentStore.readSeqStore)-1)
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
}



#endif
