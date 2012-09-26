#ifndef BENCHMARKS_READ_MAPPERS_APPS_WIT_STORE_H_
#define BENCHMARKS_READ_MAPPERS_APPS_WIT_STORE_H_

#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/sequence.h>
#include <seqan/store.h>  // For Sort* tags.

#include "intervals.h"

using namespace seqan;

// Stores closed intervals [firstPos, lastPos].
struct IntervalOfReadOnContig {
    static size_t invalidId() { return ~0ul; }
    static size_t additionalId() { return ~1ul; }
    static size_t superflousId() { return ~2ul; }
    
    size_t id;
    size_t readId;
    unsigned distance;
    size_t contigId;
    bool isForward;
    size_t firstPos;
    size_t lastPos;

    IntervalOfReadOnContig()
            : id(invalidId()), readId(0), distance(0), contigId(0), isForward(0), firstPos(0),
              lastPos(0) {}

    IntervalOfReadOnContig(size_t _readId, unsigned _distance, size_t _contigId, bool _isForward, size_t _firstPos, size_t _lastPos)
            : id(invalidId()), readId(_readId), distance(_distance), contigId(_contigId),
              isForward(_isForward), firstPos(_firstPos), lastPos(_lastPos) {}
};


template <typename TStream>
inline
TStream & operator<<(TStream & stream, IntervalOfReadOnContig const & record) {
    stream << "(id=" << record.id << ", readId=" << record.readId << ", contigId=" << record.contigId << ", distance="
           << record.distance << ", isForward=" << record.isForward << ", firstPos="
           << record.firstPos << ", lastPos=" << record.lastPos << ")";
    return stream;
}


struct WitStore {
    typedef StringSet<CharString, Owner<> > TNameSet;
    typedef String<IntervalOfReadOnContig> TIntervalStore;

    char mateSeparator;

    Holder<TNameSet> readNames;
    Holder<TNameSet> contigNames;

    // TODO(holtgrew): Rename to witRecords.
    TIntervalStore intervals;

    WitStore() {}

    WitStore(TNameSet & _readNames, TNameSet & _contigNames)
            : readNames(_readNames), contigNames(_contigNames) {}

};

template <typename TStream>
inline
TStream & operator<<(TStream & stream, WitStore const & store) {
    typedef typename Iterator<typename WitStore::TIntervalStore, Standard>::Type TIterator;
    stream << ",-- WIT Store" << std::endl;
    for (TIterator it = begin(store.intervals, Standard()); it != end(store.intervals, Standard()); ++it) {
        stream << "| " << value(store.readNames)[value(it).readId] << " read id = " << value(it).readId << "\t distance=" << value(it).distance << "\t" << value(store.contigNames)[value(it).contigId]
               << "\t" << (value(it).isForward ? "F" : "R") << "\t" << value(it).firstPos << "\t" << value(it).lastPos << std::endl;
    }
    stream << "`--" << std::endl;
    return stream;
}

void move(WitStore & target, WitStore & source) {
    target.mateSeparator = source.mateSeparator;
    move(target.readNames, source.readNames);
    move(target.contigNames, source.contigNames);
    move(target.intervals, source.intervals);
}


template <typename TSpec>
struct WitStoreLess;


template <>
struct WitStoreLess<SortReadId>
{
    WitStore const & store;

    WitStoreLess(WitStore const & _store) : store(_store) {};

    bool
    operator()(IntervalOfReadOnContig const & a, IntervalOfReadOnContig const & b)  const{
        return a.readId < b.readId;
    }
};


template <>
struct WitStoreLess<SortId>
{
    WitStore const & store;

    WitStoreLess(WitStore const & _store) : store(_store) {};

    bool
    operator()(IntervalOfReadOnContig const & a, IntervalOfReadOnContig const & b)  const{
        return a.id < b.id;
    }
};


template <>
struct WitStoreLess<SortContigId>
{
    WitStore const & store;

    WitStoreLess(WitStore const & _store) : store(_store) {};

    bool
    operator()(IntervalOfReadOnContig const & a, IntervalOfReadOnContig const & b)  const{
        return a.contigId < b.contigId;
    }
};


struct SortDistance_ {};
typedef Tag<SortDistance_> const SortDistance;

template <>
struct WitStoreLess<SortDistance>
{
    WitStore const & store;

    WitStoreLess(WitStore const & _store) : store(_store) {};

    bool
    operator()(IntervalOfReadOnContig const & a, IntervalOfReadOnContig const & b)  const{
        return a.distance < b.distance;
    }
};


struct SortFirstPos_ {};
typedef Tag<SortFirstPos_> const SortFirstPos;


template <>
struct WitStoreLess<SortFirstPos>
{
    WitStore const & store;

    WitStoreLess(WitStore const & _store) : store(_store) {};

    bool
    operator()(IntervalOfReadOnContig const & a, IntervalOfReadOnContig const & b)  const{
      if (a.firstPos < b.firstPos)
        return true;
      if (a.firstPos == b.firstPos && a.lastPos > b.lastPos)
        return true;
      if (a.firstPos == b.firstPos && a.lastPos == b.lastPos && a.distance > b.distance)
        return true;
      return false;
    }
};


struct SortLastPos_ {};
typedef Tag<SortLastPos_> const SortLastPos;


template <>
struct WitStoreLess<SortLastPos>
{
    WitStore const & store;

    WitStoreLess(WitStore const & _store) : store(_store) {};

    bool
    operator()(IntervalOfReadOnContig const & a, IntervalOfReadOnContig const & b)  const{
        return a.lastPos < b.lastPos;
    }
};


template <typename TSortTag>
void
sortWitRecords(WitStore const & store, TSortTag &) {
    std::stable_sort(begin(store.intervals, Standard()), end(store.intervals, Standard()), WitStoreLess<TSortTag>(store));
}


inline size_t
appendValue(WitStore & store, IntervalOfReadOnContig const & record) {
    IntervalOfReadOnContig tmp(record);
    tmp.id = length(store.intervals);
    appendValue(store.intervals, tmp);
    return tmp.id;
}


template <typename TFragmentStore>
void loadWitFile(WitStore & store,
                 TFragmentStore /*const*/ & fragments,
                 CharString const & fileName) {
    typedef typename Value<typename TFragmentStore::TMatePairStore>::Type TMatePairElement;

    // Assign read and contig names into wit store members.
    setValue(store.readNames, fragments.readNameStore);
    setValue(store.contigNames, fragments.contigNameStore);

    // Build mapping from existing read and contig names to their indices
    // (ids) in readNames and contigNames.
    NameStoreCache<StringSet<CharString> > readNameCache(value(store.readNames));
    refresh(readNameCache);
    NameStoreCache<StringSet<CharString> > contigNameCache(value(store.contigNames));
    refresh(contigNameCache);

    // Then, load the file.
    std::fstream fstrm(toCString(fileName), std::ios_base::in | std::ios_base::binary);
    if (!fstrm.is_open()) {
        std::cerr << "Could not open WIT file " << fileName << std::endl;
        exit(1);
    }

    // Read header.
    char c = '\0';
    // Separator char of read name and mate identifier, index of first mate id.
    char mateSeparator = '/';  // TODO(holtgrew): Un-hardcode.
    int mateStart = 0;  // TODO(holtgrew): Un-hardcode.
    readWitHeader(fstrm, c);

    // Temporary data for one record.
    CharString readName;
    int mateNo;
    size_t distance;
    CharString contigName;
    bool isForward;
    size_t firstPos;
    size_t lastPos;
    bool wasInSam = true;

    // Read WIT file.
    while (!_streamEOF(fstrm)) {
        // Skip comments.
        if (c == '#') {
            _parseSkipLine(fstrm, c);
            continue;
        }

        clear(readName);
        clear(contigName);
        mateNo = -1;

        // Read line.
        _parseReadIdentifier(fstrm, readName, c);
        //std::cout << "readName " << readName << ", " << mateNo << std::endl;
        if (readName[length(readName) - 2] == mateSeparator) {
          mateNo = readName[length(readName) - 1] - '0' - mateStart;
          resize(readName, length(readName) - 2);
        }
        _parseSkipWhitespace(fstrm, c);
        distance = _parseReadNumber(fstrm, c);
        _parseSkipWhitespace(fstrm, c);
        _parseReadIdentifier(fstrm, contigName, c);
        _parseSkipWhitespace(fstrm, c);
        isForward = (_parseReadChar(fstrm, c) == 'F');
        _parseSkipWhitespace(fstrm, c);
        firstPos = _parseReadNumber(fstrm, c);
        _parseSkipWhitespace(fstrm, c);
        lastPos = _parseReadNumber(fstrm, c);
        _parseSkipLine(fstrm, c);

        // Insert record into read store.
        //
        // We also need to insert it into the fragment store if it does not
        // exist there yet.
        IntervalOfReadOnContig record;
        if (!getIdByName(value(store.readNames), readName, record.readId, readNameCache)) {
          wasInSam = false;
          record.readId = length(value(store.readNames));
          appendName(value(store.readNames), readName, readNameCache);

          if (mateNo != -1) {
            // If read is paired, create new entry in read name store.
            TMatePairElement mateElem;
            // set the first or second read ID in the mate pair element
            size_t matePairId = length(fragments.matePairStore);
            mateElem.readId[mateNo] = record.readId;
            // get a new mate pair ID and add the new mate pair element
            appendValue(fragments.matePairStore, mateElem);
            // set the new mate pair ID in the read element
            appendRead(fragments, "", matePairId);
            SEQAN_ASSERT(getIdByName(value(store.readNames), readName, record.readId, readNameCache));
          }
        } else if (mateNo != -1) {
          // Handle case if we know this read's mate but not the read itself.
          size_t matePairId = fragments.readStore[record.readId].matePairId;
          SEQAN_ASSERT_NEQ(matePairId, TMatePairElement::INVALID_ID);
          record.readId = fragments.matePairStore[matePairId].readId[mateNo];
          if (record.readId == TMatePairElement::INVALID_ID) {
            // create new entry in read and read name store
            // set sequence and mate pair ID in new read store element
            record.readId = appendRead(fragments, "", matePairId);
            // add the identifier to the read name store
            appendName(fragments.readNameStore, readName, fragments.readNameStoreCache);
            // set the ID in the mate pair store
            fragments.matePairStore[matePairId].readId[mateNo] = record.readId;
          }
        }
        // Handle case of mate-paired read.
        //std::cerr << "-  " << readName << ", " << mateNo << std::endl;
        //std::cerr << "-  " << fragments.readNameStore[record.readId] << "/" << getMateNo(fragments, record.readId) << std::endl;
        if (mateNo != -1 && getMateNo(fragments, record.readId) != mateNo)
          record.readId = fragments.matePairStore[fragments.readStore[record.readId].matePairId].readId[mateNo];
        //std::cerr << "+  " << readName << ", " << mateNo << std::endl;
        //std::cerr << "+  " << fragments.readNameStore[record.readId] << "/" << getMateNo(fragments, record.readId) << std::endl;
        // end of "handle case of mate-paired read"
        record.distance = distance;
        if (!getIdByName((store.contigNames), contigName, record.contigId, contigNameCache)) {
          record.contigId = length(value(store.contigNames));
          appendName(value(store.contigNames), contigName, contigNameCache);
        }
        record.isForward = isForward;
        record.firstPos = firstPos;
        record.lastPos = lastPos;
        SEQAN_ASSERT_LEQ(record.firstPos, record.lastPos);
        /*int id = */appendValue(store, record);
        //std::cerr << "   record.id == " << id << std::endl;
    }
}


void
performIntervalScoreLowering(WitStore & witStore, unsigned const maxError)
{
  if (empty(witStore.intervals))
    return;
  // TODO handle different reads
  // TODO handle different contigs
  //std::cout << "performINtervalScoreLowering(witStore, " << maxError << ")" << std::endl;
  typedef WitStore::TIntervalStore TIntervalStore;
  typedef Iterator<TIntervalStore>::Type TIterator;
  typedef IntervalAndCargo<size_t, size_t> TInterval;

  // Step 1: Adjust distances.
  //std::cout << "Step 1" << std::endl;
  sortWitRecords(witStore, SortFirstPos());
  sortWitRecords(witStore, SortContigId());
  sortWitRecords(witStore, SortReadId());

  IntervalOfReadOnContig sentinel(back(witStore.intervals));
  sentinel.firstPos = -1;
  sentinel.lastPos = -1;
  sentinel.id = -1;
  appendValue(witStore.intervals, sentinel);

  //std::cout << witStore;

  String<TInterval> openIntervals;
  size_t i = 0;
  for (TIterator it = begin(witStore.intervals, Standard()), itend = end(witStore.intervals, Standard()); it != itend; ++it, ++i) {
    //std::cout << "processing " << *it << " , length(openIntervals) == " << length(openIntervals) << std::endl;
    // Remove non-overlapping intervals on top of openInterval stack.
    unsigned count = 0;
    for (unsigned j = 0; j < length(openIntervals); ++j) {
      unsigned idx = length(openIntervals) - 1 - j;
      //std::cout << rightBoundary(openIntervals[idx]) << " <= " << it->firstPos << std::endl;
      IntervalOfReadOnContig const & thisInterval = witStore.intervals[cargo(openIntervals[idx])];
      if (thisInterval.readId != it->readId || thisInterval.contigId != it->contigId || thisInterval.lastPos < it->firstPos) {
        count += 1;
        //std::cout << "popping " << witStore.intervals[cargo(openIntervals[idx])] << std::endl;
      }
    }
    resize(openIntervals, length(openIntervals) - count);

    // Perform distance lowering for containing intervals.
    for (unsigned j = 0; j < length(openIntervals); ++j) {
      unsigned idx = length(openIntervals) - 1 - j;
      unsigned id = cargo(openIntervals[idx]);
      if (witStore.intervals[id].distance <= maxError) {
        witStore.intervals[id].distance = _min(witStore.intervals[id].distance, it->distance);
        //std::cout << "witStore.intervals[" << id << "].distance = _min(" << witStore.intervals[id].distance << ", " << it->distance << ");" << std::endl;
      } else {
        //std::cout << "break;" << std::endl;
        break;  // All containing intervals must have a greater distance.
      }
    }

    // Add interval to the stack of intervals.
    appendValue(openIntervals, TInterval(it->firstPos, it->lastPos + 1, i));
    //std::cout << "pushing " << witStore.intervals[i] << std::endl;
  }
  //std::cout << witStore;

  // Step 2: Filter out intervals that are contained in intervals of lesser/equal distance.
  String<IntervalOfReadOnContig> filteredIntervals;
  clear(openIntervals);
  i = 0;
  //std::cout << "Step 2" << std::endl;
  for (TIterator it = begin(witStore.intervals, Standard()), itend = end(witStore.intervals, Standard()); it != itend; ++it, ++i) {
    //std::cout << "processing " << *it << std::endl;
    // Remove non-overlapping intervals on top of openInterval stack, appending to filtered intervals
    unsigned count = 0;
    for (unsigned j = 0; j < length(openIntervals); ++j) {
      unsigned idx = length(openIntervals) - 1 - j;
      IntervalOfReadOnContig const & thisInterval = witStore.intervals[cargo(openIntervals[idx])];
      if (thisInterval.readId != it->readId || thisInterval.contigId != it->contigId || thisInterval.lastPos < it->firstPos) {
        count += 1;
        unsigned startDistance = witStore.intervals[cargo(openIntervals[idx])].distance;
        if (!empty(filteredIntervals)) {
          if (back(filteredIntervals).lastPos >= leftBoundary(openIntervals[idx])) {
            if (back(filteredIntervals).readId == witStore.intervals[cargo(openIntervals[idx])].readId &&
                back(filteredIntervals).contigId == witStore.intervals[cargo(openIntervals[idx])].contigId) {
              // Assert current containing already written out.
              SEQAN_ASSERT_GEQ(back(filteredIntervals).firstPos, leftBoundary(openIntervals[idx]));
              SEQAN_ASSERT_LEQ(back(filteredIntervals).lastPos + 1, rightBoundary(openIntervals[idx]));
              // Get start distance.
              startDistance = back(filteredIntervals).distance + 1;
            }
          }
        }
        unsigned upperLimit = maxError;
        if (maxError < startDistance)
          upperLimit = startDistance;
        //std::cout << "startDistance = " << startDistance << ", maxError = " << upperLimit << std::endl;
        for (unsigned i = startDistance; i <= upperLimit; ++i) {
          appendValue(filteredIntervals, witStore.intervals[cargo(openIntervals[idx])]);
          back(filteredIntervals).id = length(filteredIntervals);
          back(filteredIntervals).distance = i;
          //std::cout << "appended " << back(filteredIntervals) << std::endl;
        }
      }
    }
    resize(openIntervals, length(openIntervals) - count);

    // Add interval to the stack of intervals.
    if (empty(openIntervals) || witStore.intervals[cargo(back(openIntervals))].distance > it->distance) {
      //std::cout << "appending TInterval(" << it->firstPos << ", " << it->lastPos + 1 << ", " << i << std::endl;
      appendValue(openIntervals, TInterval(it->firstPos, it->lastPos + 1, i));
    }
  }
  move(witStore.intervals, filteredIntervals);

  //std::cout << witStore;
}


template <typename TStream>
void writeWitFile(TStream & stream, WitStore const & witStore) {
    writeWitHeader(stream);
    writeWitComment(stream , WIT_COLUMN_NAMES);
    typedef typename WitStore::TIntervalStore TIntervalStore;
    typedef typename Iterator<TIntervalStore, Standard>::Type TIterator;
    for (TIterator it = begin(witStore.intervals, Standard()); it != end(witStore.intervals, Standard()); ++it) {
        for (unsigned i = 0; i < length(value(witStore.readNames)[it->readId]); ++i) {
            char c = value(witStore.readNames)[it->readId][i];
            if (c == ' ' || c == '\t' || c == '\n' || c == '\r')
                break;
            stream << c;
        }
        stream << '\t'
               << it->distance << '\t'
               << value(witStore.contigNames)[it->contigId] << '\t'
               << (it->isForward ? 'F' : 'R') << '\t'
               << it->firstPos << '\t'
               << it->lastPos << std::endl;
    }
}

#endif  // BENCHMARKS_READ_MAPPERS_APPS_WIT_STORE_H_
