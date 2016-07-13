//![bedRecord]
#include <seqan/basic.h>
#include <seqan/bed_io.h>

using namespace seqan;

class BedRecord
{
public:
   CharString ref;      // reference name
   int32_t rID;         // index in sequenceNames of BedFile
   int32_t beginPos;    // begin position of the interval
   int32_t endPos;      // end position of the interval
   CharString name;     // name of the interval
   CharString score;    // score of the interval
   char strand;         // strand of the interval

   int32_t thickBegin;  // begin position for drawing thickly
   int32_t thickEnd;    // end position for drawing thickly
   BedRgb itemRgb;      // color for the item
   int32_t blockCount;  // number of blocks/exons
   String<int32_t> blockSizes;   // block sizes
   String<int32_t> blockBegins;  // block begin positions

   CharString data;    // any data not fitting into other members

   // Constants for marking reference id and position as invalid.
   static const int32_t INVALID_REFID = -1;
   static const int32_t INVALID_POS = -1;
};
//![bedRecord]

int main()
{
	return 0;
}

