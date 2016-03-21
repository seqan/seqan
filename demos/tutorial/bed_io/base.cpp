//![bedRecord]
#include <seqan/basic.h>
#include <seqan/bed_io.h>

using namespace seqan;

class BedRecord
{
public:
   CharString ref;      // reference name
   __int32 rID;         // index in sequenceNames of BedFile
   __int32 beginPos;    // begin position of the interval
   __int32 endPos;      // end position of the interval
   CharString name;     // name of the interval
   CharString score;    // score of the interval
   char strand;         // strand of the interval

   __int32 thickBegin;  // begin position for drawing thickly
   __int32 thickEnd;    // end position for drawing thickly
   BedRgb itemRgb;      // color for the item
   __int32 blockCount;  // number of blocks/exons
   String<__int32> blockSizes;   // block sizes
   String<__int32> blockBegins;  // block begin positions

   CharString data;    // any data not fitting into other members

   // Constants for marking reference id and position as invalid.
   static const __int32 INVALID_REFID = -1;
   static const __int32 INVALID_POS = -1;
};
//![bedRecord]

int main()
{
	return 0;
}

