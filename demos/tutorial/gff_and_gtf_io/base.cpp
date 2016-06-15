#include <seqan/gff_io.h>

//![GffRecord]
using namespace seqan;

class GffRecord
{
public:
    CharString ref;      // reference name
    int32_t rID;         // index in sequenceNames of GffFile
    CharString source;   // source free text descriptor
    CharString type;     // type of the feature
    int32_t beginPos;    // begin position of the interval
    int32_t endPos;      // end position of the interval
    float score;         // score of the annotation
    char strand;         // the strand
    char phase;          // one of '0', '1', '2', and '.'

    // The key/value list, split into a list of keys and values.
    StringSet<CharString> tagNames;
    StringSet<CharString> tagValues;

    // Returns float value for an invalid score.
    static float INVALID_SCORE();

    // Constants for marking reference id and position as invalid.
    static const int32_t INVALID_IDX = -1;
    static const int32_t INVALID_POS = -1;
};
//![GffRecord]

int main()
{
	return 0;
}
