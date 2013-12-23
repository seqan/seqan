#include <iostream>
#include <seqan/file.h>
#include <seqan/bam_io.h>

using namespace seqan;

int main()
{
	CharString samStr = "AA:Z:value1\tAB:Z:value2\tAC:i:30";
	CharString bamStr;
	assignTagsSamToBam(bamStr, samStr);
	BamTagsDict tags(bamStr);
	std::cerr << length(tags) << std::endl;  // #=> "3"
	for (unsigned i = 0; i < length(tags); ++i)
	{
	    std::cerr << getTagKey(tags, i) << " -> " << getTagValue(tags, i) << std::endl;
	    if (getTagValue(tags, i)[0] == 'i')  // is 32 bit integer
	    {
	        __int32 x = 0;
	        bool res = extractTagValue(x, tags, i);
	        SEQAN_ASSERT_MSG(res, "Not a valid integer at pos %u!", i);
	        std::cerr << "     " << x << std::endl;
	    }
	}

    return 0;
}

