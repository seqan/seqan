#include <iostream>
#include <seqan/stream.h>
#include <seqan/bam_io.h>

using namespace seqan;

int main()
{
    CharString bamStr, samStr = "AA:Z:value1\tAB:Z:value2\tAC:i:30";
    assignTagsSamToBam(bamStr, samStr);
    BamTagsDict tags(bamStr);
    std::cout << length(tags) << std::endl;  // #=> "3"
    for (unsigned id = 0; id < length(tags); ++id)
    {
        std::cout << getTagKey(tags, id) << " -> ";

        if (getTagType(tags, id) == 'i')  // is 32 bit integer
        {
            int32_t x = 0;
            if (!extractTagValue(x, tags, id))
                SEQAN_ASSERT_FAIL("Not a valid integer at pos %u!", id);
            std::cout << x;
        }
        if (getTagType(tags, id) == 'Z')  // is string
        {
            CharString str;
            if (!extractTagValue(str, tags, id))
                SEQAN_ASSERT_FAIL("Not a valid string at pos %u!", id);
            std::cout << '"' << str << '"';
        }

        std::cout << std::endl;
    }

    return 0;
}
