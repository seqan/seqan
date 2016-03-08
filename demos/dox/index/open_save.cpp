#include <seqan/index.h>

using namespace seqan;

int main()
{
    typedef StringSet<String<char> >    TText;
    typedef Index<TText>                TIndex;

    TText text;
    appendValue(text, "MISSISSIPPI");
    appendValue(text, "MYMISSISAHAPPY");

    TIndex saveIndex(text);

    // Because indices are build on demand we fore the index creation here.
    indexCreate(saveIndex, FibreSA());

    const char * tempFileName = SEQAN_TEMP_FILENAME();
    std::cout << save(saveIndex, tempFileName) << std::endl;

    // In a different program
    TIndex openIndex;
    std::cout << open(openIndex, tempFileName) << std::endl;

    // clean up tmp directory (internal function)
    // only important when running seqan unit tests
    ClassTest::_deleteTempFile(ClassTest::_stripFileName(tempFileName));

    return 0;
}
