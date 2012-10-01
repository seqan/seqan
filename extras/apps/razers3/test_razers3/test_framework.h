// test_framework.h

#ifndef SEQAN_FIND_TEST_FRAMEWORK_H_
#define SEQAN_FIND_TEST_FRAMEWORK_H_

namespace seqan {
namespace ClassTest {
namespace override {

inline
void beginTestSuite(const char * testSuiteName, const char * argv0)
{
    using::seqan::ClassTest::StaticData;

    // First things first: Print the current debug level.
    printDebugLevel(std::cout);
    (void)testSuiteName;
    StaticData::testCount() = 0;
    StaticData::skippedCount() = 0;
    StaticData::errorCount() = 0;
    StaticData::totalCheckPointCount() = 0;
    StaticData::foundCheckPointCount() = 0;
    (void) argv0;

    char const * path = "../../..";
    size_t const pathlen = std::strlen(path);
    char * pathbuffer = new char[pathlen + 1];
    std::strncpy(pathbuffer, path, pathlen);
    StaticData::pathToProjects() = pathbuffer;

    // Get path to argv0.
    const char * end = argv0;
#ifdef PLATFORM_WINDOWS
    const char pathSeparator = '\\';
#else  // PLATFORM_WINDOWS
    const char pathSeparator = '/';
#endif  // PLATFORM_WINDOWS
    for (const char * ptr = strchr(argv0, pathSeparator); ptr != 0; ptr = strchr(ptr + 1, pathSeparator))
        end = ptr;
    int rpos = end - argv0;
    if (rpos <= 0)
    {
        StaticData::basePath() = new char[1];
        strcpy(StaticData::basePath(), ".");
    }
    else
    {
        int len = rpos;
        StaticData::basePath() = new char[len];
        strncpy(StaticData::basePath(), argv0, len);
    }
//        // Get path to projects.
//        const char *file = __FILE__;
//        int pos = -1;
//        for (size_t i = 0; i < strlen(file) - strlen("projects"); ++i) {
//            if (strncmp(file + i, "projects", strlen("projects")) == 0) {
//                pos = i;
//            }
//        }
//        if (pos == -1) {
//            std::cerr << "Could not extrapolate path to projects from __FILE__ == \""
//                      << __FILE__ << "\"" << std::endl;
//            exit(1);
//        }
//        StaticData::pathToProjects() = new char[pos];
//        strncpy(StaticData::pathToProjects(), file, pos);
//        StaticData::pathToProjects()[pos-1] = '\0';
}

}
}
}     // namespace seqan::ClassTest::override

#undef SEQAN_BEGIN_TESTSUITE
#define SEQAN_BEGIN_TESTSUITE(suite_name)                       \
    int main(int argc, char ** argv) {                           \
        (void) argc;                                                \
        ::seqan::ClassTest::override::beginTestSuite(# suite_name, argv[0]);

#endif // ndef SEQAN_FIND_TEST_FRAMEWORK_H_
