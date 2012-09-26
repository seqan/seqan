// FRAGMENT(header)
#include <cstdio>
#include <fstream>
#if SEQAN_HAS_ZLIB
#include <zlib.h>
#endif  // #if SEQAN_HAS_ZLIB
#if SEQAN_HAS_BZIP2
#include <bzlib.h>
#endif  // #if SEQAN_HAS_BZIP2

#include <seqan/basic.h>
#include <seqan/stream.h>

// FRAGMENT(open-gz)
int openGz(char const * filename)
{
#if SEQAN_HAS_ZLIB
    seqan::Stream<seqan::GZFile> f;
    if (!open(f, filename, "rb"))
    {
        std::cerr << "ERROR: GZip file has the wrong format!" << std::endl;
        return 1;
    }
    
    // Totally inefficient char-wise writing of characters from .gz file to stderr.
    while (!streamEof(f))
    {
        char c = '\0';
        int res = streamReadChar(c, f);
        if (res != 0)
        {
            std::cerr << "ERROR: Reading byte from GZip file." << std::endl;
            return 1;
        }
        std::cout << c;
    }
#else  // #if SEQAN_HAS_ZLIB
    (void) filename;
    std::cerr << "ZLIB not available!" << std::endl;
#endif  // #if SEQAN_HAS_ZLIB
    return 0;
}

// FRAGMENT(open-bz2)
int openBz2(char const * filename)
{
#if SEQAN_HAS_BZIP2
    seqan::Stream<seqan::BZ2File> f;
    if (!open(f, filename, "rb"))
    {
        std::cerr << "ERROR: BZ2 file has the wrong format!" << std::endl;
        return 1;
    }

    // Totally inefficient char-wise writing of characters from .bz2 file to stderr.
    while (!streamEof(f))
    {
        char c = '\0';
        int res = streamReadChar(c, f);
        if (res != 0)
        {
            std::cout << "ERROR: Reading byte from BZ2 file." << std::endl;
            return 1;
        }
        std::cerr << c;
    }
#else  // #if SEQAN_HAS_BZIP2
    (void) filename;
    std::cerr << "BZLIB not available!" << std::endl;
#endif  // #if SEQAN_HAS_BZIP2
    return 0;
}

// FRAGMENT(main)
int main(int argc, char const ** argv)
{
    if (argc != 2)
        return 1;
    openGz(argv[1]);
    openBz2(argv[1]);
    return 0;
}
