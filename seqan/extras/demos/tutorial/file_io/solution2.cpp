#include <iostream>
#include <fstream>

#include <seqan/stream.h>

#if SEQAN_HAS_ZLIB && SEQAN_HAS_BZIP2  // Guard against either not being installed.

// Copy from stream in to the stream out.

template <typename TInStream, typename TOutStream>
int copyStream(TInStream & in, TOutStream & out)
{
    seqan::CharString buffer;
    resize(buffer, 1000);

    while (!seqan::streamEof(in) && (seqan::streamError(in) == 0))
    {
        int num = seqan::streamReadBlock(&buffer[0], in, length(buffer));
        seqan::streamWriteBlock(out, &buffer[0], num);
    }
    
    return 0;
}

// The main function parses the command line, opens the files in the
// appropriate modes with the appropriate stream types and then calls either
// copyStream.

int main(int argc, char const ** argv)
{
    if (argc != 5)
    {
        std::cerr << "USAGE: " << argv[0] << " [gz|bz2] [c|x] FILE_IN FILE_OUT\n";
        return 1;
    }

    // Check first argument.
    if (seqan::CharString(argv[1]) != "gz" && seqan::CharString(argv[1]) != "bz2")
    {
        std::cerr << "ERROR: " << argv[1] << " is not a valid compression format.\n";
        return 1;
    }
    bool useGzip = (seqan::CharString(argv[1]) == "gz");

    // Check second argument.
    if (seqan::CharString(argv[2]) != "c" && seqan::CharString(argv[2]) != "x")
    {
        std::cerr << "ERROR: " << argv[2] << " is not a valid operation name.\n";
        return 1;
    }
    bool doCompress = (seqan::CharString(argv[2]) == "c");

    // Branches for stream and operation type.
    int res = 0;
    if (useGzip)
    {
        seqan::Stream<seqan::GZFile> gzFileStream;
        std::fstream fileStream;

        if (doCompress)
        {
            fileStream.open(argv[3], std::ios::binary | std::ios::in);
            if (!fileStream.good())
            {
                std::cerr << "ERROR: Could not open file " << argv[3] << "\n";
                return 1;
            }

            if (!open(gzFileStream, argv[4], "w"))
            {
                std::cerr << "ERROR: Could not open file " << argv[4] << "\n";
                return 1;
            }

            res = copyStream(fileStream, gzFileStream);
        }
        else  // extract
        {
            if (!open(gzFileStream, argv[3], "r"))
            {
                std::cerr << "ERROR: Could not open file " << argv[3] << "\n";
                return 1;
            }

            fileStream.open(argv[4], std::ios::binary | std::ios::out);
            if (!fileStream.good())
            {
                std::cerr << "ERROR: Could not open file " << argv[4] << "\n";
                return 1;
            }

            res = copyStream(gzFileStream, fileStream);
        }
    }
    else  // bz2
    {
        seqan::Stream<seqan::BZ2File> bz2FileStream;
        std::fstream fileStream;

        if (doCompress)
        {
            fileStream.open(argv[3], std::ios::binary | std::ios::in);
            if (!fileStream.good())
            {
                std::cerr << "ERROR: Could not open file " << argv[3] << "\n";
                return 1;
            }

            if (!open(bz2FileStream, argv[4], "w"))
            {
                std::cerr << "ERROR: Could not open file " << argv[4] << "\n";
                return 1;
            }

            res = copyStream(fileStream, bz2FileStream);
        }
        else  // extract
        {
            if (!open(bz2FileStream, argv[3], "r"))
            {
                std::cerr << "ERROR: Could not open file " << argv[3] << "\n";
                return 1;
            }

            fileStream.open(argv[4], std::ios::binary | std::ios::out);
            if (!fileStream.good())
            {
                std::cerr << "ERROR: Could not open file " << argv[4] << "\n";
                return 1;
            }

            res = copyStream(bz2FileStream, fileStream);
        }
    }

    if (res != 0)
        std::cerr << "ERROR: There was an error reading/writing!\n";
    return res;
}

#else  // #if SEQAN_HAS_ZLIB && SEQAN_HAS_BZIP2

int main()
{
    return 0;
}

#endif  // #if SEQAN_HAS_ZLIB && SEQAN_HAS_BZIP2
